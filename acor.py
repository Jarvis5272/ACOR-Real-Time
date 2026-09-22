#!/usr/bin/env python3
"""Minimal public runner for ACOR streaming reconstruction."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import platform
import random
import re
import secrets
import signal
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
import uuid
from collections import Counter
from collections.abc import Callable
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RUNNER = ROOT / "bin" / "acor_runner"
ED_EVALUATOR = ROOT / "bin" / "ed_pairs"
DEFAULT_KMER_SIZE = 9
ALGORITHM_MODE = "NO_STOP"
ENGINE_MODE = "LEDGER_ADAPTIVE_POOL"
MAX_DESIGN_LENGTH = 1_000_000
VERSION = "0.2.0"
READS_HEADER = [
    "dataset_id", "cluster_id", "read_id", "read_sequence",
    "original_read_sequence", "read_quality", "design_length",
]


@dataclass(frozen=True)
class FileIdentity:
    device: int
    inode: int
    size: int
    mtime_ns: int
    ctime_ns: int
    sha256: str


@dataclass(frozen=True)
class TruthSnapshot:
    path: Path
    sequences: dict[str, str]
    identity: FileIdentity
    row_count: int



class ThreadBar:
    def __init__(self, workers: int) -> None:
        self.workers = workers
        self.percent = 0.0
        self.status = "waiting"
        self.pulse: int | None = None


class Dashboard:
    """Persistent one-row-per-thread terminal dashboard."""

    BAR_WIDTH = 32
    FRAME_SECONDS = 0.08

    def __init__(self, dataset: str, thread_counts: list[int]) -> None:
        self.dataset = dataset
        self.order = thread_counts
        self.states = {workers: ThreadBar(workers) for workers in thread_counts}
        self.enabled = sys.stdout.isatty()
        self.visible = False
        self.last_draw = 0.0

    def _line(self, state: ThreadBar) -> str:
        width = self.BAR_WIDTH
        filled = min(width, int(width * state.percent / 100.0))
        cells = ["━"] * filled + ["·"] * (width - filled)
        if state.pulse is not None and filled < width:
            cells[filled + state.pulse % (width - filled)] = "◆"
        bar = "".join(cells)
        if self.enabled:
            label = f"\033[1;36mt{state.workers:<3}\033[0m"
            bar = f"\033[38;5;42m{bar[:filled]}\033[38;5;45m{bar[filled:]}\033[0m"
            percent = f"\033[1m{state.percent:6.2f}%\033[0m"
            status = f"\033[2m{state.status}\033[0m"
        else:
            label = f"t{state.workers:<3}"
            percent = f"{state.percent:6.2f}%"
            status = state.status
        return f"  {label}  [{bar}]  {percent}  {status}"

    def start(self) -> None:
        if not self.enabled:
            return
        sys.stdout.write("\033[?25l")
        for workers in self.order:
            sys.stdout.write(self._line(self.states[workers]) + "\n")
        sys.stdout.flush()
        self.visible = True

    def _render(self, force: bool = False) -> None:
        if not self.enabled or not self.visible:
            return
        now = time.monotonic()
        if not force and now - self.last_draw < self.FRAME_SECONDS:
            return
        sys.stdout.write(f"\033[{len(self.order)}A")
        for workers in self.order:
            sys.stdout.write("\r\033[2K" + self._line(self.states[workers]) + "\n")
        sys.stdout.flush()
        self.last_draw = now

    def update(
        self,
        workers: int,
        percent: float,
        status: str,
        pulse: int | None = None,
        force: bool = False,
    ) -> None:
        state = self.states[workers]
        state.percent = max(0.0, min(100.0, percent))
        state.status = status
        state.pulse = pulse
        self._render(force=force)

    def pulse(self, workers: int, status: str, step: int) -> None:
        state = self.states[workers]
        self.update(workers, state.percent, status, pulse=step)

    def log(self, lines: list[str]) -> None:
        if not self.enabled:
            for line in lines:
                print(line, flush=True)
            return
        if self.visible:
            sys.stdout.write(f"\033[{len(self.order)}A\r\033[J")
        for line in lines:
            sys.stdout.write(line + "\n")
        for workers in self.order:
            sys.stdout.write(self._line(self.states[workers]) + "\n")
        sys.stdout.flush()
        self.visible = True
        self.last_draw = time.monotonic()

    def close(self) -> None:
        if self.enabled:
            self._render(force=True)
            sys.stdout.write("\033[?25h")
            sys.stdout.flush()
        self.visible = False


def parse_list(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def parse_threads(value: str) -> list[int]:
    result = [int(part) for part in parse_list(value)]
    if not result or any(count < 1 for count in result):
        raise ValueError("--threads must contain positive integers")
    return list(dict.fromkeys(result))


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    return cleaned or "run"


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def identity_from_stat(stat_result: os.stat_result, sha256: str) -> FileIdentity:
    return FileIdentity(
        device=stat_result.st_dev,
        inode=stat_result.st_ino,
        size=stat_result.st_size,
        mtime_ns=stat_result.st_mtime_ns,
        ctime_ns=stat_result.st_ctime_ns,
        sha256=sha256,
    )


def verify_file_identity(path: Path, expected: FileIdentity, verify_sha: bool = False) -> None:
    observed_stat = path.stat()
    observed = identity_from_stat(observed_stat, expected.sha256)
    for field in ("device", "inode", "size", "mtime_ns", "ctime_ns"):
        if getattr(observed, field) != getattr(expected, field):
            raise RuntimeError(
                f"file identity changed for {path.name}: {field} "
                f"{getattr(expected, field)} -> {getattr(observed, field)}"
            )
    if verify_sha and sha256_path(path) != expected.sha256:
        raise RuntimeError(f"file content SHA changed for {path.name}")


def preflight_binary(path: Path, label: str) -> dict[str, object]:
    if not path.is_file():
        raise RuntimeError(f"missing {label}: {path.name}")
    if not os.access(path, os.X_OK):
        raise RuntimeError(f"{label} is not executable: {path.name}")
    with path.open("rb") as handle:
        header = handle.read(20)
    if len(header) < 20 or header[:4] != b"\x7fELF" or header[4] != 2:
        raise RuntimeError(f"{label} is not a 64-bit ELF executable")
    machine = int.from_bytes(header[18:20], "little")
    if machine != 62:
        raise RuntimeError(f"{label} is not x86_64 ELF (e_machine={machine})")
    return {
        "name": path.name,
        "sha256": sha256_path(path),
        "executable": True,
        "elf_class": "ELF64",
        "architecture": "x86_64",
    }


def inspect_reads(
    source: Path,
    length_override: int | None,
    effective_path: Path | None,
    minimum_design_length: int = DEFAULT_KMER_SIZE,
) -> tuple[str, list[tuple[str, int]], Path, int, dict[str, object]]:
    if length_override is not None and not (
        minimum_design_length <= length_override <= MAX_DESIGN_LENGTH
    ):
        raise ValueError(
            f"--length must be between {minimum_design_length} and {MAX_DESIGN_LENGTH}"
        )
    cluster_sizes: list[tuple[str, int]] = []
    lengths: Counter[int] = Counter()
    dataset_ids: set[str] = set()
    cluster_ids: set[str] = set()
    completed_clusters: set[str] = set()
    current_read_ids: set[str] = set()
    current = None
    current_count = 0
    current_design_length: int | None = None
    reads_count = 0
    before = source.stat()
    digest = hashlib.sha256()

    destination = None
    if effective_path is not None:
        effective_path.parent.mkdir(parents=True, exist_ok=True)
        destination = effective_path.open("w", encoding="utf-8", newline="")

    try:
        with source.open("rb") as handle:
            raw_header = handle.readline()
            digest.update(raw_header)
            header_line = raw_header.decode("utf-8")
            header = header_line.rstrip("\n\r").split("\t")
            if header != READS_HEADER:
                raise ValueError("reads.tsv has an unsupported header")
            indexes = {name: header.index(name) for name in header}
            if destination is not None:
                destination.write("\t".join(header) + "\n")
            for source_row, raw_line in enumerate(handle, 2):
                digest.update(raw_line)
                line = raw_line.decode("utf-8")
                fields = line.rstrip("\n\r").split("\t")
                if len(fields) != len(READS_HEADER):
                    raise ValueError(
                        f"reads.tsv row {source_row} must contain exactly 7 columns"
                    )
                cid = fields[indexes["cluster_id"]]
                dataset_id = fields[indexes["dataset_id"]]
                read_id = fields[indexes["read_id"]]
                read_sequence = fields[indexes["read_sequence"]]
                if not dataset_id or not cid or not read_id or not read_sequence:
                    raise ValueError(
                        f"reads.tsv row {source_row} has an empty required field"
                    )
                dataset_ids.add(dataset_id)
                try:
                    design_length = int(fields[indexes["design_length"]])
                except ValueError as error:
                    raise ValueError(
                        f"invalid design_length at row {source_row}"
                    ) from error
                if not (minimum_design_length <= design_length <= MAX_DESIGN_LENGTH):
                    raise ValueError(
                        f"design_length out of safe range at row {source_row}: {design_length}"
                    )
                lengths[design_length] += 1
                if current is None:
                    current = cid
                    current_design_length = design_length
                elif cid != current:
                    cluster_sizes.append((current, current_count))
                    completed_clusters.add(current)
                    if cid in completed_clusters:
                        raise ValueError(f"duplicate/non-contiguous input cluster_id: {cid}")
                    current = cid
                    current_count = 0
                    current_design_length = design_length
                    current_read_ids.clear()
                elif design_length != current_design_length:
                    raise ValueError(
                        f"design_length changed within cluster {cid} at row {source_row}"
                    )
                if read_id in current_read_ids:
                    raise ValueError(f"duplicate read_id in cluster {cid}: {read_id}")
                current_read_ids.add(read_id)
                cluster_ids.add(cid)
                current_count += 1
                reads_count += 1
                if destination is not None:
                    fields[indexes["design_length"]] = str(length_override)
                    destination.write("\t".join(fields) + "\n")
            inside = os.fstat(handle.fileno())
        if current is not None:
            cluster_sizes.append((current, current_count))
    finally:
        if destination is not None:
            destination.close()

    if len(dataset_ids) != 1 or not cluster_sizes:
        raise ValueError("reads.tsv must contain one non-empty dataset")
    after = source.stat()
    for name, observed in (("open", inside), ("after", after)):
        for field in ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns"):
            if getattr(observed, field) != getattr(before, field):
                raise RuntimeError(f"reads.tsv changed during scan ({name}:{field})")
    dataset = next(iter(dataset_ids))
    inferred_length = lengths.most_common(1)[0][0]
    source_identity = identity_from_stat(after, digest.hexdigest())
    runner_reads = effective_path or source
    runner_identity = (
        identity_from_stat(runner_reads.stat(), sha256_path(runner_reads))
        if effective_path is not None
        else source_identity
    )
    audit = {
        "input_rows": reads_count,
        "input_unique_clusters": len(cluster_ids),
        "input_ids": cluster_ids,
        "duplicate_input": 0,
        "source_identity": source_identity,
        "runner_identity": runner_identity,
    }
    return dataset, cluster_sizes, effective_path or source, inferred_length, audit


def labels(seeds: list[str] | None, rounds: int) -> list[tuple[str, int]]:
    if seeds is not None:
        return [(seed, index) for index, seed in enumerate(seeds, 1)]
    return [(secrets.token_hex(16), index) for index in range(1, rounds + 1)]


def deterministic_entropy(seed: str, dataset: str, cluster_id: str) -> int:
    key = f"acor|minimal|{seed}|{dataset}|{cluster_id}"
    return int.from_bytes(hashlib.sha256(key.encode("utf-8")).digest()[:16], "big")



def write_order(
    path: Path,
    dataset: str,
    seed: str,
    cluster_sizes: list[tuple[str, int]],
    on_progress: Callable[[float], None],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    raw = temporary.open("wb")
    compressed = gzip.GzipFile(filename="", fileobj=raw, mode="wb", compresslevel=1, mtime=0)
    output = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    total = len(cluster_sizes)
    try:
        for index, (cluster_id, read_count) in enumerate(cluster_sizes, 1):
            entropy = deterministic_entropy(seed, dataset, cluster_id)
            permutation = list(range(read_count))
            random.Random(entropy).shuffle(permutation)
            output.write(json.dumps({
                "cluster_id": cluster_id,
                "entropy_bits": f"{entropy:032x}",
                "permutation": permutation,
            }, separators=(",", ":"), sort_keys=True) + "\n")
            on_progress(index / total)
    finally:
        output.close()
        compressed.close()
        raw.close()
    os.replace(temporary, path)


def _physical_core_key(cpu: int) -> tuple[int, int]:
    topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
    try:
        package = int((topology / "physical_package_id").read_text().strip())
        core = int((topology / "core_id").read_text().strip())
        return package, core
    except (OSError, ValueError):
        return 0, cpu


def select_cpus(workers: int) -> tuple[list[int], int]:
    available = sorted(os.sched_getaffinity(0))
    if len(available) < workers + 1:
        raise ValueError(f"threads={workers} needs {workers + 1} available CPUs")
    physical: list[int] = []
    observed_cores: set[tuple[int, int]] = set()
    for cpu in available:
        key = _physical_core_key(cpu)
        if key not in observed_cores:
            observed_cores.add(key)
            physical.append(cpu)
    if len(physical) < workers:
        raise ValueError(
            f"threads={workers} needs {workers} distinct physical cores; "
            f"only {len(physical)} are available"
        )
    worker_cpus = physical[:workers]
    control_cpu = next(cpu for cpu in available if cpu not in set(worker_cpus))
    return worker_cpus, control_cpu


def terminate_process_group(process: subprocess.Popen, grace_seconds: float = 5.0) -> int:
    if process.poll() is not None:
        return process.wait()
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return process.wait()
    try:
        return process.wait(timeout=grace_seconds)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        return process.wait()


def log_tail(path: Path, limit: int = 16 * 1024) -> str:
    if not path.is_file():
        return ""
    with path.open("rb") as handle:
        handle.seek(0, os.SEEK_END)
        size = handle.tell()
        handle.seek(max(0, size - limit))
        return handle.read().decode("utf-8", errors="replace").strip()


def atomic_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.{uuid.uuid4().hex}.tmp")
    try:
        with source.open("rb") as src, temporary.open("wb") as dst:
            shutil.copyfileobj(src, dst, length=4 * 1024 * 1024)
            dst.flush()
            os.fsync(dst.fileno())
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()



def run_engine(
    reads: Path,
    order: Path,
    engine_dir: Path,
    dataset: str,
    workers: int,
    on_pulse: Callable[[int], None],
    config: str = "kmc_k9_b16_lognormal",
    engine_mode: str = ENGINE_MODE,
) -> float:
    engine_dir.mkdir(parents=True, exist_ok=True)
    worker_cpus, control_cpu = select_cpus(workers)
    command = [
        str(RUNNER), str(reads), str(order), str(engine_dir), "ACOR_MINIMAL",
        dataset, dataset, str(workers), ",".join(map(str, worker_cpus)),
        str(control_cpu), config, engine_mode, "0.0", "0.0", "NO_PROFILE",
    ]
    atomic_json(engine_dir / "COMMAND.json", {
        "argv": command,
        "workers": workers,
        "worker_cpus": worker_cpus,
        "control_cpu": control_cpu,
        "config": config,
        "algorithm_mode": ALGORITHM_MODE,
        "engine_mode": engine_mode,
    })
    stdout_path = engine_dir / "runner.stdout.log"
    stderr_path = engine_dir / "runner.stderr.log"
    timeout_seconds = float(os.environ.get("ACOR_ENGINE_TIMEOUT_SECONDS", "0") or "0")
    if timeout_seconds < 0:
        raise ValueError("ACOR_ENGINE_TIMEOUT_SECONDS must be non-negative")
    with stdout_path.open("wb") as stdout_handle, stderr_path.open("wb") as stderr_handle:
        process = subprocess.Popen(
            command,
            stdout=stdout_handle,
            stderr=stderr_handle,
            start_new_session=True,
            env={
                **os.environ,
                "ACOR_CONFIG": config,
                "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "NUMEXPR_NUM_THREADS": "1",
            },
        )
        frame = 0
        started = time.monotonic()
        try:
            while process.poll() is None:
                if timeout_seconds and time.monotonic() - started > timeout_seconds:
                    raise TimeoutError(f"engine exceeded {timeout_seconds:g} seconds")
                on_pulse(frame)
                frame += 1
                time.sleep(Dashboard.FRAME_SECONDS)
            return_code = process.wait()
        except BaseException:
            terminate_process_group(process)
            raise
    if return_code != 0:
        raise RuntimeError(log_tail(stderr_path) or f"engine exited {return_code}")
    marker_path = engine_dir / "RUNNER_COMPLETE.json"
    if not marker_path.is_file():
        raise RuntimeError("runner exited without completion marker")
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    if marker.get("state") != "COMPLETE" or marker.get("status") != "PASS":
        raise RuntimeError(f"runner completion marker rejected: {marker}")
    with (engine_dir / "RUNNER_METRICS.tsv").open("r", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    return float(row["total_wall_seconds"])


def load_predictions(path: Path) -> dict[str, str]:
    predictions: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cluster_id = row["cluster_id"]
            if cluster_id in predictions:
                raise ValueError(f"duplicate pred cluster_id: {cluster_id}")
            predictions[cluster_id] = row["reconstructed_sequence"]
    if not predictions:
        raise ValueError("pred.tsv is empty")
    return predictions


def load_truth(path: Path) -> TruthSnapshot:
    before = path.stat()
    digest = hashlib.sha256()
    truth: dict[str, str] = {}
    row_count = 0
    with path.open("rb") as raw_handle:
        raw_header = raw_handle.readline()
        digest.update(raw_header)
        header = raw_header.decode("utf-8").rstrip("\n\r").split("\t")
        simple = header == ["cluster_id", "sequence"]
        legacy = header == ["row_key", "dataset_id", "reference_sequence"]
        if not simple and not legacy:
            raise ValueError("truth.tsv has an unsupported header")
        for source_row, raw_line in enumerate(raw_handle, 2):
            digest.update(raw_line)
            fields = raw_line.decode("utf-8").rstrip("\n\r").split("\t")
            if len(fields) != len(header):
                raise ValueError(
                    f"truth.tsv row {source_row} has an invalid column count"
                )
            if simple:
                cluster_id, sequence = fields
            else:
                row_key, dataset_id, sequence = fields
                parts = row_key.split("|")
                if len(parts) < 2 or not dataset_id:
                    raise ValueError(f"invalid truth row_key at row {source_row}")
                cluster_id = parts[1]
            if not cluster_id:
                raise ValueError(f"empty truth cluster_id at row {source_row}")
            if cluster_id in truth:
                raise ValueError(f"duplicate truth cluster_id: {cluster_id}")
            truth[cluster_id] = sequence
            row_count += 1
        inside = os.fstat(raw_handle.fileno())
    if not truth:
        raise ValueError("truth.tsv is empty")
    after = path.stat()
    for name, observed in (("open", inside), ("after", after)):
        for field in ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns"):
            if getattr(observed, field) != getattr(before, field):
                raise RuntimeError(f"truth.tsv changed during load ({name}:{field})")
    return TruthSnapshot(
        path=path,
        sequences=truth,
        identity=identity_from_stat(after, digest.hexdigest()),
        row_count=row_count,
    )


def require_equal_ids(name: str, expected: set[str], observed: set[str]) -> None:
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing or extra:
        raise ValueError(
            f"{name} cluster set mismatch: missing={missing[:10]} extra={extra[:10]} "
            f"missing_count={len(missing)} extra_count={len(extra)}"
        )


def evaluate(
    pred_path: Path,
    truth_snapshot: TruthSnapshot,
    input_ids: set[str],
    evaluation_path: Path,
) -> dict[str, object]:
    pred = load_predictions(pred_path)
    truth = truth_snapshot.sequences
    verify_file_identity(truth_snapshot.path, truth_snapshot.identity)
    require_equal_ids("pred/input", input_ids, set(pred))
    require_equal_ids("truth/input", input_ids, set(truth))
    evaluation_path.parent.mkdir(parents=True, exist_ok=True)
    token = uuid.uuid4().hex
    input_path = evaluation_path.parent / f".evaluation_input.{token}.tsv"
    raw_path = evaluation_path.parent / f".evaluation_raw.{token}.tsv"
    stderr_path = evaluation_path.parent / "evaluator.stderr.log"
    with input_path.open("w", encoding="utf-8", newline="") as handle:
        for cluster_id, sequence in pred.items():
            handle.write(f"{cluster_id}\t{sequence}\t{truth[cluster_id]}\n")
    with input_path.open("rb") as stdin_handle, raw_path.open("wb") as stdout_handle, stderr_path.open("wb") as stderr_handle:
        completed = subprocess.run(
            [str(ED_EVALUATOR)],
            stdin=stdin_handle,
            stdout=stdout_handle,
            stderr=stderr_handle,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(log_tail(stderr_path) or "ED evaluation failed")
    ed_by_cluster: dict[str, int] = {}
    evaluation_rows: list[dict[str, object]] = []
    with raw_path.open("r", encoding="utf-8", newline="") as handle:
      for output_row, line in enumerate(handle, 1):
        if line.endswith("\r\n"):
            raise RuntimeError(f"evaluator emitted CRLF at row {output_row}")
        parts = line.split("\t")
        if len(parts) != 4:
            raise RuntimeError(
                f"evaluator row {output_row} must contain exactly four columns"
            )
        cluster_id, pred_length_text, truth_length_text, ed_text = parts
        ed_text = ed_text.rstrip("\n")
        if not cluster_id or cluster_id in ed_by_cluster:
            raise RuntimeError(f"invalid/duplicate evaluator ID at row {output_row}")
        try:
            pred_length = int(pred_length_text)
            truth_length = int(truth_length_text)
            edit_distance = int(ed_text)
        except ValueError as error:
            raise RuntimeError(f"non-integer evaluator field at row {output_row}") from error
        if pred_length != len(pred.get(cluster_id, "")):
            raise RuntimeError(f"evaluator pred length mismatch for {cluster_id}")
        if truth_length != len(truth.get(cluster_id, "")):
            raise RuntimeError(f"evaluator truth length mismatch for {cluster_id}")
        if edit_distance < 0 or edit_distance > max(pred_length, truth_length):
            raise RuntimeError(f"evaluator ED out of bounds for {cluster_id}")
        ed_by_cluster[cluster_id] = edit_distance
        evaluation_rows.append({
            "cluster_id": cluster_id,
            "pred_length": pred_length,
            "truth_length": truth_length,
            "edit_distance": edit_distance,
            "exact": int(edit_distance == 0),
            "normalized_accuracy": 1.0 - edit_distance / max(1, truth_length),
        })
    require_equal_ids("evaluator/input", input_ids, set(ed_by_cluster))
    verify_file_identity(truth_snapshot.path, truth_snapshot.identity)
    input_path.unlink()
    raw_path.unlink()
    _write_rows(evaluation_path, evaluation_rows)
    eds = [ed_by_cluster[cid] for cid in pred]
    exact = sum(ed == 0 for ed in eds) / len(eds)
    accuracy = statistics.mean(
        1.0 - ed_by_cluster[cid] / max(1, len(truth[cid])) for cid in pred
    )
    return {
        "exact": exact,
        "accuracy": accuracy,
        "mean_ed": statistics.mean(eds),
        "evaluated_clusters": len(eds),
        "exact_count": sum(ed == 0 for ed in eds),
        "total_ed": sum(eds),
        "total_truth_bases": sum(len(truth[cid]) for cid in pred),
    }


def closure_fields(
    pred_path: Path,
    truth_snapshot: TruthSnapshot | None,
    input_audit: dict[str, object],
) -> dict[str, object]:
    input_ids = set(input_audit["input_ids"])
    pred = load_predictions(pred_path)
    require_equal_ids("pred/input", input_ids, set(pred))
    truth: dict[str, str] | None = None
    if truth_snapshot is not None:
        truth = truth_snapshot.sequences
        require_equal_ids("truth/input", input_ids, set(truth))
    return {
        "input_rows": input_audit["input_rows"],
        "input_unique_clusters": input_audit["input_unique_clusters"],
        "pred_rows": len(pred),
        "pred_unique_clusters": len(pred),
        "truth_rows": truth_snapshot.row_count if truth_snapshot is not None else "NA",
        "truth_unique_clusters": len(truth) if truth is not None else "NA",
        "duplicate_input": input_audit["duplicate_input"],
        "duplicate_pred": 0,
        "duplicate_truth": 0 if truth is not None else "NA",
        "missing_pred": 0,
        "extra_pred": 0,
        "missing_truth": 0 if truth is not None else "NA",
        "extra_truth": 0 if truth is not None else "NA",
        "metric_closure": "PASS",
    }


def write_summary(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def mean_value(rows: list[dict[str, object]], field: str) -> str:
    values = [float(row[field]) for row in rows if row[field] != "NA"]
    return f"{statistics.mean(values):.12g}" if values else "NA"


def mean_row(rows: list[dict[str, object]], workers: int) -> dict[str, object]:
    row = dict(rows[0])
    row.update({
        "threads": workers,
        "round": "MEAN",
        "seed": "NA",
        "exact": mean_value(rows, "exact"),
        "accuracy": mean_value(rows, "accuracy"),
        "mean_ed": mean_value(rows, "mean_ed"),
        "engine_time_s": mean_value(rows, "engine_time_s"),
        "round_time_s": mean_value(rows, "round_time_s"),
        "pred": "NA",
        "metric_closure": "PASS_ALL_ROUNDS",
    })
    return row



def run(args: argparse.Namespace) -> int:
    data = args.data.resolve()
    reads = data / "reads.tsv"
    truth = data / "truth.tsv"
    if not reads.is_file():
        raise ValueError(f"missing {reads}")
    if args.aim and not truth.is_file():
        raise ValueError("--aim requires truth.tsv")
    effective_truth = truth if truth.is_file() else None

    if args.seeds is not None and args.rounds is not None:
        raise ValueError("--seeds and --rounds cannot be used together")
    seeds = parse_list(args.seeds) if args.seeds is not None else None
    rounds = args.rounds if args.rounds is not None else 1
    if rounds < 1:
        raise ValueError("--rounds must be at least 1")
    thread_counts = parse_threads(args.threads)
    config = f"kmc_k{args.k}_b16_lognormal"

    result_root = args.results.resolve()
    result_root.mkdir(parents=True, exist_ok=True)
    session_started = time.monotonic()
    session_id = (
        datetime.now().strftime("%Y%m%dT%H%M%S")
        + "_" + uuid.uuid4().hex
    )
    final_session_root = result_root / f"acor_{session_id}"
    if final_session_root.exists():
        raise RuntimeError(f"refuse existing session: {final_session_root}")
    session_root = Path(tempfile.mkdtemp(prefix=f".{final_session_root.name}.staging.", dir=result_root))
    status_path = session_root / "RUN_STATUS.json"
    manifest_path = session_root / "RUN_MANIFEST.json"
    status: dict[str, object] = {
        "session_id": session_id,
        "state": "RUNNING",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "result_dir": final_session_root.name,
    }
    atomic_json(status_path, status)
    manifest: dict[str, object] = {
        "version": VERSION,
        "release_mode": ALGORITHM_MODE,
        "execution_model": "CAUSAL_OFFLINE_STREAM_REPLAY",
        "execution_backend": ENGINE_MODE,
        "argv": sys.argv,
        "dataset_argument": data.name,
        "seed_source": "EXPLICIT_SEEDS" if seeds is not None else "OS_RANDOM",
        "threads": thread_counts,
        "config": config,
        "k": args.k,
        "mode": ALGORITHM_MODE,
        "state": "RUNNING",
        "jobs": [],
        "environment": {
            "python": platform.python_version(),
            "os": platform.platform(),
            "machine": platform.machine(),
            "logical_cpus": os.cpu_count(),
        },
    }
    atomic_json(manifest_path, manifest)
    dashboard: Dashboard | None = None
    try:
        manifest["runner"] = preflight_binary(RUNNER, "runner")
        manifest["evaluator"] = preflight_binary(ED_EVALUATOR, "evaluator")
        effective_path = session_root / ".input" / "reads.tsv" if args.length is not None else None
        dataset, cluster_sizes, effective_reads, inferred_length, input_audit = inspect_reads(
            reads, args.length, effective_path, args.k,
        )
        chosen_length = args.length if args.length is not None else inferred_length
        run_labels = labels(seeds, rounds)
        manifest.update({
            "dataset": dataset,
            "seeds": [seed for seed, _index in run_labels],
            "reads": asdict(input_audit["source_identity"]),
            "runner_reads": asdict(input_audit["runner_identity"]),
            "input_clusters": input_audit["input_unique_clusters"],
            "input_rows": input_audit["input_rows"],
            "truth_present": effective_truth is not None,
            "truth_sha256": "PENDING_AFTER_PRED" if effective_truth is not None else "NOT_APPLICABLE",
        })
        atomic_json(manifest_path, manifest)
        status.update({
            "dataset": dataset,
            "threads": thread_counts,
            "seeds": [seed for seed, _index in run_labels],
            "reads_identity": asdict(input_audit["source_identity"]),
            "runner_reads_identity": asdict(input_audit["runner_identity"]),
            "state": "RUNNING",
        })
        atomic_json(status_path, status)
        rows: list[dict[str, object]] = []
        session_rows: list[dict[str, object]] = []
        truth_snapshot: TruthSnapshot | None = None
        dashboard = Dashboard(dataset, thread_counts)
        print(
            f"ACOR  dataset={dataset}  rounds={len(run_labels)}  "
            f"threads={','.join(map(str, thread_counts))}",
            flush=True,
        )
        dashboard.start()

        for workers in thread_counts:
            thread_rows: list[dict[str, object]] = []
            round_span = 100.0 / len(run_labels)
            for seed, round_index in run_labels:
                round_started = time.monotonic()
                base = (round_index - 1) * round_span
                relative_run_dir = (
                    Path(f"threads_{workers}")
                    / f"round_{round_index:02d}_seed_{safe_name(seed)}"
                )
                run_dir = session_root / ".round_staging" / (
                    f"threads_{workers}_round_{round_index:02d}_{uuid.uuid4().hex}"
                )
                run_dir.mkdir(parents=True)
                internal = run_dir / ".internal"
                order = (
                    session_root / ".orders"
                    / f"round_{round_index:02d}_seed_{safe_name(seed)}.jsonl.gz"
                )

                order_status = f"round {round_index}/{len(run_labels)} · ordering"
                dashboard.update(workers, base, order_status, force=True)
                if not order.is_file():
                    write_order(
                        order, dataset, seed, cluster_sizes,
                        lambda fraction: dashboard.update(
                            workers,
                            base + round_span * 0.35 * fraction,
                            order_status,
                        ),
                    )
                else:
                    dashboard.update(
                        workers, base + round_span * 0.35,
                        f"round {round_index}/{len(run_labels)} · order reused",
                    )
                job_manifest: dict[str, object] = {
                    "dataset": dataset,
                    "threads": workers,
                    "round": round_index,
                    "seed": seed,
                    "state": "RUNNING",
                    "order_sha256": sha256_path(order),
                    "order_artifact": str(
                        Path(".orders")
                        / f"round_{round_index:02d}_seed_{safe_name(seed)}.jsonl.gz"
                    ),
                    "runner_command_artifact": str(
                        relative_run_dir / ".internal/engine/COMMAND.json"
                    ),
                }
                manifest["jobs"].append(job_manifest)
                atomic_json(manifest_path, manifest)

                engine_status = f"round {round_index}/{len(run_labels)} · reconstructing"
                dashboard.update(workers, base + round_span * 0.35, engine_status, force=True)
                verify_file_identity(reads, input_audit["source_identity"])
                verify_file_identity(effective_reads, input_audit["runner_identity"])
                engine_time = run_engine(
                    effective_reads,
                    order,
                    internal / "engine",
                    dataset,
                    workers,
                    lambda frame: dashboard.pulse(workers, engine_status, frame),
                    config,
                    ENGINE_MODE,
                )
                verify_file_identity(reads, input_audit["source_identity"])
                verify_file_identity(effective_reads, input_audit["runner_identity"])
                dashboard.update(
                    workers,
                    base + round_span * 0.82,
                    f"round {round_index}/{len(run_labels)} · evaluating",
                    force=True,
                )

                source_pred = internal / "engine" / "pred.tsv"
                pred_path = run_dir / "pred.tsv"
                atomic_copy(source_pred, pred_path)
                if effective_truth is not None and truth_snapshot is None:
                    truth_snapshot = load_truth(effective_truth)
                    manifest["truth_sha256"] = truth_snapshot.identity.sha256
                    manifest["truth_identity"] = asdict(truth_snapshot.identity)
                    atomic_json(manifest_path, manifest)
                elif truth_snapshot is not None:
                    verify_file_identity(truth_snapshot.path, truth_snapshot.identity)
                closure = closure_fields(pred_path, truth_snapshot, input_audit)
                if effective_truth is not None:
                    metrics = evaluate(
                        pred_path,
                        truth_snapshot,
                        set(input_audit["input_ids"]),
                        run_dir / "evaluation.tsv",
                    )
                else:
                    metrics = {
                        "exact": "NA",
                        "accuracy": "NA",
                        "mean_ed": "NA",
                        "evaluated_clusters": "NA",
                        "exact_count": "NA",
                        "total_ed": "NA",
                        "total_truth_bases": "NA",
                    }

                round_time = time.monotonic() - round_started
                row: dict[str, object] = {
                    "dataset": dataset,
                    "threads": workers,
                    "round": round_index,
                    "seed": seed,
                    "design_length": chosen_length,
                    **metrics,
                    "engine_time_s": engine_time,
                    "round_time_s": round_time,
                    "pred": str(relative_run_dir / "pred.tsv"),
                    **closure,
                }
                write_summary(run_dir / "summary.tsv", row)
                committed_run_dir = session_root / relative_run_dir
                committed_run_dir.parent.mkdir(parents=True, exist_ok=True)
                os.replace(run_dir, committed_run_dir)
                job_manifest.update({
                    "state": "COMPLETE",
                    "pred_artifact": str(relative_run_dir / "pred.tsv"),
                    "pred_sha256": sha256_path(committed_run_dir / "pred.tsv"),
                    "engine_time_s": engine_time,
                    "round_time_s": round_time,
                })
                atomic_json(manifest_path, manifest)
                rows.append(row)
                thread_rows.append(row)
                session_rows.append(row)
                completed = base + round_span
                progress_status = (
                    "complete"
                    if round_index == len(run_labels)
                    else f"round {round_index}/{len(run_labels)} complete"
                )
                dashboard.update(workers, completed, progress_status, force=True)

                divider = "─" * 78
                dashboard.log([
                    divider,
                    f"ROUND {round_index}/{len(run_labels)}  threads={workers}  seed={seed}",
                    (
                        f"RESULT exact={metrics['exact']}  accuracy={metrics['accuracy']}  "
                        f"mean_ed={metrics['mean_ed']}  "
                        f"time={engine_time:.3f}s"
                    ),
                    f"PRED   {final_session_root / relative_run_dir / 'pred.tsv'}",
                    divider,
                ])

            if len(thread_rows) > 1:
                aggregate = mean_row(thread_rows, workers)
                session_rows.append(aggregate)
                dashboard.log([
                    (
                        f"MEAN   threads={workers}  exact={aggregate['exact']}  "
                        f"accuracy={aggregate['accuracy']}  "
                        f"mean_ed={aggregate['mean_ed']}  "
                        f"time={aggregate['engine_time_s']}s"
                    )
                ])

        if len(session_rows) == 1:
            write_summary(session_root / "SUMMARY.tsv", session_rows[0])
        else:
            _write_rows(session_root / "SUMMARY.tsv", session_rows)
        status.update({
            "state": "COMPLETE",
            "completed_jobs": len(rows),
            "completed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
            "summary": "SUMMARY.tsv",
        })
        manifest.update({
            "state": "COMPLETE",
            "completed_jobs": len(rows),
            "session_wall_time_s": time.monotonic() - session_started,
            "summary_artifact": "SUMMARY.tsv",
        })
        atomic_json(status_path, status)
        atomic_json(manifest_path, manifest)
        if dashboard is not None:
            dashboard.close()
            dashboard = None
        os.replace(session_root, final_session_root)
        print(
            f"DONE  jobs={len(rows)}  summary={final_session_root / 'SUMMARY.tsv'}",
            flush=True,
        )
        return 0
    except KeyboardInterrupt:
        status.update({
            "state": "INTERRUPTED",
            "interrupted_utc": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
        })
        manifest.update({
            "state": "INTERRUPTED",
            "session_wall_time_s": time.monotonic() - session_started,
        })
        atomic_json(status_path, status)
        atomic_json(manifest_path, manifest)
        interrupted = result_root / f"{final_session_root.name}.interrupted"
        os.replace(session_root, interrupted)
        raise
    except BaseException as error:
        status.update({
            "state": "FAILED",
            "failed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
            "error": f"{type(error).__name__}: {error}",
        })
        manifest.update({
            "state": "FAILED",
            "session_wall_time_s": time.monotonic() - session_started,
            "error": f"{type(error).__name__}: {error}",
        })
        atomic_json(status_path, status)
        atomic_json(manifest_path, manifest)
        failed = result_root / f"{final_session_root.name}.failed"
        os.replace(session_root, failed)
        raise
    finally:
        if dashboard is not None:
            dashboard.close()


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description="Minimal ACOR streaming reconstruction")
    subparsers = parser.add_subparsers(dest="command", required=True)
    runner = subparsers.add_parser("run")
    runner.add_argument("--data", type=Path, required=True)
    mode = runner.add_mutually_exclusive_group()
    mode.add_argument("--aim", action="store_true")
    mode.add_argument("--no-aim", action="store_true")
    runner.add_argument("--length", type=int)
    runner.add_argument(
        "--k", type=int, choices=range(5, 16), default=DEFAULT_KMER_SIZE,
        help="k-mer size (5--15; default: 9)",
    )
    runner.add_argument("--threads", default="1")
    runner.add_argument("--seeds")
    runner.add_argument("--rounds", type=int)
    runner.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    try:
        return run(args)
    except (ValueError, RuntimeError, KeyError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
