#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
DATASETS = [
    "oligo0", "ncomms19_365dishes", "binned_nanopore", "microsoft_cnr",
    "trellisbma", "dnaformer_nanopore", "dnaformer_illumina",
    "zenodo10943282", "err1816980_fountain", "err1816980_official",
]
THREADS = [1, 4, 16, 20, 32]
FIXED_SEED = "6c2f91a47bd803e52a9c14f8d730be65"
COUNT_FIELDS = ("exact_count", "total_ed", "total_truth_bases")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"empty TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def prepare_data_view(
    dataset: str,
    source: Path,
    output: Path,
    expected_reads_sha256: str,
    original_truth_sha256: str,
) -> tuple[Path, str, int]:
    view = output / "INPUT_VIEWS" / dataset
    summary_path = view / "TRUTH_SCOPE_AUDIT.tsv"
    if summary_path.is_file():
        summary = read_tsv(summary_path)[0]
        if summary["reads_sha256"] != expected_reads_sha256:
            raise RuntimeError(f"{dataset} cached reads SHA does not match frozen manifest")
        if summary["original_truth_sha256"] != original_truth_sha256:
            raise RuntimeError(f"{dataset} cached truth SHA does not match frozen manifest")
        if summary["view_required"] == "1":
            reads_view = view / "reads.tsv"
            effective_truth = view / "truth.tsv"
            if not reads_view.is_file() or not os.path.samefile(reads_view, source / "reads.tsv"):
                raise RuntimeError(f"{dataset} cached reads view is not the frozen source file")
            if not effective_truth.is_file() or sha256_file(effective_truth) != summary["effective_truth_sha256"]:
                raise RuntimeError(f"{dataset} cached effective truth SHA mismatch")
        return (
            source if summary["view_required"] == "0" else view,
            summary["effective_truth_sha256"],
            int(summary["truth_only_excluded"]),
        )
    input_ids: set[str] = set()
    with (source / "reads.tsv").open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "cluster_id" not in reader.fieldnames:
            raise RuntimeError(f"{dataset} reads source lacks cluster_id")
        for line_number, row in enumerate(reader, start=2):
            cluster_id = row["cluster_id"]
            if not cluster_id:
                raise RuntimeError(f"{dataset} empty reads cluster_id at row {line_number}")
            input_ids.add(cluster_id)
    if not input_ids:
        raise RuntimeError(f"{dataset} reads source has no clusters")
    with (source / "truth.tsv").open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"empty truth source: {dataset}")
        truth_fields = reader.fieldnames
        truth_rows = list(reader)
    def truth_id(row: dict[str, str]) -> str:
        return row.get("cluster_id") or row["row_key"].split("|")[1]
    truth_ids: set[str] = set()
    for line_number, row in enumerate(truth_rows, start=2):
        cluster_id = truth_id(row)
        if not cluster_id:
            raise RuntimeError(f"{dataset} empty truth cluster_id at row {line_number}")
        if cluster_id in truth_ids:
            raise RuntimeError(f"{dataset} duplicate truth cluster_id={cluster_id}")
        truth_ids.add(cluster_id)
    missing_truth = sorted(input_ids - truth_ids)
    truth_only = sorted(truth_ids - input_ids)
    if missing_truth:
        raise RuntimeError(f"{dataset} missing truth for input clusters: {len(missing_truth)}")
    effective_truth_sha256 = original_truth_sha256
    effective_root = source
    reads_view_mode = "ORIGINAL"
    if truth_only:
        view.mkdir(parents=True, exist_ok=True)
        reads_view = view / "reads.tsv"
        if reads_view.exists():
            if not os.path.samefile(reads_view, source / "reads.tsv"):
                raise RuntimeError(f"{dataset} existing reads view points to a different file")
            reads_view_mode = "EXISTING_SAMEFILE"
        else:
            try:
                os.link(source / "reads.tsv", reads_view)
                reads_view_mode = "HARDLINK"
            except OSError:
                reads_view.symlink_to((source / "reads.tsv").resolve())
                reads_view_mode = "SYMLINK"
        filtered_truth = view / "truth.tsv"
        temporary = view / ".truth.tsv.tmp"
        effective_rows = [row for row in truth_rows if truth_id(row) in input_ids]
        if len(effective_rows) != len(input_ids):
            raise RuntimeError(
                f"{dataset} effective truth closure failed: rows={len(effective_rows)} input_clusters={len(input_ids)}"
            )
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=truth_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(effective_rows)
        temporary.replace(filtered_truth)
        write_tsv(view / "TRUTH_ONLY_CLUSTER_IDS.tsv", [
            {"dataset": dataset, "cluster_id": cluster_id, "reason": "TRUTH_ONLY_NO_READS"}
            for cluster_id in truth_only
        ])
        effective_truth_sha256 = sha256_file(filtered_truth)
        effective_root = view
    actual_reads_sha256 = sha256_file(source / "reads.tsv")
    actual_truth_sha256 = sha256_file(source / "truth.tsv")
    if actual_reads_sha256 != expected_reads_sha256:
        raise RuntimeError(f"{dataset} reads SHA does not match frozen data manifest")
    if actual_truth_sha256 != original_truth_sha256:
        raise RuntimeError(f"{dataset} truth SHA does not match frozen data manifest")
    summary = [{
        "dataset": dataset,
        "input_clusters": len(input_ids),
        "original_truth_clusters": len(truth_ids),
        "effective_truth_clusters": len(input_ids),
        "truth_only_excluded": len(truth_only),
        "missing_truth": len(missing_truth),
        "view_required": int(bool(truth_only)),
        "policy": "EXCLUDE_TRUTH_ONLY_NO_READS",
        "original_truth_sha256": original_truth_sha256,
        "effective_truth_sha256": effective_truth_sha256,
        "reads_sha256": actual_reads_sha256,
        "reads_view_mode": reads_view_mode,
    }]
    write_tsv(summary_path, summary)
    return effective_root, effective_truth_sha256, len(truth_only)


def retry_task_root(base: Path) -> Path:
    if (base / "TASK_COMPLETE.json").is_file() or not base.exists() or not any(base.iterdir()):
        return base
    if not (base / "TASK_FAILED.json").is_file():
        raise RuntimeError(f"ambiguous partial task root: {base}")
    for attempt in range(1, 100):
        candidate = base.with_name(f"{base.name}.retry_truth_view_v{attempt}")
        if (candidate / "TASK_COMPLETE.json").is_file() or not candidate.exists() or not any(candidate.iterdir()):
            return candidate
        if not (candidate / "TASK_FAILED.json").is_file():
            raise RuntimeError(f"ambiguous partial retry task root: {candidate}")
    raise RuntimeError(f"too many failed retry roots: {base}")


def json_dump(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def canonical(path: Path) -> tuple[dict[str, str], int, str]:
    mapped: dict[str, str] = {}
    duplicates = 0
    for row in read_tsv(path):
        cluster_id = row["cluster_id"]
        duplicates += int(cluster_id in mapped)
        mapped[cluster_id] = row["reconstructed_sequence"]
    payload = "".join(f"{cluster_id}\t{mapped[cluster_id]}\n" for cluster_id in sorted(mapped))
    return mapped, duplicates, hashlib.sha256(payload.encode()).hexdigest()


def locate_session(result_root: Path) -> Path:
    sessions = [path for path in result_root.iterdir() if path.is_dir() and (path / "SUMMARY.tsv").is_file()]
    if len(sessions) != 1:
        raise RuntimeError(f"expected one session under {result_root}: {sessions}")
    return sessions[0]


def locate_order(session: Path) -> Path:
    paths = list(session.glob(".orders/*.jsonl.gz")) + list(session.glob("threads_*/round_*_seed_*/.internal/order.jsonl.gz"))
    if len(paths) != 1:
        raise RuntimeError(f"expected one order in {session}: {paths}")
    return paths[0]


def locate_pred(session: Path, threads: int) -> Path:
    paths = list(session.glob(f"threads_{threads}/round_01_seed_*/pred.tsv"))
    if len(paths) != 1:
        raise RuntimeError(f"expected one pred in {session}: {paths}")
    return paths[0]


def locate_runner_metrics(session: Path, threads: int) -> Path:
    paths = list(session.glob(f"threads_{threads}/round_01_seed_*/.internal/engine/RUNNER_METRICS.tsv"))
    if len(paths) != 1:
        raise RuntimeError(f"expected one runner metrics in {session}: {paths}")
    return paths[0]


def run_one(
    package: Path,
    package_label: str,
    dataset: str,
    data_dir: Path,
    threads: int,
    result_root: Path,
    cpu_pool: str,
    seed: str | None,
    reads_sha256: str,
    truth_original_sha256: str,
    truth_effective_sha256: str,
    truth_only_excluded: int,
) -> dict[str, object]:
    marker_path = result_root / "TASK_COMPLETE.json"
    if marker_path.is_file():
        record = json.loads(marker_path.read_text(encoding="utf-8"))
        record.setdefault("truth_original_sha256", record.get("truth_sha256", truth_original_sha256))
        record.setdefault("truth_effective_sha256", record.get("truth_sha256", truth_effective_sha256))
        record.setdefault("truth_only_excluded", truth_only_excluded)
        record.setdefault(
            "metric_counts_source",
            "NATIVE_EVALUATOR"
            if all(str(record.get(field, "NA")) != "NA" for field in COUNT_FIELDS)
            else "LEGACY_SUMMARY_UNAVAILABLE",
        )
        return record
    if result_root.exists() and any(result_root.iterdir()):
        raise RuntimeError(f"ambiguous partial task root: {result_root}")
    result_root.mkdir(parents=True, exist_ok=True)
    log_path = result_root / "console.log"
    command = [
        "taskset", "-c", cpu_pool, sys.executable, str(package / "acor.py"), "run",
        "--data", str(data_dir), "--aim", "--threads", str(threads),
    ]
    if seed is not None:
        command += ["--seeds", seed]
    command += ["--results", str(result_root)]
    started = time.monotonic()
    with log_path.open("wb") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    outer_wall = time.monotonic() - started
    if completed.returncode != 0:
        json_dump(result_root / "TASK_FAILED.json", {
            "dataset": dataset, "threads": threads, "package": package_label,
            "returncode": completed.returncode, "command": command,
        })
        raise RuntimeError(f"task failed {package_label}/{dataset}/P{threads}; see {log_path}")
    session = locate_session(result_root)
    summary_rows = [row for row in read_tsv(session / "SUMMARY.tsv") if row["round"] != "MEAN"]
    if len(summary_rows) != 1:
        raise RuntimeError(f"unexpected summary rows: {session}")
    summary = summary_rows[0]
    observed_seed = summary["seed"]
    pred = locate_pred(session, threads)
    order = locate_order(session)
    runner_metrics = read_tsv(locate_runner_metrics(session, threads))[0]
    mapped, duplicates, canonical_sha = canonical(pred)
    record: dict[str, object] = {
        "package": package_label,
        "dataset": dataset,
        "threads": threads,
        "seed": observed_seed,
        "seed_source": "EXPLICIT" if seed is not None else "OS_RANDOM",
        "reads_sha256": reads_sha256,
        "truth_sha256": truth_effective_sha256,
        "truth_original_sha256": truth_original_sha256,
        "truth_effective_sha256": truth_effective_sha256,
        "truth_only_excluded": truth_only_excluded,
        "runner_sha256": sha256_file(package / "bin/acor_runner"),
        "evaluator_sha256": sha256_file(package / "bin/ed_pairs"),
        "session": str(session.relative_to(result_root)),
        "order_path": str(order.relative_to(session)),
        "order_sha256": sha256_file(order),
        "pred_path": str(pred.relative_to(session)),
        "raw_pred_sha256": sha256_file(pred),
        "canonical_pred_sha256": canonical_sha,
        "pred_clusters": len(mapped),
        "duplicate": duplicates,
        "exact": summary["exact"],
        "accuracy": summary["accuracy"],
        "mean_ed": summary["mean_ed"],
        "exact_count": summary.get("exact_count", "NA"),
        "total_ed": summary.get("total_ed", "NA"),
        "total_truth_bases": summary.get("total_truth_bases", "NA"),
        "metric_counts_source": (
            "NATIVE_EVALUATOR"
            if all(str(summary.get(field, "NA")) != "NA" for field in COUNT_FIELDS)
            else "LEGACY_SUMMARY_UNAVAILABLE"
        ),
        "input_clusters": summary["input_unique_clusters"],
        "input_reads": summary["input_rows"],
        "engine_time_s": summary["engine_time_s"],
        "round_time_s": summary["round_time_s"],
        "outer_wall_s": f"{outer_wall:.9f}",
        "runner_workers": runner_metrics["workers"],
        "runner_worker_cpus": runner_metrics["worker_cpus"],
        "future_read_violations": runner_metrics["future_read_violations"],
        "safe_batch_violations": runner_metrics["safe_batch_violations"],
        "stop_index_violations": runner_metrics["stop_index_violations"],
        "exit_status": "PASS",
    }
    json_dump(marker_path, record)
    return record


def compare_pair(old: dict[str, object], new: dict[str, object], old_root: Path, new_root: Path) -> dict[str, object]:
    old_session = old_root / str(old["session"])
    new_session = new_root / str(new["session"])
    old_pred = old_session / str(old["pred_path"])
    new_pred = new_session / str(new["pred_path"])
    old_map, old_duplicates, old_canonical = canonical(old_pred)
    new_map, new_duplicates, new_canonical = canonical(new_pred)
    old_ids = set(old_map)
    new_ids = set(new_map)
    mismatches = sorted(cluster_id for cluster_id in old_ids & new_ids if old_map[cluster_id] != new_map[cluster_id])
    row = {
        "dataset": old["dataset"],
        "threads": old["threads"],
        "seed": old["seed"],
        "old_order_sha256": old["order_sha256"],
        "new_order_sha256": new["order_sha256"],
        "reads_sha256": old["reads_sha256"],
        "truth_sha256": old["truth_sha256"],
        "truth_original_sha256": old["truth_original_sha256"],
        "truth_effective_sha256": old["truth_effective_sha256"],
        "truth_only_excluded": old["truth_only_excluded"],
        "old_runner_sha256": old["runner_sha256"],
        "new_runner_sha256": new["runner_sha256"],
        "old_runner_workers": old["runner_workers"],
        "new_runner_workers": new["runner_workers"],
        "old_runner_worker_cpus": old["runner_worker_cpus"],
        "new_runner_worker_cpus": new["runner_worker_cpus"],
        "old_raw_pred_sha256": old["raw_pred_sha256"],
        "new_raw_pred_sha256": new["raw_pred_sha256"],
        "old_canonical_pred_sha256": old_canonical,
        "new_canonical_pred_sha256": new_canonical,
        "missing": len(old_ids - new_ids),
        "extra": len(new_ids - old_ids),
        "duplicate": old_duplicates + new_duplicates,
        "mismatch_clusters": len(mismatches),
        "mismatch_ids": ",".join(mismatches[:50]) or "NONE",
        "old_exact": old["exact"],
        "new_exact": new["exact"],
        "old_accuracy": old["accuracy"],
        "new_accuracy": new["accuracy"],
        "old_mean_ed": old["mean_ed"],
        "new_mean_ed": new["mean_ed"],
        "old_engine_time_s": old["engine_time_s"],
        "new_engine_time_s": new["engine_time_s"],
        "old_round_time_s": old["round_time_s"],
        "new_round_time_s": new["round_time_s"],
        "protocol_old": "LEGACY_HARDCODED_EXCLUDED",
        "protocol_new": "NOT_APPLICABLE",
        "schema_only_raw_difference": int(old["raw_pred_sha256"] != new["raw_pred_sha256"]),
        "status": "PASS",
    }
    passed = (
        row["old_order_sha256"] == row["new_order_sha256"]
        and old["truth_original_sha256"] == new["truth_original_sha256"]
        and old["truth_effective_sha256"] == new["truth_effective_sha256"]
        and old["truth_only_excluded"] == new["truth_only_excluded"]
        and row["missing"] == row["extra"] == row["duplicate"] == row["mismatch_clusters"] == 0
        and row["old_exact"] == row["new_exact"]
        and row["old_accuracy"] == row["new_accuracy"]
        and row["old_mean_ed"] == row["new_mean_ed"]
        and int(row["old_runner_workers"]) == int(row["threads"])
        and int(row["new_runner_workers"]) == int(row["threads"])
        and len(str(row["old_runner_worker_cpus"]).split(",")) == int(row["threads"])
        and len(str(row["new_runner_worker_cpus"]).split(",")) == int(row["threads"])
    )
    if not passed:
        row["status"] = "FAIL"
    return row


def backfill_identical_pair_counts(
    old: dict[str, object],
    new: dict[str, object],
    pair: dict[str, object],
) -> None:
    if pair["status"] != "PASS":
        return
    if any(str(new.get(field, "NA")) == "NA" for field in COUNT_FIELDS):
        raise RuntimeError(
            f"new evaluator omitted metric counts: {new['dataset']}/P{new['threads']}"
        )
    for field in COUNT_FIELDS:
        if str(old.get(field, "NA")) == "NA":
            old[field] = new[field]
        if str(old[field]) != str(new[field]):
            raise RuntimeError(
                f"paired metric count differs after canonical identity: "
                f"{old['dataset']}/P{old['threads']}/{field}"
            )
        pair[f"old_{field}"] = old[field]
        pair[f"new_{field}"] = new[field]
    if old["metric_counts_source"] == "LEGACY_SUMMARY_UNAVAILABLE":
        old["metric_counts_source"] = "PAIRED_NEW_EVALUATOR_AFTER_CANONICAL_IDENTITY"
    pair["old_metric_counts_source"] = old["metric_counts_source"]
    pair["new_metric_counts_source"] = new["metric_counts_source"]


def build_final_evidence(
    output: Path,
    run_rows: list[dict[str, object]],
    pair_rows: list[dict[str, object]],
    random_rows: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    fixed_new = {
        (str(row["dataset"]), int(row["threads"])): row
        for row in run_rows
        if row["package"] == "new"
    }
    pairs = {(str(row["dataset"]), int(row["threads"])): row for row in pair_rows}
    random_by_dataset = {str(row["dataset"]): row for row in random_rows}
    expected_fixed = {(dataset, threads) for dataset in DATASETS for threads in THREADS}
    if set(fixed_new) != expected_fixed or set(pairs) != expected_fixed:
        raise RuntimeError("fixed regression matrix is incomplete")
    if set(random_by_dataset) != set(DATASETS):
        raise RuntimeError("random replay matrix is incomplete")

    data_manifest: list[dict[str, object]] = []
    expected_predictions: list[dict[str, object]] = []
    thread_audit: list[dict[str, object]] = []
    for dataset in DATASETS:
        reference = fixed_new[(dataset, 1)]
        scope = read_tsv(output / "INPUT_VIEWS" / dataset / "TRUTH_SCOPE_AUDIT.tsv")[0]
        truth_only_ids = output / "INPUT_VIEWS" / dataset / "TRUTH_ONLY_CLUSTER_IDS.tsv"
        data_manifest.append({
            "dataset": dataset,
            "input_clusters": reference["input_clusters"],
            "input_reads": reference["input_reads"],
            "reads_sha256": reference["reads_sha256"],
            "original_truth_clusters": scope["original_truth_clusters"],
            "effective_truth_clusters": scope["effective_truth_clusters"],
            "truth_only_excluded": scope["truth_only_excluded"],
            "missing_truth": scope["missing_truth"],
            "truth_scope_policy": scope["policy"],
            "original_truth_sha256": reference["truth_original_sha256"],
            "effective_truth_sha256": reference["truth_effective_sha256"],
            "truth_only_ids_sha256": sha256_file(truth_only_ids) if truth_only_ids.is_file() else "NOT_APPLICABLE",
            "fixed_seed": FIXED_SEED,
            "thread_matrix": ",".join(map(str, THREADS)),
            "status": "PASS",
        })
        source_pairs = [pairs[(dataset, threads)] for threads in THREADS]
        for pair in source_pairs:
            expected_predictions.append({
                "dataset": dataset,
                "threads": pair["threads"],
                "fixed_seed": pair["seed"],
                "order_sha256": pair["new_order_sha256"],
                "canonical_prediction_sha256": pair["new_canonical_pred_sha256"],
                "old_raw_pred_sha256": pair["old_raw_pred_sha256"],
                "new_raw_pred_sha256": pair["new_raw_pred_sha256"],
                "exact": pair["new_exact"],
                "accuracy": pair["new_accuracy"],
                "mean_ed": pair["new_mean_ed"],
                "exact_count": pair["new_exact_count"],
                "total_ed": pair["new_total_ed"],
                "total_truth_bases": pair["new_total_truth_bases"],
                "mismatch_clusters": pair["mismatch_clusters"],
                "status": pair["status"],
            })
        canonical_hashes = {str(row["new_canonical_pred_sha256"]) for row in source_pairs}
        order_hashes = {str(row["new_order_sha256"]) for row in source_pairs}
        thread_propagation = all(
            int(row["new_runner_workers"]) == int(row["threads"])
            and len(str(row["new_runner_worker_cpus"]).split(",")) == int(row["threads"])
            for row in source_pairs
        )
        random_row = random_by_dataset[dataset]
        status = (
            "PASS"
            if len(canonical_hashes) == 1
            and len(order_hashes) == 1
            and thread_propagation
            and random_row["status"] == "PASS"
            and int(random_row["order_regenerated"]) == 1
            else "FAIL"
        )
        thread_audit.append({
            "dataset": dataset,
            "thread_matrix": ",".join(map(str, THREADS)),
            "unique_canonical_prediction_hashes": len(canonical_hashes),
            "unique_fixed_order_hashes": len(order_hashes),
            "runner_threads_propagated": int(thread_propagation),
            "random_seed": random_row["seed"],
            "random_order_regenerated": random_row["order_regenerated"],
            "random_replay_order_sha_match": int(
                random_row["random_order_sha256"] == random_row["replay_order_sha256"]
            ),
            "random_replay_prediction_sha_match": int(
                random_row["random_canonical_pred_sha256"]
                == random_row["replay_canonical_pred_sha256"]
            ),
            "status": status,
        })
    return data_manifest, expected_predictions, thread_audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--old-package", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu-pool", default="0-31,33-39")
    parser.add_argument("--datasets", default=",".join(DATASETS))
    parser.add_argument("--threads", default=",".join(map(str, THREADS)))
    parser.add_argument("--skip-random", action="store_true")
    args = parser.parse_args()
    data_root = args.data_root.resolve()
    old_package = args.old_package.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    selected_datasets = [value for value in args.datasets.split(",") if value]
    selected_threads = [int(value) for value in args.threads.split(",") if value]
    if any(value not in DATASETS for value in selected_datasets):
        raise RuntimeError(f"unknown dataset selection: {selected_datasets}")
    if any(value not in THREADS for value in selected_threads):
        raise RuntimeError(f"unknown thread selection: {selected_threads}")
    checksum_path = data_root.parent / "SHA256SUMS.txt"
    if not checksum_path.is_file():
        raise RuntimeError(f"missing data checksum manifest: {checksum_path}")
    data_hashes: dict[str, str] = {}
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        sha256, relative = line.split(maxsplit=1)
        data_hashes[relative] = sha256
    run_rows: list[dict[str, object]] = []
    pair_rows: list[dict[str, object]] = []
    random_rows: list[dict[str, object]] = []
    total = 2 * len(selected_datasets) * len(selected_threads)
    if not args.skip_random:
        total += 2 * len(selected_datasets)
    completed_tasks = 0
    started = time.monotonic()

    def progress(stage: str, dataset: str, threads: int) -> None:
        elapsed = time.monotonic() - started
        eta = elapsed / completed_tasks * (total - completed_tasks) if completed_tasks else 0.0
        status_row = {
            "stage": stage, "dataset": dataset, "threads": threads,
            "completed": completed_tasks, "total": total,
            "elapsed_s": elapsed, "eta_s": eta,
            "service_name": os.environ.get("ACOR_SERVICE_NAME", "FOREGROUND"),
            "pid": os.getpid(),
            "cpu_set": args.cpu_pool,
            "log_path": os.environ.get("ACOR_SERVICE_LOG", "STDOUT"),
        }
        json_dump(output / "CURRENT_STATUS.json", status_row)
        write_tsv(output / "RUN_STATUS.tsv", [status_row])

    for dataset in selected_datasets:
        source_data_dir = data_root / dataset
        if not (source_data_dir / "reads.tsv").is_file() or not (source_data_dir / "truth.tsv").is_file():
            raise RuntimeError(f"missing data source: {dataset}")
        reads_sha256 = data_hashes[f"real/{dataset}/reads.tsv"]
        truth_original_sha256 = data_hashes[f"real/{dataset}/truth.tsv"]
        data_dir, truth_effective_sha256, truth_only_excluded = prepare_data_view(
            dataset, source_data_dir, output, reads_sha256, truth_original_sha256)
        for threads in selected_threads:
            progress("FIXED_OLD", dataset, threads)
            old_root = retry_task_root(output / "RUNS/fixed/old" / dataset / f"P{threads}")
            old = run_one(old_package, "old", dataset, data_dir, threads, old_root, args.cpu_pool, FIXED_SEED, reads_sha256, truth_original_sha256, truth_effective_sha256, truth_only_excluded)
            completed_tasks += 1
            run_rows.append(old)
            progress("FIXED_NEW", dataset, threads)
            new_root = retry_task_root(output / "RUNS/fixed/new" / dataset / f"P{threads}")
            new = run_one(ROOT, "new", dataset, data_dir, threads, new_root, args.cpu_pool, FIXED_SEED, reads_sha256, truth_original_sha256, truth_effective_sha256, truth_only_excluded)
            completed_tasks += 1
            run_rows.append(new)
            pair = compare_pair(old, new, old_root, new_root)
            backfill_identical_pair_counts(old, new, pair)
            pair_rows.append(pair)
            write_tsv(output / "REAL10_RUN_LEDGER.tsv", run_rows)
            write_tsv(output / "REAL10_REGRESSION_LEDGER.tsv", pair_rows)
            if pair["status"] != "PASS":
                raise RuntimeError(f"prediction regression: {dataset}/P{threads}")

        if args.skip_random:
            continue
        progress("RANDOM_P20", dataset, 20)
        random_root = retry_task_root(output / "RUNS/random_replay/random" / dataset / "P20")
        random_run = run_one(ROOT, "new_random", dataset, data_dir, 20, random_root, args.cpu_pool, None, reads_sha256, truth_original_sha256, truth_effective_sha256, truth_only_excluded)
        completed_tasks += 1
        run_rows.append(random_run)
        progress("REPLAY_P20", dataset, 20)
        replay_root = retry_task_root(output / "RUNS/random_replay/replay" / dataset / "P20")
        replay = run_one(ROOT, "new_replay", dataset, data_dir, 20, replay_root, args.cpu_pool, str(random_run["seed"]), reads_sha256, truth_original_sha256, truth_effective_sha256, truth_only_excluded)
        completed_tasks += 1
        run_rows.append(replay)
        random_session = random_root / str(random_run["session"])
        replay_session = replay_root / str(replay["session"])
        random_order = random_session / str(random_run["order_path"])
        replay_order = replay_session / str(replay["order_path"])
        random_row = {
            "dataset": dataset,
            "threads": 20,
            "seed": random_run["seed"],
            "random_order_sha256": random_run["order_sha256"],
            "replay_order_sha256": replay["order_sha256"],
            "random_order_inode": random_order.stat().st_ino,
            "replay_order_inode": replay_order.stat().st_ino,
            "order_regenerated": int(random_order.stat().st_ino != replay_order.stat().st_ino),
            "random_canonical_pred_sha256": random_run["canonical_pred_sha256"],
            "replay_canonical_pred_sha256": replay["canonical_pred_sha256"],
            "random_exact": random_run["exact"],
            "replay_exact": replay["exact"],
            "status": "PASS",
        }
        if not (
            random_row["random_order_sha256"] == random_row["replay_order_sha256"]
            and random_row["order_regenerated"] == 1
            and random_row["random_canonical_pred_sha256"] == random_row["replay_canonical_pred_sha256"]
            and random_run["exact"] == replay["exact"]
            and random_run["accuracy"] == replay["accuracy"]
            and random_run["mean_ed"] == replay["mean_ed"]
        ):
            random_row["status"] = "FAIL"
        random_rows.append(random_row)
        write_tsv(output / "REAL10_RUN_LEDGER.tsv", run_rows)
        write_tsv(output / "RANDOM_REPLAY_LEDGER.tsv", random_rows)
        if random_row["status"] != "PASS":
            raise RuntimeError(f"random replay regression: {dataset}")

    data_manifest, expected_predictions, thread_audit = build_final_evidence(
        output, run_rows, pair_rows, random_rows)
    write_tsv(output / "REAL10_DATA_MANIFEST.tsv", data_manifest)
    write_tsv(output / "EXPECTED_PREDICTION_HASHES.tsv", expected_predictions)
    write_tsv(output / "THREAD_IDENTITY_AUDIT.tsv", thread_audit)
    complete_markers = list((output / "RUNS").rglob("TASK_COMPLETE.json"))
    failed_markers = list((output / "RUNS").rglob("TASK_FAILED.json"))
    if completed_tasks != total or len(run_rows) != total or len(complete_markers) != total:
        raise RuntimeError("completed-run evidence count is not closed")
    if any(row["status"] != "PASS" for row in pair_rows + random_rows + thread_audit):
        raise RuntimeError("final regression evidence contains a failed row")
    json_dump(output / "FINAL_STATUS.json", {
        "status": "PASS",
        "completed": completed_tasks,
        "total": total,
        "actual_completed_task_markers": len(complete_markers),
        "retained_failed_attempt_markers": len(failed_markers),
        "fixed_pairs": len(pair_rows),
        "random_pairs": len(random_rows),
        "truth_only_clusters_excluded": sum(int(row["truth_only_excluded"]) for row in data_manifest),
        "task_outer_wall_sum_s": sum(float(row["outer_wall_s"]) for row in run_rows),
        "ledger_rebuild_wall_s": time.monotonic() - started,
    })
    print(json.dumps({"status": "PASS", "completed": completed_tasks, "total": total}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
