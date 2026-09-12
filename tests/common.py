from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


SOURCE_ROOT = Path(__file__).resolve().parent.parent
PACKAGE_ROOT = Path(os.environ.get("ACOR_PACKAGE_ROOT", SOURCE_ROOT)).resolve()
CLI = PACKAGE_ROOT / "acor.py"
FIXTURES = SOURCE_ROOT / "tests/fixtures"
GOLDEN = SOURCE_ROOT / "tests/golden"
SEED = "6c2f91a47bd803e52a9c14f8d730be65"


def run_cli(data: Path, results: Path, *options: str) -> Path:
    completed = subprocess.run(
        [sys.executable, str(CLI), "run", "--data", str(data), *options, "--results", str(results)],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(completed.stdout + completed.stderr)
    sessions = sorted(path for path in results.iterdir() if path.is_dir() and path.name.startswith("acor_"))
    if len(sessions) != 1:
        raise AssertionError(f"expected one complete session: {sessions}")
    return sessions[0]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pred_path(session: Path, threads: int = 1, round_index: int = 1) -> Path:
    return next(session.glob(f"threads_{threads}/round_{round_index:02d}_seed_*/pred.tsv"))


def order_path(session: Path, threads: int = 1, round_index: int = 1) -> Path:
    del threads
    return next(session.glob(f".orders/round_{round_index:02d}_seed_*.jsonl.gz"))


def canonical_sha(path: Path) -> str:
    rows = read_rows(path)
    payload = "".join(
        f"{row['cluster_id']}\t{row['reconstructed_sequence']}\n"
        for row in sorted(rows, key=lambda item: item["cluster_id"])
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def manifest(session: Path) -> dict:
    return json.loads((session / "RUN_MANIFEST.json").read_text(encoding="utf-8"))
