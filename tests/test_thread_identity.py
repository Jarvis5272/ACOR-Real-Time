from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path

from common import FIXTURES, SEED, canonical_sha, manifest, pred_path, run_cli


class ThreadIdentityTests(unittest.TestCase):
    def test_threads_are_passed_and_predictions_match(self) -> None:
        available_cpus = (
            len(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity")
            else (os.cpu_count() or 1)
        )
        parallel_workers = 2
        if available_cpus < parallel_workers + 1:
            self.skipTest("thread-identity test needs at least three available CPUs")
        with tempfile.TemporaryDirectory() as temp:
            session = run_cli(
                FIXTURES / "tiny_aim", Path(temp), "--aim", "--threads",
                f"1,{parallel_workers}", "--seeds", SEED,
            )
            self.assertEqual(
                canonical_sha(pred_path(session, 1)),
                canonical_sha(pred_path(session, parallel_workers)),
            )
            jobs = manifest(session)["jobs"]
            self.assertEqual([job["threads"] for job in jobs], [1, parallel_workers])
            self.assertEqual(len({job["order_artifact"] for job in jobs}), 1)
            self.assertEqual(len(list((session / ".orders").glob("*.jsonl.gz"))), 1)
            for workers in (1, parallel_workers):
                command_path = next(session.glob(f"threads_{workers}/round_01_seed_*/.internal/engine/COMMAND.json"))
                command = json.loads(command_path.read_text())
                self.assertEqual(command["workers"], workers)
                self.assertEqual(len(command["worker_cpus"]), workers)
                self.assertEqual(command["config"], "kmc_k9_b16_lognormal")
                self.assertEqual(command["algorithm_mode"], "NO_STOP")
                self.assertEqual(command["engine_mode"], "LEDGER_ADAPTIVE_POOL")

    def test_public_k_is_recorded_and_changes_engine_config(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            session = run_cli(
                FIXTURES / "tiny_aim", Path(temp), "--aim", "--k", "12",
                "--threads", "1", "--seeds", SEED,
            )
            run_manifest = manifest(session)
            self.assertEqual(run_manifest["k"], 12)
            self.assertEqual(run_manifest["config"], "kmc_k12_b16_lognormal")
            command_path = next(session.glob("threads_1/round_01_seed_*/.internal/engine/COMMAND.json"))
            command = json.loads(command_path.read_text())
            self.assertEqual(command["config"], "kmc_k12_b16_lognormal")


if __name__ == "__main__":
    unittest.main()
