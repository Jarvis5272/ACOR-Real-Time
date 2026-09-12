from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from common import FIXTURES, SEED, canonical_sha, manifest, pred_path, run_cli


class ThreadIdentityTests(unittest.TestCase):
    def test_threads_are_passed_and_predictions_match(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            session = run_cli(FIXTURES / "tiny_aim", Path(temp), "--aim", "--threads", "1,4", "--seeds", SEED)
            self.assertEqual(canonical_sha(pred_path(session, 1)), canonical_sha(pred_path(session, 4)))
            jobs = manifest(session)["jobs"]
            self.assertEqual([job["threads"] for job in jobs], [1, 4])
            self.assertEqual(len({job["order_artifact"] for job in jobs}), 1)
            self.assertEqual(len(list((session / ".orders").glob("*.jsonl.gz"))), 1)
            for workers in (1, 4):
                command_path = next(session.glob(f"threads_{workers}/round_01_seed_*/.internal/engine/COMMAND.json"))
                command = json.loads(command_path.read_text())
                self.assertEqual(command["workers"], workers)
                self.assertEqual(len(command["worker_cpus"]), workers)


if __name__ == "__main__":
    unittest.main()
