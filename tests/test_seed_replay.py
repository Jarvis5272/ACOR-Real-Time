from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from common import FIXTURES, GOLDEN, SEED, manifest, order_path, pred_path, run_cli


class SeedReplayTests(unittest.TestCase):
    def test_fixed_seed_golden(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            session = run_cli(FIXTURES / "tiny_aim", Path(temp), "--aim", "--seeds", SEED)
            self.assertEqual(order_path(session).read_bytes(), (GOLDEN / "tiny_fixed_order.jsonl.gz").read_bytes())
            self.assertEqual(pred_path(session).read_bytes(), (GOLDEN / "tiny_fixed_pred.tsv").read_bytes())

    def test_random_seed_replay_regenerates_order(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            random_session = run_cli(FIXTURES / "tiny_aim", Path(temp) / "random", "--aim")
            seed = manifest(random_session)["seeds"][0]
            replay_session = run_cli(FIXTURES / "tiny_aim", Path(temp) / "replay", "--aim", "--seeds", seed)
            first_order = order_path(random_session)
            second_order = order_path(replay_session)
            self.assertNotEqual(first_order.stat().st_ino, second_order.stat().st_ino)
            self.assertEqual(first_order.read_bytes(), second_order.read_bytes())
            self.assertEqual(pred_path(random_session).read_bytes(), pred_path(replay_session).read_bytes())


if __name__ == "__main__":
    unittest.main()
