from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from common import FIXTURES, SEED, read_rows, run_cli


class CliTests(unittest.TestCase):
    def test_aim_and_auto_truth(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            aim = run_cli(FIXTURES / "tiny_aim", Path(temp) / "aim", "--aim", "--seeds", SEED)
            automatic = run_cli(FIXTURES / "tiny_aim", Path(temp) / "auto", "--no-aim", "--seeds", SEED)
            self.assertEqual(read_rows(aim / "SUMMARY.tsv")[0]["exact"], "0.3333333333333333")
            self.assertEqual(read_rows(automatic / "SUMMARY.tsv")[0]["exact"], "0.3333333333333333")

    def test_no_truth_is_reconstruction_only(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            session = run_cli(FIXTURES / "tiny_noaim", Path(temp), "--no-aim", "--seeds", SEED)
            row = read_rows(session / "SUMMARY.tsv")[0]
            self.assertEqual((row["exact"], row["accuracy"], row["mean_ed"]), ("NA", "NA", "NA"))
            self.assertEqual(row["metric_closure"], "PASS")

    def test_seeds_and_rounds_write_mean(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            seeded = run_cli(FIXTURES / "tiny_aim", Path(temp) / "seeded", "--aim", "--seeds", "s1,s2")
            rows = read_rows(seeded / "SUMMARY.tsv")
            self.assertEqual([row["round"] for row in rows], ["1", "2", "MEAN"])
            random_rounds = run_cli(FIXTURES / "tiny_aim", Path(temp) / "rounds", "--aim", "--rounds", "2")
            rows = read_rows(random_rounds / "SUMMARY.tsv")
            self.assertEqual(len({row["seed"] for row in rows if row["round"] != "MEAN"}), 2)


if __name__ == "__main__":
    unittest.main()
