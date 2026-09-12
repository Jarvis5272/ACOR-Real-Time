from __future__ import annotations

import csv
import shutil
import tempfile
import unittest
from pathlib import Path

from common import FIXTURES, SEED, canonical_sha, pred_path, read_rows, run_cli


def rewrite_reads(source: Path, destination: Path, field: str, value: str) -> None:
    with source.open("r", encoding="utf-8", newline="") as src, destination.open("w", encoding="utf-8", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        assert reader.fieldnames is not None
        writer = csv.DictWriter(dst, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            row[field] = value
            writer.writerow(row)


class IsolationTests(unittest.TestCase):
    def test_truth_original_and_quality_do_not_change_prediction(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cases = {}
            for name in ("normal", "wrong_truth", "no_truth", "original", "quality"):
                case = root / name
                case.mkdir()
                cases[name] = case
            shutil.copyfile(FIXTURES / "tiny_aim/reads.tsv", cases["normal"] / "reads.tsv")
            shutil.copyfile(FIXTURES / "tiny_aim/truth.tsv", cases["normal"] / "truth.tsv")
            shutil.copytree(cases["normal"], cases["wrong_truth"], dirs_exist_ok=True)
            with (cases["wrong_truth"] / "truth.tsv").open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fields = reader.fieldnames
                rows = list(reader)
            assert fields is not None
            for row in rows:
                row["reference_sequence"] = "T" * len(row["reference_sequence"])
            with (cases["wrong_truth"] / "truth.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            shutil.copyfile(FIXTURES / "tiny_aim/reads.tsv", cases["no_truth"] / "reads.tsv")
            rewrite_reads(FIXTURES / "tiny_aim/reads.tsv", cases["original"] / "reads.tsv", "original_read_sequence", "N")
            shutil.copyfile(FIXTURES / "tiny_aim/truth.tsv", cases["original"] / "truth.tsv")
            rewrite_reads(FIXTURES / "tiny_aim/reads.tsv", cases["quality"] / "reads.tsv", "read_quality", "IIII")
            shutil.copyfile(FIXTURES / "tiny_aim/truth.tsv", cases["quality"] / "truth.tsv")

            sessions = {
                name: run_cli(case, root / f"result_{name}", "--no-aim", "--seeds", SEED)
                for name, case in cases.items()
            }
            hashes = {name: canonical_sha(pred_path(session)) for name, session in sessions.items()}
            self.assertEqual(len(set(hashes.values())), 1)
            normal = read_rows(sessions["normal"] / "SUMMARY.tsv")[0]
            wrong = read_rows(sessions["wrong_truth"] / "SUMMARY.tsv")[0]
            absent = read_rows(sessions["no_truth"] / "SUMMARY.tsv")[0]
            self.assertNotEqual(normal["exact"], wrong["exact"])
            self.assertEqual((absent["exact"], absent["accuracy"], absent["mean_ed"]), ("NA", "NA", "NA"))


if __name__ == "__main__":
    unittest.main()
