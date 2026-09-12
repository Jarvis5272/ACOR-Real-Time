from __future__ import annotations

import subprocess
import unittest

from common import PACKAGE_ROOT


class EvaluatorTests(unittest.TestCase):
    def test_boundaries_and_byte_semantics(self) -> None:
        binary = PACKAGE_ROOT / "bin/ed_pairs"
        cases = [
            ("empty", "", "", 0),
            ("n", "AN", "AA", 1),
            ("case", "acgt", "ACGT", 4),
            ("indel", "A" * 65, "A" * 64, 1),
            ("different", "A" * 64, "T" * 64, 64),
        ]
        payload = "".join(f"{key}\t{left}\t{right}\n" for key, left, right, _ed in cases).encode()
        result = subprocess.run([str(binary)], input=payload, capture_output=True, check=False)
        self.assertEqual(result.returncode, 0)
        observed = {line.split("\t")[0]: int(line.split("\t")[3]) for line in result.stdout.decode().splitlines()}
        self.assertEqual(observed, {key: ed for key, _left, _right, ed in cases})

    def test_malformed_rows_fail(self) -> None:
        binary = PACKAGE_ROOT / "bin/ed_pairs"
        for payload in (b"id\tA\n", b"id\tA\tA\textra\n", b"\tA\tA\n", b"id\tA\tA\r\n"):
            result = subprocess.run([str(binary)], input=payload, capture_output=True, check=False)
            self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
