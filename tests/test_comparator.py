from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path

from common import GOLDEN, SOURCE_ROOT


class ComparatorTests(unittest.TestCase):
    def compile(self, source: Path, output: Path) -> None:
        result = subprocess.run(
            ["g++", "-std=c++17", "-O2", "-DNDEBUG", str(source), "-o", str(output)],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_properties_and_frontier_golden(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            properties = root / "properties"
            probe = root / "frontier_probe"
            self.compile(SOURCE_ROOT / "tests/tools/comparator_properties.cpp", properties)
            self.compile(SOURCE_ROOT / "tests/tools/frontier_probe.cpp", probe)
            self.assertEqual(subprocess.run([str(properties)], check=False).returncode, 0)
            fixtures = SOURCE_ROOT / "tests/fixtures/tiny_tie"
            for name, input_name in (("tiny", "frontier_tiny_ordered.tsv"), ("tie", "frontier_reads.tsv")):
                frontier = root / f"frontier_{name}.tsv"
                final = root / f"frontier_{name}_final.tsv"
                result = subprocess.run([str(probe), str(fixtures / input_name), str(frontier), str(final)], check=False)
                self.assertEqual(result.returncode, 0)
                self.assertEqual(frontier.read_bytes(), (GOLDEN / f"frontier_{name}.tsv").read_bytes())
                self.assertEqual(final.read_bytes(), (GOLDEN / f"frontier_{name}_final.tsv").read_bytes())


if __name__ == "__main__":
    unittest.main()
