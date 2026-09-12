from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
import subprocess

from common import FIXTURES, GOLDEN, PACKAGE_ROOT


class CppBoundaryTests(unittest.TestCase):
    def invoke(self, output: Path, workers: str, cpus: str, control: str, config: str) -> subprocess.CompletedProcess:
        output.mkdir()
        return subprocess.run([
            str(PACKAGE_ROOT / "bin/acor_runner"),
            str(FIXTURES / "tiny_aim/reads.tsv"),
            str(GOLDEN / "tiny_fixed_order.jsonl.gz"),
            str(output), "test", "oligo0", "oligo0", workers, cpus, control,
            config, "NO_STOP", "0.0", "0.0", "NO_PROFILE",
        ], capture_output=True, text=True, check=False)

    def test_cpu_and_config_rejections(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cases = [
                ("negative", "1", "-1", "37", "kmc_k9_b16_lognormal"),
                ("cpu_tail", "1", "36x", "37", "kmc_k9_b16_lognormal"),
                ("overlap", "1", "36", "36", "kmc_k9_b16_lognormal"),
                ("config_k_tail", "1", "36", "37", "kmc_k9x_b16_lognormal"),
                ("config_b_tail", "1", "36", "37", "kmc_k9_b16x_lognormal"),
            ]
            for name, workers, cpus, control, config in cases:
                result = self.invoke(root / name, workers, cpus, control, config)
                self.assertNotEqual(result.returncode, 0, name)
                self.assertFalse((root / name / "RUNNER_COMPLETE.json").exists())


if __name__ == "__main__":
    unittest.main()
