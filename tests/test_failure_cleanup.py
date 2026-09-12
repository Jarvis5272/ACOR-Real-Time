from __future__ import annotations

import importlib.util
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path

from common import SOURCE_ROOT


def load_acor():
    spec = importlib.util.spec_from_file_location("acor_failure_test", SOURCE_ROOT / "acor.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class FailureCleanupTests(unittest.TestCase):
    def test_large_stderr_and_callback_failure(self) -> None:
        module = load_acor()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fake = root / "fake.py"
            fake.write_text(
                "#!/usr/bin/env python3\n"
                "import json, os, pathlib, subprocess, sys, time\n"
                "out=pathlib.Path(sys.argv[3]); mode=os.environ.get('FAKE_MODE')\n"
                "if mode=='large':\n"
                " os.write(2,b'X'*(11*1024*1024)); (out/'RUNNER_METRICS.tsv').write_text('total_wall_seconds\\n0.1\\n'); (out/'RUNNER_COMPLETE.json').write_text('{\"state\":\"COMPLETE\",\"status\":\"PASS\"}'); raise SystemExit(0)\n"
                "child=subprocess.Popen(['sleep','120']); (out/'pids').write_text(f'{os.getpid()}\\n{child.pid}\\n'); time.sleep(120)\n"
            )
            fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
            module.RUNNER = fake
            os.environ["FAKE_MODE"] = "large"
            self.assertEqual(module.run_engine(Path("r"), Path("o"), root / "large", "d", 1, lambda _x: None), 0.1)
            self.assertGreater((root / "large/runner.stderr.log").stat().st_size, 10 * 1024 * 1024)
            os.environ["FAKE_MODE"] = "sleep"
            calls = 0
            def fail(_frame: int) -> None:
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise RuntimeError("callback failure")
            with self.assertRaises(RuntimeError):
                module.run_engine(Path("r"), Path("o"), root / "fail", "d", 1, fail)
            os.environ.pop("FAKE_MODE", None)


if __name__ == "__main__":
    unittest.main()
