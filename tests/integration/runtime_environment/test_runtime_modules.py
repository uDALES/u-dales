#!/usr/bin/env python3
"""Contract tests for opt-in integration-test module loading."""

from __future__ import annotations

import os
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path
from typing import Optional, Tuple


REPO_ROOT = Path(__file__).resolve().parents[3]
MODULE_HELPER = REPO_ROOT / "tests" / "integration" / "common" / "runtime_modules.sh"


class TestRuntimeModuleLoading(unittest.TestCase):
    def _run_loader(
        self,
        module_spec: Optional[str],
        module_status: int = 0,
    ) -> Tuple[subprocess.CompletedProcess, Path]:
        temporary = tempfile.TemporaryDirectory(prefix="udales-runtime-modules-")
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        fake_bin = root / "bin"
        fake_bin.mkdir()
        module_log = root / "module.log"
        fake_module = fake_bin / "module"
        fake_module.write_text(
            "#!/usr/bin/env bash\n"
            "printf '%s\\n' \"$@\" >> \"$UDALES_MODULE_LOG\"\n"
            "exit \"${UDALES_FAKE_MODULE_STATUS:-0}\"\n",
            encoding="utf-8",
        )
        fake_module.chmod(0o755)

        environment = os.environ.copy()
        for name in tuple(environment):
            if name.startswith("BASH_FUNC_module"):
                environment.pop(name)
        environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
        environment["UDALES_MODULE_LOG"] = str(module_log)
        environment["UDALES_FAKE_MODULE_STATUS"] = str(module_status)
        if module_spec is None:
            environment.pop("UDALES_RUNTIME_MODULES", None)
        else:
            environment["UDALES_RUNTIME_MODULES"] = module_spec

        command = (
            f"source {shlex.quote(str(MODULE_HELPER))} && "
            "load_udales_runtime_modules"
        )
        completed = subprocess.run(
            ["/bin/bash", "--noprofile", "--norc", "-c", command],
            env=environment,
            check=False,
            capture_output=True,
            text=True,
        )
        return completed, module_log

    def test_unset_module_list_does_not_invoke_module_command(self) -> None:
        completed, module_log = self._run_loader(None)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertFalse(module_log.exists())

    def test_explicit_module_list_is_forwarded_as_separate_arguments(self) -> None:
        completed, module_log = self._run_loader(
            "intel/2021a netCDF/4.8.0-iimpi-2021a FFTW/3.3.9-intel-2021a"
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(
            module_log.read_text(encoding="utf-8").splitlines(),
            [
                "load",
                "intel/2021a",
                "netCDF/4.8.0-iimpi-2021a",
                "FFTW/3.3.9-intel-2021a",
            ],
        )

    def test_explicit_module_load_failure_is_propagated(self) -> None:
        completed, module_log = self._run_loader("missing/module", module_status=17)

        self.assertEqual(completed.returncode, 17)
        self.assertEqual(
            module_log.read_text(encoding="utf-8").splitlines(),
            ["load", "missing/module"],
        )


if __name__ == "__main__":
    unittest.main()
