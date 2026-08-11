#!/usr/bin/env python3
"""Build-wrapper tests that do not invoke a real compiler or CMake build."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from typing import Optional


REPO_ROOT = Path(__file__).resolve().parents[3]
BUILD_SCRIPT = REPO_ROOT / "tools" / "build_executable.sh"


class TestBuildExecutableLayout(unittest.TestCase):
    def _run_wrapper(
        self,
        system: str,
        build_type: str,
        expected_relative_path: str,
        *,
        override: Optional[str] = None,
    ) -> None:
        with tempfile.TemporaryDirectory(prefix="udales-build-wrapper-") as temporary:
            project = Path(temporary)
            (project / "src").mkdir()
            (project / "tools").mkdir()
            script = project / "tools" / "build_executable.sh"
            shutil.copy2(BUILD_SCRIPT, script)

            fake_bin = project / "fake-bin"
            fake_bin.mkdir()
            fake_cmake = fake_bin / "cmake"
            fake_cmake.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
            fake_cmake.chmod(0o755)
            fake_module = fake_bin / "module"
            fake_module.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
            fake_module.chmod(0o755)

            environment = os.environ.copy()
            for name in tuple(environment):
                if name.startswith("BASH_FUNC_module"):
                    environment.pop(name)
            environment["PATH"] = f"{fake_bin}{os.pathsep}{environment['PATH']}"
            environment.setdefault("NETCDF_DIR", "")
            if override is not None:
                environment["UDALES_BUILD_DIR"] = override

            completed = subprocess.run(
                ["bash", str(script), system, build_type],
                cwd=project,
                env=environment,
                check=False,
                capture_output=True,
                text=True,
            )

            expected = project / expected_relative_path
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertTrue(expected.is_dir())
            self.assertIn(f"Build directory: {expected}", completed.stdout)

    def test_all_non_gpu_systems_use_cpu_directory(self) -> None:
        for system in ("common", "icl", "archer", "cca"):
            with self.subTest(system=system):
                self._run_wrapper(system, "debug", "build/cpu/debug")

    def test_gpu_release_uses_gpu_directory(self) -> None:
        self._run_wrapper("gpu", "release", "build/gpu/release")

    def test_relative_override_takes_precedence(self) -> None:
        self._run_wrapper(
            "common",
            "release",
            "custom/cmake-cache",
            override="custom/cmake-cache",
        )


if __name__ == "__main__":
    unittest.main()
