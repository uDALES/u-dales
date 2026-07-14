from __future__ import annotations

from pathlib import Path
import os
import shutil
import sys
import tempfile

REPO_ROOT = Path(__file__).resolve().parents[3]
PYTHON_DIR = REPO_ROOT / "tools" / "python"

# Base directory for test scratch space. Defaults to the OS temp dir (robust on
# read-only or OneDrive-synced checkouts, where an in-tree .tmp can't be
# created); override with UDALES_TEST_TMPDIR to place it elsewhere.
TEST_TMP_DIR = Path(os.environ.get("UDALES_TEST_TMPDIR") or tempfile.gettempdir())

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))


class _CaseDir:
    def __init__(self, path: Path):
        self.name = str(path)
        self._path = path

    def cleanup(self) -> None:
        shutil.rmtree(self._path, ignore_errors=True)


def copy_case(source: Path) -> tuple[_CaseDir, Path]:
    TEST_TMP_DIR.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="udales_case_", dir=TEST_TMP_DIR))
    temp_dir = _CaseDir(root)
    case_dir = root / source.name
    shutil.copytree(source, case_dir)
    return temp_dir, case_dir
