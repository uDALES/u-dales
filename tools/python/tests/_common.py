from __future__ import annotations

from pathlib import Path
import shutil
import sys
import uuid

REPO_ROOT = Path(__file__).resolve().parents[3]
PYTHON_DIR = REPO_ROOT / "tools" / "python"
TEST_TMP_DIR = REPO_ROOT / "tools" / "python" / "tests" / ".tmp"

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
    root = TEST_TMP_DIR / f"case_{uuid.uuid4().hex}"
    root.mkdir(parents=True, exist_ok=False)
    temp_dir = _CaseDir(root)
    case_dir = root / source.name
    shutil.copytree(source, case_dir)
    return temp_dir, case_dir
