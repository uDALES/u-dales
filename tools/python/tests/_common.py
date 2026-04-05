from __future__ import annotations

from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory

REPO_ROOT = Path(__file__).resolve().parents[3]
PYTHON_DIR = REPO_ROOT / "tools" / "python"

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))


def copy_case(source: Path) -> tuple[TemporaryDirectory, Path]:
    temp_dir = TemporaryDirectory()
    case_dir = Path(temp_dir.name) / source.name
    shutil.copytree(source, case_dir)
    return temp_dir, case_dir
