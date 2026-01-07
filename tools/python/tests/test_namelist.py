import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
import shutil

PYTHON_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_DIR))

from udprep.namelist import update_namelist_value  # noqa: E402


class DummySim:
    def __init__(self, expnr: str, path: Path):
        self.expnr = expnr
        self.path = path


def _get_block(lines, name):
    name_lower = name.lower()
    start = None
    end = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.lower().startswith("&") and stripped[1:].strip().lower() == name_lower:
            start = i
            break
    if start is None:
        return None
    for i in range(start + 1, len(lines)):
        if lines[i].strip().startswith("/"):
            end = i
            break
    if end is None:
        return None
    return lines[start : end + 1]


class TestNamelistUpdate(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)
        source = Path(__file__).resolve().parents[3] / "tests" / "tests_tree_input" / "namoptions.525"
        shutil.copyfile(source, self.workdir / "namoptions.525")
        self.sim = DummySim("525", self.workdir)

    def _read_lines(self):
        return (self.workdir / "namoptions.525").read_text(encoding="ascii").splitlines()

    def test_update_existing_variable(self):
        update_namelist_value(self.sim, "TREES", "ntrees", 5)
        lines = self._read_lines()
        block = _get_block(lines, "TREES")
        self.assertIsNotNone(block)
        matches = [line for line in block if line.strip().lower().startswith("ntrees")]
        self.assertEqual(len(matches), 1)
        self.assertIn("= 5", matches[0])

    def test_add_missing_variable(self):
        update_namelist_value(self.sim, "TREES", "newvar", 123)
        lines = self._read_lines()
        block = _get_block(lines, "TREES")
        self.assertIsNotNone(block)
        matches = [line for line in block if line.strip().lower().startswith("newvar")]
        self.assertEqual(len(matches), 1)
        self.assertIn("= 123", matches[0])

    def test_add_missing_namelist(self):
        update_namelist_value(self.sim, "newblock", "foo", 7)
        lines = self._read_lines()
        block = _get_block(lines, "NEWBLOCK")
        self.assertIsNotNone(block)
        self.assertTrue(block[0].strip().startswith("&NEWBLOCK"))
        self.assertTrue(any(line.strip().lower().startswith("foo") for line in block))


if __name__ == "__main__":
    unittest.main()
