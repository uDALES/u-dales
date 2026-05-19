import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udprep.udprep_init import _read_iexpnr_from_namoptions, validate_expnr  # noqa: E402


def _write_namoptions(path: Path, expnr: str, iexpnr_line: str | None = None) -> None:
    """Write a minimal namoptions file to *path*."""
    if iexpnr_line is None:
        iexpnr_line = f" iexpnr = {expnr},"
    path.write_text(
        f"&RUN\n{iexpnr_line}\n/\n",
        encoding="ascii",
    )


class TestReadIexpnrFromNameoptions(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)

    def test_reads_valid_expnr(self):
        f = self.workdir / "namoptions.042"
        _write_namoptions(f, "042")
        self.assertEqual(_read_iexpnr_from_namoptions(f), "042")

    def test_reads_expnr_without_trailing_comma(self):
        f = self.workdir / "namoptions.001"
        f.write_text("&RUN\n iexpnr = 001\n/\n", encoding="ascii")
        self.assertEqual(_read_iexpnr_from_namoptions(f), "001")

    def test_reads_expnr_with_inline_comment(self):
        f = self.workdir / "namoptions.100"
        f.write_text("&RUN\n iexpnr = 100, ! experiment number\n/\n", encoding="ascii")
        self.assertEqual(_read_iexpnr_from_namoptions(f), "100")

    def test_exits_on_missing_iexpnr(self):
        f = self.workdir / "namoptions.001"
        f.write_text("&RUN\n lbuoyancy = .true.\n/\n", encoding="ascii")
        with self.assertRaises(SystemExit):
            _read_iexpnr_from_namoptions(f)

    def test_exits_on_non_three_digit_iexpnr(self):
        f = self.workdir / "namoptions.001"
        f.write_text("&RUN\n iexpnr = 1,\n/\n", encoding="ascii")
        with self.assertRaises(SystemExit):
            _read_iexpnr_from_namoptions(f)

    def test_exits_on_four_digit_iexpnr(self):
        f = self.workdir / "namoptions.001"
        f.write_text("&RUN\n iexpnr = 1000,\n/\n", encoding="ascii")
        with self.assertRaises(SystemExit):
            _read_iexpnr_from_namoptions(f)

    def test_skips_comment_lines(self):
        f = self.workdir / "namoptions.007"
        f.write_text(
            "! iexpnr = 999\n&RUN\n iexpnr = 007,\n/\n", encoding="ascii"
        )
        self.assertEqual(_read_iexpnr_from_namoptions(f), "007")


class TestValidateExpnr(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        # Simulate an experiments directory with one case subdirectory
        self.expbase = Path(self.temp_dir.name) / "experiments"
        self.expbase.mkdir()

    def _make_case(self, name: str, iexpnr: str | None = None) -> Path:
        """Create a case directory with a matching namoptions file."""
        expdir = self.expbase / name
        expdir.mkdir()
        _write_namoptions(expdir / f"namoptions.{name}", iexpnr or name)
        return expdir

    def test_happy_path_returns_expnr(self):
        expdir = self._make_case("042")
        self.assertEqual(validate_expnr(expdir), "042")

    def test_exits_if_directory_does_not_exist(self):
        with self.assertRaises(SystemExit):
            validate_expnr(self.expbase / "999")

    def test_exits_if_no_namoptions_file(self):
        expdir = self.expbase / "001"
        expdir.mkdir()
        with self.assertRaises(SystemExit):
            validate_expnr(expdir)

    def test_exits_if_namoptions_suffix_not_three_digits(self):
        expdir = self.expbase / "abc"
        expdir.mkdir()
        (expdir / "namoptions.abc").write_text("&RUN\n iexpnr = abc,\n/\n", encoding="ascii")
        with self.assertRaises(SystemExit):
            validate_expnr(expdir)

    def test_exits_if_multiple_namoptions_files(self):
        expdir = self._make_case("010")
        (expdir / "namoptions.011").write_text("&RUN\n iexpnr = 011,\n/\n", encoding="ascii")
        with self.assertRaises(SystemExit):
            validate_expnr(expdir)

    def test_exits_if_dirname_mismatches_filename_suffix(self):
        # Directory is "050" but the namoptions file says 051
        expdir = self.expbase / "050"
        expdir.mkdir()
        _write_namoptions(expdir / "namoptions.051", "051")
        with self.assertRaises(SystemExit):
            validate_expnr(expdir)

    def test_exits_if_iexpnr_mismatches_filename_suffix(self):
        # namoptions.050 exists but iexpnr inside says 099
        expdir = self.expbase / "050"
        expdir.mkdir()
        _write_namoptions(expdir / "namoptions.050", "099")
        with self.assertRaises(SystemExit):
            validate_expnr(expdir)

    def test_exits_if_all_three_mismatch(self):
        expdir = self.expbase / "001"
        expdir.mkdir()
        _write_namoptions(expdir / "namoptions.002", "003")
        with self.assertRaises(SystemExit):
            validate_expnr(expdir)


if __name__ == "__main__":
    unittest.main()
