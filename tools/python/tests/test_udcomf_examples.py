"""Parser, logging and configuration-reader tests for the thermal-comfort command-line scripts."""

from importlib.util import module_from_spec, spec_from_file_location
import json
import logging
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
import unittest

from udcomf import logs


_SCRIPTS = Path(__file__).parents[1] / "udcomf"


def _load_script(name: str):
    spec = spec_from_file_location(name, _SCRIPTS / f"{name}.py")
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class TestThermalComfortExamples(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.prepare = _load_script("udales_to_comfort")
        cls.calculate = _load_script("calculate_indices")

    def test_prepare_parser_requires_case_and_heights(self):
        with self.assertRaises(SystemExit):
            self.prepare.build_parser().parse_args(["--case-dir", "/some/case"])
        args = self.prepare.build_parser().parse_args(
            ["--case-dir", "/some/case", "--heights", "1.1", "1.5"])
        self.assertEqual(args.heights, [1.1, 1.5])
        self.assertEqual(args.statistics_source, "stats_t")
        self.assertEqual(args.data_source, [])
        args = self.prepare.build_parser().parse_args(
            ["--case-dir", "/c/4", "--heights", "1.1", "--data-source", "/c/1", "--data-source", "/c/2"])
        self.assertEqual([p.name for p in args.data_source], ["1", "2"])
        with self.assertRaises(SystemExit):
            self.prepare.build_parser().parse_args(["--case-dir", "/c", "--heights", "1", "--output-dir", "x"])

    def test_calculate_parser_defaults_to_all_indices(self):
        args = self.calculate.build_parser().parse_args(["a.nc"])
        self.assertEqual(args.indices, ("mrt", "pet", "utci", "wbgt"))
        args = self.calculate.build_parser().parse_args(["--indices", "wbgt,mrt", "a.nc", "b.nc"])
        self.assertEqual(args.indices, ("mrt", "wbgt"))
        self.assertEqual(len(args.exchange_files), 2)
        with self.assertRaises(SystemExit):
            self.calculate.build_parser().parse_args(["--indices", "pmv", "a.nc"])
        with self.assertRaises(SystemExit):
            self.calculate.build_parser().parse_args(["--indices", "", "a.nc"])

    def test_parameter_reader_uses_the_defaults_file_and_requires_complete_files(self):
        defaults = self.calculate.load_comfort_parameters(None)
        self.assertEqual(defaults.pet_met, 1.37)
        with open(self.calculate.DEFAULT_PARAMETERS_PATH, encoding="utf-8") as stream:
            raw = json.load(stream)
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "parameters.json"
            path.write_text(json.dumps(dict(raw, pet_clo=0.7)), encoding="utf-8")
            self.assertEqual(
                self.calculate.load_comfort_parameters(path).pet_clo, 0.7,
            )
            path.write_text(json.dumps({"pet_clo": 0.7}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Missing comfort parameters"):
                self.calculate.load_comfort_parameters(path)
            path.write_text(json.dumps(dict(raw, unknown=1)), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Unknown comfort parameters"):
                self.calculate.load_comfort_parameters(path)

class TestWrapperAndLaunchers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.wrapper = _load_script("ud_run_thermal_comfort")

    def test_wrapper_parser_is_the_union_of_both_steps(self):
        parser = self.wrapper.build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args(["--case-dir", "/c"])
        args = parser.parse_args([
            "--case-dir", "/c/004", "--heights", "1.1", "1.5", "--indices", "mrt,wbgt",
            "--data-source", "/c/001", "--statistics-source", "stats_tkslice",
            "--parameters", "p.json", "--tile-size", "64", "--overwrite",
        ])
        self.assertEqual(args.heights, [1.1, 1.5])
        self.assertEqual(args.indices, ("mrt", "wbgt"))
        self.assertEqual(args.statistics_source, "stats_tkslice")
        self.assertEqual([p.name for p in args.data_source], ["001"])
        self.assertEqual(args.tile_size, 64)
        self.assertTrue(args.overwrite)
        self.assertIsNone(args.log_file)
        with self.assertRaises(SystemExit):
            parser.parse_args(["--case-dir", "/c", "--heights", "1", "--output-dir", "x"])

    def test_wrapper_refuses_a_missing_case_directory(self):
        with TemporaryDirectory() as tmp:
            status = self.wrapper.main(["--case-dir", f"{tmp}/none", "--heights", "1.1"])
        self.assertEqual(status, 2)

    def test_wrapper_writes_the_chain_log_and_stops_at_preflight(self):
        with TemporaryDirectory() as tmp:
            case = Path(tmp) / "007"
            case.mkdir()
            (case / "namoptions.007").write_text(
                "&DOMAIN\n itot=1\n jtot=1\n ktot=2\n xlen=10.\n ylen=10.\n zsize=2.\n"
                " year=2023\n month=8\n day=21\n hour=0\n timezone=0.\n"
                " longitude=2.33\n latitude=48.86\n elevation=35.\n/\n"
                "&OUTPUT\n receptor_height=1.1\n tstatsdump=900.\n/\n", encoding="ascii")
            (case / "prof.inp.007").write_text("# z thl\n0.5 298.15\n1.5 298.15\n", encoding="ascii")
            status = self.wrapper.main(["--case-dir", str(case), "--heights", "1.1", "--indices", "mrt"])
            self.assertEqual(status, 1)
            log_text = (case / "thermal_comfort.007.log").read_text()
        self.assertIn("[preflight ]", log_text)
        self.assertIn("BLOCKER", log_text)
        self.assertIn("stopped before the index step", log_text)
        # headers sit at the margin, details are indented
        self.assertRegex(log_text, r"\[pipeline  \] thermal-comfort chain")
        self.assertRegex(log_text, r"\[preflight \]   receptor heights: 1.1 m")

    def test_shell_launcher_parses(self):
        script = Path(__file__).parents[2] / "ud_run_thermal_comfort.sh"
        self.assertTrue(script.is_file())
        subprocess.run(["bash", "-n", str(script)], check=True)
        result = subprocess.run(["bash", str(script)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("usage", result.stdout)


class TestStageLogging(unittest.TestCase):
    def test_stage_format_and_file_output(self):
        with TemporaryDirectory() as tmp:
            log_file = Path(tmp) / "chain.log"
            logs.configure(log_file, console=False)
            try:
                logs.header("radiation", "tracing %d tiles", 3)
                logs.stage_logger("radiation").info("tile %d/%d done", 1, 3)
                logs.stage_logger("pipeline").error("failed")
            finally:
                logs.configure(None, console=False)
            lines = log_file.read_text().splitlines()
        self.assertRegex(lines[0], r"^\d\d:\d\d:\d\d INFO    \[radiation \] tracing 3 tiles$")
        self.assertRegex(lines[1], r"^\d\d:\d\d:\d\d INFO    \[radiation \]   tile 1/3 done$")
        self.assertRegex(lines[2], r"^\d\d:\d\d:\d\d ERROR   \[pipeline  \]   failed$")
        with self.assertRaises(ValueError):
            logs.stage_logger("unknown")
        self.assertFalse(logging.getLogger("udcomf").propagate)


if __name__ == "__main__":
    unittest.main()
