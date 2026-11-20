"""Configuration-reader tests for the runnable thermal-comfort examples."""

from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest


_EXAMPLES = Path(__file__).parents[1] / "examples" / "thermal_comfort"


def _load_script(name: str):
    spec = spec_from_file_location(name, _EXAMPLES / f"{name}.py")
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

if __name__ == "__main__":
    unittest.main()
