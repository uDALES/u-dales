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

    def test_prepare_parser_needs_only_the_case_directory(self):
        args = self.prepare.build_parser().parse_args(["--case-dir", "/some/case"])
        self.assertEqual(args.statistics_source, "stats_tkslice")
        self.assertEqual(args.facet_source, [])
        self.assertEqual(args.stats_source, [])

    def test_parameter_reader_uses_defaults_and_rejects_unknown_fields(self):
        defaults = self.calculate.load_comfort_parameters(None)
        self.assertEqual(defaults.pet_met, 1.37)
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "parameters.json"
            path.write_text(json.dumps({"pet_clo": 0.7}), encoding="utf-8")
            self.assertEqual(
                self.calculate.load_comfort_parameters(path).pet_clo, 0.7,
            )
            path.write_text(json.dumps({"unknown": 1}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Unknown comfort parameters"):
                self.calculate.load_comfort_parameters(path)


if __name__ == "__main__":
    unittest.main()
