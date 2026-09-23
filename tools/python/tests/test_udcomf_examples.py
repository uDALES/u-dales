"""Configuration-reader tests for the runnable thermal-comfort examples."""

from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np


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

    def test_loads_relative_georeference_archive_without_pickle(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            np.savez(
                root / "geo.npz",
                x=np.array([1.0]),
                y=np.array([2.0]),
                x_bounds=np.array([[0.5, 1.5]]),
                y_bounds=np.array([[1.5, 2.5]]),
                longitude=np.array([[2.3]]),
                latitude=np.array([[48.8]]),
                z_ground=np.array([[77.0]]),
            )
            config = {
                "georeference_npz": "geo.npz",
                "crs_attributes": {
                    "grid_mapping_name": "transverse_mercator",
                    "crs_wkt": "TEST WKT",
                },
                "vertical_datum": "IGN69",
                "model_version": "test",
                "institution": "test",
                "building_representation": "IBM facets",
                "terrain_convention": "flat ground",
                "spinup_start_utc": "2023-08-20T18:00:00Z",
                "spinup_end_utc": "2023-08-21T00:00:00Z",
            }
            path = root / "metadata.json"
            path.write_text(json.dumps(config), encoding="utf-8")

            metadata = self.prepare.load_exchange_metadata(path)

            np.testing.assert_array_equal(metadata.x_bounds, [[0.5, 1.5]])
            self.assertEqual(metadata.vertical_datum, "IGN69")

    def test_metadata_reader_rejects_unknown_fields(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "metadata.json"
            path.write_text(json.dumps({"unexpected": 1}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "fields are invalid"):
                self.prepare.load_exchange_metadata(path)

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
