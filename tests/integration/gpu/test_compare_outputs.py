#!/usr/bin/env python3
"""Hardware-independent tests for the strict GPU parity comparator."""

from __future__ import annotations

import copy
import io
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from typing import Optional
from unittest.mock import patch

import netCDF4
import numpy as np

from compare_outputs import compare_output_directories
from run_gpu_tests import (
    CUDA_SELFTEST_FAIL,
    CUDA_SELFTEST_PASS,
    ConfigurationError,
    DEFAULT_MATRIX,
    REPO_ROOT,
    _read_matrix,
    _require_cuda_selftest,
    _validate_output_assertions,
    main,
    validate_matrix,
)


def _write_output(
    directory: Path,
    values: np.ndarray,
    *,
    include_scalar: bool = True,
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(directory / "fielddump.103.nc", "w") as dataset:
        dataset.createDimension("time", 1)
        dataset.createDimension("z", values.shape[0])
        dataset.createDimension("y", values.shape[1])
        dataset.createDimension("x", values.shape[2])
        time = dataset.createVariable("time", "f8", ("time",))
        time[:] = [1.0]
        field = dataset.createVariable("u", "f8", ("time", "z", "y", "x"))
        field[0] = values
        if include_scalar:
            scalar = dataset.createVariable("sca1", "f8", ("time", "z", "y", "x"))
            scalar[0] = values * 0.5


def _write_ranked_scalar_output(
    path: Path,
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
    peak_index: Optional[tuple[int, int, int]] = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    z_coordinates = np.array([0.5, 1.5, 2.5], dtype=np.float32)
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("time", 1)
        dataset.createDimension("zt", len(z_coordinates))
        dataset.createDimension("yt", len(y_coordinates))
        dataset.createDimension("xt", len(x_coordinates))
        dataset.createVariable("time", "f4", ("time",))[:] = [0.25]
        dataset.createVariable("zt", "f4", ("zt",))[:] = z_coordinates
        dataset.createVariable("yt", "f4", ("yt",))[:] = y_coordinates
        dataset.createVariable("xt", "f4", ("xt",))[:] = x_coordinates
        scalar = dataset.createVariable("sca1", "f4", ("time", "zt", "yt", "xt"))
        values = np.zeros((1, len(z_coordinates), len(y_coordinates), len(x_coordinates)))
        if peak_index is not None:
            values[(0, *peak_index)] = 2.0
        scalar[:] = values


class TestGpuOutputComparator(unittest.TestCase):
    def test_validate_config_output_explains_scope_and_next_step(self) -> None:
        output = io.StringIO()
        with patch.object(sys, "argv", ["run_gpu_tests.py", "full", "--validate-config"]):
            with redirect_stdout(output):
                status = main()

        message = output.getvalue()
        self.assertEqual(status, 0)
        self.assertIn("PASSED: GPU test matrix configuration is valid", message)
        self.assertIn("no solver executable was built or run and no GPU was used", message)
        self.assertIn("python tests/run_tests.py gpu-smoke", message)

    def test_committed_matrix_is_valid(self) -> None:
        validate_matrix(_read_matrix(DEFAULT_MATRIX))

    def test_matrix_rejects_unmapped_cuda_routine(self) -> None:
        matrix = copy.deepcopy(_read_matrix(DEFAULT_MATRIX))
        del matrix["ported_routines"]["cuda"]["src/modcuda.f90:initfield"]

        with self.assertRaisesRegex(ConfigurationError, "misses source routines"):
            validate_matrix(matrix)

    def test_matrix_rejects_incorrect_activation_assertion(self) -> None:
        matrix = copy.deepcopy(_read_matrix(DEFAULT_MATRIX))
        matrix["cases"][0]["namelist_assertions"]["ipoiss"] = 999

        with self.assertRaisesRegex(ConfigurationError, "namelist assertion ipoiss"):
            validate_matrix(matrix)

    def test_matrix_rejects_unsupported_gpu_option(self) -> None:
        matrix = copy.deepcopy(_read_matrix(DEFAULT_MATRIX))
        source = Path(matrix["cases"][0]["namelist"])
        source = REPO_ROOT / source
        with tempfile.TemporaryDirectory() as temporary:
            namelist = Path(temporary) / "namoptions.103"
            namelist.write_text(
                source.read_text(encoding="utf-8").replace("ipoiss       = 3", "ipoiss       = 2"),
                encoding="utf-8",
            )
            matrix["cases"][0]["namelist"] = str(namelist)

            with self.assertRaisesRegex(ConfigurationError, "unsupported GPU option ipoiss=2"):
                validate_matrix(matrix)

    def test_matrix_rejects_invalid_output_assertion(self) -> None:
        matrix = copy.deepcopy(_read_matrix(DEFAULT_MATRIX))
        matrix["cases"][0]["output_assertions"] = [
            {
                "type": "global_peak",
                "output": "missing-output",
                "variable": "u",
                "expected_coordinates": {"xt": 0.5},
            }
        ]

        with self.assertRaisesRegex(ConfigurationError, "is not in required_outputs"):
            validate_matrix(matrix)

    def test_global_peak_assertion_checks_coordinate_and_rank(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = Path(temporary)
            y_coordinates = np.array([0.5, 1.5, 2.5, 3.5, 4.5], dtype=np.float32)
            _write_ranked_scalar_output(
                run_dir / "fielddump.000.000.103.nc",
                np.array([0.5, 1.5, 2.5, 3.5, 4.5], dtype=np.float32),
                y_coordinates,
            )
            _write_ranked_scalar_output(
                run_dir / "fielddump.001.000.103.nc",
                np.array([5.5, 6.5, 7.5, 8.5], dtype=np.float32),
                y_coordinates,
                peak_index=(1, 4, 2),
            )
            case = {
                "case_id": 103,
                "output_assertions": [
                    {
                        "type": "global_peak",
                        "output": "fielddump",
                        "variable": "sca1",
                        "record": -1,
                        "minimum_value": 0.1,
                        "expected_coordinates": {"xt": 7.5, "yt": 4.5, "zt": 1.5},
                        "coordinate_atol": 1.0e-6,
                        "expected_rank": [1, 0],
                    }
                ],
            }

            self.assertEqual(_validate_output_assertions(case, run_dir, "CPU"), [])

    def test_global_peak_assertion_rejects_shifted_position(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = Path(temporary)
            _write_ranked_scalar_output(
                run_dir / "fielddump.001.000.103.nc",
                np.array([5.5, 6.5, 7.5, 8.5], dtype=np.float32),
                np.array([0.5, 1.5, 2.5, 3.5, 4.5], dtype=np.float32),
                peak_index=(1, 4, 3),
            )
            case = {
                "case_id": 103,
                "output_assertions": [
                    {
                        "type": "global_peak",
                        "output": "fielddump",
                        "variable": "sca1",
                        "minimum_value": 0.1,
                        "expected_coordinates": {"xt": 7.5},
                        "coordinate_atol": 1.0e-6,
                        "expected_rank": [1, 0],
                    }
                ],
            }

            failures = _validate_output_assertions(case, run_dir, "GPU")

            self.assertTrue(any("peak xt=8.5" in failure for failure in failures))

    def test_allows_values_inside_absolute_and_relative_tolerance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "cpu"
            candidate = root / "gpu"
            values = np.arange(8, dtype=np.float64).reshape(2, 2, 2)
            _write_output(reference, values)
            _write_output(candidate, values + 5.0e-7)

            result = compare_output_directories(
                reference,
                candidate,
                103,
                ["fielddump"],
                {"default": {"atol": 1.0e-6, "rtol": 1.0e-6}},
            )

            self.assertTrue(result.passed, result.failures)
            self.assertEqual(result.files_compared, 1)
            self.assertEqual(result.variables_compared, 3)

    def test_reports_variable_and_location_outside_tolerance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "cpu"
            candidate = root / "gpu"
            values = np.zeros((2, 2, 2), dtype=np.float64)
            changed = values.copy()
            changed[1, 0, 1] = 1.0
            _write_output(reference, values)
            _write_output(candidate, changed)

            result = compare_output_directories(
                reference,
                candidate,
                103,
                ["fielddump"],
                {"default": {"atol": 1.0e-6, "rtol": 1.0e-6}},
            )

            self.assertFalse(result.passed)
            self.assertTrue(any("fielddump.103.nc:u" in failure for failure in result.failures))
            variable = next(item for item in result.variable_results if item["variable"] == "u")
            self.assertEqual(tuple(variable["location"]), (0, 1, 0, 1))

    def test_rejects_missing_variable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "cpu"
            candidate = root / "gpu"
            values = np.zeros((2, 2, 2), dtype=np.float64)
            _write_output(reference, values, include_scalar=True)
            _write_output(candidate, values, include_scalar=False)

            result = compare_output_directories(
                reference,
                candidate,
                103,
                ["fielddump"],
            )

            self.assertFalse(result.passed)
            self.assertTrue(any("variables missing from GPU" in item for item in result.failures))

    def test_rejects_nonfinite_unmasked_value(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "cpu"
            candidate = root / "gpu"
            values = np.zeros((2, 2, 2), dtype=np.float64)
            changed = values.copy()
            changed[0, 0, 0] = np.inf
            _write_output(reference, values)
            _write_output(candidate, changed)

            result = compare_output_directories(
                reference,
                candidate,
                103,
                ["fielddump"],
            )

            self.assertFalse(result.passed)
            self.assertTrue(any("NaN or infinity" in item for item in result.failures))

    def test_cuda_selftest_requires_one_marker_per_rank(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "run.log"
            log.write_text(
                f"{CUDA_SELFTEST_PASS} rank=1\n{CUDA_SELFTEST_PASS} rank=0\n",
                encoding="utf-8",
            )

            self.assertIsNone(_require_cuda_selftest(log, expected_ranks=2))

    def test_cuda_selftest_rejects_missing_rank(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "run.log"
            log.write_text(f"{CUDA_SELFTEST_PASS} rank=0\n", encoding="utf-8")

            failure = _require_cuda_selftest(log, expected_ranks=2)

            self.assertIsNotNone(failure)
            self.assertIn("expected=[0, 1]", failure or "")

    def test_cuda_selftest_rejects_fortran_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "run.log"
            log.write_text(f"{CUDA_SELFTEST_FAIL} rank=0\n", encoding="utf-8")

            failure = _require_cuda_selftest(log, expected_ranks=1)

            self.assertEqual(
                failure,
                "CUDA device self-test suite reported failure",
            )


if __name__ == "__main__":
    unittest.main()
