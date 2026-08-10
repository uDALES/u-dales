#!/usr/bin/env python3
"""Hardware-independent tests for the strict GPU parity comparator."""

from __future__ import annotations

import copy
import tempfile
import unittest
from pathlib import Path

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


class TestGpuOutputComparator(unittest.TestCase):
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
