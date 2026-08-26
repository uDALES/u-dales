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

from compare_outputs import _output_files, _read_fortran_records, compare_output_directories
from run_restart_roundtrip import first_restart, warm_namelist
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

    # ---- restart files -------------------------------------------------------
    #
    # The writer's layout is fixed in modsave::writerestartfiles; these build
    # files of that shape directly, with the 4-byte record markers gfortran and
    # nvfortran both emit, so the reader is tested without a solver run and
    # without a binary fixture in the repository.

    @staticmethod
    def _fortran_file(path: Path, records: list[np.ndarray], marker_bytes: int = 4) -> None:
        marker_dtype = "<i4" if marker_bytes == 4 else "<i8"
        with path.open("wb") as handle:
            for record in records:
                payload = np.ascontiguousarray(record).tobytes()
                marker = np.array([len(payload)], dtype=marker_dtype).tobytes()
                handle.write(marker + payload + marker)

    @classmethod
    def _restart_pair(cls, root: Path, case_id: int, *, perturb: Optional[str] = None,
                      nudge: float = 0.0) -> tuple[Path, Path]:
        rng = np.random.default_rng(7)
        cells, halo_cells = 4 * 4 * 4, 6 * 6 * 5
        base = {
            "mindist": rng.random(cells),
            "wall": rng.integers(0, 9, size=5 * cells).astype("<i4"),
        }
        for name in ("u0", "v0", "w0", "pres0", "thl0", "e120", "ekm", "qt0", "ql0", "ql0h"):
            base[name] = rng.random(halo_cells)
        base["time"] = np.array([2.5, 0.25])
        sv = {"sv0": rng.random(halo_cells * 2), "time": np.array([2.5])}

        dirs = []
        for side in ("cpu", "gpu"):
            side_dir = root / side
            side_dir.mkdir()
            main = {k: v.copy() for k, v in base.items()}
            scal = {k: v.copy() for k, v in sv.items()}
            if side == "gpu" and perturb is not None:
                target = scal if perturb == "sv0" else main
                if target[perturb].dtype.kind == "i":
                    target[perturb][0] += 1
                else:
                    target[perturb][0] += nudge
            cls._fortran_file(side_dir / f"initd00000010_000_000.{case_id:03d}", list(main.values()))
            cls._fortran_file(side_dir / f"inits00000010_000_000.{case_id:03d}", list(scal.values()))
            dirs.append(side_dir)
        return dirs[0], dirs[1]

    def test_restart_files_compare_record_by_record(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cpu, gpu = self._restart_pair(Path(tmp), 103)
            result = compare_output_directories(cpu, gpu, 103, ["initd", "inits"], {})
            self.assertTrue(result.passed, result.failures)
            self.assertEqual(result.files_compared, 2)
            names = sorted(r["variable"] for r in result.variable_results if r["file"].startswith("initd"))
            self.assertEqual(names, sorted(["mindist", "wall", "u0", "v0", "w0", "pres0", "thl0",
                                            "e120", "ekm", "qt0", "ql0", "ql0h", "time"]))

    def test_restart_real_record_outside_tolerance_is_named(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cpu, gpu = self._restart_pair(Path(tmp), 103, perturb="u0", nudge=1.0e-3)
            result = compare_output_directories(cpu, gpu, 103, ["initd", "inits"], {})
            self.assertFalse(result.passed)
            self.assertTrue(any("initd00000010_000_000.103:u0" in f for f in result.failures), result.failures)

    def test_restart_integer_record_must_match_exactly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cpu, gpu = self._restart_pair(Path(tmp), 103, perturb="wall")
            result = compare_output_directories(cpu, gpu, 103, ["initd"], {})
            self.assertFalse(result.passed)
            self.assertTrue(any(":wall" in f for f in result.failures), result.failures)

    def test_restart_scalar_file_is_compared_too(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cpu, gpu = self._restart_pair(Path(tmp), 103, perturb="sv0", nudge=1.0e-3)
            result = compare_output_directories(cpu, gpu, 103, ["initd", "inits"], {})
            self.assertFalse(result.passed)
            self.assertTrue(any("inits00000010_000_000.103:sv0" in f for f in result.failures), result.failures)

    def test_restart_time_record_uses_the_time_tolerance(self) -> None:
        # The matrix pins time at 1e-12. A 1e-9 shift in (timee, dt) is inside
        # the default tolerance and must still fail under the override.
        with tempfile.TemporaryDirectory() as tmp:
            cpu, gpu = self._restart_pair(Path(tmp), 103, perturb="time", nudge=1.0e-9)
            loose = compare_output_directories(cpu, gpu, 103, ["initd"], {})
            self.assertTrue(loose.passed, loose.failures)
            tight = compare_output_directories(
                cpu, gpu, 103, ["initd"],
                {"variables": {"time": {"atol": 1e-12, "rtol": 1e-12}}},
            )
            self.assertFalse(tight.passed)
            self.assertTrue(any(":time" in f for f in tight.failures), tight.failures)

    def test_restart_reader_rejects_foreign_record_markers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "initd00000010_000_000.103"
            self._fortran_file(path, [np.arange(4.0)], marker_bytes=8)
            with self.assertRaises(ValueError):
                _read_fortran_records(path)
            truncated = Path(tmp) / "initd00000011_000_000.103"
            self._fortran_file(truncated, [np.arange(4.0)])
            truncated.write_bytes(truncated.read_bytes()[:-3])
            with self.assertRaises(ValueError):
                _read_fortran_records(truncated)

    # ---- restart round trip ---------------------------------------------------

    def test_warm_namelist_substitutes_rather_than_appends(self) -> None:
        cold = (
            "&RUN\niexpnr = 103\nlwarmstart = .false.\nstartfile = ''\n"
            "runtime = 2.0\ntrestart = 1.0\n/\n"
        )
        warm = warm_namelist(cold, "initd00000004_000_000.103", 1.0)
        self.assertIn("lwarmstart = .true.", warm)
        self.assertIn("startfile = 'initd00000004_000_000.103'", warm)
        # runtime is relative to the restart time on a warm start
        self.assertIn("runtime = 1", warm)
        self.assertEqual(warm.count("lwarmstart"), 1)
        self.assertEqual(warm.count("runtime"), 1)
        with self.assertRaises(ValueError):
            warm_namelist("&RUN\nruntime = 2.0\n/\n", "x", 1.0)

    def test_first_restart_picks_the_earliest_step_and_its_companions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            for name in ("initd00000008_000_000.103", "inits00000008_000_000.103",
                         "initd00000004_000_000.103", "inits00000004_000_000.103",
                         "fielddump.000.000.103.nc"):
                (run_dir / name).write_bytes(b"")
            earliest, companions = first_restart(run_dir, 103)
            self.assertEqual(earliest.name, "initd00000004_000_000.103")
            self.assertEqual(sorted(p.name for p in companions),
                             ["initd00000004_000_000.103", "inits00000004_000_000.103"])
            with self.assertRaises(FileNotFoundError):
                first_restart(run_dir / "empty", 103)

    def test_output_prefix_prefers_exact_file_over_glob(self) -> None:
        # fac.064.nc, facEB.064.nc and facT.064.nc are three different outputs
        # that share a prefix. "fac" must name only its own file, while a
        # prefix with no exact file, like fielddump, still finds its
        # rank-suffixed one through the glob.
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            for name in ("fac.064.nc", "facEB.064.nc", "facT.064.nc",
                         "fielddump.000.000.064.nc"):
                (run_dir / name).write_bytes(b"")
            self.assertEqual(sorted(_output_files(run_dir, "fac", 64)), ["fac.064.nc"])
            self.assertEqual(sorted(_output_files(run_dir, "facEB", 64)), ["facEB.064.nc"])
            self.assertEqual(
                sorted(_output_files(run_dir, "fielddump", 64)),
                ["fielddump.000.000.064.nc"],
            )
            self.assertEqual(_output_files(run_dir, "initd", 64), {})

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
