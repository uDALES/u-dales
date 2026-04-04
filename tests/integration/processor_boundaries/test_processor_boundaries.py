#!/usr/bin/env python3

import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import netCDF4 as nc
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
TEST_DIR = Path(__file__).resolve().parent
TREE_CASE_ID = 526
NO_TREE_CASE_ID = 100
UDALES_BUILD = Path(os.environ.get("UDALES_BUILD", REPO_ROOT / "build" / "release" / "u-dales"))
RUNTIME_MODULES = os.environ.get(
    "UDALES_RUNTIME_MODULES",
    "intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a "
    "FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0",
)
TREE_FIELDS = ("tr_u", "tr_v", "tr_w")
NO_TREE_FIELDS = ("ut", "vt", "wt")
ABS_TOL = 1.0e-9
# The no-tree Xie/Castro parity check compares tdump output, which is written
# to NetCDF as NF90_FLOAT in modstat_nc.f90. The residual y-split/xy-split wt
# mismatch observed on Ubuntu is therefore quantized at float32 precision even
# though the solver runs in r8. Keep this tolerance tight and scoped only to
# that global no-tree comparison.
NO_TREE_GLOBAL_ABS_TOL = 2.0e-8
CONFIGS = {
    "serial": (1, 1),
    "x_split": (2, 1),
    "y_split": (1, 2),
    "xy_split": (2, 2),
}


def _read_int_setting(namelist: Path, key: str) -> int:
    match = re.search(
        rf"(?m)^\s*{re.escape(key)}\s*=\s*(\d+)\s*$",
        namelist.read_text(encoding="utf-8"),
    )
    if not match:
        raise RuntimeError(f"Missing integer setting '{key}' in {namelist}")
    return int(match.group(1))


def _patch_namelist(run_dir: Path, case_id: int, nprocx: int, nprocy: int, vegetation_enabled: bool) -> None:
    namelist = run_dir / f"namoptions.{case_id}"
    text = namelist.read_text(encoding="utf-8")
    text = re.sub(r"(?m)^(\s*nprocx\s*=\s*)\d+(\s*)$", rf"\g<1>{nprocx}\2", text)
    text = re.sub(r"(?m)^(\s*nprocy\s*=\s*)\d+(\s*)$", rf"\g<1>{nprocy}\2", text)
    text = re.sub(r"(?m)^(\s*runtime\s*=\s*)[0-9.]+(\s*)$", r"\g<1>0.01\2", text)
    text = re.sub(r"(?m)^(\s*dtmax\s*=\s*)[0-9.]+(\s*)$", r"\g<1>0.01\2", text)
    text = re.sub(r"(?m)^(\s*tstatstart\s*=\s*)[0-9.]+(\s*)$", r"\g<1>0.0\2", text)
    text = re.sub(r"(?m)^(\s*tstatsdump\s*=\s*)[0-9.]+(\s*)$", r"\g<1>0.005\2", text)
    text = re.sub(r"(?m)^(\s*tsample\s*=\s*)[0-9.]+(\s*)$", r"\g<1>0.005\2", text)
    text = re.sub(r"(?m)^(\s*lxytdump\s*=\s*)\.true\.(\s*)$", r"\g<1>.false.\2", text)
    if re.search(r"(?m)^\s*randu\s*=", text):
        text = re.sub(r"(?m)^(\s*randu\s*=\s*)[0-9.]+(\s*)$", r"\g<1>1.0\2", text)
    else:
        text = re.sub(r"(?m)^(&RUN\s*)$", r"\1\nrandu        = 1.0", text, count=1)
    if not vegetation_enabled:
        text = re.sub(r"(?m)^(\s*ltreesfile\s*=\s*)\.true\.(\s*)$", r"\g<1>.false.\2", text)
        text = re.sub(r"(?m)^(\s*ltrees\s*=\s*)\.true\.(\s*)$", r"\g<1>.false.\2", text)
        text = re.sub(r"(?m)^(\s*ltreedump\s*=\s*)\.true\.(\s*)$", r"\g<1>.false.\2", text)
    namelist.write_text(text, encoding="utf-8")


def _load_veg_mask(case_dir: Path, case_id: int) -> np.ndarray:
    namelist = case_dir / f"namoptions.{case_id}"
    itot = _read_int_setting(namelist, "itot")
    jtot = _read_int_setting(namelist, "jtot")
    ktot = _read_int_setting(namelist, "ktot")
    mask = np.zeros((ktot, jtot, itot), dtype=bool)
    with (case_dir / f"veg.inp.{case_id}").open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            i, j, k = (int(part) for part in stripped.split()[:3])
            mask[k - 1, j - 1, i - 1] = True
    return mask


def _support_mask(field: str, veg_mask: np.ndarray) -> np.ndarray:
    support = np.zeros_like(veg_mask, dtype=bool)
    if field == "tr_u":
        support[:, :, 1:] = veg_mask[:, :, :-1] | veg_mask[:, :, 1:]
    elif field == "tr_v":
        support[:, 1:, :] = veg_mask[:, :-1, :] | veg_mask[:, 1:, :]
    elif field == "tr_w":
        support[1:, :, :] = veg_mask[:-1, :, :] | veg_mask[1:, :, :]
    else:
        raise ValueError(f"Unknown field: {field}")
    return support


def _processor_boundary_band(field: str, shape: Tuple[int, int, int], nprocx: int, nprocy: int) -> np.ndarray:
    band = np.zeros(shape, dtype=bool)
    if field == "tr_u" and nprocx > 1:
        tile_width = shape[2] // nprocx
        for boundary in range(tile_width, shape[2], tile_width):
            lo = max(boundary - 1, 0)
            hi = min(boundary + 1, shape[2] - 1)
            band[:, :, lo : hi + 1] = True
    elif field == "tr_v" and nprocy > 1:
        tile_height = shape[1] // nprocy
        for boundary in range(tile_height, shape[1], tile_height):
            lo = max(boundary - 1, 0)
            hi = min(boundary + 1, shape[1] - 1)
            band[:, lo : hi + 1, :] = True
    return band


def _parse_tile_indices(path: Path, prefix: str, case_id: int) -> Tuple[int, int]:
    match = re.match(rf"{re.escape(prefix)}\.(\d{{3}})\.(\d{{3}})\.{case_id}\.nc$", path.name)
    if not match:
        raise RuntimeError(f"Unexpected treedump file name: {path.name}")
    return int(match.group(1)), int(match.group(2))


def _load_global_fields(run_dir: Path, case_id: int, prefix: str, fields: Iterable[str]) -> Dict[str, np.ndarray]:
    namelist = run_dir / f"namoptions.{case_id}"
    itot = _read_int_setting(namelist, "itot")
    jtot = _read_int_setting(namelist, "jtot")
    ktot = _read_int_setting(namelist, "ktot")
    nprocx = _read_int_setting(namelist, "nprocx")
    nprocy = _read_int_setting(namelist, "nprocy")

    field_names = tuple(fields)
    files = sorted(run_dir.glob(f"{prefix}.*.*.{case_id}.nc"))
    if not files:
        raise RuntimeError(f"No {prefix} outputs found in {run_dir}")

    fields = {
        field: np.zeros((ktot, jtot, itot), dtype=np.float64) for field in field_names
    }
    for path in files:
        px, py = _parse_tile_indices(path, prefix, case_id)
        with nc.Dataset(path) as ds:
            sample = np.asarray(ds.variables[field_names[0]][0], dtype=np.float64)
            nz, ny, nx = sample.shape
            istart = px * (itot // nprocx)
            jstart = py * (jtot // nprocy)
            for field in field_names:
                arr = np.asarray(ds.variables[field][0], dtype=np.float64)
                fields[field][:nz, jstart : jstart + ny, istart : istart + nx] = arr
    return fields


def _mpi_exec_and_args() -> Tuple[str, str]:
    mpiexec = os.environ.get("MPIEXEC")
    if not mpiexec:
        mpiifort = shutil.which("mpiifort")
        if mpiifort:
            mpiexec = str(Path(mpiifort).parent / "mpiexec")
        else:
            mpiexec = "mpiexec"

    extra_args = os.environ.get("MPI_LAUNCH_EXTRA_ARGS", "").strip()
    try:
        version = subprocess.run(
            [mpiexec, "--version"],
            check=False,
            capture_output=True,
            text=True,
        ).stdout
    except OSError:
        version = ""
    if re.search(r"Open MPI|OpenRTE", version, flags=re.IGNORECASE) and "--oversubscribe" not in extra_args:
        extra_args = f"--oversubscribe {extra_args}".strip()
    return mpiexec, extra_args


def _run_case(executable: Path, run_dir: Path, case_id: int, nprocs: int) -> None:
    mpiexec, extra_args = _mpi_exec_and_args()
    diag_log = run_dir / "fortran_diagnostics.log"
    for stale in run_dir.glob("monitor*.txt"):
        stale.unlink()
    for stale in run_dir.glob(f"treedump.*.*.{case_id}.nc"):
        stale.unlink()
    for stale in run_dir.glob(f"tdump.*.*.{case_id}.nc"):
        stale.unlink()
    if diag_log.exists():
        diag_log.unlink()

    module_setup = ""
    if Path("/etc/profile.d/modules.sh").is_file():
        module_setup = "source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; "
    shell_command = (
        f"{module_setup}"
        f"if command -v module >/dev/null 2>&1; then module load {RUNTIME_MODULES}; fi && "
        f"export HDF5_USE_FILE_LOCKING=FALSE && "
        f"export FOR_DISABLE_DIAGNOSTIC_DISPLAY=TRUE && "
        f"export FOR_DIAGNOSTIC_LOG_FILE='{diag_log}' && "
        f"cd '{run_dir}' && "
        f"'{mpiexec}' {extra_args} -n {nprocs} '{executable}' namoptions.{case_id}"
    )
    completed = subprocess.run(
        ["bash", "-lc", shell_command],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode == 0:
        return

    def _tail(label: str, text: str, limit: int = 40) -> str:
        stripped = text.strip()
        if not stripped:
            return f"{label}: <empty>"
        lines = stripped.splitlines()
        snippet = "\n".join(lines[-limit:])
        return f"{label} (last {min(len(lines), limit)} lines):\n{snippet}"

    details = [
        f"u-dales run failed for case {case_id} in {run_dir}",
        f"command: {shell_command}",
        f"exit code: {completed.returncode}",
        _tail("stdout", completed.stdout),
        _tail("stderr", completed.stderr),
    ]
    if diag_log.exists():
        details.append(
            _tail(
                "fortran diagnostics",
                diag_log.read_text(encoding="utf-8", errors="replace"),
            )
        )
    raise RuntimeError("\n\n".join(details))


def _launcher_unavailable_reason() -> Optional[str]:
    mpiexec, extra_args = _mpi_exec_and_args()
    module_setup = ""
    if Path("/etc/profile.d/modules.sh").is_file():
        module_setup = "source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; "
    shell_command = (
        f"{module_setup}"
        f"if command -v module >/dev/null 2>&1; then module load {RUNTIME_MODULES}; fi && "
        f"'{mpiexec}' {extra_args} -n 1 /bin/true"
    )
    probe = subprocess.run(
        ["bash", "-lc", shell_command],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if probe.returncode == 0:
        return None
    details = (probe.stderr or "").strip() or (probe.stdout or "").strip()
    if not details:
        details = f"mpiexec probe failed with exit code {probe.returncode}"
    return details


def _max_diff_with_location(diff: np.ndarray, mask: np.ndarray) -> Tuple[float, Optional[Tuple[int, int, int]]]:
    if not np.any(mask):
        return 0.0, None
    masked = np.where(mask, diff, -1.0)
    flat = int(np.argmax(masked))
    idx = tuple(int(v) for v in np.unravel_index(flat, diff.shape))
    return float(diff[idx]), idx


class _BaseProcessorBoundaryParity(unittest.TestCase):
    CASE_ID = 0
    PREFIX = ""
    FIELDS = ()
    VEGETATION_ENABLED = True

    @classmethod
    def setUpClass(cls) -> None:
        if not UDALES_BUILD.is_file():
            raise RuntimeError(f"u-dales executable not found at {UDALES_BUILD}")
        cls.case_source = REPO_ROOT / "tests" / "cases" / str(cls.CASE_ID)
        if not cls.case_source.is_dir():
            raise RuntimeError(f"Case source not found: {cls.case_source}")
        launcher_issue = _launcher_unavailable_reason()
        if launcher_issue is not None:
            raise unittest.SkipTest(
                f"MPI launcher is not usable in this environment: {launcher_issue}"
            )

        cls.veg_mask = None
        if cls.VEGETATION_ENABLED:
            cls.veg_mask = _load_veg_mask(cls.case_source, cls.CASE_ID)
        cls.workdir = tempfile.TemporaryDirectory(prefix="udales-veg-proc-boundary-", dir="/tmp")
        cls.run_dirs: Dict[str, Path] = {}
        cls.outputs: Dict[str, Dict[str, np.ndarray]] = {}

        for label, (nprocx, nprocy) in CONFIGS.items():
            run_dir = Path(cls.workdir.name) / label
            shutil.copytree(cls.case_source, run_dir)
            _patch_namelist(run_dir, cls.CASE_ID, nprocx, nprocy, cls.VEGETATION_ENABLED)
            _run_case(UDALES_BUILD, run_dir, cls.CASE_ID, nprocx * nprocy)
            cls.run_dirs[label] = run_dir
            cls.outputs[label] = _load_global_fields(run_dir, cls.CASE_ID, cls.PREFIX, cls.FIELDS)

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(cls, "workdir"):
            cls.workdir.cleanup()

    def _assert_fields_match(
        self,
        reference_label: str,
        candidate_label: str,
        fields: Iterable[str],
        mask_builder,
        atol: float = ABS_TOL,
    ) -> None:
        reference = self.outputs[reference_label]
        candidate = self.outputs[candidate_label]
        nprocx, nprocy = CONFIGS[candidate_label]
        failures = []

        for field in fields:
            mask = mask_builder(field, reference[field].shape, nprocx, nprocy)
            diff = np.abs(candidate[field] - reference[field])
            max_diff, idx = _max_diff_with_location(diff, mask)
            if max_diff > atol:
                failures.append(
                    f"{field} max abs diff {max_diff:.3e} at {idx} for {candidate_label}"
                )

        if failures:
            self.fail(
                f"Processor-boundary parity failed against {reference_label}:\n- "
                + "\n- ".join(failures)
            )


class TestNoTreeProcessorBoundaryParity(_BaseProcessorBoundaryParity):
    CASE_ID = NO_TREE_CASE_ID
    PREFIX = "tdump"
    FIELDS = NO_TREE_FIELDS
    VEGETATION_ENABLED = False

    def test_x_split_matches_serial_on_u_boundary_band(self) -> None:
        def mask_builder(field: str, shape: Tuple[int, int, int], nprocx: int, nprocy: int) -> np.ndarray:
            if field == "ut":
                return _processor_boundary_band("tr_u", shape, nprocx, nprocy)
            return np.zeros(shape, dtype=bool)

        self._assert_fields_match("serial", "x_split", ("ut",), mask_builder)

    def test_y_split_matches_serial_on_v_boundary_band(self) -> None:
        def mask_builder(field: str, shape: Tuple[int, int, int], nprocx: int, nprocy: int) -> np.ndarray:
            if field == "vt":
                return _processor_boundary_band("tr_v", shape, nprocx, nprocy)
            return np.zeros(shape, dtype=bool)

        self._assert_fields_match("serial", "y_split", ("vt",), mask_builder)

    def test_xy_split_matches_serial_on_boundary_bands(self) -> None:
        def mask_builder(field: str, shape: Tuple[int, int, int], nprocx: int, nprocy: int) -> np.ndarray:
            if field == "ut":
                return _processor_boundary_band("tr_u", shape, nprocx, nprocy)
            if field == "vt":
                return _processor_boundary_band("tr_v", shape, nprocx, nprocy)
            return np.zeros(shape, dtype=bool)

        self._assert_fields_match("serial", "xy_split", ("ut", "vt"), mask_builder)

    def test_xy_split_matches_serial_globally(self) -> None:
        def mask_builder(field: str, shape: Tuple[int, int, int], nprocx: int, nprocy: int) -> np.ndarray:
            return np.ones(shape, dtype=bool)

        self._assert_fields_match(
            "serial",
            "xy_split",
            NO_TREE_FIELDS,
            mask_builder,
            atol=NO_TREE_GLOBAL_ABS_TOL,
        )

class TestTreeProcessorBoundaryParity(_BaseProcessorBoundaryParity):
    CASE_ID = TREE_CASE_ID
    PREFIX = "treedump"
    FIELDS = TREE_FIELDS
    VEGETATION_ENABLED = True

    def test_fixture_places_vegetation_on_y_processor_boundary(self) -> None:
        y_boundary = self.veg_mask.shape[1] // 2
        self.assertTrue(
            np.any(self.veg_mask[:, y_boundary - 1, :]) or np.any(self.veg_mask[:, y_boundary, :]),
            f"Case {self.CASE_ID} should place vegetation on the y processor boundary for 1x2/2x2 splits",
        )

    def test_fixture_does_not_place_vegetation_on_x_processor_boundary(self) -> None:
        x_boundary = self.veg_mask.shape[2] // 2
        self.assertFalse(
            np.any(self.veg_mask[:, :, x_boundary - 1]) or np.any(self.veg_mask[:, :, x_boundary]),
            f"Case {self.CASE_ID} unexpectedly places vegetation on the x processor boundary",
        )

    def test_x_split_matches_serial_on_u_support(self) -> None:
        def mask_builder(field: str, shape: Tuple[int, int, int], nprocx: int, nprocy: int) -> np.ndarray:
            support = _support_mask(field, self.veg_mask)
            if field == "tr_u":
                return support
            return np.zeros(shape, dtype=bool)

        self._assert_fields_match("serial", "x_split", ("tr_u",), mask_builder)

    def test_y_split_matches_serial_on_v_support(self) -> None:
        def mask_builder(field: str, shape: Tuple[int, int, int], nprocx: int, nprocy: int) -> np.ndarray:
            support = _support_mask(field, self.veg_mask)
            if field == "tr_v":
                return support
            return np.zeros(shape, dtype=bool)

        self._assert_fields_match("serial", "y_split", ("tr_v",), mask_builder)

    def test_xy_split_matches_serial_on_boundary_support(self) -> None:
        def mask_builder(field: str, shape: Tuple[int, int, int], nprocx: int, nprocy: int) -> np.ndarray:
            support = _support_mask(field, self.veg_mask)
            band = _processor_boundary_band(field, shape, nprocx, nprocy)
            return support & band

        self._assert_fields_match("serial", "xy_split", ("tr_u", "tr_v"), mask_builder)

    def test_xy_split_matches_serial_globally(self) -> None:
        def mask_builder(field: str, shape: Tuple[int, int, int], nprocx: int, nprocy: int) -> np.ndarray:
            return _support_mask(field, self.veg_mask)

        self._assert_fields_match("serial", "xy_split", TREE_FIELDS, mask_builder)


if __name__ == "__main__":
    unittest.main()
