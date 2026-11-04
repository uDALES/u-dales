"""View3D interface for computing view factors between facets.

Wraps the external View3D executable: converts STL geometry to the
View3D input format, runs the solver, and reads/writes view-factor
and sky-view-factor files used by the radiation preprocessing.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
import os
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Iterable, Optional

import numpy as np
from scipy import sparse

try:
    import netCDF4
except ImportError:  # pragma: no cover - optional runtime dependency
    netCDF4 = None

from exceptions import DependencyError
try:
    import trimesh
except ImportError as exc:  # pragma: no cover - required for STL conversion
    raise DependencyError("trimesh is required for View3D geometry export") from exc


@dataclass(frozen=True)
class View3DRunStats:
    """Runtime diagnostics from one external View3D process."""

    elapsed_seconds: float
    peak_rss_kb: int | None
    returncode: int


@dataclass(frozen=True)
class ViewFactorValidationReport:
    """Physical diagnostics for a facet/sky view-factor set."""

    nfacets: int
    nnz: int
    nonfinite_matrix: int
    nonfinite_sky: int
    negative_matrix: int
    above_one_matrix: int
    sky_outside: int
    bad_closure_rows: int
    matrix_min: float
    matrix_max: float
    sky_min: float
    sky_max: float
    closure_min: float
    closure_max: float
    reciprocity_l1_relative: float | None = None
    bad_reciprocity: bool = False

    @property
    def is_valid(self) -> bool:
        return not self.problem_messages()

    @property
    def closure_only(self) -> bool:
        return bool(self.bad_closure_rows) and not any(
            (
                self.nonfinite_matrix,
                self.nonfinite_sky,
                self.negative_matrix,
                self.above_one_matrix,
                self.sky_outside,
                self.bad_reciprocity,
            )
        )

    def problem_messages(self, tolerance: float = 1.0e-3) -> list[str]:
        problems: list[str] = []
        if self.nonfinite_matrix or self.nonfinite_sky:
            problems.append(
                "non-finite values "
                f"(matrix={self.nonfinite_matrix}, sky={self.nonfinite_sky})"
            )
        if self.negative_matrix:
            problems.append(
                f"{self.negative_matrix} matrix entries below {-tolerance:g} "
                f"(minimum={self.matrix_min:.6g})"
            )
        if self.above_one_matrix:
            problems.append(
                f"{self.above_one_matrix} matrix entries above {1.0 + tolerance:g} "
                f"(maximum={self.matrix_max:.6g})"
            )
        if self.sky_outside:
            problems.append(
                f"{self.sky_outside} sky factors outside physical bounds "
                f"(range={self.sky_min:.6g} to {self.sky_max:.6g})"
            )
        if self.bad_closure_rows:
            problems.append(
                f"{self.bad_closure_rows} rows fail sum(Fij) + SVF = 1 within "
                f"{tolerance:g} (range={self.closure_min:.6g} to "
                f"{self.closure_max:.6g})"
            )
        if self.bad_reciprocity and self.reciprocity_l1_relative is not None:
            problems.append(
                "reciprocity L1 error "
                f"{self.reciprocity_l1_relative:.6g} exceeds the allowed limit"
            )
        return problems


class ViewFactorValidationError(ValueError):
    """Raised when view-factor validation finds non-physical values."""

    def __init__(self, report: ViewFactorValidationReport, tolerance: float):
        self.report = report
        self.tolerance = tolerance
        super().__init__(
            "Invalid View3D factors: "
            + "; ".join(report.problem_messages(tolerance))
            + ". Repair or regenerate the view factors before radiation preprocessing."
        )


@dataclass(frozen=True)
class ViewFactorRepairLimits:
    """Safety limits that prevent automatic repair from masking bad geometry."""

    max_row_sum: float = 2.0
    max_overfull_area_fraction: float = 0.05
    max_exchange_area_reduction_fraction: float = 0.01
    max_reciprocity_l1_relative: float = 0.01

    def __post_init__(self) -> None:
        values = (
            self.max_row_sum,
            self.max_overfull_area_fraction,
            self.max_exchange_area_reduction_fraction,
            self.max_reciprocity_l1_relative,
        )
        if not all(np.isfinite(value) for value in values):
            raise ValueError("View-factor repair limits must be finite")
        if self.max_row_sum <= 1.0:
            raise ValueError("max_row_sum must be greater than one")
        for name, value in (
            ("max_overfull_area_fraction", self.max_overfull_area_fraction),
            (
                "max_exchange_area_reduction_fraction",
                self.max_exchange_area_reduction_fraction,
            ),
        ):
            if value < 0.0 or value > 1.0:
                raise ValueError(f"{name} must lie between zero and one")
        if self.max_reciprocity_l1_relative <= 0.0:
            raise ValueError("max_reciprocity_l1_relative must be positive")


@dataclass(frozen=True)
class ViewFactorRepairReport:
    """Summary of one open-domain view-factor repair."""

    algorithm: str
    repaired: bool
    nfacets: int
    nnz_before: int
    nnz_after: int
    overfull_rows: int
    materially_overfull_rows: int
    clamped_negative_entries: int
    scaled_entries: int
    materially_overfull_area_fraction: float
    exchange_area_reduction_fraction: float
    maximum_local_reduction_fraction: float
    max_row_sum_before: float
    max_row_sum_after: float
    reciprocity_l1_relative_before: float
    reciprocity_l1_relative_after: float


class ViewFactorRepairError(ValueError):
    """Raised when an automatic view-factor repair would be too extensive."""


_VIEW_FACTOR_REPAIR_ALGORITHM = "reciprocal-open-domain-v1"


_VIEW3D_CONTROL_ENV_KEYS = (
    "VIEW3D_EXE",
    "VIEW3D_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VIEW3D_MAX_DENSE_MATRIX_GIB",
    "VIEW3D_DISABLE_OPENMP",
    "VIEW3D_DISABLE_SPARSE_DIRECT",
    "VIEW3D_DISABLE_DENSE_MEMORY_GUARD",
)


def default_view3d_config_path() -> Path:
    """Return the repo-level View3D runtime configuration path."""

    return Path(__file__).resolve().parents[2] / "view3d_config.sh"


def _parse_env0(data: bytes) -> dict[str, str]:
    env: dict[str, str] = {}
    for item in data.split(b"\0"):
        if not item:
            continue
        key, sep, value = item.partition(b"=")
        if not sep:
            continue
        env[key.decode(errors="surrogateescape")] = value.decode(errors="surrogateescape")
    return env


def _source_view3d_config(config_path: Path, base_env: dict[str, str]) -> dict[str, str]:
    bash = shutil.which("bash")
    if bash is None:
        print(
            f"[view3d] config skipped because bash is unavailable: {config_path}",
            file=sys.stderr,
            flush=True,
        )
        return dict(base_env)

    env = dict(base_env)
    env["VIEW3D_CONFIG"] = str(config_path)
    command = [
        bash,
        "-c",
        (
            "set -a; "
            'cd "$(dirname "$1")" || exit 1; '
            'source "$1" >/dev/null && '
            "env -0"
        ),
        "bash",
        str(config_path),
    ]
    result = subprocess.run(command, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        stderr = result.stderr.decode(errors="replace").strip()
        msg = f"Failed to source View3D config {config_path}"
        if stderr:
            msg = f"{msg}: {stderr}"
        raise RuntimeError(msg)
    return _parse_env0(result.stdout)


def load_view3d_runtime_env(
    base_env: dict[str, str] | None = None,
    config_path: str | Path | None = None,
) -> tuple[dict[str, str], Path | None]:
    """
    Return the environment used for View3D, after sourcing the config file.
    """

    env = dict(os.environ if base_env is None else base_env)
    explicit_config = config_path is not None or bool(env.get("VIEW3D_CONFIG"))
    if config_path is None:
        config_value = env.get("VIEW3D_CONFIG")
        path = Path(config_value).expanduser() if config_value else default_view3d_config_path()
    else:
        path = Path(config_path).expanduser()

    path = path.resolve()
    if not path.exists():
        if explicit_config:
            raise FileNotFoundError(f"View3D config file not found: {path}")
        return env, None

    return _source_view3d_config(path, env), path


def _format_view3d_controls(env: dict[str, str]) -> str:
    controls = [f"{key}={env[key]}" for key in _VIEW3D_CONTROL_ENV_KEYS if env.get(key)]
    return ", ".join(controls)


def _read_proc_memory_kb(pid: int) -> int | None:
    """
    Read the best available RSS-like memory value for a running Linux process.
    """
    status_path = Path("/proc") / str(pid) / "status"
    values: dict[str, int] = {}
    try:
        with status_path.open("r", encoding="ascii") as f:
            for line in f:
                name, _, rest = line.partition(":")
                if name in {"VmHWM", "VmRSS"}:
                    parts = rest.strip().split()
                    if parts:
                        values[name] = int(parts[0])
    except (FileNotFoundError, PermissionError, ProcessLookupError, ValueError):
        return None

    return values.get("VmHWM") or values.get("VmRSS")


def _format_peak_rss(peak_rss_kb: int | None) -> str:
    if peak_rss_kb is None:
        return "unavailable"
    return f"{peak_rss_kb} kB ({peak_rss_kb / 1024.0**2:.3f} GiB)"


def _rusage_maxrss_kb(ru_maxrss: int) -> int:
    if sys.platform == "darwin":
        return int((ru_maxrss + 1023) // 1024)
    return int(ru_maxrss)


def resolve_view3d_exe(override: str | Path | None = None) -> Path:
    """
    Resolve the View3D executable path.

    Priority:
    1) explicit override
    2) VIEW3D_EXE environment variable
    3) tools/View3D defaults (Windows vs. others)
    """
    if override is not None:
        return Path(override).expanduser().resolve()

    env, _ = load_view3d_runtime_env()
    env_path = env.get("VIEW3D_EXE")
    if env_path:
        return Path(env_path).expanduser().resolve()

    tools_dir = Path(__file__).resolve().parents[2]
    candidates = []
    if os.name == "nt":
        candidates.append(tools_dir / "View3D" / "src" / "View3D.exe")
    else:
        candidates.append(tools_dir / "preprocessing" / "build" / "bin" / "view3d")
        candidates.append(tools_dir / "View3D" / "build" / "src" / "view3d")
        candidates.append(tools_dir / "View3D" / "src" / "View3D")

    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def stl_to_view3d(
    stl_path: str | Path,
    out_path: str | Path,
    outformat: int,
    maxD: float = np.inf,
    row: int = 0,
    col: int = 0,
) -> Path:
    """
    Convert STL geometry to View3D .vs3 format (STLtoView3D.m equivalent).
    """
    stl_path = Path(stl_path)
    out_path = Path(out_path)

    mesh = trimesh.load_mesh(str(stl_path))
    if isinstance(mesh, trimesh.Scene):
        mesh = trimesh.util.concatenate(
            tuple(trimesh.Trimesh(vertices=g.vertices, faces=g.faces) for g in mesh.geometry.values())
        )

    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=int) + 1  # View3D expects 1-based indices

    # Match tools/SEB/STLtoView3D.m byte-for-byte as closely as practical:
    # CRLF line endings, 6-decimal vertex formatting, and the historical
    # trailing literal "f" on each surface line. View3D is sensitive enough
    # that keeping the export contract aligned avoids tiny sparse-VF deltas on
    # larger cases such as 065.
    with out_path.open("w", encoding="ascii", newline="") as f:
        f.write("T\r\n")
        if maxD < np.inf:
            f.write(f"C out={outformat} maxD={maxD} row={row} col={col}\r\n")
        else:
            f.write(f"C out={outformat} row={row} col={col}\r\n")
        f.write("F 3\r\n")
        f.write("!    #      x      y      z\r\n")
        for idx, (x, y, z) in enumerate(vertices, start=1):
            f.write(f"V {idx:4d} {x:6.6f} {y:6.6f} {z:6.6f}\r\n")
        f.write("!    #     v1     v2     v3     v4   base    cmb   emit   name\r\n")
        for idx, (v1, v2, v3) in enumerate(faces, start=1):
            f.write(
                f"S {idx:4d} {v1:6d} {v2:6d} {v3:6d} {0:6d} {0:6d} {0:6d} {0:6d} {idx:6d}f\r\n"
            )
        f.write("End of Data\r\n")

    return out_path


def run_view3d(
    view3d_exe: str | Path,
    vs3_path: str | Path,
    out_path: str | Path,
    check: bool = True,
    nfacets: int | None = None,
    memory_poll_interval: float = 0.10,
) -> View3DRunStats:
    """
    Execute View3D as an external process.
    """
    cmd = [str(view3d_exe), str(vs3_path), str(out_path)]
    if nfacets is not None:
        print(f"[view3d] facets: {nfacets}", flush=True)
    env, config_path = load_view3d_runtime_env()
    if config_path is not None:
        print(f"[view3d] config: {config_path}", flush=True)
    controls = _format_view3d_controls(env)
    if controls:
        print(f"[view3d] controls: {controls}", flush=True)

    start = time.perf_counter()
    proc = subprocess.Popen(cmd, env=env)
    peak_rss_kb = _read_proc_memory_kb(proc.pid)
    wait4_available = hasattr(os, "wait4") and hasattr(os, "waitstatus_to_exitcode")

    try:
        if wait4_available:
            while True:
                pid, status, rusage = os.wait4(proc.pid, os.WNOHANG)
                sample_kb = _read_proc_memory_kb(proc.pid)
                if sample_kb is not None and (peak_rss_kb is None or sample_kb > peak_rss_kb):
                    peak_rss_kb = sample_kb
                if pid == proc.pid:
                    proc.returncode = os.waitstatus_to_exitcode(status)
                    rusage_peak_kb = _rusage_maxrss_kb(rusage.ru_maxrss)
                    if peak_rss_kb is None or rusage_peak_kb > peak_rss_kb:
                        peak_rss_kb = rusage_peak_kb
                    break
                time.sleep(memory_poll_interval)
        else:
            while proc.poll() is None:
                sample_kb = _read_proc_memory_kb(proc.pid)
                if sample_kb is not None and (peak_rss_kb is None or sample_kb > peak_rss_kb):
                    peak_rss_kb = sample_kb
                time.sleep(memory_poll_interval)
    except BaseException:
        try:
            proc.kill()
        except ProcessLookupError:
            pass
        if wait4_available:
            try:
                os.wait4(proc.pid, 0)
            except ChildProcessError:
                pass
        else:
            proc.wait()
        raise

    elapsed_seconds = time.perf_counter() - start
    stats = View3DRunStats(
        elapsed_seconds=elapsed_seconds,
        peak_rss_kb=peak_rss_kb,
        returncode=proc.returncode,
    )
    print(f"[view3d] runtime: {stats.elapsed_seconds:.3f} s", flush=True)
    print(f"[view3d] peak memory: {_format_peak_rss(stats.peak_rss_kb)}", flush=True)

    if check and proc.returncode:
        raise subprocess.CalledProcessError(proc.returncode, cmd)
    return stats


def count_sparse_entries(path: str | Path) -> int:
    """
    Count non-empty lines in a sparse View3D output file.
    """
    path = Path(path)
    with path.open("rb") as f:
        return sum(1 for line in f if line.strip())


def read_view3d_output(
    out_path: str | Path,
    nfacets: int,
    outformat: int,
    one_based: bool = True,
) -> sparse.csr_matrix:
    """
    Read View3D output into a sparse matrix.

    outformat:
      0: text
      1: binary
      2: sparse text
    """
    out_path = Path(out_path)

    if outformat == 0:
        data = np.loadtxt(out_path, skiprows=2)
        if data.size == 0:
            return sparse.csr_matrix((nfacets, nfacets))
        data = data[:-1, :]  # drop trailing line
        return sparse.csr_matrix(data, shape=(nfacets, nfacets))

    if outformat == 1:
        with out_path.open("rb") as f:
            _ = np.fromfile(f, dtype=np.float32, count=8 + nfacets)
            raw = np.fromfile(f, dtype=np.float32, count=nfacets * nfacets)
        # MATLAB reads this as reshape(raw, [n, n])' after fread. Since
        # MATLAB reshape is column-major, NumPy's default row-major reshape
        # is already the equivalent final orientation. fread(..., 'single')
        # returns doubles by default, so promote the decoded single-precision
        # bytes before downstream sums/reflections.
        vf = raw.reshape((nfacets, nfacets)).astype(float)
        return sparse.csr_matrix(vf)

    if outformat == 2:
        ijs = np.loadtxt(out_path)
        if ijs.size == 0:
            return sparse.csr_matrix((nfacets, nfacets))
        if ijs.ndim == 1:
            ijs = ijs.reshape(1, -1)
        rows = ijs[:, 0].astype(int)
        cols = ijs[:, 1].astype(int)
        vals = ijs[:, 2].astype(float)
        if one_based:
            rows = rows - 1
            cols = cols - 1
        return sparse.csr_matrix((vals, (rows, cols)), shape=(nfacets, nfacets))

    raise ValueError(f"Unsupported view3d output format: {outformat}")


def compute_svf(vf: sparse.spmatrix) -> np.ndarray:
    """
    Compute sky view factors from view factor matrix.
    """
    row_sum = np.asarray(vf.sum(axis=1)).reshape(-1)
    return np.maximum(1.0 - row_sum, 0.0)


def _reciprocity_l1_relative(
    vf: sparse.spmatrix,
    areas: np.ndarray,
) -> float:
    exchange = sparse.diags(areas) @ sparse.csr_matrix(vf)
    denominator = float(np.sum(np.abs(exchange.data)))
    if denominator == 0.0:
        return 0.0
    residual = exchange - exchange.T
    return float(np.sum(np.abs(residual.data))) / denominator


def inspect_view_factors(
    vf: sparse.spmatrix,
    svf: np.ndarray,
    *,
    areas: np.ndarray | None = None,
    tolerance: float = 1.0e-3,
    reciprocity_tolerance: float = 1.0e-2,
) -> ViewFactorValidationReport:
    """Return structured physical diagnostics without raising."""
    if tolerance <= 0.0:
        raise ValueError("view-factor validation tolerance must be positive")
    if reciprocity_tolerance <= 0.0:
        raise ValueError("view-factor reciprocity tolerance must be positive")

    matrix = sparse.csr_matrix(vf)
    sky = np.asarray(svf, dtype=float).reshape(-1)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"View-factor matrix must be square; got {matrix.shape}")
    if sky.shape != (matrix.shape[0],):
        raise ValueError(
            f"Sky view factors must have shape ({matrix.shape[0]},); got {sky.shape}"
        )

    data = matrix.data
    finite_data = data[np.isfinite(data)]
    finite_sky = sky[np.isfinite(sky)]
    nonfinite_matrix = int(np.count_nonzero(~np.isfinite(data)))
    nonfinite_sky = int(np.count_nonzero(~np.isfinite(sky)))
    matrix_min = float(finite_data.min()) if finite_data.size else 0.0
    matrix_max = float(finite_data.max()) if finite_data.size else 0.0
    sky_min = float(finite_sky.min()) if finite_sky.size else 0.0
    sky_max = float(finite_sky.max()) if finite_sky.size else 0.0
    negative_matrix = int(np.count_nonzero(finite_data < -tolerance))
    above_one_matrix = int(np.count_nonzero(finite_data > 1.0 + tolerance))
    sky_outside = int(
        np.count_nonzero(
            (finite_sky < -tolerance) | (finite_sky > 1.0 + tolerance)
        )
    )

    bad_closure_rows = 0
    closure_min = float("nan")
    closure_max = float("nan")
    if not nonfinite_matrix and not nonfinite_sky:
        closure = np.asarray(matrix.sum(axis=1)).reshape(-1) + sky
        if closure.size:
            closure_min = float(closure.min())
            closure_max = float(closure.max())
        else:
            closure_min = closure_max = 1.0
        bad_closure_rows = int(
            np.count_nonzero(np.abs(closure - 1.0) > tolerance)
        )

    reciprocity = None
    if areas is not None:
        area_values = np.asarray(areas, dtype=float).reshape(-1)
        if area_values.shape != (matrix.shape[0],):
            raise ValueError(
                f"Facet areas must have shape ({matrix.shape[0]},); "
                f"got {area_values.shape}"
            )
        if np.any(~np.isfinite(area_values)) or np.any(area_values <= 0.0):
            raise ValueError("Facet areas must be finite and strictly positive")
        if not nonfinite_matrix:
            reciprocity = _reciprocity_l1_relative(matrix, area_values)

    return ViewFactorValidationReport(
        nfacets=matrix.shape[0],
        nnz=int(matrix.nnz),
        nonfinite_matrix=nonfinite_matrix,
        nonfinite_sky=nonfinite_sky,
        negative_matrix=negative_matrix,
        above_one_matrix=above_one_matrix,
        sky_outside=sky_outside,
        bad_closure_rows=bad_closure_rows,
        matrix_min=matrix_min,
        matrix_max=matrix_max,
        sky_min=sky_min,
        sky_max=sky_max,
        closure_min=closure_min,
        closure_max=closure_max,
        reciprocity_l1_relative=reciprocity,
        bad_reciprocity=(
            reciprocity is not None and reciprocity > reciprocity_tolerance
        ),
    )


def validate_view_factors(
    vf: sparse.spmatrix,
    svf: np.ndarray,
    *,
    areas: np.ndarray | None = None,
    tolerance: float = 1.0e-3,
    reciprocity_tolerance: float = 1.0e-2,
) -> ViewFactorValidationReport:
    """Reject non-physical facet and sky view-factor combinations."""
    report = inspect_view_factors(
        vf,
        svf,
        areas=areas,
        tolerance=tolerance,
        reciprocity_tolerance=reciprocity_tolerance,
    )
    if not report.is_valid:
        raise ViewFactorValidationError(report, tolerance)
    return report


def repair_view_factors(
    vf: sparse.spmatrix,
    areas: np.ndarray,
    *,
    tolerance: float = 1.0e-3,
    closure_margin: float = 5.0e-4,
    limits: ViewFactorRepairLimits | None = None,
) -> tuple[sparse.csr_matrix, np.ndarray, ViewFactorRepairReport]:
    """Conservatively repair overfull open-domain view-factor rows.

    The correction is performed on reciprocal interchange areas ``A_i F_ij``.
    Pair estimates are symmetrized, then every edge incident on an overfull row
    is reduced by the stricter endpoint scale. The residual is assigned to sky.
    """
    if tolerance <= 0.0:
        raise ValueError("view-factor repair tolerance must be positive")
    if closure_margin <= 0.0 or closure_margin >= 1.0:
        raise ValueError("view-factor closure margin must lie between zero and one")
    limits = limits or ViewFactorRepairLimits()

    matrix = sparse.csr_matrix(vf, dtype=float).copy()
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"View-factor matrix must be square; got {matrix.shape}")
    matrix.sum_duplicates()
    matrix.sort_indices()

    area_values = np.asarray(areas, dtype=float).reshape(-1)
    if area_values.shape != (matrix.shape[0],):
        raise ValueError(
            f"Facet areas must have shape ({matrix.shape[0]},); got {area_values.shape}"
        )
    if np.any(~np.isfinite(area_values)) or np.any(area_values <= 0.0):
        raise ValueError("Facet areas must be finite and strictly positive")

    data = matrix.data
    if np.any(~np.isfinite(data)):
        raise ViewFactorRepairError("Cannot repair non-finite view factors")
    materially_negative = int(np.count_nonzero(data < -tolerance))
    above_one = int(np.count_nonzero(data > 1.0 + tolerance))
    if materially_negative:
        raise ViewFactorRepairError(
            f"Cannot repair {materially_negative} view factors below {-tolerance:g}"
        )
    if above_one:
        raise ViewFactorRepairError(
            f"Cannot repair {above_one} individual view factors above "
            f"{1.0 + tolerance:g}"
        )
    tiny_negative = data < 0.0
    n_tiny_negative = int(np.count_nonzero(tiny_negative))
    if n_tiny_negative:
        data[tiny_negative] = 0.0
        matrix.eliminate_zeros()

    reciprocity_before = _reciprocity_l1_relative(matrix, area_values)
    if reciprocity_before > limits.max_reciprocity_l1_relative:
        raise ViewFactorRepairError(
            "Cannot safely repair view factors with reciprocity L1 error "
            f"{reciprocity_before:.6g}; limit is "
            f"{limits.max_reciprocity_l1_relative:.6g}"
        )

    row_sum_raw = np.asarray(matrix.sum(axis=1)).reshape(-1)
    max_row_sum_raw = float(row_sum_raw.max()) if row_sum_raw.size else 0.0
    if max_row_sum_raw > limits.max_row_sum:
        raise ViewFactorRepairError(
            f"Maximum view-factor row sum {max_row_sum_raw:.6g} exceeds "
            f"automatic-repair limit {limits.max_row_sum:.6g}"
        )

    raw_materially_overfull = row_sum_raw > 1.0 + tolerance
    if not np.any(raw_materially_overfull):
        sky = np.clip(1.0 - row_sum_raw, 0.0, 1.0)
        report = ViewFactorRepairReport(
            algorithm=_VIEW_FACTOR_REPAIR_ALGORITHM,
            repaired=bool(n_tiny_negative),
            nfacets=matrix.shape[0],
            nnz_before=int(sparse.csr_matrix(vf).nnz),
            nnz_after=int(matrix.nnz),
            overfull_rows=int(np.count_nonzero(row_sum_raw > 1.0)),
            materially_overfull_rows=0,
            clamped_negative_entries=n_tiny_negative,
            scaled_entries=0,
            materially_overfull_area_fraction=0.0,
            exchange_area_reduction_fraction=0.0,
            maximum_local_reduction_fraction=0.0,
            max_row_sum_before=max_row_sum_raw,
            max_row_sum_after=max_row_sum_raw,
            reciprocity_l1_relative_before=reciprocity_before,
            reciprocity_l1_relative_after=reciprocity_before,
        )
        validate_view_factors(
            matrix,
            sky,
            areas=area_values,
            tolerance=tolerance,
            reciprocity_tolerance=limits.max_reciprocity_l1_relative,
        )
        return matrix, sky, report

    exchange = (sparse.diags(area_values) @ matrix).tocsr()
    presence = exchange.copy()
    presence.data = np.ones_like(presence.data)
    counts = (presence + presence.T).tocsr()
    inverse_counts = counts.copy()
    inverse_counts.data = 1.0 / inverse_counts.data
    symmetric_exchange = ((exchange + exchange.T).multiply(inverse_counts)).tocsr()
    symmetric_exchange.eliminate_zeros()

    row_sum = (
        np.asarray(symmetric_exchange.sum(axis=1)).reshape(-1) / area_values
    )
    overfull = row_sum > 1.0
    materially_overfull = row_sum > 1.0 + tolerance
    max_row_sum_before = float(row_sum.max()) if row_sum.size else 0.0
    if max_row_sum_before > limits.max_row_sum:
        raise ViewFactorRepairError(
            f"Reciprocity adjustment gives maximum row sum {max_row_sum_before:.6g}, "
            f"above automatic-repair limit {limits.max_row_sum:.6g}"
        )

    total_area = float(np.sum(area_values))
    materially_overfull_area_fraction = (
        float(np.sum(area_values[materially_overfull])) / total_area
        if total_area
        else 0.0
    )
    if materially_overfull_area_fraction > limits.max_overfull_area_fraction:
        raise ViewFactorRepairError(
            "Materially overfull facets represent "
            f"{materially_overfull_area_fraction:.3%} of facet area; "
            f"automatic-repair limit is {limits.max_overfull_area_fraction:.3%}"
        )

    scale = np.ones(matrix.shape[0], dtype=float)
    scale[overfull] = (1.0 - closure_margin) / row_sum[overfull]
    coo = symmetric_exchange.tocoo(copy=True)
    edge_scale = np.minimum(scale[coo.row], scale[coo.col])
    changed = edge_scale < 1.0
    coo.data *= edge_scale
    repaired_exchange = coo.tocsr()

    exchange_before = float(np.sum(symmetric_exchange.data))
    exchange_removed = exchange_before - float(np.sum(repaired_exchange.data))
    reduction_fraction = (
        max(exchange_removed, 0.0) / exchange_before if exchange_before else 0.0
    )
    if reduction_fraction > limits.max_exchange_area_reduction_fraction:
        raise ViewFactorRepairError(
            f"Repair would remove {reduction_fraction:.3%} of surface exchange area; "
            "automatic-repair limit is "
            f"{limits.max_exchange_area_reduction_fraction:.3%}"
        )

    repaired = (sparse.diags(1.0 / area_values) @ repaired_exchange).tocsr()
    repaired.eliminate_zeros()
    row_sum_after = np.asarray(repaired.sum(axis=1)).reshape(-1)
    sky = np.clip(1.0 - row_sum_after, 0.0, 1.0)
    reciprocity_after = _reciprocity_l1_relative(repaired, area_values)

    report = ViewFactorRepairReport(
        algorithm=_VIEW_FACTOR_REPAIR_ALGORITHM,
        repaired=True,
        nfacets=matrix.shape[0],
        nnz_before=int(sparse.csr_matrix(vf).nnz),
        nnz_after=int(repaired.nnz),
        overfull_rows=int(np.count_nonzero(overfull)),
        materially_overfull_rows=int(np.count_nonzero(materially_overfull)),
        clamped_negative_entries=n_tiny_negative,
        scaled_entries=int(np.count_nonzero(changed)),
        materially_overfull_area_fraction=materially_overfull_area_fraction,
        exchange_area_reduction_fraction=reduction_fraction,
        maximum_local_reduction_fraction=float(1.0 - scale.min()),
        max_row_sum_before=max_row_sum_before,
        max_row_sum_after=(
            float(row_sum_after.max()) if row_sum_after.size else 0.0
        ),
        reciprocity_l1_relative_before=reciprocity_before,
        reciprocity_l1_relative_after=reciprocity_after,
    )
    validate_view_factors(
        repaired,
        sky,
        areas=area_values,
        tolerance=tolerance,
        reciprocity_tolerance=limits.max_reciprocity_l1_relative,
    )
    return repaired, sky, report


def _atomic_output(path: Path, writer) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        writer(tmp_path)
        os.replace(tmp_path, path)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise
    return path


def write_view_factor_repair_report(
    path: str | Path,
    report: ViewFactorRepairReport,
) -> Path:
    """Write machine-readable repair diagnostics atomically."""
    path = Path(path)

    def _write(tmp_path: Path) -> None:
        tmp_path.write_text(
            json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
            encoding="ascii",
        )

    return _atomic_output(path, _write)


def write_svf(path: str | Path, svf: np.ndarray) -> Path:
    """
    Write sky view factors to svf.inp.* format.
    """
    path = Path(path)

    def _write(tmp_path: Path) -> None:
        with tmp_path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# sky view factors\n")
            np.savetxt(f, svf, fmt="%.8f")

    return _atomic_output(path, _write)


def write_vfsparse(path: str | Path, vf: sparse.spmatrix, threshold: float = 5e-7) -> Path:
    """
    Write sparse view factors (vfsparse.inp.*).

    Values are rounded to 8 decimal places and sorted by (row, col).
    """
    path = Path(path)
    vf = vf.tocoo()
    mask = vf.data >= threshold
    rows = vf.row[mask] + 1
    cols = vf.col[mask] + 1
    vals = vf.data[mask]

    def _write(tmp_path: Path) -> None:
        if rows.size == 0:
            tmp_path.write_text("", encoding="ascii")
            return
        stacked = np.column_stack([rows, cols, vals])
        stacked = stacked[np.lexsort((stacked[:, 1], stacked[:, 0]))]
        np.savetxt(tmp_path, stacked, fmt="%d %d %.8f")

    return _atomic_output(path, _write)


def write_vf(path: str | Path, vf: np.ndarray) -> Path:
    """
    Write full view factor matrix to NetCDF (vf.nc.inp.*).

    MATLAB's low-level netcdf.putVar writes the in-memory matrix with the
    opposite dimension order as seen by Python netCDF readers, so store the
    transpose here to keep the generated file byte-contract compatible with
    the MATLAB preprocessing route.
    """
    if netCDF4 is None:
        raise DependencyError("netCDF4 is required to write full view factor matrices")

    path = Path(path)

    def _write(tmp_path: Path) -> None:
        with netCDF4.Dataset(tmp_path, "w", format="NETCDF4") as ds:
            n = vf.shape[0]
            ds.createDimension("rows", n)
            ds.createDimension("columns", n)
            var = ds.createVariable("view factor", "f4", ("rows", "columns"))
            var[:, :] = np.asarray(vf).T

    return _atomic_output(path, _write)
