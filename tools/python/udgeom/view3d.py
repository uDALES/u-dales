"""View3D interface for computing view factors between facets.

Wraps the external View3D executable: converts STL geometry to the
View3D input format, runs the solver, and reads/writes view-factor
and sky-view-factor files used by the radiation preprocessing.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import shutil
import subprocess
import sys
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


_VIEW3D_CONTROL_ENV_KEYS = (
    "VIEW3D_EXE",
    "VIEW3D_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VIEW3D_MAX_DENSE_MATRIX_GIB",
    "VIEW3D_MAX_DENSE_MATRIX_BYTES",
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
            'VIEW3D_CONFIG_DIR="$(cd "$(dirname "$1")" && pwd)"; '
            'source "$1" >/dev/null; '
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
        proc.kill()
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


def write_svf(path: str | Path, svf: np.ndarray) -> Path:
    """
    Write sky view factors to svf.inp.* format.
    """
    path = Path(path)
    with path.open("w", encoding="ascii", newline="\n") as f:
        f.write("# sky view factors\n")
    with path.open("a", encoding="ascii", newline="\n") as f:
        np.savetxt(f, svf, fmt="%4f")
    return path


def write_vfsparse(path: str | Path, vf: sparse.spmatrix, threshold: float = 5e-7) -> Path:
    """
    Write sparse view factors (vfsparse.inp.*).

    Values are rounded to 6 decimal places and sorted by (row, col).
    """
    path = Path(path)
    vf = vf.tocoo()
    mask = vf.data >= threshold
    rows = vf.row[mask] + 1
    cols = vf.col[mask] + 1
    vals = vf.data[mask]

    if rows.size == 0:
        path.write_text("", encoding="ascii")
        return path

    stacked = np.column_stack([rows, cols, vals])
    stacked = stacked[np.lexsort((stacked[:, 1], stacked[:, 0]))]
    np.savetxt(path, stacked, fmt="%d %d %.6f")
    return path


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
    with netCDF4.Dataset(path, "w", format="NETCDF4") as ds:
        n = vf.shape[0]
        ds.createDimension("rows", n)
        ds.createDimension("columns", n)
        var = ds.createVariable("view factor", "f4", ("rows", "columns"))
        var[:, :] = np.asarray(vf).T
    return path
