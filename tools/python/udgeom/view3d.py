from __future__ import annotations

from pathlib import Path
import os
import subprocess
from typing import Iterable, Optional

import numpy as np
from scipy import sparse

try:
    import netCDF4
except ImportError:  # pragma: no cover - optional runtime dependency
    netCDF4 = None

try:
    import trimesh
except ImportError as exc:  # pragma: no cover - required for STL conversion
    raise ImportError("trimesh is required for View3D geometry export") from exc


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

    env_path = os.environ.get("VIEW3D_EXE")
    if env_path:
        return Path(env_path).expanduser().resolve()

    tools_dir = Path(__file__).resolve().parents[2]
    if os.name == "nt":
        return tools_dir / "View3D" / "src" / "View3D.exe"
    return tools_dir / "View3D" / "build" / "src" / "view3d"


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

    with out_path.open("w", encoding="ascii", newline="\n") as f:
        f.write("T\n")
        if maxD < np.inf:
            f.write(f"C out={outformat} maxD={maxD} row={row} col={col}\n")
        else:
            f.write(f"C out={outformat} row={row} col={col}\n")
        f.write("F 3\n")
        f.write("!    #      x      y      z\n")
        for idx, (x, y, z) in enumerate(vertices, start=1):
            f.write(f"V {idx:4d} {x:6.6f} {y:6.6f} {z:6.6f}\n")
        f.write("!    #     v1     v2     v3     v4   base    cmb   emit   name\n")
        for idx, (v1, v2, v3) in enumerate(faces, start=1):
            f.write(f"S {idx:4d} {v1:6d} {v2:6d} {v3:6d} {0:6d} {0:6d} {0:6d} {0:6d} {idx:6d}\n")
        f.write("End of Data\n")

    return out_path


def run_view3d(
    view3d_exe: str | Path,
    vs3_path: str | Path,
    out_path: str | Path,
    check: bool = True,
) -> None:
    """
    Execute View3D as an external process.
    """
    cmd = [str(view3d_exe), str(vs3_path), str(out_path)]
    subprocess.run(cmd, check=check)


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
        vf = raw.reshape((nfacets, nfacets)).T
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
    """
    if netCDF4 is None:
        raise ImportError("netCDF4 is required to write full view factor matrices")

    path = Path(path)
    with netCDF4.Dataset(path, "w", format="NETCDF4") as ds:
        n = vf.shape[0]
        ds.createDimension("rows", n)
        ds.createDimension("columns", n)
        var = ds.createVariable("view factor", "f4", ("rows", "columns"))
        var[:, :] = vf
    return path

