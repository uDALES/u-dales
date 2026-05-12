import argparse
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
PYTHON_DIR = REPO_ROOT / "tools" / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udbase import UDBase
from udprep import UDPrep
from udprep.udprep_radiation import read_view3d_output
from udprep.solar import nsun_from_angles


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare direct shortwave f2py and standalone Fortran scanline backends."
    )
    parser.add_argument("case", help="Case number or path to a case directory")
    parser.add_argument(
        "--ray-density",
        type=float,
        default=12.0,
        help="Ray density used to derive scanline resolution from the grid",
    )
    parser.add_argument(
        "--azimuth-deg",
        type=float,
        default=None,
        help="Solar azimuth in degrees for the comparison; defaults to the case settings when omitted",
    )
    parser.add_argument(
        "--elevation-deg",
        type=float,
        default=None,
        help="Solar elevation in degrees for the comparison; defaults to the case settings when omitted",
    )
    parser.add_argument(
        "--irradiance",
        type=float,
        default=800.0,
        help="Direct normal irradiance in W/m^2",
    )
    parser.add_argument(
        "--reference-sdir",
        type=Path,
        default=None,
        help="Optional reference Sdir.txt to compare against both backends",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=None,
        help="Explicit scanline resolution; defaults to psc_res for ishortwave=1 cases or the ray-density-derived spacing otherwise",
    )
    parser.add_argument(
        "--include-knet",
        action="store_true",
        help="Also compare full shortwave preprocessing output (knet/netsw) using the two scanline backends.",
    )
    return parser.parse_args()


def _resolve_case(case_arg: str) -> tuple[str, Path]:
    case_path = Path(case_arg)
    if case_path.is_dir():
        expnr = case_path.name
        return expnr, case_path.resolve()
    case_dir = REPO_ROOT / "tests" / "cases" / case_arg
    if case_dir.is_dir():
        return case_arg, case_dir
    raise FileNotFoundError(f"Could not resolve case {case_arg!r} to a case directory")


def _sun_vector(azimuth_deg: float, elevation_deg: float) -> np.ndarray:
    az = np.deg2rad(azimuth_deg)
    el = np.deg2rad(elevation_deg)
    return np.array(
        [np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)],
        dtype=float,
    )


def _scanline_resolution(sim: UDBase, ray_density: float) -> float:
    cell_min = min(sim.dx, sim.dy, float(np.min(sim.dzt)))
    return 0.25 * cell_min / ray_density


def _case_solar_state(sim: UDBase) -> tuple[np.ndarray, float, float, float]:
    isolar = int(getattr(sim, "isolar", 1))
    if isolar != 1:
        raise ValueError(f"Only isolar=1 is supported by this comparison utility, got {isolar}")
    zenith_deg = float(getattr(sim, "solarzenith"))
    azimuth_deg = float(getattr(sim, "solarazimuth")) - float(getattr(sim, "xazimuth", 0.0))
    nsun = nsun_from_angles(zenith_deg, azimuth_deg)
    elevation_deg = 90.0 - zenith_deg
    return nsun, azimuth_deg, elevation_deg, float(getattr(sim, "I"))


def _run_fortran_scanline(
    repo_root: Path,
    faces: np.ndarray,
    centers: np.ndarray,
    normals: np.ndarray,
    vertices: np.ndarray,
    nsun: np.ndarray,
    irradiance: float,
    resolution: float,
) -> tuple[np.ndarray, float]:
    compiler = shutil.which("gfortran")
    if compiler is None:
        raise RuntimeError("gfortran not found")

    with tempfile.TemporaryDirectory(prefix="udales-scanline-fortran-") as td:
        td_path = Path(td)
        exe_path = td_path / "DS.exe"
        source = repo_root / "tools" / "SEB" / "directShortwave.f90"
        subprocess.run(
            [compiler, "-O3", str(source), "-o", str(exe_path)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        np.savetxt(td_path / "vertices.txt", vertices, fmt="%15.10f %15.10f %15.10f")
        face_rows = np.hstack([faces, centers, normals])
        np.savetxt(
            td_path / "faces.txt",
            face_rows,
            fmt="%8d %8d %8d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f",
        )
        with (td_path / "info_directShortwave.txt").open("w", encoding="ascii") as f:
            f.write(f"{len(faces):8d} {len(vertices):8d}\n")
            f.write(f"{nsun[0]:15.10f} {nsun[1]:15.10f} {nsun[2]:15.10f}\n")
            f.write(f"{irradiance:15.10f}\n")
            f.write(f"{resolution:15.10f}\n")

        t0 = time.perf_counter()
        subprocess.run(
            [str(exe_path)],
            cwd=td_path,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        elapsed = time.perf_counter() - t0
        return np.loadtxt(td_path / "Sdir.txt"), elapsed


def _run_f2py_scanline(
    faces: np.ndarray,
    centers: np.ndarray,
    normals: np.ndarray,
    vertices: np.ndarray,
    nsun: np.ndarray,
    irradiance: float,
    resolution: float,
) -> tuple[np.ndarray, float]:
    import udprep.directshortwave_f2py as dsroot

    dsmod = getattr(dsroot, "directshortwave_mod", dsroot)
    t0 = time.perf_counter()
    sdir = dsmod.calculate_direct_shortwave_f2py(
        np.asfortranarray(np.asarray(faces, dtype=np.int32)),
        np.asfortranarray(np.asarray(centers, dtype=np.float32)),
        np.asfortranarray(np.asarray(normals, dtype=np.float32)),
        np.asfortranarray(np.asarray(vertices, dtype=np.float32)),
        np.asarray(nsun, dtype=float),
        float(irradiance),
        float(resolution),
    )
    elapsed = time.perf_counter() - t0
    return np.asarray(sdir, dtype=float), elapsed


def _summarize(candidate: np.ndarray, reference: np.ndarray, areas: np.ndarray) -> dict[str, float]:
    diff = candidate - reference
    active = (reference > 1.0e-6) & (candidate > 1.0e-6)
    if np.any(active):
        corr = float(np.corrcoef(reference[active], candidate[active])[0, 1])
        p95_rel = float(np.quantile(np.abs(diff[active]) / np.maximum(reference[active], 1.0), 0.95))
    else:
        corr = 1.0
        p95_rel = 0.0
    return {
        "max_abs": float(np.max(np.abs(diff))),
        "mean_abs": float(np.mean(np.abs(diff))),
        "corr": corr,
        "p95_rel": p95_rel,
        "energy_ratio": float(np.sum(candidate * areas) / np.sum(reference * areas)),
    }


def _load_reference_vf_svf(prep: UDPrep):
    sim = prep.sim
    rad = prep.radiation
    expnr = getattr(sim, "expnr", "")
    out_dir = Path(sim.path)

    svf_path = out_dir / f"svf.inp.{expnr}"
    svf = np.loadtxt(svf_path)

    view3d_out = int(rad.view3d_out)
    if view3d_out == 0:
        vf_path = out_dir / "vf.txt"
    elif view3d_out == 1:
        vf_path = out_dir / "vf.bin"
    else:
        vf_path = out_dir / f"vfsparse.inp.{expnr}"
    vf = read_view3d_output(vf_path, nfacets=sim.geom.stl.faces.shape[0], outformat=view3d_out)
    return vf, svf


def main() -> int:
    args = _parse_args()
    expnr, case_dir = _resolve_case(args.case)

    prep = UDPrep(expnr, case_dir, load_geometry=True)
    sim = prep.sim
    mesh = sim.geom.stl
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=np.int32) + 1
    normals = np.asarray(mesh.face_normals, dtype=float)
    centers = np.asarray(sim.geom.face_incenters, dtype=float)
    case_nsun, case_azimuth_deg, case_elevation_deg, case_irradiance = _case_solar_state(sim)
    azimuth_deg = case_azimuth_deg if args.azimuth_deg is None else args.azimuth_deg
    elevation_deg = case_elevation_deg if args.elevation_deg is None else args.elevation_deg
    irradiance = case_irradiance if args.irradiance == 800.0 and args.azimuth_deg is None and args.elevation_deg is None else args.irradiance
    nsun = case_nsun if args.azimuth_deg is None and args.elevation_deg is None else _sun_vector(azimuth_deg, elevation_deg)
    if args.resolution is not None:
        resolution = float(args.resolution)
    elif int(getattr(sim, "ishortwave", 1)) == 1:
        resolution = float(getattr(sim, "psc_res", 0.1))
    else:
        resolution = _scanline_resolution(sim, args.ray_density)
    areas = np.asarray(mesh.area_faces, dtype=float)

    print(
        {
            "case": expnr,
            "case_dir": str(case_dir),
            "ray_density": args.ray_density,
            "resolution": resolution,
            "azimuth_deg": azimuth_deg,
            "elevation_deg": elevation_deg,
            "irradiance": irradiance,
        }
    )

    reference = None
    if args.reference_sdir is not None:
        reference = np.loadtxt(args.reference_sdir)
        print({"reference_sdir": str(args.reference_sdir), "reference_size": int(reference.size)})
    else:
        default_sdir = case_dir / "Sdir.txt"
        if default_sdir.exists():
            reference = np.loadtxt(default_sdir)
            print({"reference_sdir": str(default_sdir), "reference_size": int(reference.size)})

    try:
        sdir_fortran, fortran_seconds = _run_fortran_scanline(
            REPO_ROOT,
            faces,
            centers,
            normals,
            vertices,
            nsun,
            irradiance,
            resolution,
        )
    except subprocess.CalledProcessError as exc:
        print(
            {
                "fortran_error": f"standalone executable failed with return code {exc.returncode}",
            }
        )
        return 0

    sdir_f2py, f2py_seconds = _run_f2py_scanline(
        faces,
        centers,
        normals,
        vertices,
        nsun,
        irradiance,
        resolution,
    )
    print(
        {
            "fortran_seconds": round(fortran_seconds, 3),
            "f2py_seconds": round(f2py_seconds, 3),
            "fortran_vs_f2py": _summarize(sdir_f2py, sdir_fortran, areas),
        }
    )
    if reference is not None:
        print(
            {
                "fortran_seconds": round(fortran_seconds, 3),
                "f2py_seconds": round(f2py_seconds, 3),
                "fortran_vs_reference": _summarize(sdir_fortran, reference, areas),
                "f2py_vs_reference": _summarize(sdir_f2py, reference, areas),
            }
        )

    if args.include_knet:
        rad = prep.radiation
        nsun_case, _, _, irradiance_case, dsky = rad._solar_state_time(
            datetime(
                int(rad.year),
                int(rad.month),
                int(rad.day),
                int(rad.hour),
                int(rad.minute),
                int(rad.second),
            )
        )
        if not np.allclose(nsun, nsun_case, rtol=0.0, atol=1.0e-12):
            raise ValueError("knet comparison requires the case-default solar state")
        if abs(float(irradiance) - float(irradiance_case)) > 1.0e-12:
            raise ValueError("knet comparison requires the case-default irradiance")

        lscatter = bool(getattr(sim, "lEB", False))
        albedo = sim.assign_prop_to_fac("al")
        face_normals = sim.geom.stl.face_normals
        fss = (1.0 + face_normals[:, 2]) * 0.5 if not lscatter else None
        vf = svf = None
        if lscatter:
            vf, svf = _load_reference_vf_svf(prep)

        _, knet_fortran, _ = rad._compute_knet(
            nsun,
            irradiance,
            dsky,
            "scanline_legacy",
            resolution,
            lscatter,
            albedo,
            vf,
            svf,
            fss,
        )
        _, knet_f2py, _ = rad._compute_knet(
            nsun,
            irradiance,
            dsky,
            "scanline",
            resolution,
            lscatter,
            albedo,
            vf,
            svf,
            fss,
        )
        print({"knet_fortran_vs_f2py": _summarize(knet_f2py, knet_fortran, areas)})

        netsw_ref = case_dir / f"netsw.inp.{expnr}"
        if netsw_ref.exists():
            knet_reference = np.loadtxt(netsw_ref, comments="#")
            print({"reference_netsw": str(netsw_ref), "reference_netsw_size": int(knet_reference.size)})
            print({"fortran_knet_vs_reference": _summarize(knet_fortran, knet_reference, areas)})
            print({"f2py_knet_vs_reference": _summarize(knet_f2py, knet_reference, areas)})

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
