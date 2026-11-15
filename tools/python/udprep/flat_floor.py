"""Immersed-boundary input files for a flat floor made of facets.

A flat floor at z = 0 covers every column of the domain. The solver input it needs is
small and regular, so it is written directly here instead of through the STL pipeline
(``udprep_ibm``): the w points at z = 0 are solid, the first u, v and cell-centre points
(k = 1) are fluid boundary points, and each of them has one facet section of area dx*dy
at the distance z1 = zf(1) from the floor. The w points at k = 2 are listed as fluid
boundary points so the subgrid diffusion across the solid w plane is cancelled
(``modibm:diffw_corr``), as for an STL floor.

Every column belongs to facet 1 unless ``facet_of`` assigns another one. The facet
normals are taken as given, so a test can declare a facet tilted while the grid floor
stays flat (only the normal, the type and the section geometry enter the wall functions).

Usage (prints the &WALLS counts to copy into namoptions)::

    python -m udprep.flat_floor <run_dir> <expnr> --itot 64 --jtot 64 --dx 100 --dy 100 \
        --z1 20 --z0 2e-4 --z0h 2e-4 [--no-lgr] [--tfacinit 302]
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, Sequence

NFACLYRS = 3  # initfac default: factypes has 6 + 4*nfaclyrs + 1 columns

FLOOR_TYPE = -1  # facet type of the floor; negative = ground, not a building (initfac)


def _factypes_row(wallid: int, lgr: bool, z0: float, z0h: float) -> str:
    # wallid lGR z0 z0h albedo emissivity d1..d3 C1..C3 l1..l3 k1..k4 (the layer
    # properties only enter the energy balance, lEB)
    layers = [0.1] * NFACLYRS + [1.875e6] * NFACLYRS + [0.75] * NFACLYRS + [4.0e-7] * (NFACLYRS + 1)
    cols = [f"{wallid:8d}", f"{int(lgr):4d}", f"{z0:.6e}", f"{z0h:.6e}", "0.50", "0.85"]
    cols += [f"{x:g}" for x in layers]
    return " ".join(cols)


def write_floor_factypes(run_dir: Path, expnr: str, z0: float, z0h: float, lgr: bool = True) -> None:
    """Write factypes.inp with the single floor type (also used for STL floors in tests)."""
    with open(Path(run_dir) / f"factypes.inp.{expnr}", "w") as f:
        f.write("# facet types of a floor (udprep.flat_floor)\n"
                "# lGR = 1: the facet exchanges moisture (iwallmoist)\n"
                "# wallid lGR z0[m] z0h[m] al[-] em[-] d1..d3[m] C1..C3[J/(K m^3)] l1..l3[W/(m K)] k1..k4[W/(m K)]\n")
        f.write(_factypes_row(FLOOR_TYPE, lgr, z0, z0h) + "\n")


def write_flat_floor(run_dir: Path, expnr: str, itot: int, jtot: int, dx: float, dy: float,
                     z1: float, z0: float, z0h: float, lgr: bool = True,
                     normals: Sequence[Sequence[float]] | None = None,
                     facet_of: Callable[[int, int], int] | None = None,
                     tfacinit: float | None = None) -> dict:
    """Write the IBM input files of a flat facet floor into run_dir.

    Parameters
    ----------
    z1 : height of the first full level zf(1) [m] (the section distance).
    z0, z0h : roughness lengths of the floor facet type [m].
    lgr : green flag (factypes column 2) of the floor type. In the solver it means
        "this facet exchanges moisture" (modibm:wallfunheat); set it for a moist run.
    normals : facet normals, one per facet (default: one horizontal facet).
    facet_of : maps the 1-based column (i, j) to its 1-based facet (default: 1).
    tfacinit : if given, also writes Tfacinit.inp with this facet temperature [K]
        (needed only by the temperature-dependent wall functions and the energy balance).

    Returns the &WALLS counts {name: value}.
    """
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    normals = [(0.0, 0.0, 1.0)] if normals is None else [tuple(n) for n in normals]
    nfcts = len(normals)
    facet_of = facet_of or (lambda i, j: 1)

    cols = [(i, j) for j in range(1, jtot + 1) for i in range(1, itot + 1)]
    fac = [facet_of(i, j) for i, j in cols]
    if min(fac) < 1 or max(fac) > nfcts:
        raise ValueError(f"facet_of returns facets outside 1..{nfcts}")
    area = dx * dy

    def write(name: str, header: str, rows):
        with open(run_dir / name, "w") as f:
            f.write(header + "\n")
            for r in rows:
                f.write(r + "\n")

    hdr_pts = "# position (i,j,k)"
    hdr_sec = " # facet      area flux point distance"
    pts_k1 = [f"{i:5d} {j:5d} {1:5d}" for i, j in cols]
    pts_k2 = [f"{i:5d} {j:5d} {2:5d}" for i, j in cols]
    secs = [f"{fac[n]:8d} {area:.6e} {n + 1:10d} {z1:.6e}" for n in range(len(cols))]

    for g in ("u", "v", "c"):
        write(f"solid_{g}.txt", hdr_pts, [])
        write(f"fluid_boundary_{g}.txt", hdr_pts, pts_k1)
        write(f"facet_sections_{g}.txt", hdr_sec, secs)
    write("solid_w.txt", hdr_pts, pts_k1)
    write("fluid_boundary_w.txt", hdr_pts, pts_k2)
    write("facet_sections_w.txt", hdr_sec, [])

    write(f"facets.inp.{expnr}", "# type, normal",
          [f"{FLOOR_TYPE:5d} {n[0]:.6f} {n[1]:.6f} {n[2]:.6f}" for n in normals])
    write(f"facetarea.inp.{expnr}", "# area of facets",
          [f"{fac.count(n + 1) * area:.6f}" for n in range(nfcts)])
    write_floor_factypes(run_dir, expnr, z0, z0h, lgr)
    if tfacinit is not None:
        write(f"Tfacinit.inp.{expnr}", "# initial facet temperature [K]",
              [f"{tfacinit:.6f}"] * nfcts)

    n = itot * jtot
    return {"nfcts": nfcts,
            "nsolpts_u": 0, "nsolpts_v": 0, "nsolpts_w": n, "nsolpts_c": 0,
            "nbndpts_u": n, "nbndpts_v": n, "nbndpts_w": n, "nbndpts_c": n,
            "nfctsecs_u": n, "nfctsecs_v": n, "nfctsecs_w": 0, "nfctsecs_c": n}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("run_dir", type=Path)
    p.add_argument("expnr")
    p.add_argument("--itot", type=int, required=True)
    p.add_argument("--jtot", type=int, required=True)
    p.add_argument("--dx", type=float, required=True)
    p.add_argument("--dy", type=float, required=True)
    p.add_argument("--z1", type=float, required=True, help="zf(1), height of the first full level [m]")
    p.add_argument("--z0", type=float, required=True)
    p.add_argument("--z0h", type=float, required=True)
    p.add_argument("--no-lgr", action="store_true", help="floor does not exchange moisture")
    p.add_argument("--tfacinit", type=float, default=None,
                   help="also write Tfacinit.inp with this facet temperature [K] (iwalltemp = 2, lEB)")
    a = p.parse_args(argv)
    counts = write_flat_floor(a.run_dir, a.expnr, a.itot, a.jtot, a.dx, a.dy, a.z1, a.z0, a.z0h,
                              lgr=not a.no_lgr, tfacinit=a.tfacinit)
    for k, v in counts.items():
        print(f"{k:<12} = {v}")


if __name__ == "__main__":
    main()
