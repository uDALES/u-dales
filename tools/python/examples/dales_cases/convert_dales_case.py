#!/usr/bin/env python3
"""Convert the profile inputs of a DALES case to uDALES input files.

Usage: python convert_dales_case.py <dales_case_dir> <expnr> <out_dir>

Reads  prof.inp.001   (z thl qt u v tke)
       lscale.inp.001 (z ug vg wfls dqtdxls dqtdyls dqtdtls dthlrad)
       scalar.inp.001 (z sv1 ... svn, optional)
       namoptions.001 (only printed, not converted)
and writes prof.inp.<expnr>, lscale.inp.<expnr> (10 columns:
z ug vg pgx pgy wfls dqtdxls dqtdyls dqtdtls thlpcar, with pgx = pgy = 0 and
thlpcar = dthlrad) and scalar.inp.<expnr>. Every output file gets exactly two
header lines, as the uDALES reader skips two lines.

DALES files have one or two header lines; header lines are detected as lines
that do not parse as floats. The data columns are copied unchanged (no unit
change: z [m], thl [K], qt [kg/kg], u, v [m/s], tke column = sqrt(e) [m/s] in
both models, wfls [m/s], dqtdxls etc. [kg/kg/s], dthlrad = thlpcar [K/s]).

Namelist values of the DALES case are printed so that they can be transcribed
into namoptions.<expnr>. Surface fluxes are printed already negated:
DALES wtsurf/wqsurf are positive upward, uDALES prescribed fluxes are positive
in -z (decision D1); in uDALES they are the facet fluxes bctfz/bcqfz of a facet
floor (tools/python/udprep/flat_floor.py). isurf 4 -> iwallmom = 4 (u* from
Monin-Obukhov), isurf 3 -> iwallmom = 4 with ustin in &WALLS (prescribed u*).

After writing, the output is re-read and compared with the input (row 1 and
all rows, column for column); the script exits non-zero on a mismatch.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

ISURF_TO_UDALES = {3: "iwallmom = 4 with ustin (&WALLS)", 4: "iwallmom = 4"}
# Namelist entries to report (DALES name, uDALES name or None)
REPORT = ["ps", "thls", "wtsurf", "wqsurf", "z0", "z0mav", "z0hav", "ustin", "isurf",
          "itot", "jtot", "kmax", "xsize", "ysize", "runtime", "lcoriol", "lmoist",
          "xlat", "cu", "cv", "iadv_mom", "iadv_thl", "iadv_qt", "iadv_sv", "nsv",
          "ladaptive", "courant", "igrw_damp", "iradiation", "rad_longw", "rad_smoke",
          "isvsmoke", "rka", "dlwtop", "dlwbot", "timerad", "randthl", "randqt", "randu"]


def read_table(path: Path, ncol: int | None = None) -> np.ndarray:
    """Return the numeric rows of a DALES profile file (header lines skipped)."""
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        try:
            vals = [float(x.replace("d", "e").replace("D", "E")) for x in parts]
        except ValueError:
            continue  # header line
        rows.append(vals)
    if not rows:
        raise ValueError(f"{path}: no numeric rows")
    n = {len(r) for r in rows}
    if len(n) != 1:
        raise ValueError(f"{path}: rows have different column counts {sorted(n)}")
    arr = np.array(rows, dtype=float)
    if ncol is not None and arr.shape[1] != ncol:
        raise ValueError(f"{path}: expected {ncol} columns, found {arr.shape[1]}")
    return arr


def write_table(path: Path, header: str, colnames: str, arr: np.ndarray) -> None:
    with open(path, "w") as f:
        f.write(f"# {header}\n# {colnames}\n")
        for row in arr:
            f.write(" ".join(f"{x:.9e}" for x in row) + "\n")


def read_namelist(path: Path) -> dict[str, str]:
    """Flat key -> raw value dict of a Fortran namelist file (groups ignored)."""
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.split("!")[0].strip()
        m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+?)\s*,?$", line)
        if m:
            out[m.group(1).lower()] = m.group(2)
    return out


def fnum(s: str) -> float:
    return float(s.replace("d", "e").replace("D", "E"))


def convert(case: Path, expnr: str, out: Path, dales_expnr: str = "001") -> dict:
    out.mkdir(parents=True, exist_ok=True)
    prof = read_table(case / f"prof.inp.{dales_expnr}", 6)
    write_table(out / f"prof.inp.{expnr}", f"converted from DALES {case.name} prof.inp.{dales_expnr}",
                "z thl qt u v e12", prof)

    ls = read_table(case / f"lscale.inp.{dales_expnr}", 8)
    if len(ls) > len(prof) and np.allclose(ls[:len(prof), 0], prof[:, 0]):
        # DALES lscale.inp may carry extra rows above the model top (BOMEX has one,
        # at the k1 level); uDALES reads exactly ktot rows.
        print(f"note: dropping {len(ls) - len(prof)} lscale.inp row(s) above the model top "
              f"(z = {ls[len(prof):, 0].tolist()})")
        ls = ls[:len(prof)]
    zeros = np.zeros(len(ls))
    ls10 = np.column_stack([ls[:, 0], ls[:, 1], ls[:, 2], zeros, zeros,
                            ls[:, 3], ls[:, 4], ls[:, 5], ls[:, 6], ls[:, 7]])
    write_table(out / f"lscale.inp.{expnr}", f"converted from DALES {case.name} lscale.inp.{dales_expnr}",
                "z ug vg pgx pgy wfls dqtdxls dqtdyls dqtdtls thlpcar", ls10)

    result = {"prof": (prof, out / f"prof.inp.{expnr}"), "lscale": (ls10, out / f"lscale.inp.{expnr}")}
    if len(prof) != len(ls):
        raise ValueError("prof.inp and lscale.inp have different numbers of levels")
    if not np.allclose(prof[:, 0], ls[:, 0]):
        raise ValueError("prof.inp and lscale.inp heights differ")

    sc_in = case / f"scalar.inp.{dales_expnr}"
    if sc_in.exists():
        sc = read_table(sc_in)
        # uDALES reader: z, then (svprof(k,n), n = 1..nsv), two header lines
        if len(sc) != len(prof) or not np.allclose(sc[:, 0], prof[:, 0]):
            raise ValueError("scalar.inp heights differ from prof.inp")
        write_table(out / f"scalar.inp.{expnr}", f"converted from DALES {case.name} scalar.inp.{dales_expnr}",
                    "z " + " ".join(f"sv{n + 1}" for n in range(sc.shape[1] - 1)), sc)
        result["scalar"] = (sc, out / f"scalar.inp.{expnr}")
    return result


def unit_check(result: dict) -> bool:
    """Re-read every output and compare it with the input: row 1 and all rows."""
    ok = True
    for name, (src, path) in result.items():
        back = read_table(path)
        exp = src
        same = back.shape == exp.shape and np.allclose(back, exp, rtol=1e-8, atol=0)
        print(f"unit check {path.name}: row 1 in  = {np.array2string(exp[0], precision=6)}")
        print(f"unit check {path.name}: row 1 out = {np.array2string(back[0], precision=6)}"
              f"  all {len(exp)} rows match: {same}")
        ok &= same
    return ok


def report_namoptions(case: Path, dales_expnr: str = "001") -> None:
    nml = read_namelist(case / f"namoptions.{dales_expnr}")
    print("\nDALES namoptions values to transcribe (uDALES names/signs in brackets):")
    for key in REPORT:
        if key in nml:
            print(f"  {key:<11} = {nml[key]}")
    if "isurf" in nml:
        isurf = int(fnum(nml["isurf"]))
        if isurf in ISURF_TO_UDALES:
            print(f"  -> {ISURF_TO_UDALES[isurf]} on a facet floor  (DALES isurf = {isurf})")
        else:
            print(f"  -> DALES isurf = {isurf} has no uDALES MO equivalent")
    for key, ukey in (("wtsurf", "bctfz"), ("wqsurf", "bcqfz")):
        if key in nml:
            print(f"  -> uDALES {ukey} = {-fnum(nml[key]):+.6g}  (DALES {key} = {nml[key]} is upward; "
                  "uDALES positive = downward, decision D1)")
    if "xlat" in nml:
        f = 2 * 7.292e-5 * np.sin(np.radians(fnum(nml["xlat"])))
        print(f"  -> Coriolis parameter f = {f:.4e} 1/s at xlat = {nml['xlat']}")
    if all(k in nml for k in ("xsize", "itot")):
        print(f"  -> dx = {fnum(nml['xsize']) / fnum(nml['itot']):g} m")
    if "cu" in nml or "cv" in nml:
        print("  note: DALES cu/cv (Galilean shift) have no uDALES counterpart; "
              "prof.inp u,v and lscale ug,vg are physical winds and are copied unchanged")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("dales_case_dir", type=Path)
    ap.add_argument("expnr")
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--dales-expnr", default="001")
    a = ap.parse_args(argv)
    res = convert(a.dales_case_dir, a.expnr, a.out_dir, a.dales_expnr)
    ok = unit_check(res)
    report_namoptions(a.dales_case_dir, a.dales_expnr)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
