#!/usr/bin/env python3
"""Build the parent ("big brother") case directory for the V1 validation.

The parent is a periodic, neutral urban LES over an aligned array of cubes,
driven by a fixed mean pressure gradient.  It is run in two phases out of one
directory:

``namoptions_spinup.<nr>``
    cold start, no field dumps, one restart file written at ``t = spinup``.
``namoptions.<nr>``
    warm start from that restart, field dumps of ``u, v, w`` every ``dtdump``
    seconds over the production window.  This is the canonical namelist -- it
    is the one ``UDPrep`` reads, and the one whose ``iexpnr`` names the outputs.

The two-phase split exists so the dumps cover the production window only: a
single run dumping from ``t = 0`` would triple the output volume and the I/O
time for data that is thrown away.

Usage
-----
    python make_parent_case.py <outdir> [--preset tiny|production]

``<outdir>`` gets a subdirectory named after the experiment number, because
``UDPrep`` requires the directory name, the ``namoptions`` suffix and ``iexpnr``
to agree.
"""

from __future__ import annotations

import argparse
import json
from collections import OrderedDict
from pathlib import Path

import numpy as np

from caselib import (
    cube_geometry,
    run_preprocessing,
    write_lscale,
    write_namoptions,
    write_prof,
)
from config import Preset, get_preset


def parent_sections(preset: Preset, *, warmstart: bool, startfile: str) -> "OrderedDict":
    """Namelist content for one phase of the parent run."""
    nr = preset.parent_expnr
    production = warmstart
    return OrderedDict([
        ("RUN", OrderedDict([
            ("iexpnr", int(nr)),
            ("runmode", 1),
            ("runtime", preset.production if production else preset.spinup),
            ("dtmax", preset.dtmax),
            ("ladaptive", True),
            ("trestart", 1.0e9 if production else preset.spinup),
            ("lwarmstart", warmstart),
            ("startfile", startfile),
            ("irandom", 43),
            ("randu", 0.1 * preset.u0),
            ("krand", preset.ktot),
            ("lrandomize", not warmstart),
            ("libm", True),
            ("lles", True),
            ("nprocx", preset.nprocx),
            ("nprocy", preset.nprocy),
        ])),
        ("DOMAIN", OrderedDict([
            ("itot", preset.itot),
            ("jtot", preset.jtot),
            ("ktot", preset.ktot),
            ("xlen", preset.xlen),
            ("ylen", preset.ylen),
        ])),
        ("PHYSICS", OrderedDict([
            ("ltempeq", False),
            ("lmoist", False),
            ("lbuoyancy", False),
            ("lcoriol", False),
            ("lprofforc", False),
            ("luvolflowr", False),
            ("lvvolflowr", False),
            ("igrw_damp", 0),
            # The whole momentum source.  lscale.inp carries pgx = 0 so that
            # dpdxl = -pgx - dpdx is exactly -dpdx (src/modstartup.f90:2236).
            ("dpdx", preset.dpdx),
        ])),
        ("DYNAMICS", OrderedDict([
            ("ipoiss", 0),
        ])),
        ("BC", OrderedDict([
            ("BCxm", 1),   # periodic
            ("BCym", 1),   # periodic
            ("BCbotm", 3),  # neutral wall function at the floor
            ("BCtopm", 1),  # free slip lid
        ])),
        ("NAMSUBGRID", OrderedDict([
            ("lvreman", True),
        ])),
        ("WALLS", OrderedDict([
            ("iwallmom", 2),
            ("iwalltemp", 1),
            ("lwritefac", False),
        ])),
        ("NAMCHECKSIM", OrderedDict([
            ("tcheck", max(1.0, preset.production / 20.0)),
        ])),
        ("OUTPUT", OrderedDict([
            ("lfielddump", production),
            ("tfielddump", preset.dtdump if production else 1.0e9),
            ("fieldvars", "u0,v0,w0"),
        ])),
        ("INP", OrderedDict([
            ("zsize", preset.zsize),
            ("stl_file", f"geom.{nr}.stl"),
            ("u0", preset.u0),
            ("v0", 0.0),
            ("thl0", 288.0),
        ])),
    ])


def build(outdir: Path, preset: Preset, ibm_backend: str = "auto") -> Path:
    """Create and preprocess the parent case; return its directory."""
    nr = preset.parent_expnr
    casedir = Path(outdir) / nr
    casedir.mkdir(parents=True, exist_ok=True)

    zf = (np.arange(preset.ktot) + 0.5) * preset.dz
    write_prof(casedir / f"prof.inp.{nr}", zf,
               u=preset.u0, v=0.0, e12=preset.tke0,
               comment=f"V1 big-brother parent ({preset.name}): uniform u0")
    # pgx = 0 on purpose; the forcing lives in &PHYSICS dpdx.  See the
    # docstring of caselib.write_lscale before changing this.
    write_lscale(casedir / f"lscale.inp.{nr}", zf,
                 comment="V1 big-brother parent: forcing is dpdx in &PHYSICS")

    cube_geometry(preset, 0.0, 0.0, preset.xlen, preset.ylen,
                  casedir / f"geom.{nr}.stl")

    # The canonical namelist is the production one: UDPrep reads namoptions.<nr>
    # and writes the &WALLS point counts back into it.  The spin-up variant is a
    # copy with the phase keys changed, written afterwards from the same dict.
    write_namoptions(
        casedir / f"namoptions.{nr}",
        parent_sections(preset, warmstart=True, startfile=f"initd00000000_000_000.{nr}"),
        header=[
            f"V1 Big Brother parent, preset '{preset.name}' -- PRODUCTION phase",
            "generated by tests/validation/nesting/make_parent_case.py; do not hand-edit",
            "run namoptions_spinup.%s first, then patch startfile here" % nr,
        ],
    )
    write_namoptions(
        casedir / f"namoptions_spinup.{nr}",
        parent_sections(preset, warmstart=False, startfile=f"initd00000000_000_000.{nr}"),
        header=[
            f"V1 Big Brother parent, preset '{preset.name}' -- SPIN-UP phase",
            "generated by tests/validation/nesting/make_parent_case.py; do not hand-edit",
        ],
    )

    run_preprocessing(casedir, ibm_backend=ibm_backend)

    # UDPrep only touches namoptions.<nr>; mirror the &WALLS counts it added into
    # the spin-up namelist so the two phases see identical geometry input.
    _mirror_walls(casedir / f"namoptions.{nr}", casedir / f"namoptions_spinup.{nr}")

    (casedir / "preset.json").write_text(
        json.dumps({"preset": preset.name, "role": "parent",
                    "expnr": nr, "dpdx": preset.dpdx,
                    "ustar": preset.ustar,
                    "t_start": preset.t_start, "t_end": preset.t_end},
                   indent=2) + "\n",
        encoding="ascii",
    )
    return casedir


def _mirror_walls(src: Path, dst: Path) -> None:
    """Copy the &WALLS block of ``src`` over the one in ``dst``."""
    def walls_block(text: str) -> str:
        lines = text.splitlines()
        start = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("&WALLS"))
        end = next(i for i in range(start, len(lines)) if lines[i].strip() == "/")
        return "\n".join(lines[start:end + 1])

    src_text = src.read_text(encoding="ascii")
    dst_text = dst.read_text(encoding="ascii")
    old = walls_block(dst_text)
    dst.write_text(dst_text.replace(old, walls_block(src_text)), encoding="ascii")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--preset", default="production")
    parser.add_argument("--ibm-backend", default="auto")
    args = parser.parse_args()

    preset = get_preset(args.preset)
    print(preset.summary())
    casedir = build(args.outdir, preset, ibm_backend=args.ibm_backend)
    print(f"\nparent case written to {casedir}")


if __name__ == "__main__":
    main()
