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

A parent can instead be **warm-started from another run's restart files**
(``--restart-dir``): the C0b fine-cadence parent (``config.C0_FINE``) picks up
the converged parent's end-of-spin-up state and only runs the production phase,
dumping at 0.5 s.  The restart files are symlinked into the new case under the
new experiment number -- ``readrestartfiles`` builds each rank's file name from
``startfile`` by overwriting the rank fields (``modstartup.f90``, ``name(15:17)
= cmyidx``), so the extension has to be whatever ``startfile`` says and the
namelist is self-consistent when it is the new ``iexpnr`` -- and no spin-up
namelist is written.  The new case must share the source's grid, rank layout and
geometry, which is what makes the restart's IBM arrays (``mindist``, ``wall``)
valid for it; ``restart_source`` checks the first two and the preset guarantees
the third.

Usage
-----
    python make_parent_case.py <outdir> [--preset tiny|production]
    python make_parent_case.py <outdir> --preset c0-fine \\
        --restart-dir $EPHEMERAL/nesting-v1-converged/903

``<outdir>`` gets a subdirectory named after the experiment number, because
``UDPrep`` requires the directory name, the ``namoptions`` suffix and ``iexpnr``
to agree.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from caselib import (
    cube_geometry,
    run_preprocessing,
    write_lscale,
    write_namoptions,
    write_prof,
)
from config import Preset, get_preset


def parent_sections(preset: Preset, *, warmstart: bool, startfile: str,
                    nestparent_child: Optional[Preset] = None) -> "OrderedDict":
    """Namelist content for one phase of the parent run.

    What the production phase writes is ``preset.parent_output``: full field
    dumps (``&OUTPUT``), the child's band only (``&NESTPARENT``,
    ``src/nesting_parent.f90``), or both.  ``nestparent_child`` sizes the band for a
    child other than the preset's own -- a V0 driver dumping for the fine
    child it will drive -- see :meth:`config.Preset.nestparent_sections`.
    """
    nr = preset.parent_expnr
    production = warmstart
    nestparent = OrderedDict([("lnestparent", False)])
    if production and preset.writes_nestparent:
        nestparent = preset.nestparent_sections(nestparent_child)
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
            ("lfielddump", production and preset.writes_fielddump),
            ("tfielddump", preset.fielddump_interval if production else 1.0e9),
            ("fieldvars", "u0,v0,w0"),
        ])),
        ("NESTPARENT", nestparent),
        ("INP", OrderedDict([
            ("zsize", preset.zsize),
            ("stl_file", f"geom.{nr}.stl"),
            ("u0", preset.u0),
            ("v0", 0.0),
            ("thl0", 288.0),
        ])),
    ])


#: ``initd<ntrun:08d>_<x:03d>_<y:03d>.<expnr>`` -- modsave's restart file name.
_RESTART = re.compile(r"^initd(\d{8})_(\d{3})_(\d{3})\.(\d+)$")


def restart_source(restart_dir: Path, preset: Preset,
                   source_expnr: Optional[str] = None) -> Tuple[str, List[Path]]:
    """The complete restart set in ``restart_dir`` that this preset can start from.

    Returns ``(ntrun, files)``: the step counter the set was written at and
    one file per rank, sorted.  Picks the latest set (largest ``ntrun``) of the
    given ``source_expnr`` -- or of the only expnr present -- and refuses one
    that does not have exactly ``nprocx * nprocy`` ranks laid out as
    ``000..nprocx-1`` x ``000..nprocy-1``, because a restart is a per-rank
    dump of the decomposed field and cannot be re-decomposed here.
    """
    restart_dir = Path(restart_dir)
    found: Dict[Tuple[str, str], Dict[Tuple[int, int], Path]] = {}
    for p in restart_dir.iterdir():
        m = _RESTART.match(p.name)
        if not m:
            continue
        ntrun, ix, iy, nr = m.groups()
        if source_expnr is not None and nr != source_expnr:
            continue
        found.setdefault((nr, ntrun), {})[(int(ix), int(iy))] = p
    if not found:
        raise FileNotFoundError(
            f"no initd????????_???_???.{source_expnr or '*'} restart files in {restart_dir}")
    expnrs = {nr for nr, _ in found}
    if len(expnrs) > 1:
        raise ValueError(f"{restart_dir} holds restart files of several experiments "
                         f"({', '.join(sorted(expnrs))}); pass source_expnr")
    nr, ntrun = max(found, key=lambda k: int(k[1]))
    ranks = found[(nr, ntrun)]
    want = {(i, j) for i in range(preset.nprocx) for j in range(preset.nprocy)}
    if set(ranks) != want:
        raise ValueError(
            f"restart set initd{ntrun}_*.{nr} in {restart_dir} covers ranks "
            f"{sorted(ranks)[:3]}... ({len(ranks)} files), but preset "
            f"'{preset.name}' runs {preset.nprocx} x {preset.nprocy} ranks; a "
            "restart cannot be re-decomposed")
    return ntrun, [ranks[k] for k in sorted(ranks)]


def link_restart(casedir: Path, preset: Preset, restart_dir: Path,
                 source_expnr: Optional[str] = None) -> str:
    """Symlink a restart set into ``casedir`` under this preset's expnr.

    Returns the ``startfile`` name (rank 0,0) to put in the namelist.  Symlinks
    rather than copies: the set is 64 x 7.8 MB for the converged parent and is
    read once, at start-up.  An existing link or file of the same name is
    replaced, so re-building a case is idempotent.
    """
    nr = preset.parent_expnr
    ntrun, files = restart_source(restart_dir, preset, source_expnr)
    for src in files:
        m = _RESTART.match(src.name)
        assert m is not None
        dst = casedir / f"initd{m.group(1)}_{m.group(2)}_{m.group(3)}.{nr}"
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(src.resolve())
    return f"initd{ntrun}_000_000.{nr}"


def build(outdir: Path, preset: Preset, ibm_backend: str = "auto",
          restart_dir: Optional[Path] = None,
          restart_expnr: Optional[str] = None,
          nestparent_child: Optional[Preset] = None) -> Path:
    """Create and preprocess the parent case; return its directory.

    ``nestparent_child`` is passed on to :func:`parent_sections`: the child the
    ``&NESTPARENT`` band is sized for when it is not the preset's own.

    With ``restart_dir`` the case is warm-started from the restart files found
    there (see :func:`link_restart`): only the production namelist is written,
    with ``lwarmstart = .true.`` and ``startfile`` already set, and there is no
    spin-up phase to run.  The clock continues from the restart's ``timee``,
    so the preset's ``spinup`` must be the time the source's restart was
    written at for ``t_start``/``t_end`` to mean what they say.
    """
    nr = preset.parent_expnr
    casedir = Path(outdir) / nr
    casedir.mkdir(parents=True, exist_ok=True)
    startfile = f"initd00000000_000_000.{nr}"
    if restart_dir is not None:
        startfile = link_restart(casedir, preset, Path(restart_dir), restart_expnr)

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
        parent_sections(preset, warmstart=True, startfile=startfile,
                        nestparent_child=nestparent_child),
        header=[
            f"V1 Big Brother parent, preset '{preset.name}' -- PRODUCTION phase",
            f"production output: {preset.parent_output}"
            + (f" (OUTPUT every {preset.fielddump_interval:g} s, NESTPARENT every "
               f"{preset.dtdump:g} s)" if preset.parent_output == "both" else ""),
            "generated by tests/validation/nesting/make_parent_case.py; do not hand-edit",
        ] + ([f"warm-started from the restart files of {Path(restart_dir).resolve()}"]
             if restart_dir is not None else
             ["run namoptions_spinup.%s first, then patch startfile here" % nr]),
    )
    if restart_dir is None:
        write_namoptions(
            casedir / f"namoptions_spinup.{nr}",
            parent_sections(preset, warmstart=False, startfile=startfile),
            header=[
                f"V1 Big Brother parent, preset '{preset.name}' -- SPIN-UP phase",
                "generated by tests/validation/nesting/make_parent_case.py; do not hand-edit",
            ],
        )

    run_preprocessing(casedir, ibm_backend=ibm_backend)

    # UDPrep only touches namoptions.<nr>; mirror the &WALLS counts it added into
    # the spin-up namelist so the two phases see identical geometry input.
    if restart_dir is None:
        _mirror_walls(casedir / f"namoptions.{nr}", casedir / f"namoptions_spinup.{nr}")

    (casedir / "preset.json").write_text(
        json.dumps({"preset": preset.name, "role": "parent",
                    "expnr": nr, "dpdx": preset.dpdx,
                    "ustar": preset.ustar,
                    "t_start": preset.t_start, "t_end": preset.t_end,
                    "dtdump": preset.dtdump,
                    "fielddump_interval": preset.fielddump_interval,
                    "parent_output": preset.parent_output,
                    "warmstart_from": (None if restart_dir is None
                                       else str(Path(restart_dir).resolve())),
                    "startfile": startfile if restart_dir is not None else None},
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
    parser.add_argument("--restart-dir", type=Path, default=None,
                        help="warm-start from the restart files in this directory "
                             "(no spin-up phase is written)")
    parser.add_argument("--restart-expnr", default=None,
                        help="which experiment's restart set to take from "
                             "--restart-dir when it holds several")
    args = parser.parse_args()

    preset = get_preset(args.preset)
    print(preset.summary())
    casedir = build(args.outdir, preset, ibm_backend=args.ibm_backend,
                    restart_dir=args.restart_dir, restart_expnr=args.restart_expnr)
    print(f"\nparent case written to {casedir}")
    if args.restart_dir is not None:
        print("  warm start: run namoptions.%s directly, there is no spin-up phase"
              % preset.parent_expnr)


if __name__ == "__main__":
    main()
