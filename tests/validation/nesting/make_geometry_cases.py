#!/usr/bin/env python3
"""Case builders for the V3 and V4 geometry-mismatch experiments.

Only the *periodic* runs need anything new.  V1's ``make_parent_case`` builds a
periodic urban LES on one aligned cube array driven by a fixed ``dpdx``; V3 and
V4 additionally need

* a **flat** parent -- no buildings at all (V3);
* a **staggered** cube array (V4);
* a **volume-flow-rate** forced parent, so that a flat parent can be run at the
  canopy's own bulk velocity rather than at its own much faster one (V3; see the
  ``V3_PARENT`` docstring in ``presets_geometry`` for why that matters).

The namelist itself is still ``make_parent_case.parent_sections`` -- there is
one description of what a periodic parent run is, and this module patches three
keys of it rather than restating it.  The nested children need nothing new at
all: ``make_child_case.build`` already takes the child's layout from
``Preset.child_cube_centres`` and its zone rule from
``Preset.building_free_zone``, both of which ``GeoPreset`` overrides.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import numpy as np

import caselib
import make_parent_case
from caselib import (
    cube_geometry,
    run_preprocessing,
    write_lscale,
    write_namoptions,
    write_prof,
)
from presets_geometry import GeoPreset, get_experiment


def periodic_geometry(preset: GeoPreset, out_stl: Path) -> Path:
    """Write the STL for a periodic run's whole domain.

    ``"none"``       a flat ground plane -- ``udgeom.create_flat_surface``.
    ``"aligned"``    the unbroken aligned array, through the same
                     ``caselib.cube_mesh`` path V1 uses, so parent and child
                     tessellate identically.
    ``"staggered"``  ``udgeom.create_cubes(..., 'SC')`` directly, because the
                     staggered array's displaced columns put a half cube on each
                     spanwise periodic face and ``create_cubes`` is what owns
                     the clipping that produces them.  ``caselib.cube_mesh``
                     places whole cubes only and would put geometry outside the
                     domain.
    """
    from udgeom.geometry_generation import create_cubes, create_flat_surface

    out_stl.parent.mkdir(parents=True, exist_ok=True)
    if preset.parent_layout == "none":
        geom = create_flat_surface(preset.xlen, preset.ylen, preset.edgelength)
        geom.save(str(out_stl))
        return out_stl
    if preset.parent_layout == "staggered":
        geom = create_cubes(preset.xlen, preset.ylen,
                            preset.building_width, preset.building_width,
                            preset.building_height,
                            preset.street_width, preset.street_width,
                            "SC", preset.edgelength)
        geom.save(str(out_stl))
        return out_stl
    return cube_geometry(preset, 0.0, 0.0, preset.xlen, preset.ylen, out_stl)


def _apply_forcing(sections, preset: GeoPreset, uflowrate: Optional[float]):
    """Swap ``dpdx`` for a volume-flow-rate controller when one is asked for.

    ``lscale.inp`` carries ``pgx = 0`` in either case (see
    ``caselib.write_lscale``), so ``dpdx = 0`` really does mean no mean pressure
    gradient forcing and the controller is the only momentum source.
    """
    if uflowrate is None:
        return sections
    phys = sections["PHYSICS"]
    phys["dpdx"] = 0.0
    phys["luvolflowr"] = True
    phys["uflowrate"] = float(uflowrate)
    return sections


def build_periodic(outdir: Path, preset: GeoPreset,
                   uflowrate: Optional[float] = None,
                   ibm_backend: str = "auto") -> Path:
    """Create and preprocess a periodic run -- a reference or a parent.

    Two namelists as in V1: ``namoptions_spinup.<nr>`` (cold start, one restart
    file, no dumps) and ``namoptions.<nr>`` (warm start, dumps over the
    production window).  ``uflowrate`` overrides ``preset.uflowrate``; that is
    how V3's flat parent gets the bulk velocity its reference run measured.
    """
    nr = preset.parent_expnr
    casedir = Path(outdir) / nr
    casedir.mkdir(parents=True, exist_ok=True)
    rate = preset.uflowrate if uflowrate is None else float(uflowrate)

    zf = (np.arange(preset.ktot) + 0.5) * preset.dz
    u_init = preset.u0 if rate is None else rate
    write_prof(casedir / f"prof.inp.{nr}", zf, u=u_init, v=0.0, e12=preset.tke0,
               comment=f"{preset.name} ({preset.role}): uniform u = {u_init:g} m/s")
    # pgx = 0 on purpose; see caselib.write_lscale before changing this.
    write_lscale(casedir / f"lscale.inp.{nr}", zf,
                 comment=f"{preset.name}: forcing is dpdx or uflowrate in &PHYSICS")

    periodic_geometry(preset, casedir / f"geom.{nr}.stl")

    startfile = f"initd00000000_000_000.{nr}"
    write_namoptions(
        casedir / f"namoptions.{nr}",
        _apply_forcing(
            make_parent_case.parent_sections(preset, warmstart=True,
                                             startfile=startfile),
            preset, rate),
        header=[
            f"{preset.name} ({preset.role}), layout '{preset.parent_layout}' "
            "-- PRODUCTION phase",
            "generated by tests/validation/nesting/make_geometry_cases.py; "
            "do not hand-edit",
            "run namoptions_spinup.%s first, then patch startfile here" % nr,
        ],
    )
    write_namoptions(
        casedir / f"namoptions_spinup.{nr}",
        _apply_forcing(
            make_parent_case.parent_sections(preset, warmstart=False,
                                             startfile=startfile),
            preset, rate),
        header=[
            f"{preset.name} ({preset.role}), layout '{preset.parent_layout}' "
            "-- SPIN-UP phase",
            "generated by tests/validation/nesting/make_geometry_cases.py; "
            "do not hand-edit",
        ],
    )

    run_preprocessing(casedir, ibm_backend=ibm_backend)
    make_parent_case._mirror_walls(casedir / f"namoptions.{nr}",
                                   casedir / f"namoptions_spinup.{nr}")

    (casedir / "preset.json").write_text(
        json.dumps({
            "preset": preset.name,
            "role": preset.role,
            "expnr": nr,
            "parent_layout": preset.parent_layout,
            "n_cubes": int(len(preset.cube_centres())),
            "forcing": ("uflowrate" if rate is not None else "dpdx"),
            "uflowrate": None if rate is None else float(rate),
            "dpdx": 0.0 if rate is not None else preset.dpdx,
            "ustar_nominal": preset.ustar,
            "t_start": preset.t_start,
            "t_end": preset.t_end,
            "dtdump": preset.dtdump,
            "fielddump_interval": preset.fielddump_interval,
            "parent_output": preset.parent_output,
        }, indent=2) + "\n", encoding="ascii")
    return casedir


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--experiment", default="v3")
    ap.add_argument("--run", default="parent", help="which periodic run to build")
    ap.add_argument("--uflowrate", type=float, default=None)
    ap.add_argument("--ibm-backend", default="auto")
    ns = ap.parse_args()

    exp = get_experiment(ns.experiment)
    run = exp.periodic_run(ns.run)
    print(run.preset.summary())
    casedir = build_periodic(ns.outdir, run.preset, uflowrate=ns.uflowrate,
                             ibm_backend=ns.ibm_backend)
    print(f"\ncase written to {casedir}")


if __name__ == "__main__":
    main()
