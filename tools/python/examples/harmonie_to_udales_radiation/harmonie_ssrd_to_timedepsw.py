#!/usr/bin/env python3
"""Generate uDALES timedepsw.inp from HARMONIE SSRD.

Usage:

    /data/simulateParis/NWP_demo_data/.venv/bin/python \
        tools/python/examples/harmonie_to_udales_radiation/harmonie_ssrd_to_timedepsw.py \
        --case-dir /home/dipanjan/simulation/udtest/experiments/300 \
        --overwrite

Cheap atmospheric-only check, without loading the STL or tracing facets:

    /data/simulateParis/NWP_demo_data/.venv/bin/python \
        tools/python/examples/harmonie_to_udales_radiation/harmonie_ssrd_to_timedepsw.py \
        --case-dir /home/dipanjan/simulation/udtest/experiments/300 \
        --atmos-only

The script reads HARMONIE ``ssrd`` accumulated J/m2, averages it over the same
central Paris Lambert-93 box used by the NWP nudging workflow, converts it to
GHI in W/m2, splits GHI into DNI and diffuse sky irradiance with an empirical
Erbs clearness-index model, and maps it to uDALES facets using the existing
udprep direct-shortwave and net-shortwave routines.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


UDALES_ROOT = Path(__file__).resolve().parents[4]
PYTOOLS = UDALES_ROOT / "tools" / "python"
if str(PYTOOLS) not in sys.path:
    sys.path.insert(0, str(PYTOOLS))

from udprep.harmonie_radiation import (  # noqa: E402
    DEFAULT_NWP_ROOT,
    DEFAULT_VERSION,
    generate_timedepsw_from_harmonie,
)


DEFAULT_CASE_DIR = Path.home() / "simulation/udtest/experiments/300"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate timedepsw.inp.<expnr> from HARMONIE ssrd."
    )
    parser.add_argument("--case-dir", type=Path, default=DEFAULT_CASE_DIR)
    parser.add_argument("--nwp-root", type=Path, default=DEFAULT_NWP_ROOT)
    parser.add_argument("--version", default=DEFAULT_VERSION)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing the downloaded URL tree, or the day directory "
            "containing GRIBPFDEOD+ files."
        ),
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Run the NWP demo downloader if the local data directory is missing.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Override the timedepsw output path.",
    )
    parser.add_argument(
        "--sdir-nc",
        type=Path,
        default=None,
        help="Override the Sdir.nc output path.",
    )
    parser.add_argument(
        "--no-sdir-nc",
        action="store_true",
        help="Do not write Sdir.nc alongside timedepsw.",
    )
    parser.add_argument(
        "--method",
        default=None,
        help="Override direct-shortwave method, e.g. scanline_f2py, facsec, moller.",
    )
    parser.add_argument(
        "--atmos-only",
        action="store_true",
        help="Only compute/print GHI, DNI, and Dsky; do not load facets or write files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing timedepsw/Sdir outputs.",
    )
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        result = generate_timedepsw_from_harmonie(
            case_dir=args.case_dir,
            nwp_root=args.nwp_root,
            version=args.version,
            data_dir=args.data_dir,
            download=args.download,
            output=args.output,
            sdir_nc=args.sdir_nc,
            write_sdir_nc=not args.no_sdir_nc,
            overwrite=args.overwrite,
            atmos_only=args.atmos_only,
            method=args.method,
            verbose=not args.quiet,
        )
    except (RuntimeError, FileExistsError, FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if args.atmos_only:
        print(
            "Atmospheric conversion complete; no facet files written.\n"
            f"  nt     : {result.times.size}\n"
            f"  GHI    : {result.ghi.min():.3f} to {result.ghi.max():.3f} W/m2\n"
            f"  DNI    : {result.dni.min():.3f} to {result.dni.max():.3f} W/m2\n"
            f"  Dsky   : {result.dsky.min():.3f} to {result.dsky.max():.3f} W/m2"
        )
        return 0

    print(f"Wrote {result.timedepsw_path}")
    if result.sdir_nc_path is not None:
        print(f"Wrote {result.sdir_nc_path}")
    if result.timedepsveg_path is not None:
        print(f"Wrote {result.timedepsveg_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
