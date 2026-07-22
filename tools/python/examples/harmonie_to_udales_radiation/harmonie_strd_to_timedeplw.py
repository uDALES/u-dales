#!/usr/bin/env python3
"""Generate uDALES timedeplw.inp from HARMONIE STRD.

Usage:

    /data/simulateParis/NWP_demo_data/.venv/bin/python \
        tools/python/examples/harmonie_to_udales_radiation/harmonie_strd_to_timedeplw.py \
        --case-dir /home/dipanjan/simulation/udtest/experiments/300 \
        --overwrite

Dry run:

    /data/simulateParis/NWP_demo_data/.venv/bin/python \
        tools/python/examples/harmonie_to_udales_radiation/harmonie_strd_to_timedeplw.py \
        --case-dir /home/dipanjan/simulation/udtest/experiments/300 \
        --dry-run

The script reads HARMONIE ``strd`` accumulated J/m2, averages it over the same
central Paris Lambert-93 box used by the NWP nudging workflow, converts it to
downward sky longwave flux in W/m2, and writes the two-column uDALES
``timedeplw.inp.<expnr>`` format.
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
    generate_timedeplw_from_harmonie,
    longwave_config_from_case_dir,
    prepare_harmonie_strd_longwave,
)


DEFAULT_CASE_DIR = Path.home() / "simulation/udtest/experiments/300"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate timedeplw.inp.<expnr> from HARMONIE strd."
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
        help="Override the timedeplw output path.",
    )
    parser.add_argument("--output-interval-seconds", type=int, default=None)
    parser.add_argument(
        "--difference-interval-seconds",
        type=int,
        default=None,
        help="Accumulation differencing interval; defaults to the native NWP frequency.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the conversion summary without writing the uDALES input file.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing timedeplw file.",
    )
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    case_dir = args.case_dir.expanduser().resolve()
    output = (
        args.output.expanduser().resolve()
        if args.output is not None
        else case_dir / f"timedeplw.inp.{case_dir.name}"
    )

    try:
        if args.dry_run:
            config = longwave_config_from_case_dir(case_dir)
            series = prepare_harmonie_strd_longwave(
                config=config,
                nwp_root=args.nwp_root,
                version=args.version,
                data_dir=args.data_dir,
                download=args.download,
                output_interval_seconds=args.output_interval_seconds,
                difference_interval_seconds=args.difference_interval_seconds,
                verbose=not args.quiet,
            )
            print(
                "Longwave conversion complete; no file written.\n"
                f"  nt       : {series.times.size}\n"
                f"  mask     : {series.mask_points} NWP grid cells\n"
                f"  LWsky    : {series.lwsky.min():.6f} to {series.lwsky.max():.6f} W/m2\n"
                f"  target   : {output}"
            )
            return 0

        result = generate_timedeplw_from_harmonie(
            case_dir=case_dir,
            nwp_root=args.nwp_root,
            version=args.version,
            data_dir=args.data_dir,
            download=args.download,
            output=output,
            overwrite=args.overwrite,
            output_interval_seconds=args.output_interval_seconds,
            difference_interval_seconds=args.difference_interval_seconds,
            verbose=not args.quiet,
        )
    except (RuntimeError, FileExistsError, FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {result.timedeplw_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
