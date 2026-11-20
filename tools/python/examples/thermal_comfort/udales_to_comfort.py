#!/usr/bin/env python3
"""Prepare thermal-comfort exchange files from uDALES outputs.

Usage
-----
    python udales_to_comfort.py --case-dir CASE_DIR --heights H [H ...] [options]

    --case-dir DIR       The uDALES case directory (required). Its name is the
                         experiment number. It must hold namoptions.<expnr>, the
                         preprocessing files, the gathered facet energy-balance
                         output facEB.<expnr>.nc, and the gathered time-averaged
                         statistics. The site longitude, latitude
                         and elevation, the model axis rotation (xazimuth), the
                         UTC start date, the saved receptor heights and the
                         statistics schedule (tstatstart, tstatsdump, tstatsgap)
                         are all read from namoptions. Every output is written
                         into this directory.
    --heights H [H ...]  Receptor heights in metres to process (required). Each
                         must be one of the receptor heights saved by the
                         simulation (receptor_heights in namoptions).

Options
-------
    --statistics-source {stats_t,stats_tkslice}
                         Gathered statistics file to read. Default: stats_t.
                         Without --data-source the run reads this file of the
                         case. With --data-source this file is read from every
                         segment (it must exist in each; no other file is
                         tried) and merged.
    --data-source DIR    Case directory of one warm-start segment. Repeat for
                         every segment, including --case-dir itself if its own
                         records belong to the chain. From each directory the
                         files facEB.<case>.nc and <statistics-source>.<case>.nc
                         are read following the uDALES naming rule. Only what
                         the indices need is merged:
                         facet times t and emitted longwave LWout into
                         merged_facEB.<expnr>.nc, and tha, pabs, qt, ql at the
                         bracketing scalar levels plus ws_local and ws_10 into
                         merged_stats_tkslice.<expnr>.nc, both written into
                         --case-dir and rebuilt on every run. The original
                         facEB and statistics files are never modified; facT is
                         not needed and not merged. The rest of the run then
                         reads the two merged files.
    --skip-atmosphere    Reuse existing pedestrian_atmosphere files.
    --skip-radiation     Reuse existing pedestrian_shortwave/longwave files.
    --tile-size N        Receptors traced per radiation tile (checkpoint
                         granularity). Default: 256.
    --n-mu N             Zenith quadrature points for the sky/facet radiation.
                         Default: 8.
    --n-azimuth N        Azimuth quadrature points. Default: 32.
    --max-gap S          Largest accepted interval (s) between facet records
                         inside one averaging window. Default: 1.1 times the
                         median interval.
    --checkpoint-root DIR
                         Radiation checkpoint directory. Default: inside the case.
    --no-resume          Ignore existing radiation tile checkpoints.
    --overwrite          Replace existing atmosphere, radiation and exchange files.

Outputs (per receptor height, in the case directory)
----------------------------------------------------
    pedestrian_atmosphere.h<height>.<expnr>.nc
    pedestrian_shortwave.h<height>.<expnr>.nc
    pedestrian_longwave.h<height>.<expnr>.nc
    thermal_comfort_inputs.h<height>.<expnr>.nc   (input of calculate_indices.py)

Examples
--------
    python udales_to_comfort.py --case-dir /data/cases/001 --heights 1.1 1.5
    python udales_to_comfort.py --case-dir /data/cases/004 --heights 1.1 \\
        --data-source /data/cases/001 --data-source /data/cases/002 \\
        --data-source /data/cases/003

This script does not run uDALES or repeat the surface energy balance. Every
saved statistics record becomes one output record: the mean over the
``tstatsdump`` seconds preceding its timestamp, with solar angles evaluated
at the midpoint of that window.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from udbase import UDBase
from udcomf.udcomf_export import check_comfort_preflight
from udcomf.sources import MERGED_FACETS, MERGED_STATISTICS
from udcomf.udcomf_io import merge_facet_outputs, merge_statistics_outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate height-specific atmosphere/radiation files and exchange "
            "input files for the comfort indices from existing uDALES outputs."
        )
    )
    parser.add_argument("--case-dir", type=Path, required=True,
                        help="uDALES case directory; all outputs are written here.")
    parser.add_argument(
        "--heights", type=float, nargs="+", required=True,
        help="Receptor heights (m) to process; each must be saved by the simulation.",
    )
    parser.add_argument(
        "--statistics-source", choices=("stats_t", "stats_tkslice"), default="stats_t",
        help=(
            "Gathered statistics file to read, from the case or from every "
            "--data-source segment; default: stats_t."
        ),
    )
    parser.add_argument(
        "--data-source", type=Path, action="append", default=[],
        help=(
            "Warm-start segment case directory; repeat for every segment. Its facEB and "
            "--statistics-source files are found by the uDALES naming rule; t and LWout "
            "are merged into merged_facEB and the comfort subset into merged_stats_tkslice."
        ),
    )
    parser.add_argument(
        "--skip-atmosphere", action="store_true",
        help="Reuse existing pedestrian_atmosphere files from the output directory.",
    )
    parser.add_argument(
        "--skip-radiation", action="store_true",
        help="Reuse existing pedestrian_shortwave/longwave files from the output directory.",
    )
    parser.add_argument("--tile-size", type=int, default=256)
    parser.add_argument("--n-mu", type=int, default=8)
    parser.add_argument("--n-azimuth", type=int, default=32)
    parser.add_argument(
        "--max-gap", type=float, default=None,
        help="Maximum accepted source interval in seconds for a complete radiation window.",
    )
    parser.add_argument(
        "--checkpoint-root", type=Path, default=None,
        help="Radiation checkpoint directory; default: case checkpoint directory.",
    )
    parser.add_argument(
        "--no-resume", action="store_true",
        help="Do not reuse matching radiation tile checkpoints.",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Replace generated atmospheric, radiation, and exchange files.",
    )
    return parser


def _print_preflight(report) -> None:
    print("Comfort preflight:")
    print(
        "  receptor heights: "
        + ", ".join(f"{value:g} m" for value in report.receptor_heights_m)
    )
    schedule = report.schedule
    print(
        f"  statistics      : tstatsdump={schedule.tstatsdump:g} s, "
        f"tstatsgap={schedule.tstatsgap:g} s, tstatstart={schedule.tstatstart:g} s, "
        f"first dump at {report.first_dump_s:g} s"
    )
    for caution in report.cautions:
        print(f"  CAUTION: {caution}")
    for blocker in report.blockers:
        print(f"  BLOCKER: {blocker}")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    case = args.case_dir.expanduser().resolve()
    checkpoint_root = (
        None if args.checkpoint_root is None
        else args.checkpoint_root.expanduser().resolve()
    )
    heights = tuple(args.heights)
    data_sources = [path.expanduser().resolve() for path in args.data_source]
    gathered_source = args.statistics_source
    statistics_source = MERGED_STATISTICS if data_sources else gathered_source
    facet_source = MERGED_FACETS if data_sources else "facEB"

    try:
        report = check_comfort_preflight(case)
        _print_preflight(report)
        if not report.ready:
            raise ValueError("Comfort preflight has blockers; no postprocessing was started")

        if data_sources:
            for path in data_sources:
                if not path.is_dir():
                    raise ValueError(f"--data-source must be a case directory: {path}")
            merged = merge_facet_outputs(case, data_sources)
            print(
                f"Merged facet t and LWout of {len(merged.source_cases)} segments into "
                f"{merged.path.name} ({merged.record_count} records)"
            )
            merged = merge_statistics_outputs(case, data_sources, source=gathered_source)
            print(
                f"Merged the comfort subset of {len(merged.source_cases)} {gathered_source} "
                f"segments into {merged.path.name} ({merged.record_count} records, scalar levels "
                f"{list(k + 1 for k in merged.selected_levels)})"
            )
            print(f"This run reads {statistics_source} and {facet_source}")

        sim = UDBase(case.name, case, load_geometry=not args.skip_radiation)
        times = sim.comf.statistics_times(statistics_source)
        schedule = sim.comf.schedule
        irregular = schedule.off_schedule(times)
        print(
            f"Statistics records: {times.size} from t = {times[0]:g} s to {times[-1]:g} s, "
            f"each a {schedule.tstatsdump:g} s mean"
        )
        if irregular.any():
            print(
                f"  NOTE: {int(irregular.sum())} record(s) are not on the "
                "tstatstart + n*(tstatsdump + tstatsgap) schedule, as happens after a warm "
                "start inside a window; they remain full-length means"
            )
        if not args.skip_atmosphere:
            outputs = sim.comf.atmosphere.write_heights(
                source=statistics_source,
                heights=heights,
                overwrite=args.overwrite,
            )
            for height, path in outputs.items():
                print(f"Atmosphere {height:g} m: {path}")

        if not args.skip_radiation:
            outputs = sim.comf.radiation.write_heights(
                times,
                duration=schedule.tstatsdump,
                heights=heights,
                checkpoint_root=checkpoint_root,
                tile_size=args.tile_size,
                n_mu=args.n_mu,
                n_azimuth=args.n_azimuth,
                max_gap=args.max_gap,
                resume=not args.no_resume,
                overwrite=args.overwrite,
                facet_source=facet_source,
            )
            for height, kinds in outputs.items():
                for kind, path in kinds.items():
                    print(f"{kind.capitalize()} {height:g} m: {path}")

        exchange = sim.comf.export_exchange(heights=heights, overwrite=args.overwrite)
        print("Exchange files written to:")
        for height, path in exchange.items():
            print(f"  {height:g} m: {path}")
    except (
        FileNotFoundError, FileExistsError, OSError, RuntimeError, TypeError, ValueError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
