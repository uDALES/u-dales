#!/usr/bin/env python3
"""Prepare thermal-comfort exchange files from uDALES outputs.

Usage
-----
    python udales_to_comfort.py --case-dir CASE_DIR [options]

    --case-dir DIR       The uDALES case directory (required). It must hold
                         namoptions.<expnr>, the preprocessing files, the
                         gathered facet outputs facEB/facT and the gathered
                         time-averaged statistics. The site longitude, latitude
                         and elevation, the model axis rotation (xazimuth), the
                         UTC start date, the receptor heights and the statistics
                         schedule (tstatstart, tstatsdump, tstatsgap) are all
                         read from namoptions.

Options
-------
    --statistics-source {stats_tkslice,stats_t}
                         Statistics file to read the receptor-level means from.
                         Default: stats_tkslice.
    --heights H [H ...]  Receptor heights (m) to process; they must be saved
                         receptor heights of the case. Default: all.
    --output-dir DIR     Where the intermediate planes and exchange files go.
                         Default: the case directory.
    --facet-source PATH  Warm-start segment (case directory or facEB.<case>.nc)
                         whose facet records are merged into this case first.
                         Repeat for every segment.
    --stats-source PATH  Warm-start segment (case directory or
                         stats_tkslice/stats_t.<case>.nc) whose comfort subset is
                         merged into stats_tkslice of this case first. Repeat
                         for every segment.
    --overwrite-merge    Replace an existing merged facEB/facT or stats_tkslice.
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

Outputs (per receptor height, in the output directory)
-----------------------------------------------------
    pedestrian_atmosphere.h<height>.<expnr>.nc
    pedestrian_shortwave.h<height>.<expnr>.nc
    pedestrian_longwave.h<height>.<expnr>.nc
    thermal_comfort_inputs.h<height>.<expnr>.nc   (input of calculate_indices.py)

Example
-------
    python udales_to_comfort.py --case-dir /data/cases/001 --heights 1.1 1.5

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
from udcomf.udcomf_io import merge_facet_outputs, merge_statistics_outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate height-specific atmosphere/radiation files and exchange "
            "input files for the comfort indices from existing uDALES outputs."
        )
    )
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument(
        "--statistics-source", choices=("stats_tkslice", "stats_t"),
        default="stats_tkslice",
    )
    parser.add_argument(
        "--heights", type=float, nargs="+", default=None,
        help="Configured receptor heights to process; default: all configured heights.",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="Directory for intermediate planes and exchange files; default: case directory.",
    )
    parser.add_argument(
        "--facet-source", type=Path, action="append", default=[],
        help=(
            "Warm-start case directory or facEB.<case>.nc file. Repeat as needed; "
            "records are sorted and checked before merging into the case."
        ),
    )
    parser.add_argument(
        "--stats-source", type=Path, action="append", default=[],
        help=(
            "Warm-start case directory or stats_tkslice/stats_t.<case>.nc file. Repeat as "
            "needed; only the comfort subset is merged into stats_tkslice of the case."
        ),
    )
    parser.add_argument(
        "--overwrite-merge", action="store_true",
        help="Replace existing merged facEB/facT or stats_tkslice files in the case.",
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
    output_dir = (
        case if args.output_dir is None else args.output_dir.expanduser().resolve()
    )
    checkpoint_root = (
        None if args.checkpoint_root is None
        else args.checkpoint_root.expanduser().resolve()
    )
    heights = None if args.heights is None else tuple(args.heights)

    try:
        report = check_comfort_preflight(case)
        _print_preflight(report)
        if not report.ready:
            raise ValueError("Comfort preflight has blockers; no postprocessing was started")

        if args.facet_source:
            merged = merge_facet_outputs(
                case, args.facet_source, overwrite=args.overwrite_merge,
            )
            print(
                f"Merged {len(merged.source_cases)} facet segments into "
                f"{merged.eb_path.name} and {merged.temperature_path.name} "
                f"({merged.record_count} records)"
            )
        if args.stats_source:
            if args.statistics_source != "stats_tkslice":
                raise ValueError("--stats-source writes stats_tkslice; use --statistics-source stats_tkslice")
            merged = merge_statistics_outputs(
                case, args.stats_source, overwrite=args.overwrite_merge,
            )
            print(
                f"Merged {len(merged.source_cases)} statistics segments into "
                f"{merged.path.name} ({merged.record_count} records, scalar levels "
                f"{list(k + 1 for k in merged.selected_levels)})"
            )

        sim = UDBase(case.name, case, load_geometry=not args.skip_radiation)
        times = sim.comf.statistics_times(args.statistics_source)
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
                source=args.statistics_source,
                heights=heights,
                output_dir=output_dir,
                overwrite=args.overwrite,
            )
            for height, path in outputs.items():
                print(f"Atmosphere {height:g} m: {path}")

        if not args.skip_radiation:
            outputs = sim.comf.radiation.write_heights(
                times,
                duration=schedule.tstatsdump,
                heights=heights,
                output_dir=output_dir,
                checkpoint_root=checkpoint_root,
                tile_size=args.tile_size,
                n_mu=args.n_mu,
                n_azimuth=args.n_azimuth,
                max_gap=args.max_gap,
                resume=not args.no_resume,
                overwrite=args.overwrite,
            )
            for height, kinds in outputs.items():
                for kind, path in kinds.items():
                    print(f"{kind.capitalize()} {height:g} m: {path}")

        exchange = sim.comf.export_exchange(
            heights=heights,
            input_dir=output_dir,
            output_dir=output_dir,
            overwrite=args.overwrite,
        )
        for height, path in exchange.items():
            print(f"Exchange {height:g} m: {path}")
    except (
        FileNotFoundError, FileExistsError, OSError, RuntimeError, TypeError, ValueError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
