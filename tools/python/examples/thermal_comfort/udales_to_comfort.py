#!/usr/bin/env python3
"""Prepare model-neutral thermal-comfort exchange files from uDALES outputs.

This script does not run uDALES or repeat the surface energy balance. It reads
saved statistics, facet radiation outputs, and preprocessing files from one
case directory. Verified georeferencing is supplied separately because it
cannot be inferred safely from a local or rotated LES grid.
"""

from __future__ import annotations

import argparse
from datetime import date
import json
from pathlib import Path
import sys

import numpy as np

from udbase import UDBase
from udcomf.udcomf_export import ExchangeMetadata, check_comfort_preflight
from udcomf.udcomf_io import merge_facet_outputs


_ARRAY_FIELDS = (
    "x", "y", "x_bounds", "y_bounds", "longitude", "latitude", "z_ground",
)
_TEXT_FIELDS = (
    "vertical_datum", "model_version", "institution",
    "building_representation", "terrain_convention",
    "spinup_start_utc", "spinup_end_utc",
)


def _analysis_day(value: str) -> date:
    try:
        return date.fromisoformat(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected YYYY-MM-DD") from exc


def load_exchange_metadata(path: Path) -> ExchangeMetadata:
    """Load checked scalar metadata and georeference arrays without pickle."""
    config_path = Path(path).expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)
    if not isinstance(config, dict):
        raise ValueError("Metadata JSON must contain one object")

    required = {"georeference_npz", "crs_attributes", *_TEXT_FIELDS}
    missing = sorted(required - config.keys())
    unexpected = sorted(config.keys() - required)
    if missing or unexpected:
        raise ValueError(
            f"Metadata JSON fields are invalid; missing={missing}, unexpected={unexpected}"
        )
    if not isinstance(config["georeference_npz"], str):
        raise TypeError("georeference_npz must be a path string")
    crs = config["crs_attributes"]
    if not isinstance(crs, dict) or not all(
        isinstance(key, str) and isinstance(value, (str, int, float))
        and not isinstance(value, bool)
        for key, value in crs.items()
    ):
        raise TypeError("crs_attributes must map strings to string or numeric values")
    for field in _TEXT_FIELDS:
        if not isinstance(config[field], str) or not config[field].strip():
            raise ValueError(f"Metadata field {field} must be a nonempty string")

    arrays_path = Path(config["georeference_npz"]).expanduser()
    if not arrays_path.is_absolute():
        arrays_path = config_path.parent / arrays_path
    arrays_path = arrays_path.resolve()
    with np.load(arrays_path, allow_pickle=False) as archive:
        missing_arrays = sorted(set(_ARRAY_FIELDS) - set(archive.files))
        if missing_arrays:
            raise ValueError(f"Georeference NPZ lacks arrays: {missing_arrays}")
        arrays = {
            name: np.asarray(archive[name], dtype=float).copy()
            for name in _ARRAY_FIELDS
        }

    return ExchangeMetadata(
        **arrays,
        crs_attributes=dict(crs),
        **{field: config[field] for field in _TEXT_FIELDS},
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate height-specific atmosphere/radiation files and strict "
            "model-neutral exchange files from existing uDALES outputs."
        )
    )
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--analysis-day", type=_analysis_day, required=True)
    parser.add_argument(
        "--metadata", type=Path, required=True,
        help="JSON metadata file pointing to the verified georeference NPZ.",
    )
    parser.add_argument(
        "--statistics-source", choices=("stats_kslice", "stats_t"),
        default="stats_kslice",
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
            "records are sorted and checked before merging."
        ),
    )
    parser.add_argument(
        "--overwrite-facet-merge", action="store_true",
        help="Replace existing merged facEB/facT files in the destination case.",
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
    print(
        f"  requested hours  : {report.analysis_times_s[0]:g} to "
        f"{report.analysis_times_s[-1]:g} s"
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
        report = check_comfort_preflight(case, args.analysis_day)
        _print_preflight(report)
        if not report.ready:
            raise ValueError("Comfort preflight has blockers; no postprocessing was started")
        metadata = load_exchange_metadata(args.metadata)

        if args.facet_source:
            merged = merge_facet_outputs(
                case, args.facet_source, overwrite=args.overwrite_facet_merge,
            )
            print(
                f"Merged {len(merged.source_cases)} facet segments into "
                f"{merged.eb_path.name} and {merged.temperature_path.name} "
                f"({merged.record_count} records)"
            )

        sim = UDBase(case.name, case, load_geometry=not args.skip_radiation)
        target_times = np.asarray(report.analysis_times_s, dtype=float)
        if not args.skip_atmosphere:
            outputs = sim.comf.atmosphere.write_hourly_heights(
                target_times,
                source=args.statistics_source,
                heights=heights,
                output_dir=output_dir,
                overwrite=args.overwrite,
            )
            for height, path in outputs.items():
                print(f"Atmosphere {height:g} m: {path}")

        if not args.skip_radiation:
            outputs = sim.comf.radiation.write_hourly_heights(
                target_times,
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
            args.analysis_day,
            metadata,
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
