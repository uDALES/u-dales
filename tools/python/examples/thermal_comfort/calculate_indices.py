#!/usr/bin/env python3
"""Calculate model-neutral MRT, PET, UTCI, and WBGT maps.

Accepts exchange NetCDF files with any number of time records; optional
metadata (bounds, georeference, CRS, provenance, pedestrian mask) is used
when present. The physical fields, their units and per-record solar angles
are required. The separate intercomparison-protocol tool, maintained
outside this repository, enforces the strict 24-hour exchange contract.
"""

from __future__ import annotations

import argparse
from dataclasses import fields
import json
from pathlib import Path
import sys

from udcomf.thermalcomfort import ComfortParameters, calculate_indices


def load_comfort_parameters(path: Path | None) -> ComfortParameters:
    """Read shared comfort parameters, or use the documented defaults."""
    if path is None:
        return ComfortParameters()
    config_path = Path(path).expanduser().resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)
    if not isinstance(config, dict):
        raise ValueError("Comfort-parameter JSON must contain one object")
    allowed = {field.name for field in fields(ComfortParameters)}
    unexpected = sorted(config.keys() - allowed)
    if unexpected:
        raise ValueError(f"Unknown comfort parameters: {unexpected}")
    return ComfortParameters(**config)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate native-grid MRT, PET, UTCI, and WBGT from one or more "
            "model-neutral thermal-comfort exchange files."
        )
    )
    parser.add_argument("exchange_files", type=Path, nargs="+")
    parser.add_argument(
        "--parameters", type=Path, default=None,
        help="Shared ComfortParameters JSON; default: documented protocol values.",
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help="Output NetCDF path; valid only with one exchange file.",
    )
    parser.add_argument("--tile-x", type=int, default=64)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.output is not None and len(args.exchange_files) != 1:
        print("ERROR: --output can be used only with one exchange file", file=sys.stderr)
        return 2

    try:
        parameters = load_comfort_parameters(args.parameters)
        seen: set[Path] = set()
        for source in args.exchange_files:
            source = source.expanduser().resolve()
            if source in seen:
                raise ValueError(f"Duplicate exchange input: {source}")
            seen.add(source)
            output = (
                None if args.output is None
                else args.output.expanduser().resolve()
            )
            result = calculate_indices(
                source,
                output_path=output,
                parameters=parameters,
                overwrite=args.overwrite,
                tile_x=args.tile_x,
            )
            print(f"Thermal-comfort indices: {result}")
    except (
        FileNotFoundError, FileExistsError, ImportError, OSError,
        RuntimeError, TypeError, ValueError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
