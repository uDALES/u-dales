#!/usr/bin/env python3
"""Calculate MRT, PET, UTCI and WBGT maps from thermal-comfort exchange files.

Usage
-----
    python calculate_indices.py [options] EXCHANGE_FILE [EXCHANGE_FILE ...]

    EXCHANGE_FILE        One or more exchange NetCDF files, for example the
                         thermal_comfort_inputs.h<height>.<expnr>.nc files
                         written by udales_to_comfort.py. Each file is processed
                         on its own and its result is written beside it as
                         thermal_comfort_indices.h<height>.<expnr>.nc (the
                         "thermal_comfort_inputs" prefix replaced). A file whose
                         name does not start with that prefix needs --output.

Options
-------
    --indices LIST       Comma-separated indices to calculate: any of mrt,
                         pet, utci, wbgt (for example --indices mrt,wbgt).
                         Default: all four. Only the exchange fields the chosen
                         indices need are required: mrt reads the radiation
                         planes and solar zenith; utci adds ta, pabs, qv and
                         ws_10; pet adds ta, pabs, qv and ws_local; wbgt reads
                         ta, pabs, qv, ws_local and the radiation planes. Globe
                         and natural wet-bulb temperatures are written with wbgt.
    --parameters FILE    JSON file with every person, sensor and domain
                         parameter. Default: udcomf/thermalcomfort/
                         comfort_defaults.json. To change a value, copy that
                         file, edit it and pass the copy here. The file must
                         contain all keys and no others.
    --output FILE        Output path; allowed only with one exchange file.
    --overwrite          Replace an existing output file.

Examples
--------
    python calculate_indices.py case/thermal_comfort_inputs.h1p1.001.nc
    python calculate_indices.py --indices mrt,wbgt case/thermal_comfort_inputs.h*.001.nc
    python calculate_indices.py --parameters my_parameters.json \\
        --output mrt_only.nc --indices mrt case/thermal_comfort_inputs.h1p1.001.nc

Each exchange file must carry the model grid with cell bounds, the pedestrian
mask, the receptor height, time bounds, per-record solar angles and the
physical fields the selected indices need, with their units. Inputs outside an
index's validated domain are flagged and written as missing, never clipped.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from udcomf.thermalcomfort import (
    DEFAULT_PARAMETERS_PATH, INDEX_NAMES, ComfortParameters, calculate_indices,
    selected_indices,
)


def load_comfort_parameters(path: Path | None) -> ComfortParameters:
    """Read a complete parameter file, or the packaged comfort_defaults.json."""
    source = DEFAULT_PARAMETERS_PATH if path is None else Path(path).expanduser().resolve()
    return ComfortParameters.from_file(source)


def parse_indices(text: str) -> tuple[str, ...]:
    """Comma-separated index names in canonical order; unknown names are rejected."""
    try:
        return selected_indices(name.strip() for name in text.split(",") if name.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate native-grid MRT, PET, UTCI and WBGT from one or more "
            "thermal-comfort exchange files."
        )
    )
    parser.add_argument("exchange_files", type=Path, nargs="+")
    parser.add_argument(
        "--indices", type=parse_indices, default=INDEX_NAMES, metavar="LIST",
        help=(
            "Comma-separated indices to calculate, any of "
            + ",".join(INDEX_NAMES) + "; default: all four."
        ),
    )
    parser.add_argument(
        "--parameters", type=Path, default=None,
        help=f"Complete comfort-parameter JSON; default: {DEFAULT_PARAMETERS_PATH.name}.",
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help="Output NetCDF path; valid only with one exchange file.",
    )
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
                indices=args.indices,
                overwrite=args.overwrite,
            )
            print(f"Thermal-comfort indices ({', '.join(args.indices)}): {result}")
    except (
        FileNotFoundError, FileExistsError, ImportError, OSError,
        RuntimeError, TypeError, ValueError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
