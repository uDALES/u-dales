#!/usr/bin/env python3
"""Calculate MRT, PET, UTCI and WBGT maps from thermal-comfort exchange files.

Usage
-----
    python calculate_indices.py [options] EXCHANGE_FILE [EXCHANGE_FILE ...]

    EXCHANGE_FILE        One or more exchange NetCDF files, for example the
                         thermal_comfort_inputs.h<height>.<expnr>.nc files
                         written by udales_to_comfort.py. Each file is processed
                         on its own and its result is written in the same
                         directory as thermal_comfort_indices.h<height>.<expnr>.nc
                         (the "thermal_comfort_inputs" prefix replaced). File
                         names must start with that prefix. The full path of
                         every result is printed at the end.

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
    --overwrite          Replace an existing output file.

Examples
--------
    python calculate_indices.py case/thermal_comfort_inputs.h1p1.001.nc
    python calculate_indices.py --indices mrt,wbgt case/thermal_comfort_inputs.h*.001.nc
    python calculate_indices.py --parameters my_parameters.json \\
        --indices mrt case/thermal_comfort_inputs.h1p1.001.nc

Each exchange file must carry the model grid with cell bounds, the pedestrian
mask, the receptor height, time bounds, per-record solar angles and the
physical fields the selected indices need, with their units. Inputs outside an
index's validated domain are flagged and written as missing, never clipped.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

# Run as a file from inside the udcomf package directory, Python would put that
# directory first on sys.path, where udcomf.py shadows the udcomf package.
_HERE = Path(__file__).resolve().parent
if sys.path and Path(sys.path[0]).resolve() == _HERE:
    sys.path[0] = str(_HERE.parent)

from udcomf.logs import configure as configure_logging, header, stage_logger
from udcomf.thermalcomfort import (
    DEFAULT_PARAMETERS_PATH, INDEX_NAMES, ComfortParameters, calculate_indices,
    selected_indices,
)

INDEX_ERRORS = (
    FileNotFoundError, FileExistsError, ImportError, OSError, RuntimeError, TypeError,
    ValueError,
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


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Options of the index step; shared with the one-step wrapper."""
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate native-grid MRT, PET, UTCI and WBGT from one or more "
            "thermal-comfort exchange files."
        )
    )
    parser.add_argument("exchange_files", type=Path, nargs="+")
    add_arguments(parser)
    parser.add_argument("--overwrite", action="store_true",
                        help="Replace existing index files.")
    return parser


def run(exchange_files, *, indices, parameters: Path | None, overwrite: bool) -> list[Path]:
    """Calculate the selected indices for every exchange file; return the index files.

    Raises one of ``INDEX_ERRORS`` on any failure.
    """
    header("indices", "calculating %s for %d exchange file(s) with parameters %s",
           ", ".join(indices), len(exchange_files),
           DEFAULT_PARAMETERS_PATH.name if parameters is None else parameters)
    config = load_comfort_parameters(parameters)
    seen: set[Path] = set()
    results: list[Path] = []
    for source in exchange_files:
        source = Path(source).expanduser().resolve()
        if source in seen:
            raise ValueError(f"Duplicate exchange input: {source}")
        seen.add(source)
        results.append(calculate_indices(
            source, parameters=config, indices=indices, overwrite=overwrite,
        ))
    return results


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(None)
    log = stage_logger("pipeline")
    try:
        results = run(args.exchange_files, indices=args.indices,
                      parameters=args.parameters, overwrite=args.overwrite)
    except INDEX_ERRORS as exc:
        log.error("%s", exc)
        return 1
    header("pipeline", "thermal-comfort indices (%s) written to:", ", ".join(args.indices))
    for result in results:
        log.info("%s", result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
