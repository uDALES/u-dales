#!/usr/bin/env python3
"""Run the whole thermal-comfort chain of a gathered uDALES case in one go.

Usage
-----
    python ud_run_thermal_comfort.py --case-dir CASE_DIR --heights H [H ...] [options]

After the simulation has run and its outputs have been gathered, this one
command produces thermal_comfort_indices.h<height>.<expnr>.nc in the case
directory for every requested receptor height. It performs, in order:

    1. preflight   check namoptions and the gathered files
    2. merge       (only with --data-source) merge warm-start segments into
                   merged_facEB.<expnr>.nc and merged_stats_tkslice.<expnr>.nc
    3. atmosphere  receptor-level means -> pedestrian_atmosphere.h<h>.<expnr>.nc
    4. radiation   receptor shortwave/longwave -> pedestrian_shortwave/longwave.h<h>.<expnr>.nc
    5. exchange    thermal_comfort_inputs.h<h>.<expnr>.nc
    6. indices     thermal_comfort_indices.h<h>.<expnr>.nc

Required
--------
    --case-dir DIR       The uDALES case directory (name = experiment number).
    --heights H [H ...]  Receptor heights in metres, each saved by the simulation.

Options
-------
    --indices LIST       Comma-separated subset of mrt,pet,utci,wbgt. Default: all.
    --parameters FILE    Complete comfort-parameter JSON (copy and edit
                         udcomf/thermalcomfort/comfort_defaults.json, or the
                         copy in examples/thermal_comfort/comfort_parameters.json).
                         Default: the packaged defaults.
    --statistics-source {stats_t,stats_tkslice}
                         Gathered statistics file to read; default stats_t.
    --data-source DIR    Warm-start segment case directory; repeat per segment.
    --skip-atmosphere, --skip-radiation
                         Reuse the existing per-height files of the case.
    --tile-size, --n-mu, --n-azimuth, --max-gap, --checkpoint-root, --no-resume
                         Radiation settings, as in udales_to_comfort.py.
    --overwrite          Replace existing atmosphere, radiation, exchange and
                         index files. Without it an existing file is an error.
    --log-file FILE      Log file; default CASE_DIR/thermal_comfort.<expnr>.log
                         (appended). Messages also go to the console.

The radiation step is the expensive one; it checkpoints every tile and resumes
by default, so re-running the same command after an interruption continues
where it stopped. To recompute only the indices with another parameter file,
run calculate_indices.py on the existing exchange files instead.

Example
-------
    python ud_run_thermal_comfort.py --case-dir /data/cases/001 --heights 1.1 1.5 \\
        --parameters my_comfort_parameters.json
"""

from __future__ import annotations

import argparse
from importlib.util import find_spec
from pathlib import Path
import sys

# Run as a file from inside the udcomf package directory, Python would put that
# directory first on sys.path, where udcomf.py shadows the udcomf package.
_HERE = Path(__file__).resolve().parent
if sys.path and Path(sys.path[0]).resolve() == _HERE:
    sys.path[0] = str(_HERE.parent)

from udcomf import calculate_indices as index_step
from udcomf import udales_to_comfort as exchange_step
from udcomf.heights import height_tag
from udcomf.logs import configure as configure_logging, header, stage_logger


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Thermal-comfort indices from a gathered uDALES case in one command: "
            "preflight, warm-start merge, atmosphere, radiation, exchange and indices."
        )
    )
    exchange_step.add_case_arguments(parser)
    exchange_step.add_arguments(parser)
    index_step.add_arguments(parser)
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Replace existing atmosphere, radiation, exchange and index files.",
    )
    parser.add_argument(
        "--log-file", type=Path, default=None,
        help="Log file; default: CASE_DIR/thermal_comfort.<expnr>.log (appended).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    case = args.case_dir.expanduser().resolve()
    if not case.is_dir():
        print(f"ERROR: case directory does not exist: {case}", file=sys.stderr)
        return 2
    log_file = (case / f"thermal_comfort.{case.name}.log" if args.log_file is None
                else args.log_file.expanduser().resolve())
    configure_logging(log_file)
    log = stage_logger("pipeline")
    header("pipeline", "thermal-comfort chain for case %s, heights %s, indices %s (log %s)",
           case.name, ", ".join(f"{h:g} m" for h in args.heights), ", ".join(args.indices),
           log_file)
    if {"pet", "utci"} & set(args.indices) and find_spec("pythermalcomfort") is None:
        log.error("PET and UTCI need pythermalcomfort; install it with "
                  "pip install -e 'tools/python[comfort]' or choose --indices mrt,wbgt")
        return 1

    try:
        exchange = exchange_step.run(args)
    except exchange_step.PIPELINE_ERRORS as exc:
        log.error("%s", exc)
        log.error("stopped before the index step")
        return 1
    expected = {h: case / f"thermal_comfort_inputs.{height_tag(h)}.{case.name}.nc" for h in args.heights}
    for height, path in expected.items():
        if exchange.get(height) != path:
            log.error("exchange file for %g m not where expected: %s", height, path)
            return 1

    try:
        results = index_step.run(
            [expected[h] for h in args.heights], indices=args.indices,
            parameters=args.parameters, overwrite=args.overwrite,
        )
    except index_step.INDEX_ERRORS as exc:
        log.error("%s", exc)
        log.error("exchange files are complete; rerun calculate_indices.py on them after fixing the problem")
        return 1

    header("pipeline", "finished; thermal-comfort indices written to:")
    for result in results:
        log.info("%s", result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
