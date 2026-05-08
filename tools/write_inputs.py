"""Python preprocessing entry point for uDALES case directories.

This script mirrors the role of the legacy MATLAB `write_inputs.m` workflow,
but drives the Python `UDPrep` stack instead. It is intended to be run on a
case directory containing `namoptions.<expnr>`.

Usage
-----
    All commands below assume the working directory is the repository root.

    python tools/write_inputs.py [case_dir] [--force]

Arguments
---------
    case_dir   Path to a case directory containing namoptions.<expnr>.
               The experiment number is inferred from the directory name
               (e.g. a directory named "101" sets expnr=101).

    --force    Force regeneration of radiation outputs even when they
               already exist. Only has an effect when surface energy
               balance is enabled (lEB = True in namoptions). Bypasses
               the output-file existence checks in run_short_wave (which
               would otherwise skip if Sdir.txt and netsw.inp.<expnr>
               are present) and run_short_wave_timedep (which would
               otherwise skip if Sdir.nc and timedepsw.inp.<expnr> are
               present). All other preprocessing sections (grid, IBM, IC,
               scalars, vegetation, SEB) are unaffected by this flag.

Default behaviour
-----------------
    When invoked with no arguments the script prints help and exits.

Environment
-----------
    The script requires the Python virtual environment created by
    tools/python/setup_venv.sh. Following commands can be run
    from the repository root to set up and activate it:

        bash tools/python/setup_venv.sh
        source tools/python/.venv/bin/activate
        python tools/write_inputs.py

    If the virtual environment is not active, a clear error message is
    printed with the exact commands needed to set it up.

Examples
--------
    Example commands from the repository root.

    # Process the default example case (examples/999)
    python tools/write_inputs.py

    # Process a specific case directory
    python tools/write_inputs.py examples/101

    # Force regeneration of all outputs
    python tools/write_inputs.py examples/101 --force
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


# Import the in-repo Python tooling without requiring an installed package.
script_dir = Path(__file__).resolve().parent
tools_python = script_dir / "python"
if str(tools_python) not in sys.path:
    sys.path.insert(0, str(tools_python))


try:
    from udprep import UDPrep  # noqa: E402
except ImportError:
    _setup_script = script_dir / "python" / "setup_venv.sh"
    print("Required Python packages are not available.")
    print("Set up the virtual environment by running:")
    print(f"   {_setup_script}")
    print("Then activate it and re-run this script as below:")
    print(f"   source {script_dir / 'python' / '.venv' / 'bin' / 'activate'}")
    print(f"   python {script_dir / Path(__file__).name}")
    sys.exit(1)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Python preprocessing for a uDALES case directory.")
    parser.add_argument("case_dir", help="Path to case directory containing namoptions.<expnr>")
    parser.add_argument("--force", action="store_true", help="Force regeneration where supported.")
    args = parser.parse_args(argv)

    case_dir = Path(args.case_dir).resolve()

    print("Initializing UDPrep...")
    sys.stdout.flush()
    prep = UDPrep(case_dir, load_geometry=True, suppress_load_warnings=True)
    print("-------------------------------------------------------------------")
    print("Summary of the derived preprocessing configuration")
    print(prep)
    print("-------------------------------------------------------------------")
    # Run the configured preprocessing sections.
    # Namelist writeback (write_changed_params) is intentionally NOT called
    # here — only sections that derive new values (radiation, vegetation)
    # call it internally inside their own run_all/save methods.
    prep.run_all(force=args.force)
    
    ### Optional inspection plots (uncomment to enable)
    # prep.sim.vis.plot_profiles()
    # prep.sim.vis.plot_lscale()
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
