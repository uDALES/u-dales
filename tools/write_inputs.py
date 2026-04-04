"""Python preprocessing entry point for uDALES case directories.

This script mirrors the role of the legacy MATLAB `write_inputs.m` workflow,
but drives the Python `UDPrep` stack instead. It is intended to be run on a
case directory containing `namoptions.<expnr>`.
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


from udprep import UDPrep  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Python preprocessing for a uDALES case directory.")
    parser.add_argument("case_dir", nargs="?", help="Path to case directory containing namoptions.<expnr>")
    parser.add_argument("--force", action="store_true", help="Force regeneration where supported.")
    args = parser.parse_args(argv)

    if args.case_dir:
        case_dir = Path(args.case_dir).resolve()
    else:
        # Keep the historical default example behavior when no case is given.
        udales_root = script_dir.parent
        case_dir = (udales_root / "examples" / "999").resolve()

    expnr = case_dir.name
    print("Initializing UDPrep...")
    sys.stdout.flush()
    prep = UDPrep(expnr, case_dir, load_geometry=True)
    print("-------------------------------------------------------------------")
    print("Summary of the derived preprocessing configuration")
    print(prep)
    print("-------------------------------------------------------------------")
    # Run the configured preprocessing sections and persist any derived
    # namelist updates needed by downstream solver workflows.
    prep.run_all(force=args.force)
    prep.write_changed_params()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
