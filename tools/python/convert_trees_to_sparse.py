#!/usr/bin/env python3
"""Convert legacy trees.inp vegetation to sparse vegetation inputs."""

import argparse
import sys
from pathlib import Path


def main() -> int:
    repo_python_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(repo_python_dir))

    from udprep import UDPrep

    parser = argparse.ArgumentParser()
    parser.add_argument("expnr")
    parser.add_argument("path")
    args = parser.parse_args()

    prep = UDPrep(expnr=args.expnr, path=args.path, load_geometry=False)
    out = prep.vegetation.block_to_veg()
    print(f"Conversion complete: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
