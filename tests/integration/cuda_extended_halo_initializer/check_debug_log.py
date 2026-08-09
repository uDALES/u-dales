#!/usr/bin/env python3
"""Verify that a Debug GPU run completed the extended-halo initializer test."""

from __future__ import annotations

import argparse
from pathlib import Path


PASS_MARKER = "CUDA extended-halo initfield self-test passed."
FAIL_MARKER = "CUDA extended-halo initfield self-test failed."


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path, help="uDALES Debug GPU output log")
    args = parser.parse_args()

    contents = args.log.read_text(encoding="utf-8", errors="replace")
    if FAIL_MARKER in contents:
        raise SystemExit(f"initializer self-test failed in {args.log}")
    if PASS_MARKER not in contents:
        raise SystemExit(f"initializer self-test did not run in {args.log}")

    print(f"initializer self-test passed in {args.log}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
