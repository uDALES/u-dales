#!/usr/bin/env python3
"""Refresh committed MATLAB preprocessing outputs for a test case.

This is a developer tool. It runs the MATLAB preprocessing path in a temporary
copy of a case directory and copies the generated reference outputs back into
the source case directory. Normal integration tests do not require MATLAB.
"""

import argparse
import shutil
import tempfile
from pathlib import Path

from compare_preprocessing import (
    REPO_ROOT,
    _cleanup_targets,
    _copy_case_tree,
    _discover_outputs,
    _run_matlab,
    _temp_parent,
)


def harvest_case(case_dir: Path, keep_temp: bool = False) -> int:
    case_dir = case_dir.expanduser().resolve()
    case = case_dir.name
    outputs = _discover_outputs(case_dir, case)
    if not outputs:
        raise SystemExit(f"No supported preprocessing outputs found in {case_dir}")

    temp_ctx = tempfile.TemporaryDirectory(
        prefix=f"udales-preproc-harvest-{case}-",
        dir=str(_temp_parent()),
    )
    temp_root = Path(temp_ctx.name)
    if keep_temp:
        print(f"Keeping temp directory: {temp_root}")

    try:
        matlab_case = temp_root / case
        _copy_case_tree(case_dir, matlab_case)
        for target in _cleanup_targets(matlab_case, case):
            target.unlink()

        print(f"Refreshing MATLAB preprocessing references for case {case}")
        _run_matlab(matlab_case, case, outputs)

        for relpath in outputs:
            shutil.copy2(matlab_case / relpath, case_dir / relpath)
            print(f"  updated {relpath}")
    finally:
        if keep_temp:
            temp_ctx._finalizer.detach()
        else:
            temp_ctx.cleanup()

    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case-dir",
        required=True,
        help="Case directory whose committed MATLAB reference outputs should be refreshed",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep the temporary MATLAB run directory for inspection",
    )
    args = parser.parse_args()
    return harvest_case(Path(args.case_dir), keep_temp=args.keep_temp)


if __name__ == "__main__":
    raise SystemExit(main())
