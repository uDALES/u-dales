#!/usr/bin/env python3
"""
Compare u-DALES input files across multiple experiment cases (all pairwise combinations).

Internally uses ud_compare_inputs.py for all file comparisons — no code duplication.

Usage:
  ud_compare_multiple_inputs.py <exppath1> <exp1> <exppath2> <exp2> [exppath3 exp3 ...] [tolerance]

  exppathN   Parent experiments directory for case N (the exp_num subdirectory is appended)
  expN       Experiment number for case N (integer 1-999)
  tolerance  Max absolute error for numeric/NetCDF comparisons (default: 1e-10)
             Detected automatically as the last argument if it parses as a float.

  At least two cases (4 arguments) are required.
  All pairwise combinations are compared.
  A timestamped log file is written to the current working directory.

Examples:
  ud_compare_multiple_inputs.py tests/system/experiments/ 224 ref_data/ 224
  ud_compare_multiple_inputs.py path_a/experiments/ 100 path_b/experiments/ 100 path_c/ 100 1e-8
"""

import importlib.util
import io
import itertools
import os
import sys
from contextlib import redirect_stdout
from datetime import datetime


# ---------------------------------------------------------------------------
# Load ud_compare_inputs from the same tools/ directory
# ---------------------------------------------------------------------------

def _load_ud_compare_inputs():
    here = os.path.dirname(os.path.abspath(__file__))
    module_path = os.path.join(here, 'ud_compare_inputs.py')
    spec = importlib.util.spec_from_file_location('ud_compare_inputs', module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _usage():
    print("Usage: ud_compare_multiple_inputs.py <exppath1> <exp1> <exppath2> <exp2> [exppath3 exp3 ...] [tolerance]")
    print("")
    print("  exppathN   Parent experiments directory for case N")
    print("  expN       Experiment number for case N (integer 1-999)")
    print("  tolerance  Max absolute error (default: 1e-10); detected as last arg if it parses as a float")
    print("")
    print("  At least two cases (4 arguments) are required.")


def _parse_args(argv):
    """Return list of (exppath, exp_num) pairs and tolerance."""
    args = argv[1:]

    if len(args) < 4:
        _usage()
        sys.exit(1)

    # Detect optional trailing tolerance: last arg is a float and count is odd
    tolerance = 1e-10
    if len(args) % 2 == 1:
        try:
            tolerance = float(args[-1])
            args = args[:-1]
        except ValueError:
            print(f"[ERROR] Odd number of arguments and last argument '{args[-1]}' is not a float tolerance.")
            _usage()
            sys.exit(1)

    if len(args) % 2 != 0:
        print("[ERROR] Arguments must come in pairs: <exppath> <exp_num>.")
        _usage()
        sys.exit(1)

    cases = []
    for i in range(0, len(args), 2):
        exppath = args[i]
        raw = args[i + 1]
        try:
            exp_num = int(raw)
        except ValueError:
            print(f"[ERROR] exp_num must be an integer, got: '{raw}'")
            sys.exit(1)
        if not (1 <= exp_num <= 999):
            print(f"[ERROR] exp_num must be between 1 and 999, got: {exp_num}")
            sys.exit(1)
        cases.append((exppath, exp_num))

    return cases, tolerance


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    cases, tolerance = _parse_args(sys.argv)

    # Load the comparison module and initialise NetCDF libs (exits with help if missing)
    uci = _load_ud_compare_inputs()
    uci.try_import_netcdf()

    # Resolve and validate all case directories
    resolved = []  # list of (case_dir, exp_str)
    for exppath, exp_num in cases:
        exp_str = f"{exp_num:03d}"
        case_dir = os.path.join(exppath, exp_str)
        if not os.path.isdir(case_dir):
            print(f"[ERROR] Directory not found: {case_dir}")
            sys.exit(1)
        resolved.append((case_dir, exp_str))

    n_cases = len(resolved)
    pairs = list(itertools.combinations(range(n_cases), 2))

    # Timestamped log file in the repo-level logdir/
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    _repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    _logdir = os.path.join(_repo_root, "logdir")
    os.makedirs(_logdir, exist_ok=True)
    log_path = os.path.join(_logdir, f"ud_compare_multiple_inputs_{timestamp}.log")

    # Print run header
    sep = "=" * 70
    print(sep)
    print(f"ud_compare_multiple_inputs")
    print(sep)
    print(f"  Cases     : {n_cases}")
    print(f"  Pairs     : {len(pairs)}")
    print(f"  Tolerance : {tolerance}")
    print(f"  Log file  : {log_path}")
    print("")
    for idx, (case_dir, exp_str) in enumerate(resolved):
        print(f"  Case {idx + 1}: {case_dir}  (exp {exp_str})")
    print("")

    pair_results = {}  # (i, j) -> (all_passed: bool, passed: int, total: int, skipped: int)

    with open(log_path, 'w') as log_fh:
        # Log header
        log_fh.write(f"ud_compare_multiple_inputs — {datetime.now().isoformat()}\n")
        log_fh.write(f"Tolerance: {tolerance}\n")
        log_fh.write("Cases:\n")
        for idx, (case_dir, exp_str) in enumerate(resolved):
            log_fh.write(f"  [{idx + 1}] {case_dir}  (exp {exp_str})\n")
        log_fh.write(f"\n{sep}\n\n")

        for i, j in pairs:
            dir1, exp_str1 = resolved[i]
            dir2, exp_str2 = resolved[j]

            pair_header = f"Pair ({i + 1} vs {j + 1}):  exp {exp_str1}  vs  exp {exp_str2}"
            print(f"  Running {pair_header} ...", end='', flush=True)

            log_fh.write(f"{sep}\n{pair_header}\n")
            log_fh.write(f"  Dir 1: {dir1}\n")
            log_fh.write(f"  Dir 2: {dir2}\n")
            log_fh.write(f"{sep}\n")

            # Capture full comparison output — log only
            buf = io.StringIO()
            with redirect_stdout(buf):
                counters = uci.run_comparison(dir1, exp_str1, dir2, exp_str2, tolerance)

            output = buf.getvalue()
            log_fh.write(output)

            passed = counters['pass']
            total  = counters['test']
            skipped = counters['skip']
            all_passed = (passed == total)

            result_line = (f"[{'PASS' if all_passed else 'FAIL'}] "
                           f"{passed}/{total} tests passed"
                           + (f", {skipped} skipped" if skipped else ""))
            print(f" {result_line}")
            log_fh.write(f"\n{result_line}\n\n")

            pair_results[(i, j)] = (all_passed, passed, total, skipped)

        # Write overall summary to log
        log_fh.write(f"\n{sep}\nOVERALL SUMMARY\n{sep}\n")
        for (i, j), (ok, passed, total, skipped) in pair_results.items():
            _, exp_str1 = resolved[i]
            _, exp_str2 = resolved[j]
            skip_note = f", {skipped} skipped" if skipped else ""
            line = (f"  [{'PASS' if ok else 'FAIL'}] "
                    f"Case {i + 1} (exp {exp_str1}) vs Case {j + 1} (exp {exp_str2}): "
                    f"{passed}/{total}{skip_note}")
            log_fh.write(line + "\n")
        overall_ok = all(ok for ok, *_ in pair_results.values())
        log_fh.write("\n" + ("[PASS] All pairwise comparisons PASSED.\n"
                              if overall_ok else "[FAIL] Some pairwise comparisons FAILED.\n"))

    # Terminal: overall summary only
    print(f"\n{sep}")
    print("OVERALL SUMMARY")
    print(sep)
    for (i, j), (ok, passed, total, skipped) in pair_results.items():
        _, exp_str1 = resolved[i]
        _, exp_str2 = resolved[j]
        skip_note = f", {skipped} skipped" if skipped else ""
        print(f"  [{'PASS' if ok else 'FAIL'}] "
              f"Case {i + 1} (exp {exp_str1}) vs Case {j + 1} (exp {exp_str2}): "
              f"{passed}/{total}{skip_note}")
    overall_ok = all(ok for ok, *_ in pair_results.values())
    print("")
    print("[PASS] All pairwise comparisons PASSED." if overall_ok
          else "[FAIL] Some pairwise comparisons FAILED.")
    print(f"\nFull details written to: {log_path}")
    sys.exit(0 if overall_ok else 1)


if __name__ == "__main__":
    main()
