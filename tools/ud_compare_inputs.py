#!/usr/bin/env python3
"""
Compare u-DALES input files between two experiment directories.

Text comparison (exact, unified diff on mismatch):
  namoptions.<exp>

Numerical comparison (max absolute error, # comments skipped):
  prof.inp.<exp>          lscale.inp.<exp>       probe.inp.<exp>
  facetarea.inp.<exp>     facets.inp.<exp>
  facets_unused.<exp>     factypes.inp.<exp>
  netsw.inp.<exp>         sveg.inp.<exp>         svf.inp.<exp>      vfsparse.inp.<exp>
  heatpump.inp.<exp>
  trees.inp.<exp>         veg.inp.<exp>          veg_params.inp.<exp>
  scalar.inp.<exp>        scalarsourcep.inp.N.<exp>  scalarsourcel.inp.N.<exp>  (N = 1, 2, ...)
  Tfacinit.inp.<exp>      timedeplw.inp.<exp>      timedepsw.inp.<exp>      Sdir.txt
  facet_sections_c/u/v/w.txt
  fluid_boundary_c/u/v/w.txt
  solid_c/u/v/w.txt

NetCDF comparison (all variables, max absolute error, requires netCDF4):
  vf.nc.inp.<exp>

Usage:
  ud_compare_inputs.py <exp_num> <exppath> <ref_path> [tolerance]

  exp_num   Experiment number (e.g. 100)
  exppath   Path to the experiments directory (e.g. tests/system/experiments/); the case subdirectory <exp_str> is appended automatically
  ref_path  Parent directory containing reference input cases (e.g. /path/to/ref_data)
  tolerance Max absolute error for numeric/NetCDF files (default: 1e-10)
"""

import difflib
import os
import sys

# Module-level NetCDF library handles (populated by try_import_netcdf)
nc = None
np = None

def try_import_netcdf():
    global nc, np

    try:
        import numpy as np_module
        import netCDF4 as nc_module
        np = np_module
        nc = nc_module
        return
    except ImportError:
        pass

    _tools_dir = os.path.dirname(os.path.abspath(__file__))
    _setup_script = os.path.join(_tools_dir, 'python', 'setup_venv.sh')
    print("netCDF4 or numpy not available.")
    print("Run the following command to set up the required environment:")
    print(f"   bash {_setup_script} <common|icl>")
    sys.exit(1)


def compare_text_file(filename: str, paths: list, counters: dict) -> None:
    """Exact text comparison with unified diff on mismatch."""
    exists = [os.path.exists(p) for p in paths]

    if not any(exists):
        counters['skip'] += 1
        print(f"[SKIP] {filename} not found in either directory - skipping")
        return

    if not all(exists):
        counters['test'] += 1
        for present, path in zip(exists, paths):
            if not present:
                print(f"[FAIL] {filename}: missing in one directory: {path}")
        return

    counters['test'] += 1
    contents = []
    for path in paths:
        with open(path, 'r') as f:
            contents.append(f.read())

    if contents[0] == contents[1]:
        print(f"[PASS] {filename}: files are identical")
        counters['pass'] += 1
    else:
        diff = list(difflib.unified_diff(
            contents[0].splitlines(keepends=True),
            contents[1].splitlines(keepends=True),
            fromfile=paths[0],
            tofile=paths[1],
        ))
        print(f"[FAIL] {filename}: files differ")
        for line in diff:
            print(line, end='')


def compare_numeric_file(
    filename: str,
    paths: list,
    tolerance: float,
    counters: dict,
    skiprows: int = 0,
) -> None:
    """Numerical comparison using numpy; comment lines starting with # are skipped.
    Files that are empty (or contain only comments/header rows) are treated as identical if both are empty."""
    if np is None:
        print(f"[ERROR] numpy is required for numeric comparison of {filename}")
        sys.exit(1)

    exists = [os.path.exists(p) for p in paths]

    if not any(exists):
        counters['skip'] += 1
        print(f"[SKIP] {filename} not found in either directory - skipping")
        return

    if not all(exists):
        counters['test'] += 1
        for present, path in zip(exists, paths):
            if not present:
                print(f"[FAIL] {filename}: missing in one directory: {path}")
        return

    counters['test'] += 1
    arrays = []
    for path in paths:
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                data = np.loadtxt(path, comments='#', skiprows=skiprows)
            arrays.append(data)
        except Exception as e:
            print(f"[ERROR] {filename}: Failed to read {path}: {e}")
            return

    a, b = arrays

    # Both empty
    if a.size == 0 and b.size == 0:
        print(f"[PASS] {filename}: both files are empty")
        counters['pass'] += 1
        return

    # One empty, one not
    if a.size == 0 or b.size == 0:
        empty_idx = 0 if a.size == 0 else 1
        print(f"[FAIL] {filename}: file {empty_idx + 1} is empty, file {2 - empty_idx} has {max(a.size, b.size)} elements")
        return

    if a.shape != b.shape:
        print(f"[FAIL] {filename}: shape mismatch {a.shape} vs {b.shape}")
        return

    max_error = np.max(np.abs(a - b))
    if max_error <= tolerance:
        print(f"[PASS] {filename}: max error = {max_error:.2e}")
        counters['pass'] += 1
    else:
        print(f"[FAIL] {filename}: max error = {max_error:.2e}")


def compare_netcdf_file(filename: str, paths: list, tolerance: float, counters: dict) -> None:
    """Numerical comparison of a NetCDF file; all numeric variables are compared."""
    if nc is None:
        print(f"[WARN] {filename}: NetCDF libraries not available - skipping")
        counters['skip'] += 1
        return

    exists = [os.path.exists(p) for p in paths]

    if not any(exists):
        counters['skip'] += 1
        print(f"[SKIP] {filename} not found in either directory - skipping")
        return

    if not all(exists):
        counters['test'] += 1
        for present, path in zip(exists, paths):
            if not present:
                print(f"[FAIL] {filename}: missing in one directory: {path}")
        return

    def read_all_vars(path):
        """Return dict of {var_name: flattened float array} for all numeric variables."""
        result = {}
        try:
            with nc.Dataset(path, 'r') as ds:
                for name, var in ds.variables.items():
                    if np.issubdtype(var.dtype, np.number):
                        result[name] = var[:].flatten().astype(float)
        except Exception as e:
            print(f"[ERROR] {filename}: Failed to read {path}: {e}")
        return result

    vars1 = read_all_vars(paths[0])
    vars2 = read_all_vars(paths[1])

    if not vars1 and not vars2:
        print(f"[WARN] {filename}: no numeric variables found in either file - skipping")
        counters['skip'] += 1
        return

    all_vars = sorted(set(vars1) | set(vars2))
    overall_pass = True
    max_overall = 0.0

    for var in all_vars:
        if var not in vars1 or var not in vars2:
            print(f"[WARN] {filename}: variable '{var}' missing in one file")
            overall_pass = False
            continue
        a, b = vars1[var], vars2[var]
        if a.shape != b.shape:
            print(f"[FAIL] {filename} [{var}]: shape mismatch {a.shape} vs {b.shape}")
            overall_pass = False
            continue
        err = np.max(np.abs(a - b))
        max_overall = max(max_overall, err)
        if err > tolerance:
            overall_pass = False

    counters['test'] += 1
    if overall_pass:
        print(f"[PASS] {filename}: max error = {max_overall:.2e}")
        counters['pass'] += 1
    else:
        print(f"[FAIL] {filename}: max error = {max_overall:.2e}")


def run_comparison(dir1: str, exp_str1: str, dir2: str, exp_str2: str, tolerance: float) -> dict:
    """Run all input file comparisons between two case directories.

    dir1, dir2       Absolute paths to the two case directories.
    exp_str1/str2    Zero-padded experiment number strings (e.g. '224').
                     May differ when comparing cases with different experiment numbers.
    tolerance        Max absolute error for numeric/NetCDF files.

    Returns counters dict with keys 'pass', 'test', 'skip'.
    """
    counters = {'pass': 0, 'test': 0, 'skip': 0}

    def pp(name1, name2=None):
        """Return [path_in_dir1, path_in_dir2]; name2 defaults to name1 (fixed-name files)."""
        return [os.path.join(dir1, name1), os.path.join(dir2, name2 if name2 is not None else name1)]

    # namoptions: exact text comparison
    print("\n=== Comparing namoptions ===")
    compare_text_file(
        f"namoptions.{exp_str1}",
        pp(f"namoptions.{exp_str1}", f"namoptions.{exp_str2}"),
        counters)

    # Vertical profile input files
    print("\n=== Comparing other input files ===")
    for tag in ("prof.inp", "lscale.inp", "probe.inp", "scalar.inp"):
        compare_numeric_file(
            f"{tag}.{exp_str1}",
            pp(f"{tag}.{exp_str1}", f"{tag}.{exp_str2}"),
            tolerance, counters)

    # Scalar source files: scalarsourcep/scalarsourcel.inp.N.exp (scalar specie for N = 1..9)
    for source_kind in ("scalarsourcep", "scalarsourcel"):
        for n in range(1, 10):
            b1 = f"{source_kind}.inp.{n}.{exp_str1}"
            b2 = f"{source_kind}.inp.{n}.{exp_str2}"
            paths_n = pp(b1, b2)
            if any(os.path.exists(p) for p in paths_n):
                compare_numeric_file(b1, paths_n, tolerance, counters)

    # Facet-related input files
    for tag in (
        "facetarea.inp", "facets.inp", "facets_unused",
        "factypes.inp", "svf.inp", "vfsparse.inp", "netsw.inp", "sveg.inp",
        "Tfacinit.inp", "heatpump.inp",
        "trees.inp", "veg.inp", "veg_params.inp",
    ):
        compare_numeric_file(
            f"{tag}.{exp_str1}",
            pp(f"{tag}.{exp_str1}", f"{tag}.{exp_str2}"),
            tolerance, counters)

    # timedeplw.inp has two text header rows followed by numeric data.
    compare_numeric_file(
        f"timedeplw.inp.{exp_str1}",
        pp(f"timedeplw.inp.{exp_str1}", f"timedeplw.inp.{exp_str2}"),
        tolerance, counters, skiprows=2)

    # timedepsw.inp has one text header row followed by numeric time and facet data rows.
    compare_numeric_file(
        f"timedepsw.inp.{exp_str1}",
        pp(f"timedepsw.inp.{exp_str1}", f"timedepsw.inp.{exp_str2}"),
        tolerance, counters, skiprows=1)

    # vf.nc.inp: binary NetCDF file
    compare_netcdf_file(
        f"vf.nc.inp.{exp_str1}",
        pp(f"vf.nc.inp.{exp_str1}", f"vf.nc.inp.{exp_str2}"),
        tolerance, counters)

    # Fixed-name text files (no exp suffix)
    print("\n=== Comparing fixed-name input files ===")
    for basename in (
        "Sdir.txt",
        "facet_sections_c.txt", "facet_sections_u.txt",
        "facet_sections_v.txt", "facet_sections_w.txt",
        "fluid_boundary_c.txt", "fluid_boundary_u.txt",
        "fluid_boundary_v.txt", "fluid_boundary_w.txt",
        "solid_c.txt", "solid_u.txt", "solid_v.txt", "solid_w.txt",
    ):
        compare_numeric_file(basename, pp(basename), tolerance, counters)

    return counters


def main():
    # Ensure NetCDF libraries are available (exits with help message if not)
    try_import_netcdf()

    if len(sys.argv) < 4:
        print("Usage: ud_compare_inputs.py <exp_num> <exppath> <ref_path> [tolerance]")
        print("")
        print("  exp_num   Experiment number (e.g. 100)")
        print("  exppath   Path to the local experiment directory for this case")
        print("  ref_path  Parent directory containing reference input cases (e.g. /path/to/ref_data)")
        print("  tolerance Max absolute error for numeric files (default: 1e-10)")
        sys.exit(1)

    try:
        exp_num = int(sys.argv[1])
    except ValueError:
        print(f"[ERROR] exp_num must be an integer, got: '{sys.argv[1]}'")
        sys.exit(1)
    if not (1 <= exp_num <= 999):
        print(f"[ERROR] exp_num must be between 1 and 999, got: {exp_num}")
        sys.exit(1)

    exppath = sys.argv[2]
    ref_data_path = sys.argv[3]
    tolerance = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-10

    exp_str = f"{exp_num:03d}"
    dirs = [
        os.path.join(exppath, exp_str),
        os.path.join(ref_data_path, exp_str),
    ]

    missing = [d for d in dirs if not os.path.isdir(d)]
    if missing:
        for d in missing:
            print(f"[ERROR] Directory not found: {d}")
        sys.exit(1)

    print(f"Comparing input files for experiment {exp_str}:")
    print(f"  Dir 1: {dirs[0]}")
    print(f"  Dir 2: {dirs[1]}")
    print(f"  Numeric tolerance: {tolerance}")

    counters = run_comparison(dirs[0], exp_str, dirs[1], exp_str, tolerance)

    # Summary
    print(f"\n=== SUMMARY ===")
    print(f"{counters['pass']} out of {counters['test']} tests passed.")
    if counters['skip'] > 0:
        print(f"{counters['skip']} tests skipped (missing files).")

    if counters['test'] == 0:
        print("[FAIL] No tests were performed!")
        sys.exit(1)
    elif counters['pass'] == counters['test']:
        print("[PASS] All tests PASSED!" if counters['skip'] == 0
              else "[PASS] All attempted tests PASSED! (some files were skipped)")
        sys.exit(0)
    else:
        print("[FAIL] Some tests FAILED!")
        sys.exit(1)


if __name__ == "__main__":
    main()
