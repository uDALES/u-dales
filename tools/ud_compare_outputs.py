#!/usr/bin/env python3
"""
NetCDF-based comparison tool for u-dales simulation outputs.
Requires netCDF4 and numpy.
"""

import os
import sys
from typing import List, Dict, Any

nc = None
np = None

DEFAULT_TOLERANCE = 1e-6
DEFAULT_TOL_THL = 1e-6

def try_import_netcdf():
    global nc, np

    try:
        import numpy as np_module
        import netCDF4 as nc_module
        np = np_module
        nc = nc_module
        return
    except ImportError:
        _setup_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ud_set_nc_venv.sh')
        print("netCDF4 or numpy not available.")
        print("Run the following command to set up the required environment:")
        print(f"   {_setup_script}")
        sys.exit(1)

class SimulationComparison:
    def __init__(self, output_dirs: List[str], exp_num: int, tolerance: float = DEFAULT_TOLERANCE, tol_thl: float = DEFAULT_TOL_THL):
        self.output_dirs = output_dirs
        self.exp_num = exp_num
        self.exp_str = f"{exp_num:03d}"
        self.tolerance = tolerance
        self.tol_thl = tol_thl
        self.result_count = 0
        self.test_count = 0
        self.skip_count = 0
        
        # NetCDF file paths for each directory
        self.nc_files = {}
        for i, dir_path in enumerate(output_dirs):
            self.nc_files[i] = {
                'xytdump': os.path.join(dir_path, f'xytdump.{self.exp_str}.nc'),
                'tdump': os.path.join(dir_path, f'tdump.{self.exp_str}.nc'),
                'fielddump': os.path.join(dir_path, f'fielddump.{self.exp_str}.nc'),
                'treedump': os.path.join(dir_path, f'treedump.{self.exp_str}.nc'),
            }
    
    def ncread(self, file_path: str, var_name: str):
        """Read a variable from a NetCDF file using netCDF4."""
        try:
            with nc.Dataset(file_path, 'r') as dataset:
                if var_name in dataset.variables:
                    return dataset.variables[var_name][:].flatten()
                return None
        except Exception as e:
            print(f"Warning: Failed to read {var_name} from {file_path}: {e}")
            return None
    
    def load_variable(self, file_type: str, var_name: str) -> Dict[int, Any]:
        """Load a variable from both NetCDF files. Assumes both files exist."""
        data = {}
        for i in range(len(self.output_dirs)):
            file_path = self.nc_files[i][file_type]
            var_data = self.ncread(file_path, var_name)
            if var_data is not None:
                data[i] = var_data
        return data

    def _check_file_pair(self, file_type: str) -> str:
        """Check existence of a dump file in both directories.

        Returns:
          'skip'    -- absent from both directories (not counted)
          'fail'    -- present in one directory only (counted as failure)
          'compare' -- present in both directories (proceed with comparison)
        """
        paths = [self.nc_files[i][file_type] for i in range(len(self.output_dirs))]
        exists = [os.path.exists(p) for p in paths]

        if not any(exists):
            return 'skip'

        if not all(exists):
            self.test_count += 1
            for present, path in zip(exists, paths):
                if not present:
                    print(f"[FAIL] {os.path.basename(path)}: file missing in one directory: {path}")
            return 'fail'

        return 'compare'

    def _compare_variable(self, file_type: str, var_name: str, description: str, tol: float) -> None:
        """Load a variable and apply the three-way skip/fail/compare logic.

        - Not in either file : skip (increment skip_count, no test recorded)
        - In one file only   : fail (increment test_count, no result recorded)
        - In both files      : numerical comparison via error_test
        """
        data = self.load_variable(file_type, var_name)
        label = f"{var_name} ({description})"

        if len(data) == 0:
            self.skip_count += 1
            print(f"[SKIP] {label}: not present in either {file_type} file - skipping")
        elif len(data) < 2:
            self.test_count += 1
            missing_idx = next(i for i in range(len(self.output_dirs)) if i not in data)
            print(f"[FAIL] {label}: variable missing in {self.output_dirs[missing_idx]}")
        else:
            passed = self.error_test(data, label, tol)
            if not passed:
                return

    def error_test(self, data: Dict[int, Any], name: str, tolerance: float) -> bool:
        """Perform numerical comparison between datasets"""
        self.test_count += 1

        try:
            data1 = data[0]
            data2 = data[1]
            
            # Ensure same shape
            if data1.shape != data2.shape:
                print(f"[FAIL] {name}: shape mismatch {data1.shape} vs {data2.shape}")
                return False
            
            # Calculate maximum absolute difference
            diff = np.abs(data1 - data2)
            max_error = np.max(diff)
            
            # Check if within tolerance
            if max_error <= tolerance:
                print(f"[PASS] {name}: max error = {max_error:.2e} (tolerance = {tolerance:.2e})")
                self.result_count += 1
                return True
            else:
                print(f"[FAIL] {name}: max error = {max_error:.2e} (tolerance = {tolerance:.2e})")
                return False
                
        except Exception as e:
            print(f"[ERROR] {name}: Failed to compare - {e}")
            return False
    
    def compare_xyt_variables(self) -> None:
        """Compare variables from xytdump files"""
        print("\n=== Comparing XYT dump variables ===")

        status = self._check_file_pair('xytdump')
        if status == 'skip':
            print("[SKIP] xytdump files not present in either directory - skipping")
            return
        if status == 'fail':
            return
        
        xyt_variables = [
            # Basic velocity components
            ('uxyt', 'Streamwise velocity'),
            ('vxyt', 'Spanwise velocity'), 
            ('wxyt', 'Vertical velocity'),
            
            # Velocity variances
            ('upuptxyc', 'Variance of streamwise velocity'),
            ('vpvptxyc', 'Variance of spanwise velocity'),
            ('wpwptxyc', 'Variance of vertical velocity'),
            ('wwxyt', 'Kinematic mom. flux ww'),
            
            # TKE
            ('tketxyc', 'TKE'),
            
            # Reynolds stresses and momentum fluxes
            ('upwpxyt', 'Reynolds stress u\'w\''),
            ('uwxyt', 'Kinematic mom. flux uw'),
            ('vpwpxyt', 'Reynolds stress v\'w\''),
            ('vwxyt', 'Kinematic mom. flux vw'),
            ('upvpxyt', 'Reynolds stress u\'v\''),
            ('uvxyt', 'Kinematic mom. flux uv'),
            
            # Temperature related
            ('thlxyt', 'Liquid potential temperature'),
            ('wthlxyt', 'Kinematic heat flux'),
            ('wpthlpxyt', 'Turbulent heat flux'),
            ('thlpthlptxy', 'Variance of liquid potential temperature'),
            
            # Moisture
            ('qtxyt', 'Moisture'),
            
            # Pressure
            ('pxyt', 'Pressure'),
            
            # SGS fluxes
            ('usgsxyt', 'SGS momentum flux u'),
            ('vsgsxyt', 'SGS momentum flux v'),
            ('thlsgsxyt', 'SGS temperature flux'),
        ]
        
        for var_name, description in xyt_variables:
            tol = self.tol_thl if 'thl' in var_name.lower() else self.tolerance
            self._compare_variable('xytdump', var_name, description, tol)

    def compare_t_variables(self) -> None:
        """Compare variables from tdump files"""
        print("\n=== Comparing T dump variables ===")

        status = self._check_file_pair('tdump')
        if status == 'skip':
            print("[SKIP] tdump files not present in either directory - skipping")
            return
        if status == 'fail':
            return
        
        t_variables = [
            # Basic velocity components
            ('ut', 'Streamwise velocity'),
            ('vt', 'Spanwise velocity'),
            ('wt', 'Vertical velocity'),
            
            # Velocity variances
            ('upuptc', 'Variance of streamwise velocity'),
            ('vpvptc', 'Variance of spanwise velocity'),
            ('wpwptc', 'Variance of vertical velocity'),
            
            # TKE
            ('tketc', 'TKE'),
            
            # Reynolds stresses
            ('upwpt', 'Reynolds stress u\'w\''),
            
            # Temperature related
            ('thlt', 'Liquid potential temperature'),
            ('wpthlpt', 'Turbulent heat flux w\'thl\''),
            ('thlpthlpt', 'Variance of liquid potential temperature'),
            
            # Moisture
            ('qtt', 'Moisture'),
            
            # Pressure
            ('pt', 'Pressure'),
            
            # Scalar concentrations
            ('sca1t', 'Scalar concentration 1'),
            ('sca2t', 'Scalar concentration 2'),
            ('wpsca1pt', 'Turbulent scalar flux 1'),
            ('wpsca2pt', 'Turbulent scalar flux 2'),
            ('sca1psca1pt', 'Variance of scalar concentration 1'),
            ('sca2psca2pt', 'Variance of scalar concentration 2'),
            
            # SGS variables
            ('sv1sgs', 'sv1sgs'),
            ('sv2sgs', 'sv2sgs'),
        ]
        
        for var_name, description in t_variables:
            tol = self.tol_thl if 'thl' in var_name.lower() else self.tolerance
            self._compare_variable('tdump', var_name, description, tol)

    def compare_field_variables(self) -> None:
        """Compare variables from fielddump files"""
        print("\n=== Comparing Field dump variables ===")

        status = self._check_file_pair('fielddump')
        if status == 'skip':
            print("[SKIP] fielddump files not present in either directory - skipping")
            return
        if status == 'fail':
            return
        
        field_variables = [
            ('u',    'Streamwise velocity'),
            ('v',    'Spanwise velocity'),
            ('w',    'Vertical velocity'),
            ('thl',  'Liquid potential temperature'),
            ('qt',   'Moisture'),
            ('sca1', 'Scalar concentration 1'),
        ]
        
        for var_name, description in field_variables:
            tol = self.tol_thl if 'thl' in var_name.lower() else self.tolerance
            self._compare_variable('fielddump', var_name, description, tol)

    def compare_tree_variables(self) -> None:
        """Compare variables from treedump files"""
        print("\n=== Comparing Tree dump variables ===")

        status = self._check_file_pair('treedump')
        if status == 'skip':
            print("[SKIP] treedump files not present in either directory - skipping")
            return
        if status == 'fail':
            return

        tree_variables = [
            ('tr_u',     'Streamwise velocity in tree'),
            ('tr_v',     'Spanwise velocity in tree'),
            ('tr_w',     'Vertical velocity in tree'),
            ('tr_thl',   'Liquid potential temperature in tree'),
            ('tr_qt',    'Moisture in tree'),
            ('tr_qtR',   'Rain water in tree'),
            ('tr_qtA',   'Total water in tree'),
            ('tr_sv1',   'Scalar concentration 1 in tree'),
            ('tr_sv2',   'Scalar concentration 2 in tree'),
            ('tr_omega', 'Vorticity in tree'),
        ]

        for var_name, description in tree_variables:
            tol = self.tol_thl if 'thl' in var_name.lower() else self.tolerance
            self._compare_variable('treedump', var_name, description, tol)

    def run_comparison(self) -> bool:
        """Run the complete numerical comparison"""
        print(f"Comparing simulation outputs:")
        for i, dir_path in enumerate(self.output_dirs):
            print(f"  Directory {i+1}: {dir_path}")
        print(f"Experiment number: {self.exp_num}")
        print(f"Tolerance: {self.tolerance}")
        print(f"Temperature tolerance: {self.tol_thl}")

        self.compare_xyt_variables()
        self.compare_t_variables()
        self.compare_field_variables()
        self.compare_tree_variables()
        
        # Final summary
        print(f"\n=== SUMMARY ===")
        print(f"{self.result_count} out of {self.test_count} tests passed.")
        if self.skip_count > 0:
            print(f"{self.skip_count} tests skipped (missing variables/files).")
        
        total_attempted = self.test_count + self.skip_count
        if total_attempted == 0:
            print("[FAIL] No comparison tests were attempted!")
            return False
        elif self.test_count == 0:
            print("[FAIL] No comparison tests were performed (all skipped)!")
            return False
        elif self.result_count == self.test_count:
            if self.skip_count > 0:
                print("[PASS] All attempted tests PASSED! (some variables were skipped)")
            else:
                print("[PASS] All tests PASSED!")
            return True
        else:
            print("[FAIL] Some tests FAILED!")
            return False


def run_comparison(dir1, exp_str1, dir2, exp_str2, tolerance=DEFAULT_TOLERANCE, tol_thl=DEFAULT_TOL_THL):
    """Compare two output directories pairwise.

    Returns a dict with keys 'pass', 'test', 'skip' so callers can aggregate
    results uniformly (mirrors the interface of ud_compare_inputs.run_comparison).
    """
    exp_num = int(exp_str1)
    comparison = SimulationComparison([dir1, dir2], exp_num, tolerance, tol_thl)
    # If the two cases carry different exp_strs, fix the file paths for dir2
    if exp_str1 != exp_str2:
        for ftype in ('xytdump', 'tdump', 'fielddump', 'treedump'):
            comparison.nc_files[1][ftype] = os.path.join(dir2, f'{ftype}.{exp_str2}.nc')
    comparison.run_comparison()
    return {
        'pass': comparison.result_count,
        'test': comparison.test_count,
        'skip': comparison.skip_count,
    }


def main():
    """Main function to run the comparison"""
    # Try to import NetCDF libraries (exits if not available)
    try_import_netcdf()

    if len(sys.argv) < 4:
        print("Usage: ud_compare_outputs.py <exp_num> <exppath> <ref_data_path> [tolerance] [tol_thl]")
        print("")
        print("  exp_num        Experiment number (e.g. 100)")
        print("  exppath        Path to the outputs directory (e.g. tests/system/outputs/); the case subdirectory <exp_str> is appended automatically")
        print("  ref_data_path  Parent directory containing reference output cases (e.g. /path/to/ref_data)")
        print("  tolerance      Max absolute error tolerance (default: 1e-6)")
        print("  tol_thl        Tolerance for temperature variables (default: 1e-6)")
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
    tolerance = float(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_TOLERANCE
    tol_thl = float(sys.argv[5]) if len(sys.argv) > 5 else DEFAULT_TOL_THL

    # Set up output directories
    exp_str = f"{exp_num:03d}"
    output_dirs = [
        os.path.join(exppath, exp_str),
        os.path.join(ref_data_path, exp_str),
    ]

    missing = [d for d in output_dirs if not os.path.isdir(d)]
    if missing:
        for d in missing:
            print(f"[ERROR] Directory not found: {d}")
        sys.exit(1)
    
    # Create comparison object and run
    comparison = SimulationComparison(output_dirs, exp_num, tolerance, tol_thl)
    
    try:
        success = comparison.run_comparison()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"Error during comparison: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
