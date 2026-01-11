#!/usr/bin/env python3
"""
gather_out1d.py - Concatenate x-slab output files using Python/netCDF4

Usage: python gather_out1d.py exp_number [input_dir] [output_dir]

Example: python gather_out1d.py 106 ./ ./

Processes multiple file types:
  X-direction concatenation (xt/xm):
    - test_t_out.xxx.expnr.nc  -> test_t_out.expnr.nc
    - ins_kslice.xxx.expnr.nc      -> ins_kslice.expnr.nc
    - stats_t_out.xxx.expnr.nc -> stats_t_out.expnr.nc
    - ins_islice.xxx.expnr.nc      -> ins_islice.expnr.nc
    - field_out.xxx.expnr.nc   -> field_out.expnr.nc
    - test_islice.xxx.expnr.nc -> test_islice.expnr.nc
    - test_kslice.xxx.expnr.nc -> test_kslice.expnr.nc
  
  Y-direction concatenation (yt/ym):
    - ins_jslice.xxx.expnr.nc      -> ins_jslice.expnr.nc
    - test_jslice.xxx.expnr.nc -> test_jslice.expnr.nc
"""

import os
import sys
import glob
import re
import time
import multiprocessing
import numpy as np
from netCDF4 import Dataset

# --- Configuration ---
# Set to False to disable compression (Faster writing, larger files)
USE_COMPRESSION = False
# ---------------------

def parse_args():
    # Check for required exp_num argument
    if len(sys.argv) < 2:
        print("missing input parameters, finish.")
        print("e.g.: python gather_out1d.py exp_number [input_dir] [output_dir]")
        sys.exit(1)
    
    exp_num = sys.argv[1]
    input_dir = sys.argv[2] if len(sys.argv) > 2 else "./"
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "./"
    
    return exp_num, input_dir, output_dir

def parse_x_index(filename, file_prefix):
    """Extract x-index from filename like prefix.XXX.expnr.nc"""
    base = os.path.basename(filename)
    # Look for pattern: prefix.XXX.YYY.nc where XXX is numeric
    # Escape dots in prefix for regex
    escaped_prefix = re.escape(file_prefix)
    pattern = f'{escaped_prefix}\.(\d+)\.(\d+)\.nc$'
    match = re.search(pattern, base)
    if match:
        return int(match.group(1))  # Return first number (x-index)
    return None

def copy_attributes(src_var, dst_var):
    """Copy all attributes from source to destination variable (except _FillValue)"""
    for attr in src_var.ncattrs():
        if attr != '_FillValue':  # Skip _FillValue - it's handled during creation
            try:
                setattr(dst_var, attr, getattr(src_var, attr))
            except Exception as e:
                print(f"    Warning: Could not copy attribute '{attr}': {e}")

def create_variable_with_fillvalue(out_ds, var_name, src_var, dimensions):
    """Create a variable with proper _FillValue handling"""
    # Get fill value if it exists
    fill_value = getattr(src_var, '_FillValue', None)
    
    # Create variable with or without fill_value
    # PERFORMANCE FIX: We intentionally ignore fill_value during creation to prevent 
    # HDF5 from pre-filling the entire file with this value, which causes massive delays.
    # We assume the concatenation will cover the entire variable range.
    
    # Note: This means the output variable will NOT have a _FillValue attribute,
    # and unwritten regions (if any) will contain garbage or zeros instead of fill_value.
    out_var = out_ds.createVariable(var_name, src_var.datatype, dimensions, zlib=USE_COMPRESSION)
    
    return out_var

def process_jslice_file_type(file_prefix, exp_num, input_dir, output_dir):
    """Process jslice type files that concatenate along y-direction (yt/ym)"""
    
    start_time = time.time()
    
    print(f"\n{'='*60}")
    print(f"Processing (Y-direction): {file_prefix}.*.{exp_num}.nc")
    print(f"{'='*60}")
    
    # Find all matching files
    pattern = os.path.join(input_dir, f"{file_prefix}.*.{exp_num}.nc")
    all_files = glob.glob(pattern)
    
    # Filter out files that don't have numeric y-indices
    files = []
    for f in all_files:
        y_idx = parse_x_index(f, file_prefix)  # Reuse same function, just interpret as y-index
        if y_idx is not None:
            files.append(f)
    
    if not files:
        print(f"  No files found matching pattern: {pattern}")
        print(f"  Skipping {file_prefix}")
        return False
    
    # Sort files by y-index
    files.sort(key=lambda f: parse_x_index(f, file_prefix))
    
    print(f"  Found {len(files)} files to concatenate.")
    
    # Output filename
    output_file = os.path.join(output_dir, f"{file_prefix}.{exp_num}.nc")
    
    # Remove existing output file
    if os.path.exists(output_file):
        print(f"  Removing existing output file: {output_file}")
        os.remove(output_file)
    
    # Verify output file doesn't exist in our input list
    if output_file in files:
        print(f"  Error: Output file {output_file} is in input file list!")
        return False
    
    print(f"  Creating output file: {output_file}")
    
    # Open first file to get metadata
    try:
        with Dataset(files[0], 'r') as first_ds:
            print("  Step 1: Auto-detecting variables and dimensions...")
            
            # Get all variables and separate coordinates from data variables
            coord_vars = {'time', 'xt', 'xm', 'yt', 'ym', 'zt', 'zm'}
            all_vars = set(first_ds.variables.keys())
            data_vars = sorted(list(all_vars - coord_vars))
            
            # Get dimensions from first file
            dims = {}
            for dim_name in first_ds.dimensions:
                dims[dim_name] = len(first_ds.dimensions[dim_name])
            
            # Calculate global y dimensions (concatenate along y)
            global_yt = dims.get('yt', 0) * len(files)
            global_ym = dims.get('ym', 0) * len(files)
            
            print(f"    Global dimensions: yt={global_yt}, ym={global_ym}")
            
            print("  Step 2: Creating output file structure...")
            
            # Create output dataset
            os.makedirs(output_dir, exist_ok=True)
            # Use 'w' mode with clobber=True
            with Dataset(output_file, 'w', format='NETCDF4') as out_ds:
                
                # Copy global attributes
                for attr in first_ds.ncattrs():
                    setattr(out_ds, attr, getattr(first_ds, attr))
                
                # Create dimensions
                for dim_name, size in dims.items():
                    if dim_name == 'yt':
                        out_ds.createDimension(dim_name, global_yt)
                    elif dim_name == 'ym':
                        out_ds.createDimension(dim_name, global_ym)
                    elif dim_name == 'time' and first_ds.dimensions[dim_name].isunlimited():
                        out_ds.createDimension(dim_name, None)
                    else:
                        out_ds.createDimension(dim_name, size)
                
                # Create all variables upfront
                # Store variable names that depend on Y for later processing
                y_dependent_vars = []
                
                # 1. Create Coordinate Variables
                for var_name in sorted(coord_vars & all_vars):
                    src_var = first_ds.variables[var_name]
                    out_var = create_variable_with_fillvalue(out_ds, var_name, src_var, src_var.dimensions)
                    copy_attributes(src_var, out_var)
                    
                    if var_name in ['yt', 'ym']:
                        y_dependent_vars.append(var_name)
                    else:
                        # Static coordinates (xt, zt, time etc) - copy from first file immediately
                        out_var[:] = src_var[:]
                
                # 2. Create Data Variables
                for var_name in data_vars:
                    src_var = first_ds.variables[var_name]
                    out_var = create_variable_with_fillvalue(out_ds, var_name, src_var, src_var.dimensions)
                    copy_attributes(src_var, out_var)
                    
                    if 'yt' in src_var.dimensions or 'ym' in src_var.dimensions:
                        y_dependent_vars.append(var_name)
                    else:
                        # Variable doesn't depend on Y, copy from first file
                        out_var[:] = src_var[:]

                print("  Step 3: Processing files sequentially...")
                
                # Process files one by one to save memory and file handles
                current_offset = 0
                
                for i, f in enumerate(files):
                    if i % 10 == 0:
                        print(f"    Processing file {i+1}/{len(files)}: {os.path.basename(f)}")
                        
                    with Dataset(f, 'r') as ds:
                        # Determine local size for offset update (using yt or ym)
                        local_size = 0
                        if 'yt' in ds.dimensions:
                            local_size = len(ds.dimensions['yt'])
                        elif 'ym' in ds.dimensions:
                            local_size = len(ds.dimensions['ym'])
                        
                        # Process all Y-dependent variables for this file
                        for var_name in y_dependent_vars:
                            if var_name not in ds.variables:
                                continue
                                
                            data = ds.variables[var_name][:]
                            out_var = out_ds.variables[var_name]
                            var_dims = out_var.dimensions
                            
                            # Find the Y-axis index
                            y_axis = -1
                            if 'yt' in var_dims:
                                y_axis = var_dims.index('yt')
                            elif 'ym' in var_dims:
                                y_axis = var_dims.index('ym')
                            
                            if y_axis != -1:
                                # Build slice
                                slices = [slice(None)] * len(var_dims)
                                slices[y_axis] = slice(current_offset, current_offset + data.shape[y_axis])
                                out_var[tuple(slices)] = data
                        
                        current_offset += local_size

    except Exception as e:
        print(f"  Error during processing: {e}")
        if os.path.exists(output_file):
            os.remove(output_file)
        return False
    
    elapsed_time = time.time() - start_time
    print(f"  ✓ Concatenation complete: {output_file}")
    print(f"  ⏱ Time elapsed: {elapsed_time:.2f} seconds")
    return True

def process_file_type(file_prefix, exp_num, input_dir, output_dir):
    """Process one type of file (e.g., test_t_out, kslice, etc.)"""
    
    start_time = time.time()
    
    print(f"\n{'='*60}")
    print(f"Processing: {file_prefix}.*.{exp_num}.nc")
    print(f"{'='*60}")
    
    # Find all matching files
    pattern = os.path.join(input_dir, f"{file_prefix}.*.{exp_num}.nc")
    all_files = glob.glob(pattern)
    
    # Filter out files that don't have numeric x-indices
    files = []
    for f in all_files:
        x_idx = parse_x_index(f, file_prefix)
        if x_idx is not None:
            files.append(f)
    
    if not files:
        print(f"  No files found matching pattern: {pattern}")
        print(f"  Skipping {file_prefix}")
        return False
    
    # Sort files by x-index
    files.sort(key=lambda f: parse_x_index(f, file_prefix))
    
    print(f"  Found {len(files)} files to concatenate.")
    
    # Output filename
    output_file = os.path.join(output_dir, f"{file_prefix}.{exp_num}.nc")
    
    # Remove existing output file
    if os.path.exists(output_file):
        print(f"  Removing existing output file: {output_file}")
        os.remove(output_file)
    
    # Verify output file doesn't exist in our input list
    if output_file in files:
        print(f"  Error: Output file {output_file} is in input file list!")
        return False
    
    print(f"  Creating output file: {output_file}")
    
    try:
        # Open first file to get metadata
        with Dataset(files[0], 'r') as first_ds:
            print("  Step 1: Auto-detecting variables and dimensions...")
            
            # Get all variables and separate coordinates from data variables
            coord_vars = {'time', 'xt', 'xm', 'yt', 'ym', 'zt', 'zm'}
            all_vars = set(first_ds.variables.keys())
            data_vars = sorted(list(all_vars - coord_vars))
            
            # Get dimensions from first file
            dims = {}
            for dim_name in first_ds.dimensions:
                dims[dim_name] = len(first_ds.dimensions[dim_name])
            
            # Calculate global x dimensions
            global_xt = dims.get('xt', 0) * len(files)
            global_xm = dims.get('xm', 0) * len(files)
            
            print(f"    Global dimensions: xt={global_xt}, xm={global_xm}")
            
            print("  Step 2: Creating output file structure...")
            
            # Create output dataset
            os.makedirs(output_dir, exist_ok=True)
            with Dataset(output_file, 'w', format='NETCDF4') as out_ds:
                
                # Copy global attributes
                for attr in first_ds.ncattrs():
                    setattr(out_ds, attr, getattr(first_ds, attr))
                
                # Create dimensions
                for dim_name, size in dims.items():
                    if dim_name == 'xt':
                        out_ds.createDimension(dim_name, global_xt)
                    elif dim_name == 'xm':
                        out_ds.createDimension(dim_name, global_xm)
                    elif dim_name == 'time' and first_ds.dimensions[dim_name].isunlimited():
                        out_ds.createDimension(dim_name, None)
                    else:
                        out_ds.createDimension(dim_name, size)
                
                # Create all variables upfront
                x_dependent_vars = []
                
                # 1. Create Coordinate Variables
                for var_name in sorted(coord_vars & all_vars):
                    src_var = first_ds.variables[var_name]
                    out_var = create_variable_with_fillvalue(out_ds, var_name, src_var, src_var.dimensions)
                    copy_attributes(src_var, out_var)
                    
                    if var_name in ['xt', 'xm']:
                        x_dependent_vars.append(var_name)
                    else:
                        # Static coordinates - copy from first file
                        out_var[:] = src_var[:]
                
                # 2. Create Data Variables
                for var_name in data_vars:
                    src_var = first_ds.variables[var_name]
                    out_var = create_variable_with_fillvalue(out_ds, var_name, src_var, src_var.dimensions)
                    copy_attributes(src_var, out_var)
                    
                    if 'xt' in src_var.dimensions or 'xm' in src_var.dimensions:
                        x_dependent_vars.append(var_name)
                    else:
                        # Variable doesn't depend on X, copy from first file
                        out_var[:] = src_var[:]
                
                print("  Step 3: Processing files sequentially...")
                
                # Process files one by one
                current_offset = 0
                
                for i, f in enumerate(files):
                    if i % 10 == 0:
                        print(f"    Processing file {i+1}/{len(files)}: {os.path.basename(f)}")
                        
                    with Dataset(f, 'r') as ds:
                        # Determine local size for offset update
                        local_size = 0
                        if 'xt' in ds.dimensions:
                            local_size = len(ds.dimensions['xt'])
                        elif 'xm' in ds.dimensions:
                            local_size = len(ds.dimensions['xm'])
                        
                        # Process all X-dependent variables for this file
                        for var_name in x_dependent_vars:
                            if var_name not in ds.variables:
                                continue
                                
                            data = ds.variables[var_name][:]
                            out_var = out_ds.variables[var_name]
                            var_dims = out_var.dimensions
                            
                            # Find the X-axis index
                            x_axis = -1
                            if 'xt' in var_dims:
                                x_axis = var_dims.index('xt')
                            elif 'xm' in var_dims:
                                x_axis = var_dims.index('xm')
                            
                            if x_axis != -1:
                                # Build slice
                                slices = [slice(None)] * len(var_dims)
                                slices[x_axis] = slice(current_offset, current_offset + data.shape[x_axis])
                                out_var[tuple(slices)] = data
                        
                        current_offset += local_size

    except Exception as e:
        print(f"  Error during processing: {e}")
        if os.path.exists(output_file):
            os.remove(output_file)
        return False
    
    elapsed_time = time.time() - start_time
    print(f"  ✓ Concatenation complete: {output_file}")
    print(f"  ⏱ Time elapsed: {elapsed_time:.2f} seconds")
    return True


def worker_task(args):
    """Worker function for multiprocessing"""
    func, file_prefix, exp_num, input_dir, output_dir = args
    # Set OMP_NUM_THREADS to 1 in child process to avoid oversubscription
    os.environ['OMP_NUM_THREADS'] = '1'
    return func(file_prefix, exp_num, input_dir, output_dir)

def main():
    exp_num, input_dir, output_dir = parse_args()
    
    total_start_time = time.time()
    
    print(f"\n{'='*60}")
    print(f"gather_out1d.py - Experiment {exp_num}")
    print(f"{'='*60}")
    print(f"Input directory:  {input_dir}")
    print(f"Output directory: {output_dir}")
    
    # Define file types to process (X-direction concatenation)
    file_types_x = ['stats_t_out','test_t_out', 'ins_islice', 'ins_kslice','field_out', 
                    'test_islice', 'test_kslice']
    
    # Define file types to process (Y-direction concatenation)
    file_types_y = ['ins_jslice', 'test_jslice']
    
    # Prepare tasks for multiprocessing
    tasks = []
    for ft in file_types_x:
        tasks.append((process_file_type, ft, exp_num, input_dir, output_dir))
    for ft in file_types_y:
        tasks.append((process_jslice_file_type, ft, exp_num, input_dir, output_dir))
    
    total_types = len(tasks)
    
    # Determine number of processes
    # Use min(cpu_count, len(tasks)) but cap at reasonable number if needed
    try:
        # We use min(32, ...) because user mentioned OMP_NUM_THREADS=32
        num_procs = min(multiprocessing.cpu_count(), len(tasks), 32)
    except NotImplementedError:
        num_procs = 1
        
    print(f"Running with {num_procs} parallel processes...")
    
    success_count = 0
    failed_types = []
    
    # Run tasks in parallel
    # We use a context manager for the pool
    with multiprocessing.Pool(processes=num_procs) as pool:
        results = pool.map(worker_task, tasks)
        
    # Process results
    for i, success in enumerate(results):
        # tasks[i][1] is the file_prefix
        file_type = tasks[i][1]
        if success:
            success_count += 1
        else:
            failed_types.append(file_type)
    
    # Calculate total elapsed time
    total_elapsed_time = time.time() - total_start_time
    
    # Final summary
    print(f"\n{'='*60}")
    print(f"Summary")
    print(f"{'='*60}")
    print(f"Successfully processed: {success_count}/{total_types} file types")
    if success_count > 0:
        completed = [t[1] for i, t in enumerate(tasks) if results[i]]
        print(f"✓ Completed file types: {completed}")
    if failed_types:
        print(f"✗ Skipped file types: {failed_types}")
    print(f"\n⏱ Total time elapsed: {total_elapsed_time:.2f} seconds ({total_elapsed_time/60:.2f} minutes)")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()
