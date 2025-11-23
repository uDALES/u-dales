#!/usr/bin/env python3
"""
gather_out1d.py - Concatenate x-slab output files using Python/netCDF4

Usage: python gather_out1d.py exp_number [input_dir] [output_dir]

Example: python gather_out1d.py 106 ./ ./

Processes multiple file types:
  X-direction concatenation (xt/xm):
    - test_t_out.xxx.expnr.nc  -> test_t_out.expnr.nc
    - kslice.xxx.expnr.nc      -> kslice.expnr.nc
    - stats_t_out.xxx.expnr.nc -> stats_t_out.expnr.nc
    - islice.xxx.expnr.nc      -> islice.expnr.nc
    - field_out.xxx.expnr.nc   -> field_out.expnr.nc
    - test_islice.xxx.expnr.nc -> test_islice.expnr.nc
    - test_kslice.xxx.expnr.nc -> test_kslice.expnr.nc
  
  Y-direction concatenation (yt/ym):
    - jslice.xxx.expnr.nc      -> jslice.expnr.nc
    - test_jslice.xxx.expnr.nc -> test_jslice.expnr.nc
"""

import os
import sys
import glob
import re
import time
import numpy as np
from netCDF4 import Dataset

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
    if fill_value is not None:
        out_var = out_ds.createVariable(var_name, src_var.datatype, dimensions, 
                                       fill_value=fill_value, zlib=True)
    else:
        out_var = out_ds.createVariable(var_name, src_var.datatype, dimensions, zlib=True)
    
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
    
    print(f"  Found {len(files)} files to concatenate:")
    for f in files:
        y_idx = parse_x_index(f, file_prefix)
        print(f"    {os.path.basename(f)} (y-index: {y_idx})")
    
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
    
    # Open all input files
    try:
        datasets = [Dataset(f, 'r') for f in files]
    except Exception as e:
        print(f"  Error opening input files: {e}")
        print("  Input files:")
        for f in files:
            print(f"    {f} - exists: {os.path.exists(f)}")
        return False
    
    first_ds = datasets[0]
    
    print("  Step 1: Auto-detecting variables and dimensions...")
    
    # Get all variables and separate coordinates from data variables
    coord_vars = {'time', 'xt', 'xm', 'yt', 'ym', 'zt', 'zm'}
    all_vars = set(first_ds.variables.keys())
    data_vars = all_vars - coord_vars
    
    print(f"    Coordinate variables: {sorted(coord_vars & all_vars)}")
    print(f"    Data variables: {sorted(data_vars)}")
    
    # Get dimensions from first file
    dims = {}
    for dim_name in first_ds.dimensions:
        dims[dim_name] = len(first_ds.dimensions[dim_name])
    
    print(f"    Per-file dimensions: {dims}")
    
    # Calculate global y dimensions (concatenate along y)
    global_yt = dims.get('yt', 0) * len(files)
    global_ym = dims.get('ym', 0) * len(files)
    
    print(f"    Global dimensions: yt={global_yt}, ym={global_ym}")
    
    print("  Step 2: Creating output file...")
    
    # Create output dataset
    os.makedirs(output_dir, exist_ok=True)
    out_ds = Dataset(output_file, 'w', format='NETCDF4')
    
    try:
        # Copy global attributes
        for attr in first_ds.ncattrs():
            setattr(out_ds, attr, getattr(first_ds, attr))
        
        # Create dimensions
        print("  Step 3: Creating dimensions...")
        for dim_name, size in dims.items():
            if dim_name == 'yt':
                out_ds.createDimension(dim_name, global_yt)
                print(f"    {dim_name}: {global_yt} (concatenated)")
            elif dim_name == 'ym':
                out_ds.createDimension(dim_name, global_ym)
                print(f"    {dim_name}: {global_ym} (concatenated)")
            elif dim_name == 'time' and first_ds.dimensions[dim_name].isunlimited():
                out_ds.createDimension(dim_name, None)  # Keep unlimited
                print(f"    {dim_name}: {size} (unlimited)")
            else:
                out_ds.createDimension(dim_name, size)
                print(f"    {dim_name}: {size}")
        
        print("  Step 4: Processing coordinate variables...")
        
        # Process coordinate variables
        for var_name in sorted(coord_vars & all_vars):
            print(f"    Processing coordinate: {var_name}")
            
            src_var = first_ds.variables[var_name]
            
            if var_name == 'yt':
                # Concatenate yt coordinates
                out_var = create_variable_with_fillvalue(out_ds, var_name, src_var, src_var.dimensions)
                copy_attributes(src_var, out_var)
                
                # Concatenate data from all files
                offset = 0
                for ds in datasets:
                    data = ds.variables[var_name][:]
                    out_var[offset:offset+len(data)] = data
                    offset += len(data)
                print(f"      Final size: {len(out_var)}")
                    
            elif var_name == 'ym':
                # Concatenate ym coordinates
                out_var = create_variable_with_fillvalue(out_ds, var_name, src_var, src_var.dimensions)
                copy_attributes(src_var, out_var)
                
                # Concatenate data from all files
                offset = 0
                for ds in datasets:
                    data = ds.variables[var_name][:]
                    out_var[offset:offset+len(data)] = data
                    offset += len(data)
                print(f"      Final size: {len(out_var)}")
                    
            else:
                # Copy other coordinates from first file (they should be identical)
                out_var = create_variable_with_fillvalue(out_ds, var_name, src_var, src_var.dimensions)
                copy_attributes(src_var, out_var)
                out_var[:] = src_var[:]
                print(f"      Copied from first file: {out_var.shape}")
        
        print("  Step 5: Processing data variables...")
        
        # Process data variables
        for var_name in sorted(data_vars):
            if var_name not in first_ds.variables:
                continue
                
            src_var = first_ds.variables[var_name]
            var_dims = src_var.dimensions
            
            print(f"    Processing data variable: {var_name} {var_dims}")
            
            # Create output variable with proper fill_value handling
            out_var = create_variable_with_fillvalue(out_ds, var_name, src_var, var_dims)
            copy_attributes(src_var, out_var)
            
            # Determine how to concatenate based on dimensions
            if 'yt' in var_dims:
                print(f"      Concatenating along yt dimension")
                yt_axis = var_dims.index('yt')
                
                # Concatenate along yt axis
                offset = 0
                for i, ds in enumerate(datasets):
                    data = ds.variables[var_name][:]
                    
                    # Build slice for assignment
                    slices = [slice(None)] * len(var_dims)
                    slices[yt_axis] = slice(offset, offset + data.shape[yt_axis])
                    
                    out_var[tuple(slices)] = data
                    offset += data.shape[yt_axis]
                    
                print(f"      Final shape: {out_var.shape}")
                    
            elif 'ym' in var_dims:
                print(f"      Concatenating along ym dimension")
                ym_axis = var_dims.index('ym')
                
                # Concatenate along ym axis
                offset = 0
                for i, ds in enumerate(datasets):
                    data = ds.variables[var_name][:]
                    
                    # Build slice for assignment
                    slices = [slice(None)] * len(var_dims)
                    slices[ym_axis] = slice(offset, offset + data.shape[ym_axis])
                    
                    out_var[tuple(slices)] = data
                    offset += data.shape[ym_axis]
                    
                print(f"      Final shape: {out_var.shape}")
                    
            else:
                print(f"      No y-dimension, copying from first file")
                # Variable doesn't depend on y, just copy from first file
                out_var[:] = src_var[:]
                print(f"      Shape: {out_var.shape}")
        
        print("  Step 6: Finalizing output file...")
        
    except Exception as e:
        print(f"  Error during processing: {e}")
        for ds in datasets:
            ds.close()
        out_ds.close()
        return False
    finally:
        # Close all datasets
        out_ds.close()
        for ds in datasets:
            ds.close()
    
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
        if x_idx is not None:  # Only include files with numeric x-indices
            files.append(f)
    
    if not files:
        print(f"  No files found matching pattern: {pattern}")
        print(f"  Skipping {file_prefix}")
        return False
    
    # Sort files by x-index
    files.sort(key=lambda f: parse_x_index(f, file_prefix))
    
    print(f"  Found {len(files)} files to concatenate:")
    for f in files:
        x_idx = parse_x_index(f, file_prefix)
        print(f"    {os.path.basename(f)} (x-index: {x_idx})")
    
    # Output filename (without .GLOBAL, just prefix.expnr.nc)
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
    
    # Open all input files
    try:
        datasets = [Dataset(f, 'r') for f in files]
    except Exception as e:
        print(f"  Error opening input files: {e}")
        print("  Input files:")
        for f in files:
            print(f"    {f} - exists: {os.path.exists(f)}")
        return False
    
    first_ds = datasets[0]
    
    print("  Step 1: Auto-detecting variables and dimensions...")
    
    # Get all variables and separate coordinates from data variables
    coord_vars = {'time', 'xt', 'xm', 'yt', 'ym', 'zt', 'zm'}
    all_vars = set(first_ds.variables.keys())
    data_vars = all_vars - coord_vars
    
    print(f"    Coordinate variables: {sorted(coord_vars & all_vars)}")
    print(f"    Data variables: {sorted(data_vars)}")
    
    # Get dimensions from first file
    dims = {}
    for dim_name in first_ds.dimensions:
        dims[dim_name] = len(first_ds.dimensions[dim_name])
    
    print(f"    Per-file dimensions: {dims}")
    
    # Calculate global x dimensions
    global_xt = dims.get('xt', 0) * len(files)
    global_xm = dims.get('xm', 0) * len(files)
    
    print(f"    Global dimensions: xt={global_xt}, xm={global_xm}")
    
    print("  Step 2: Creating output file...")
    
    # Create output dataset
    os.makedirs(output_dir, exist_ok=True)
    out_ds = Dataset(output_file, 'w', format='NETCDF4')
    
    try:
        # Copy global attributes
        for attr in first_ds.ncattrs():
            setattr(out_ds, attr, getattr(first_ds, attr))
        
        # Create dimensions
        print("  Step 3: Creating dimensions...")
        for dim_name, size in dims.items():
            if dim_name == 'xt':
                out_ds.createDimension(dim_name, global_xt)
                print(f"    {dim_name}: {global_xt} (concatenated)")
            elif dim_name == 'xm':
                out_ds.createDimension(dim_name, global_xm)
                print(f"    {dim_name}: {global_xm} (concatenated)")
            elif dim_name == 'time' and first_ds.dimensions[dim_name].isunlimited():
                out_ds.createDimension(dim_name, None)  # Keep unlimited
                print(f"    {dim_name}: {size} (unlimited)")
            else:
                out_ds.createDimension(dim_name, size)
                print(f"    {dim_name}: {size}")
        
        print("  Step 4: Processing coordinate variables...")
        
        # Process coordinate variables
        for var_name in sorted(coord_vars & all_vars):
            print(f"    Processing coordinate: {var_name}")
            
            src_var = first_ds.variables[var_name]
            
            if var_name == 'xt':
                # Concatenate xt coordinates
                out_var = create_variable_with_fillvalue(out_ds, var_name, src_var, src_var.dimensions)
                copy_attributes(src_var, out_var)
                
                # Concatenate data from all files
                offset = 0
                for ds in datasets:
                    data = ds.variables[var_name][:]
                    out_var[offset:offset+len(data)] = data
                    offset += len(data)
                print(f"      Final size: {len(out_var)}")
                    
            elif var_name == 'xm':
                # Concatenate xm coordinates
                out_var = create_variable_with_fillvalue(out_ds, var_name, src_var, src_var.dimensions)
                copy_attributes(src_var, out_var)
                
                # Concatenate data from all files
                offset = 0
                for ds in datasets:
                    data = ds.variables[var_name][:]
                    out_var[offset:offset+len(data)] = data
                    offset += len(data)
                print(f"      Final size: {len(out_var)}")
                    
            else:
                # Copy other coordinates from first file (they should be identical)
                out_var = create_variable_with_fillvalue(out_ds, var_name, src_var, src_var.dimensions)
                copy_attributes(src_var, out_var)
                out_var[:] = src_var[:]
                print(f"      Copied from first file: {out_var.shape}")
        
        print("  Step 5: Processing data variables...")
        
        # Process data variables
        for var_name in sorted(data_vars):
            if var_name not in first_ds.variables:
                continue
                
            src_var = first_ds.variables[var_name]
            var_dims = src_var.dimensions
            
            print(f"    Processing data variable: {var_name} {var_dims}")
            
            # Create output variable with proper fill_value handling
            out_var = create_variable_with_fillvalue(out_ds, var_name, src_var, var_dims)
            copy_attributes(src_var, out_var)
            
            # Determine how to concatenate based on dimensions
            if 'xt' in var_dims:
                print(f"      Concatenating along xt dimension")
                xt_axis = var_dims.index('xt')
                
                # Concatenate along xt axis
                offset = 0
                for i, ds in enumerate(datasets):
                    data = ds.variables[var_name][:]
                    
                    # Build slice for assignment
                    slices = [slice(None)] * len(var_dims)
                    slices[xt_axis] = slice(offset, offset + data.shape[xt_axis])
                    
                    out_var[tuple(slices)] = data
                    offset += data.shape[xt_axis]
                    
                print(f"      Final shape: {out_var.shape}")
                    
            elif 'xm' in var_dims:
                print(f"      Concatenating along xm dimension")
                xm_axis = var_dims.index('xm')
                
                # Concatenate along xm axis
                offset = 0
                for i, ds in enumerate(datasets):
                    data = ds.variables[var_name][:]
                    
                    # Build slice for assignment
                    slices = [slice(None)] * len(var_dims)
                    slices[xm_axis] = slice(offset, offset + data.shape[xm_axis])
                    
                    out_var[tuple(slices)] = data
                    offset += data.shape[xm_axis]
                    
                print(f"      Final shape: {out_var.shape}")
                    
            else:
                print(f"      No x-dimension, copying from first file")
                # Variable doesn't depend on x, just copy from first file
                out_var[:] = src_var[:]
                print(f"      Shape: {out_var.shape}")
        
        print("  Step 6: Finalizing output file...")
        
    except Exception as e:
        print(f"  Error during processing: {e}")
        for ds in datasets:
            ds.close()
        out_ds.close()
        return False
    finally:
        # Close all datasets
        out_ds.close()
        for ds in datasets:
            ds.close()
    
    elapsed_time = time.time() - start_time
    print(f"  ✓ Concatenation complete: {output_file}")
    print(f"  ⏱ Time elapsed: {elapsed_time:.2f} seconds")
    return True


def main():
    exp_num, input_dir, output_dir = parse_args()
    
    total_start_time = time.time()
    
    print(f"\n{'='*60}")
    print(f"gather_out1d.py - Experiment {exp_num}")
    print(f"{'='*60}")
    print(f"Input directory:  {input_dir}")
    print(f"Output directory: {output_dir}")
    
    # Define file types to process (X-direction concatenation)
    file_types_x = ['test_t_out', 'stats_t_out', 'islice', 'kslice','field_out', 
                    'test_islice', 'test_kslice']
    
    # Define file types to process (Y-direction concatenation)
    file_types_y = ['jslice', 'test_jslice']
    
    # Track successful processing
    success_count = 0
    failed_types = []
    total_types = len(file_types_x) + len(file_types_y)
    
    # Process X-direction file types
    for file_type in file_types_x:
        if process_file_type(file_type, exp_num, input_dir, output_dir):
            success_count += 1
        else:
            failed_types.append(file_type)
    
    # Process Y-direction file types (jslice)
    for file_type in file_types_y:
        if process_jslice_file_type(file_type, exp_num, input_dir, output_dir):
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
        completed = [ft for ft in (file_types_x + file_types_y) if ft not in failed_types]
        print(f"✓ Completed file types: {completed}")
    if failed_types:
        print(f"✗ Skipped file types: {failed_types}")
    print(f"\n⏱ Total time elapsed: {total_elapsed_time:.2f} seconds ({total_elapsed_time/60:.2f} minutes)")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()
