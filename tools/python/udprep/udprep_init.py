"""
Helper functions for reading configuration and namoptions files.

This module provides utilities to parse shell configuration files and
Fortran namoptions files for initializing UDPrep preprocessing.
"""

from __future__ import annotations

from pathlib import Path
from typing import Tuple
import subprocess
import re
import sys
import warnings


def parse_shell_config(config_path: Path) -> dict:
    """
    Parse shell config file and extract variable assignments.
    
    Parameters
    ----------
    config_path : Path
        Path to the shell configuration file (e.g., config.sh)
    
    Returns
    -------
    dict
        Dictionary of DA_* variables and their values
    """
    variables = {}
    
    # Execute shell script to handle command substitutions like $(pwd)
    try:
        result = subprocess.run(
            ["bash", "-c", f"source '{config_path}' && env"],
            capture_output=True,
            text=True,
            check=True,
            cwd=config_path.parent,
        )
        
        # Extract DA_* variables from environment
        for line in result.stdout.splitlines():
            if "=" in line:
                key, _, value = line.partition("=")
                if key.startswith("DA_"):
                    variables[key] = value
        
        return variables
        
    except subprocess.CalledProcessError:
        # Fallback: simple text parsing (won't handle command substitutions)
        with config_path.open("r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                    
                # Remove 'export ' prefix if present
                line = line.replace("export ", "", 1).strip()
                key, _, value = line.partition("=")
                key = key.strip()
                value = value.strip().strip("\"'")  # Remove quotes
                
                if key.startswith("DA_"):
                    variables[key] = value
        
        return variables


def read_iexpnr_from_namoptions(namoptions_path: Path) -> str:
    """
    Read and validate iexpnr value from namoptions file.
    
    Parameters
    ----------
    namoptions_path : Path
        Path to the namoptions file
    
    Returns
    -------
    str
        Experiment number (3-digit format: 001, 056, 999, etc.)
    """
    with namoptions_path.open("r") as f:
        for line in f:
            line = line.strip()
            # Skip comments, namelist headers, and empty lines
            if line.startswith(('&', '!')) or not line:
                continue
                
            # Look for iexpnr = value
            if "iexpnr" in line.lower() and "=" in line:
                # Extract value after =, removing inline comments and commas
                value = line.split('=', 1)[1].split('!')[0].strip().rstrip(',')
                
                # Validate 3-digit format
                if not re.match(r'^\d{3}$', value):
                    print(f"ERROR: iexpnr must be 3-digit format (001, 056, 999)", file=sys.stderr)
                    print(f"Found: '{value}' in {namoptions_path}", file=sys.stderr)
                    sys.exit(1)
                
                return value
    
    print(f"ERROR: iexpnr not found in {namoptions_path}", file=sys.stderr)
    sys.exit(1)


def set_expnr(expdir: Path) -> str:
    """
    Extract and validate experiment number from directory structure.
    
    This is the essential validation that ensures consistency between:
    - Directory name (e.g., "999")
    - Namoptions filename suffix (e.g., "namoptions.999")
    - iexpnr value inside namoptions file
    
    Parameters
    ----------
    expdir : Path
        Experiment directory containing namoptions.XXX
    
    Returns
    -------
    str
        3-digit experiment number
        
    Raises
    ------
    SystemExit
        If directory doesn't exist, namoptions not found, or inconsistent numbering
    """
    # Check experiment directory exists
    if not expdir.exists():
        print(f"ERROR: Experiment directory does not exist: {expdir}", file=sys.stderr)
        sys.exit(1)
    
    # Find and validate namoptions file
    namoptions_files = list(expdir.glob("namoptions.*"))
    
    if not namoptions_files:
        print(f"ERROR: No namoptions file in {expdir}", file=sys.stderr)
        print(f"Expected: namoptions.XXX (e.g., namoptions.001)", file=sys.stderr)
        sys.exit(1)
    
    # Validate filename pattern (namoptions.XXX where XXX is 3 digits)
    valid_files = [(f, f.name.split('.', 1)[1]) for f in namoptions_files 
                   if '.' in f.name and re.match(r'^\d{3}$', f.name.split('.', 1)[1])]
    
    if not valid_files:
        found = ', '.join(f.name for f in namoptions_files)
        print(f"ERROR: Invalid namoptions filename in {expdir}", file=sys.stderr)
        print(f"Found: {found}", file=sys.stderr)
        print(f"Expected: namoptions.XXX (e.g., namoptions.001)", file=sys.stderr)
        sys.exit(1)
      
    if len(valid_files) > 1:
        found = ', '.join(f[0].name for f in valid_files)
        print(f"ERROR: Multiple namoptions files: {found}", file=sys.stderr)
        sys.exit(1)
    
    namoptions_file, namoptions_suffix = valid_files[0]
    
    # Read and validate iexpnr
    expnr = read_iexpnr_from_namoptions(namoptions_file)
    
    # Validate consistency: directory name, filename suffix, and iexpnr must all match
    if not (namoptions_suffix == expdir.name == expnr):
        print(f"ERROR: Experiment number mismatch", file=sys.stderr)
        print(f"  Directory name: {expdir.name}", file=sys.stderr)
        print(f"  Filename: {namoptions_file.name} → {namoptions_suffix}", file=sys.stderr)
        print(f"  iexpnr in namoptions: {expnr}", file=sys.stderr)
        print(f"  All three must match. Set iexpnr under &RUN in namoptions to match the experiment case directory name and namoptions filename suffix.", file=sys.stderr)
        sys.exit(1)
    
    # Optional: validate config.sh paths (warnings only, non-blocking)
    validate_config_paths(expdir)
    
    return expnr


def validate_config_paths(expdir: Path) -> None:
    """
    Validate config.sh paths (optional validation).
    
    Checks if config.sh exists and validates that DA_EXPDIR and DA_TOOLSDIR
    point to the expected locations. This is optional for preprocessing but
    required for execution scripts (archer_execute.sh, etc.).
    
    Parameters
    ----------
    expdir : Path
        Experiment directory that should contain config.sh
        
    Notes
    -----
    If config.sh is missing, only a warning is printed.
    If config.sh exists but has incorrect paths, warnings are printed.
    """
    config_file = expdir / "config.sh"
    
    # Check config.sh exists
    if not config_file.exists():
        print("="*67, file=sys.stderr)
        print(f"WARNING: Config file not found: {config_file}", file=sys.stderr)
        print(f"  config.sh is optional for preprocessing but required for execution.", file=sys.stderr)
        print(f"  Expected variables: DA_EXPDIR, DA_TOOLSDIR, DA_WORKDIR, DA_BUILD", file=sys.stderr)
        print("="*67, file=sys.stderr)
        return
    
    # Parse config.sh
    print(f"Reading: {config_file}")
    config = parse_shell_config(config_file)
    
    # Check required variables
    missing = [var for var in ["DA_EXPDIR", "DA_TOOLSDIR"] if var not in config]
    if missing:
        print("="*67, file=sys.stderr)
        print(f"WARNING: Missing variables in config.sh: {', '.join(missing)}", file=sys.stderr)
        print(f"  These are required for execution scripts but optional for preprocessing.", file=sys.stderr)
        print("="*67, file=sys.stderr)
        return
    
    # Validate DA_EXPDIR matches the experiment directory
    config_expdir = Path(config["DA_EXPDIR"]).resolve()
    expected_expdir = expdir.parent.resolve()
    
    if config_expdir != expected_expdir:
        print("="*67, file=sys.stderr)
        print(f"WARNING: DA_EXPDIR mismatch in {config_file}", file=sys.stderr)
        print(f"  Expected: {expected_expdir}", file=sys.stderr)
        print(f"  Current: {config_expdir}", file=sys.stderr)
        print(f"  DA_EXPDIR should point to the experiments directory.", file=sys.stderr)
        print("="*67, file=sys.stderr)
    
    # Validate DA_TOOLSDIR matches script location
    expected_tools = Path(__file__).resolve().parent.parent
    config_tools = (Path(config["DA_TOOLSDIR"]) / "python").resolve()
    
    if expected_tools != config_tools:
        print("="*67, file=sys.stderr)
        print(f"WARNING: DA_TOOLSDIR mismatch in {config_file}", file=sys.stderr)
        print(f"  Expected: {expected_tools}", file=sys.stderr)
        print(f"  Current: {config_tools}", file=sys.stderr)
        print(f"  DA_TOOLSDIR should point to u-dales/tools directory.", file=sys.stderr)
        print("="*67, file=sys.stderr)



