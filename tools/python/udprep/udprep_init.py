"""
Helper functions for reading configuration and namoptions files.

This module provides utilities to parse shell configuration files and
Fortran namoptions files for initializing UDPrep preprocessing.
"""

from __future__ import annotations

from pathlib import Path
from typing import Tuple
import subprocess
import platform
import re
import sys


def parse_shell_config(config_path: Path) -> dict:
    """
    Parse shell config file and extract variable assignments.
    
    This function sources the shell script and extracts exported variables
    to properly handle shell command substitutions like $(pwd).
    
    Parameters
    ----------
    config_path : Path
        Path to the shell configuration file (e.g., config.sh)
    
    Returns
    -------
    dict
        Dictionary mapping variable names to their values
    
    Raises
    ------
    RuntimeError
        If the config file cannot be parsed
    
    Examples
    --------
    >>> config = parse_shell_config(Path("config.sh"))
    >>> expdir = config["DA_EXPDIR"]
    """
    variables = {}
    
    # Try to execute the shell script and extract environment variables
    # This handles shell command substitutions like $(pwd)
    try:
        # Detect if we're on Windows or Unix-like system
        is_windows = platform.system() == "Windows"
        
        if is_windows:
            # On Windows, try to use bash (WSL or Git Bash)
            shell_cmd = ["bash", "-c", f"source '{config_path}' && env"]
        else:
            # On Unix-like systems, use bash directly
            shell_cmd = ["bash", "-c", f"source '{config_path}' && env"]
        
        result = subprocess.run(
            shell_cmd,
            capture_output=True,
            text=True,
            check=True,
            cwd=config_path.parent,
        )
        
        # Parse the environment variables output
        for line in result.stdout.splitlines():
            if "=" in line:
                key, _, value = line.partition("=")
                # Only keep variables that start with DA_ or other relevant prefixes
                if key.startswith("DA_") or key.startswith("NCPU"):
                    variables[key] = value
        
        return variables
        
    except subprocess.CalledProcessError as e:
        # Fall back to simple text parsing if shell execution fails
        try:
            with config_path.open("r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    # Skip comments and empty lines
                    if not line or line.startswith("#"):
                        continue
                    # Handle export VAR=value or VAR=value
                    if "=" in line:
                        # Remove 'export ' prefix if present
                        line = line.replace("export ", "", 1).strip()
                        # Split on first =
                        key, _, value = line.partition("=")
                        key = key.strip()
                        value = value.strip()
                        # Remove quotes if present
                        if value.startswith('"') and value.endswith('"'):
                            value = value[1:-1]
                        elif value.startswith("'") and value.endswith("'"):
                            value = value[1:-1]
                        variables[key] = value
        except Exception as parse_error:
            raise RuntimeError(
                f"Failed to parse config file {config_path}. "
                f"Shell execution failed: {e}\n"
                f"Text parsing also failed: {parse_error}"
            )
    
    return variables


def read_iexpnr_from_namoptions(namoptions_path: Path) -> str:
    """
    Read iexpnr value from namoptions file.
    
    Parameters
    ----------
    namoptions_path : Path
        Path to the namoptions file
    
    Returns
    -------
    str
        Experiment number as written in the file (must be 3-digit format)
    
    Raises
    ------
    RuntimeError
        If the file cannot be read
    ValueError
        If iexpnr is not found or not in 3-digit format
    
    Examples
    --------
    >>> expnr = read_iexpnr_from_namoptions(Path("namoptions"))
    >>> print(expnr)  # e.g., "105"
    """
    try:
        with namoptions_path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                # Skip comments and namelist headers
                if line.startswith('&') or line.startswith('!') or not line:
                    continue
                # Look for iexpnr = value
                if "iexpnr" in line.lower() and "=" in line:
                    # Extract value after =
                    parts = line.split('=', 1)
                    value = parts[1].split('!')[0].strip()  # Remove inline comments
                    # Remove comma if present
                    value = value.rstrip(',').strip()
                    
                    # Validate that it's a 3-digit number string
                    if not re.match(r'^\d{3}$', value):
                        raise ValueError(
                            f"iexpnr value in {namoptions_path} must be a 3-digit number.\n"
                            f"Found: '{value}'\n"
                            f"Expected format: 001, 056, 105, 999, etc.\n"
                            f"Please update iexpnr in the namoptions file to use 3 digits with leading zeros."
                        )
                    
                    return value
    except ValueError:
        # Re-raise ValueError as-is (from our validation)
        raise
    except Exception as e:
        raise RuntimeError(f"Failed to read iexpnr from {namoptions_path}: {e}")
    raise ValueError(f"iexpnr not found in {namoptions_path}")


def setup_paths_from_config(config_dir: Path) -> Tuple[Path, str]:
    """
    Set up paths from a configuration directory.
    
    Reads config.sh and namoptions files from the given directory and
    extracts experiment directory and number. Validates that DA_TOOLSDIR
    matches the actual tools directory. Handles all errors internally and
    exits with helpful messages if configuration is invalid.
    
    Parameters
    ----------
    config_dir : Path
        Directory containing config.sh and namoptions files
    
    Returns
    -------
    tuple of (Path, str)
        - expdir: Path to experiment directory
        - expnr: Experiment number as zero-padded 3-digit string
    
    Examples
    --------
    >>> expdir, expnr = setup_paths_from_config(Path("/path/to/config"))
    """
    config_file = config_dir / "config.sh"
    
    # Handle missing config.sh
    if not config_file.exists():
        print(f"ERROR: Config file not found: {config_file}", file=sys.stderr)
        print(f"\nExpected files in {config_dir}:", file=sys.stderr)
        print(f"  - config.sh (containing DA_EXPDIR and DA_TOOLSDIR)", file=sys.stderr)
        print(f"  - namoptions.XXX (where XXX is a 3-digit experiment number)", file=sys.stderr)
        sys.exit(1)
    
    # Try to find namoptions file - must be in format namoptions.XXX where XXX is 3-digit integer
    namoptions_files = list(config_dir.glob("namoptions.*"))
    
    if not namoptions_files:
        print(f"ERROR: Namoptions file not found in {config_dir}", file=sys.stderr)
        print(f"\nExpected file matching pattern 'namoptions.XXX' where XXX is a 3-digit experiment number.", file=sys.stderr)
        print(f"\nExample: namoptions.001, namoptions.105, namoptions.999", file=sys.stderr)
        sys.exit(1)
    
    # Validate namoptions filename pattern
    valid_namoptions = []
    for nf in namoptions_files:
        # Extract the suffix after 'namoptions.'
        suffix = nf.name.split('.', 1)[1] if '.' in nf.name else None
        if suffix and re.match(r'^\d{3}$', suffix):
            valid_namoptions.append((nf, suffix))
    
    if not valid_namoptions:
        found_files = ', '.join([nf.name for nf in namoptions_files])
        print(f"ERROR: No valid namoptions file found in {config_dir}", file=sys.stderr)
        print(f"Found: {found_files}", file=sys.stderr)
        print(f"\nExpected format: 'namoptions.XXX' where XXX is a 3-digit number.", file=sys.stderr)
        print(f"Examples: namoptions.001, namoptions.105, namoptions.999", file=sys.stderr)
        sys.exit(1)
      
    if len(valid_namoptions) > 1:
        found_files = ', '.join([nf[0].name for nf in valid_namoptions])
        print(f"ERROR: Multiple namoptions files found in {config_dir}: {found_files}", file=sys.stderr)
        print(f"\nPlease keep only one namoptions file.", file=sys.stderr)
        sys.exit(1)
    
    namoptions_file, namoptions_suffix = valid_namoptions[0]
    
    # Read iexpnr from namoptions file
    try:
        expnr = read_iexpnr_from_namoptions(namoptions_file)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Unexpected error reading namoptions file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Validate consistency: namoptions filename suffix must match iexpnr
    if namoptions_suffix != expnr:
        print(f"ERROR: Inconsistent experiment numbers", file=sys.stderr)
        print(f"  - Namoptions filename: {namoptions_file.name} (implies expnr={namoptions_suffix})", file=sys.stderr)
        print(f"  - iexpnr in namoptions: {expnr}", file=sys.stderr)
        print(f"\nThese must match. Please rename the file or update iexpnr value.", file=sys.stderr)
        sys.exit(1)
    
    # Validate consistency: directory name must match expnr
    dir_name = config_dir.name
    if dir_name != expnr:
        print(f"ERROR: Inconsistent experiment numbers", file=sys.stderr)
        print(f"  - Directory name: {dir_name}", file=sys.stderr)
        print(f"  - Experiment number (from namoptions): {expnr}", file=sys.stderr)
        print(f"\nThe config directory name must match the experiment number.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Reading configuration from: {config_file}")
    try:
        config = parse_shell_config(config_file)
    except Exception as e:
        print(f"ERROR: Failed to parse config.sh: {e}", file=sys.stderr)
        sys.exit(1)
    
    if "DA_EXPDIR" not in config:
        print(f"ERROR: DA_EXPDIR not found in config.sh", file=sys.stderr)
        print(f"\nPlease ensure config.sh contains:", file=sys.stderr)
        print(f"  DA_EXPDIR=<path to experiment directory>", file=sys.stderr)
        print(f"  DA_TOOLSDIR=<path to tools directory>", file=sys.stderr)
        sys.exit(1)
    if "DA_TOOLSDIR" not in config:
        print(f"ERROR: DA_TOOLSDIR not found in config.sh", file=sys.stderr)
        print(f"\nPlease ensure config.sh contains:", file=sys.stderr)
        print(f"  DA_EXPDIR=<path to experiment directory>", file=sys.stderr)
        print(f"  DA_TOOLSDIR=<path to tools directory>", file=sys.stderr)
        sys.exit(1)
    
    # DA_EXPDIR might be the parent experiments directory or the specific experiment directory
    # Ensure we use the specific experiment directory
    expdir_raw = Path(config["DA_EXPDIR"]).resolve()
    if expdir_raw.name == expnr:
        # Already pointing to the specific experiment directory
        expdir = expdir_raw
    else:
        # Pointing to parent experiments directory, append experiment number
        expdir = expdir_raw / expnr
    
    # Validate that DA_TOOLSDIR is consistent with the actual tools directory
    # The script should be in tools/python, so we can derive the expected tools path
    expected_tools_path = (Path(__file__).resolve().parent.parent).resolve()
    config_tools_path = (Path(config["DA_TOOLSDIR"]) / "python").resolve()
    
    if expected_tools_path != config_tools_path:
        import warnings
        warnings.warn(
            f"DA_TOOLSDIR mismatch:\n"
            f"  Expected (from script location): {expected_tools_path}\n"
            f"  From config.sh: {config_tools_path}\n"
            f"  Using script location.",
            UserWarning
        )
    
    print(f"Using paths from config:")
    print(f"  expdir: {expdir}")
    print(f"  expnr: {expnr}")
    
    return expdir, expnr
