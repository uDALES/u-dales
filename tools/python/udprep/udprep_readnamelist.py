"""
Fortran namelist parser for uDALES namoptions files.

This module provides a simple parser for Fortran namelist format used in
uDALES namoptions files. It reads key=value pairs from namelists like &RUN,
&DOMAIN, &INPS, etc., and returns them as a Python dictionary.

Functions
---------
read_namoptions : Parse namoptions file and return dictionary of parameters
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict


def read_namoptions(namoptions_path: Path) -> Dict[str, Any]:
    """
    Read namoptions file and return all parameters as a dictionary.
    
    Parses Fortran namelist format and converts values:
    - .true./.false. -> bool
    - Numbers -> int or float
    - Strings -> str (with quotes removed)
    
    Parameters
    ----------
    namoptions_path : Path
        Path to the namoptions file
    
    Returns
    -------
    Dict[str, Any]
        Dictionary mapping parameter names to their values
        
    Raises
    ------
    FileNotFoundError
        If the namoptions file does not exist
        
    Examples
    --------
    >>> from pathlib import Path
    >>> params = read_namoptions(Path('experiments/001/namoptions.001'))
    >>> params['itot']
    128
    >>> params['lzstretch']
    True
    
    Notes
    -----
    This parser handles standard Fortran namelist format:
    - Namelist sections: &RUN, &DOMAIN, &INPS, etc.
    - Comments: ! or lines starting with &
    - Values: integers, floats, booleans (.true./.false.), strings
    - Inline comments are removed
    """
    if not namoptions_path.exists():
        raise FileNotFoundError(f"namoptions file not found: {namoptions_path}")
    
    params = {}
    
    with open(namoptions_path, 'r') as f:
        for line in f:
            line_stripped = line.strip()
            
            # Skip comments, namelist headers/footers, and empty lines
            if (line_stripped.startswith(('&', '!', '/')) or 
                not line_stripped or 
                '=' not in line_stripped):
                continue
            
            # Parse key=value pairs
            parts = line_stripped.split('=', 1)
            key = parts[0].strip()
            val_str = parts[1].split('!')[0].strip()  # Remove inline comments
            
            # Remove trailing comma if present
            val_str = val_str.rstrip(',')
            
            # Convert value types
            val = _parse_value(val_str)
            
            # Store in dictionary
            params[key] = val
    
    return params


def _parse_value(val_str: str) -> Any:
    """
    Parse a Fortran value string into appropriate Python type.
    
    Parameters
    ----------
    val_str : str
        String representation of the value
        
    Returns
    -------
    Any
        Parsed value (bool, int, float, or str)
    """
    val_lower = val_str.lower()
    
    # Boolean values
    if val_lower == '.true.':
        return True
    elif val_lower == '.false.':
        return False
    
    # Numeric values
    try:
        # Float (has decimal point or scientific notation)
        if '.' in val_str or 'e' in val_lower or 'd' in val_lower:
            # Replace Fortran double precision notation
            val_clean = val_str.lower().replace('d', 'e')
            return float(val_clean)
        # Integer
        else:
            return int(val_str)
    except ValueError:
        # String value - remove quotes
        return val_str.strip("'\"")
