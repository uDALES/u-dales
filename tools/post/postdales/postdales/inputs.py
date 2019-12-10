import f90nml
import numpy as np

def parameters(path):
    """Reads namoptions parameter file from uDales simulations."""
    parameters = f90nml.read(path)

    return parameters


def extract_f90namelist(f90namelist):
    """Extracts nested dictionaries from namelist and saves parameters to a single dictionary."""
    newdict = {}
    # first dict list, e.g. 'RUN', 'DOMAIN', etc.
    for key, val in f90namelist.items():
        # second dict list, e.g. 'runtime'
        for name, value in val.items():
            newdict[name] = value
    return newdict


def blocks(path):
    """Reads uDALES blocks configuration."""    
    # get the blocks data
    bl = np.loadtxt(path, skiprows=2)
    
    if np.size(bl) == 11:  # when only one line of blocks
        tmpblocks = [list(bl[0:6].astype(int))]
    else:
        tmpblocks = [list(b[0:6].astype(int)) for b in bl]
    
    blocks = []
    # remove blocks that are only flat surfaces
    for block in tmpblocks:
        # blocks with height 0
        if not block[4] == block[5] == 0:
            blocks.append(block)
    
    return blocks
