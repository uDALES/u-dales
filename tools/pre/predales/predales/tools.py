import os
import numpy as np
import f90nml


def environment():
    """Sets the data paths"""
    HOME = os.environ.get('HOME')
    # on hpc
    if "/rds/general/" in HOME:
        projectpath = HOME + '/dales-urban/'
        exppath = HOME + '/dales-urban/exp'
        workpath = HOME + '/dales-urban/work'
        ephemeral = '/rds/general/user/bss116/ephemeral'
    # on mac
    elif HOME == "/Users/bss116":
        PhD = "/OneDrive - Imperial College London/PhD/PhDproject"
        projectpath = HOME + PhD + '/dales-urban/'
        exppath = HOME + PhD + '/dales-urban/exp'
        workpath = HOME + PhD + '/dales-urban/work/'
        ephemeral = workpath
    else:
        print("Error: HOME undefined.")

    return projectpath, exppath, workpath, ephemeral


def colorschemes():
    light = 'dales-light-notebook'
    dark = 'dales-dark-notebook'
    emphcolors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 
           'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    
    return light, dark, emphcolors


def plotstyles():
    # solid, dotted, dashed, dashdotted, dashdotdotted
    # loosely dotted, dashed, dashdotted, dashdotdotted
    # densely dotted, dashed, dashdotted, dashdotdotted
    linestyles = [(0, ()), (0, (1, 5)), (0, (5, 5)), (0, (3, 5, 1, 5)), (0, (3, 5, 1, 5, 1, 5)),
                  (0, (1, 10)), (0, (5, 10)), (0, (3, 10, 1, 10)), (0, (3, 10, 1, 10, 1, 10)),
                  (0, (1, 1)), (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1))]
    
    return linestyles


def getinputs(datapath, inpspath=None):
    
    expnr = datapath[-4:-1]  # expnr is in last 3 characters before '/' of path

    if inpspath is None:
        # check for namoptions file
        if os.path.isfile(datapath + 'namoptions.' + expnr):
            namoptions = datapath + 'namoptions.' + expnr
#            print("Namoptions file read from %s." % (namoptions))
        else:
            namoptions = ""
            print("Namoptions file not available!")
        # check for blocks file
        if os.path.isfile(datapath + 'blocks.inp.' + expnr):
            blocksfile = datapath + 'blocks.inp.' + expnr
#            print("Blocks read from %s." % (blocksfile))
        else:
            blocksfile = ""
            print("Blocks inp file not available!")
        
    else:
        # check for namoptions file
        if os.path.isfile(datapath + 'namoptions.' + expnr):
            namoptions = datapath + 'namoptions.' + expnr
#            print("Namoptions file read from %s." % (namoptions))
        elif os.path.isfile(inpspath + 'namoptions.' + expnr):
            namoptions = inpspath + 'namoptions.' + expnr
#            print("Namoptions file read from %s." % (namoptions))
        else:
            namoptions = ""
            print("Namoptions file not available!")
        # check for blocks file
        if os.path.isfile(datapath + 'blocks.inp.' + expnr):
            blocksfile = datapath + 'blocks.inp.' + expnr
#            print("Blocks read from %s." % (blocksfile))
        elif os.path.isfile(inpspath + 'blocks.inp.' + expnr):
            blocksfile = inpspath + 'blocks.inp.' + expnr
#            print("Blocks read from %s." % (blocksfile))
        else:
            blocksfile = ""
            print("Blocks inp file not available!")
            
    return namoptions, blocksfile


def getparameters(path):
    """Reads namoptions parameter file from Dales-Urban simulations.
    Parameters saved as python dictionary: e.g. parameter['run']['runtime']"""

    parameters = f90nml.read(path)

    return parameters


def extractparameters(dicts):
    newdict = {}
    # first dict list, e.g. 'RUN', 'DOMAIN', etc.
    for key, val in dicts.items():
        if isinstance(val, dict):
        # second dict list, e.g. 'runtime', 'fieldvars', etc.
            for name, value in val.items():
                newdict[name] = value
        # if normal dict
        else:
            newdict[key] = val
    return newdict


def getblocks(blocksfile):
    
    # get the blocks data
    bl = np.loadtxt(blocksfile, skiprows=2)
    
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


def getdimblocks(blocks, xm, ym, zm):
    
    dimblocks = []
    for block in blocks:
        # only zero block, i.e. no blocks
        if all(b == 0 for b in block):
            dimblocks.append(block)
        else:
            # xm, ym zm at end + 1, index shift from 1 to 0 for x and y
            dimblock = [xm[int(block[0]) - 1], xm[int(block[1])],
                        ym[int(block[2]) - 1], ym[int(block[3])],
                        zm[int(block[4])], zm[int(block[5]) + 1]]
            dimblocks.append(dimblock)
            
    return dimblocks


def blockshift(blocks, limits, xshift, yshift):
    # shift final block layout for periodic boundaries:
    # shift all block elements by intersection (xshift, yshift)

    xmin = limits[0]
    xmax = limits[1]
    ymin = limits[2]
    ymax = limits[3]

    shiftblocks = []
    for block in blocks:
        newx = [(b - xshift) % xmax for b in block[0:2]]
        newy = [(a - yshift) % ymax for a in block[2:4]]
        z = block[4:6]
        shiftblocks.append(newx + newy + z)

    # split up blocks that go over boundaries
    tmpblocks = []
    for block in shiftblocks:
        # if imax < imin
        if block[1] <= block[0]:
            # create new block coordinates
            x1 = [block[0], xmax]
            x2 = [xmin, block[1]]
            y = block[2:4]
            z = block[4:6]
            newblock1 = x1 + y + z
            newblock2 = x2 + y + z
            tmpblocks.extend([newblock1, newblock2])
        else:
            tmpblocks.extend([block])

    newblocks = []
    for block in tmpblocks:
        # if jmax < jmin
        if block[3] <= block[2]:
            # create new block coordinates
            x = block[0:2]
            y1 = [block[2], ymax]
            y2 = [ymin, block[3]]
            z = block[4:6]
            newblock3 = x + y1 + z
            newblock4 = x + y2 + z
            newblocks.extend([newblock3, newblock4])
        else:
            newblocks.extend([block])

    return newblocks
