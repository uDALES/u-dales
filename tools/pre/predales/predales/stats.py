import math
import numpy as np


def area(x, y):
    """Calculates the area of two given intervals x = [xmin xmax] and y = [ymin ymax]"""
    a = (x[1] - x[0]) * (y[1] - y[0])
    return a


def blockplan(blocks):
    """Returns sum of all block plan areas.
    Block entries stored as block = [xmin, xmax, ymin, ymax]"""
    areas = 0
    for block in blocks:
        x = block[0:2]
        y = block[2:4]
        areas += area(x, y)
    return areas


def blockfront(blocks):
    """Returns sum of all block frontal areas facing the wind.
    Wind is U wind only, the front is the y-z surface of blocks.
    Block entries stored as block = [xmin, xmax, ymin, ymax, zmin, zmax]"""
    areas = 0
    for block in blocks:
        y = block[2:4]
        z = block[4:6]
        areas += area(y, z)
    return areas

def blockareas(blocks, z):
    """Returns a function of height that gives the area occupied by blocks.
    Block entries stored as block = [xmin, xmax, ymin, ymax, zmin, zmax]"""
    blockcover = np.zeros(len(z))
    for block in blocks:
        x = block[0:2]
        y = block[2:4]
        height = block[5]
        # should we only start from block[4],
        # since the first layer is street surface?
        for i, zi in enumerate(z):
            if zi <= height:
                blockcover[i] += area(x, y)
            
    return blockcover


def blockfronts(blocks, z):
    """Returns a function of height that gives the frontal area occupied by blocks. 
    It returns the widths of the block integrated for the height where the block is present.
    Block entries stored as block = [xmin, xmax, ymin, ymax, zmin, zmax]"""
    # we need to integrate up to z_UCL, not z_max, because at z_max we still want a value > 0.
    blockcover = np.zeros(len(z))
    for block in blocks:
        blockwidth = block[3] - block[2]
        frontalarea = block[5] - block[4]
        blockheight = block[5]
        lfblock = blockwidth*frontalarea
        zmax_index = np.argmin([abs(zi - blockheight) for zi in z])
        zucl = z[zmax_index + 1]
        for i, zi in enumerate(z):
            if zi <= zucl:
                blockcover[i] += lfblock*(zucl - zi)/blockheight
#         blockcover[0] -= blockwidth*(blockheight - frontalarea)  # correction for street surface?
            
    return blockcover


def blockfronts_discrete(blocks, z):
    """Returns a function of height that gives the frontal area occupied by blocks. It returns full frontal area (using full height of block) for       any height where the block is present.
    Block entries stored as block = [xmin, xmax, ymin, ymax, zmin, zmax]"""
    blockcover = np.zeros(len(z))
    for block in blocks:
        y = block[2:4]
        zs = block[4:6]
        height = block[5]
        # should we only start from block[4],
        # since the first layer is street surface?
        for i, zt in enumerate(z):
            if zt <= height:
                blockcover[i] += area(y, zs)
            
    return blockcover


def blockmask(blocks, z):
    """Returns a function of height that indicates whether there are any blocks at that height.
    Block entries stored as block = [xmin, xmax, ymin, ymax, zmin, zmax]"""
    blockcover = np.zeros(len(z))
    for block in blocks:
        x = block[0:2]
        y = block[2:4]
        height = block[5]
        # should we only start from block[4],
        # since the first layer is street surface?
        for i, zi in enumerate(z):
            if zi <= height:
                blockcover[i] = 1
            
    return blockcover


def blockstats(dimblocks, a0=None):
    blockstatistics = {}
    precision = 3

    # number of blocks
    nblocks = len(dimblocks)
    blockstatistics['nblocks'] = nblocks

    # blocks plan area
    afinal = blockplan(dimblocks)
    blockstatistics['blocksarea'] = afinal

    # blocks frontal area
    hfinal = blockfront(dimblocks)
    blockstatistics['blocksarea_f'] = hfinal

    # blockheights
    blockheights = [b[5] for b in dimblocks]
    blockstatistics['blockheights'] = blockheights

    # block maximum zmax
    if blockheights:  # if not empty
        zmax = np.max(blockheights)
    else:
        zmax = 0
    blockstatistics['zmax'] = round(zmax, precision)

    # block average height unweighted
    if blockheights:
        zmean = np.mean(blockheights)
    else:
        zmean = 0
    blockstatistics['zmean'] = round(zmean, precision)

    # blockfronts
    blockfronts = [blockfront([f]) for f in dimblocks]
    blockstatistics['blockfronts'] = blockfronts

    # blockplans
    blockplans = [blockplan([p]) for p in dimblocks]
    blockstatistics['blockplans'] = blockplans

    # block weighted average height zh
    if hfinal > 0:
        # weights from blockfronts
        weights = [b/hfinal for b in blockfronts]
        weightedheights = [f * b for f, b in zip(weights, blockheights)]
        zh = np.sum(weightedheights)
    else: # sets weighted heights to zero in case of no blocks
        zh = 0
    blockstatistics['zh'] = round(zh, precision)

    # block weighted average height zh_alt weighted by block area
    # not block fronts
    if afinal > 0:
        # weights from blockplan areas
        weights2 = [b/afinal for b in blockplans]
        weightedheights2 = [f * b for f, b in zip(weights2, blockheights)]
        zh_alt = np.sum(weightedheights2)
    else: # sets weighted heights to zero in case of no blocks
        zh_alt = 0
    blockstatistics['zh_alt'] = round(zh_alt, precision)

    # block height standard deviation unweighted
    if blockheights:
        sigmamean = np.sqrt(np.mean((blockheights-zmean)**2))
    else:
        sigmamean = 0
    blockstatistics['sigma_mean'] = round(sigmamean, precision)

    # block height standard deviation weighted
    if hfinal > 0:
        sigmah = np.sqrt(np.sum(weights*((blockheights-zh)**2)))
    else:
        sigmah = 0
    blockstatistics['sigma_h'] = round(sigmah, precision)

    # if a0 defined
    try:
        # building area density
        lp = afinal / a0
        blockstatistics['lp'] = round(lp, precision)
        
        # building frontal aspect ratio
        lf = hfinal / a0
        blockstatistics['lf'] = round(lf, precision)
        
    except Exception:
        print("No a0 defined.")
        pass
        
    return blockstatistics
