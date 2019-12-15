def getdimblocks(blocks, xm, ym, zm):
    """Convert uDALES blocks input file and mid cell
    coordinate axes to dimensionalised blocks coordinates."""
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
