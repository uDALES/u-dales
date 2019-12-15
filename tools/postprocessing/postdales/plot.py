import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import copy
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def vertical_profile(var, zaxis, timestep=-1):
    """Plot the vertical profile of a 1D-spatial
    netcdf variable."""
    plt.figure()
    plt.plot(var[timestep, :], zaxis[:])
    plt.ylabel(zaxis.name)
    plt.xlabel(var.name)
    plt.title(var.longname)
    plt.show()


def field_plot2D(var, xaxis, yaxis, timestep=-1, blocks=None):
    """Plot the contour field of a 2D-spatial netcdf
    variable."""
    cbarmap = plt.cm.get_cmap('RdYlBu_r')
    clevels = contour_levels(var[timestep, :, :])
    
    plt.figure()
    cs = plt.contourf(xaxis[:], yaxis[:], var[timestep, :, :], cmap=cbarmap, levels=clevels)

    if blocks is not None:
        add_2dblocks(blocks)
    plt.axis('scaled')
    plt.xlabel(xaxis.name)
    plt.ylabel(yaxis.name)
    plt.title(var.longname)
    cbar = plt.colorbar(cs)
    plt.tight_layout()
    plt.show()


def field_plot3D(var, xaxis, yaxis, timestep=-1, spacepos=0, blocks=None):
    """Plot the contour field of a 3D-spatial netcdf
    variable."""
    dataslice = slices_from_3D(var, xaxis, yaxis, timestep, spacepos)
    cbarmap = plt.cm.get_cmap('RdYlBu_r')
    clevels = contour_levels(var[dataslice])
    
    plt.figure()
    cs = plt.contourf(xaxis[:], yaxis[:], var[dataslice], cmap=cbarmap, levels=clevels)
    if blocks is not None:
        add_2dblocks(blocks)
    plt.axis('scaled')
    plt.xlabel(xaxis.name)
    plt.ylabel(yaxis.name)
    plt.title(var.longname)
    cbar = plt.colorbar(cs)
    plt.tight_layout()
    plt.show()


def add_2dblocks(blocks):
    """Function that adds blocks to 2D contour plot"""
    ax = plt.gca()
    for block in blocks:
        start = tuple([block[0], block[2]])  # xmin, ymin
        width = block[1] - block[0]  # xmax - xmin
        height = block[3] - block[2]  # ymax - ymin
        p = patches.Rectangle(start, width, height, 
        edgecolor='k', facecolor='xkcd:light grey')
        ax.add_patch(p)

    return


def blocks3d(blocks, limits=None):

    plt.figure()
    ax = plt.axes(projection='3d')
        
    fas = faces(blocks)

    for sides in fas:
        collection = Poly3DCollection(sides, edgecolor='k', facecolor='xkcd:light grey')
        ax.add_collection3d(collection)
    
    if limits is not None:
        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])
        ax.set_zlim(limits[4], limits[5])

    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.05))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.05))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.05))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_aspect('equal')
    ax.dist=11

    plt.show()


# SUPPORT FUNCTIONS # ------------
def contour_levels(fielddata):
    """Increase to maximum number of levels for 
    contour plot."""
    vmax = np.nanmax(fielddata)
    vmin = np.nanmin(fielddata)
    levels = np.linspace(vmin, vmax, 256)
    return levels


def slices_from_3D(var, xaxis, yaxis, timepos, spacepos):
    """Define the data slice of 3-D spatial netcdf variable
    for 2-D contour plot."""
    full = slice(0, None)
    if 'x' in xaxis.name and 'z' in yaxis.name:
        dataslice = [timepos, full, spacepos, full]
    elif 'x' in xaxis.name and 'y' in yaxis.name:
        dataslice = [timepos, spacepos, full, full]
    elif 'y' in xaxis.name and 'z' in yaxis.name:
        dataslice = [timepos, full, full, spacepos]
    else:
        raise ValueError('Coordinate axes must be one of the following: combinations: (x, y), (x, z), (y, z).')
    return dataslice


def get_2dblocks(blocks, xaxis, yaxis, position):
    """Get the blocks at spatial position."""    
    # deepcopy to not manipulate original blocks
    inblocks = copy.deepcopy(blocks)
    blocks2d = []
    if 'x' in xaxis.name and 'z' in yaxis.name:
        for block in inblocks:
            if block[2] <= position <= block[3]:
                # imin, imax, kmin, kmax
                newblock = block[0:2] + block[4:6]
                blocks2d.append(newblock)
                    
    elif 'x' in xaxis.name and 'y' in yaxis.name:
        for block in inblocks:
            if block[4] <= position <= block[5]:
                # imin, imax, jmin, jmax
                newblock = block[0:4]
                blocks2d.append(newblock)

    elif 'y' in xaxis.name and 'z' in yaxis.name:
        for block in inblocks:
            if block[0] <= position <= block[1]:
                # jmin, jmax, kmin, kmax
                newblock = block[2:6]
                blocks2d.append(newblock)
    else:
        raise ValueError('Coordinate axes must be one of the following: combinations: (x, y), (x, z), (y, z).')
        
    return blocks2d


def vertices(blocks):
    vertices = []
    for block in blocks:
        x = block[0:2]
        y = block[2:4]
        z = block[4:6]
        verts = [[x[0], y[0], z[0]], 
                [x[0], y[0], z[1]],
                [x[0], y[1], z[0]],
                [x[0], y[1], z[1]], 
                [x[1], y[0], z[0]],
                [x[1], y[0], z[1]],
                [x[1], y[1], z[0]], 
                [x[1], y[1], z[1]]]
        vertices.append(verts)

    return vertices


def faces(blocks):

    verts = vertices(blocks)
    faces = []
    for v in verts:
        fac = [[v[0], v[1], v[3], v[2]],
               [v[2], v[3], v[7], v[6]],
               [v[6], v[7], v[5], v[4]],
               [v[4], v[5], v[1], v[0]],
               [v[1], v[3], v[7], v[5]],
               [v[0], v[2], v[6], v[4]]]
        faces.append(fac)
    return faces