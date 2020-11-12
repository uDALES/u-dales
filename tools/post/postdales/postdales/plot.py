import numpy as np
import math
import itertools
import postdales as dap
import operator
import datetime
import time
import os
import copy
import IPython.display
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# -----------
# PLOTTING SUPPORT FUNCTIONS

def getaxis(ncax, limits=None, sparse=None):
    
    if limits is None:
        ax1 = slice(0, None, sparse)
    else:
        ax1 = slice(limits[0], limits[1], sparse)

    # axis data
    axis = ncax[ax1]
    
    return axis


def dataslice(xvar, yvar, spacepos=None, timepos=None, limits=None, sparse=None):
    """Function that transforms a 2D, 3D or 4D data array into a 2D data array. 
    Data is assumed to be in netCDF4 coordinate shape: (time, zxis, yaxis, xaxis)."""
    
    if limits is None:
        range1 = slice(0, None, sparse) 
        range2 = slice(0, None, sparse)
    else:
        range1 = slice(limits[0], limits[1], sparse) 
        range2 = slice(limits[2], limits[3], sparse)       
    
    # reduce field to two axes by fixing position in space and time
    if 'x' in xvar and 'z' in yvar:
        if spacepos is not None and timepos is not None:
            dataslice = [timepos, range2, spacepos, range1]
        elif spacepos is not None and timepos is None:
            dataslice = [range2, spacepos, range1]
        elif timepos is not None and spacepos is None:
            dataslice = [timepos, range2, range1]
        elif timepos is None and spacepos is None:
            dataslice = range1
        else:
            raise ValueError('Dimensions missmatch')

    elif 'x' in xvar and 'y' in yvar:
        if spacepos is not None and timepos is not None:
            dataslice = [timepos, spacepos, range2, range1]
        elif spacepos is not None and timepos is None:
            dataslice = [spacepos, range2, range1]
        elif timepos is not None and spacepos is None:
            dataslice = [timepos, range2, range1]
        elif timepos is None and spacepos is None:
            dataslice = range1
        else:
            raise ValueError('Dimensions missmatch')

    elif 'y' in xvar and 'z' in yvar:
        if spacepos is not None and timepos is not None:
            dataslice = [timepos, range2, range1, spacepos]
        elif spacepos is not None and timepos is None:
            dataslice = [range2, range1, spacepos]
        elif timepos is not None and spacepos is None:
            dataslice = [timepos, range2, range1]
        elif timepos is None and spacepos is None:
            dataslice = range1
        else:
            raise ValueError('Dimensions missmatch')
    
    else:
        raise KeyError('Variables xvar must be either "x" or "y", '
                       'and  "y" or "z" for yvar')
        
    return dataslice


def get_2dblocks(blocks, xvar, yvar, position, limits=None):
    
    # limits on xvar and yvar axes
    if limits is None:
        limits = [0, math.inf, 0, math.inf]
    
    inblocks = copy.deepcopy(blocks)  # so we do not manipulate original blocks
    blocks2d = []
    if 'x' in xvar and 'z' in yvar:
        for block in inblocks:
            if limits[0] <= block[0] <= limits[1] and limits[2] <= block[4] <= limits[3]:
                if block[2] <= position <= block[3]:
                    newblock = block[0:2] + block[4:6]  # imin, imax, kmin, kmax
                    # make blocks shorter if they exceed limit
                    if newblock[1] > limits[1]:
                        newblock[1] = limits[1]
                    if newblock[3] > limits[3]:
                        newblock[3] = limits[3]
                    blocks2d.append(newblock)
                    
    elif 'x' in xvar and 'y' in yvar:
        for block in inblocks:
            if limits[0] <= block[0] <= limits[1] and limits[2] <= block[2] <= limits[3]:
                if block[4] <= position <= block[5]:
                    newblock = block[0:4]  # imin, imax, jmin, jmax
                    # make blocks shorter if they exceed limit
                    if newblock[1] > limits[1]:
                        newblock[1] = limits[1]
                    if newblock[3] > limits[3]:
                        newblock[3] = limits[3]
                    blocks2d.append(newblock)

    elif 'y' in xvar and 'z' in yvar:
        for block in inblocks:
            if limits[0] <= block[2] <= limits[1]:
                if block[0] <= position <= block[1]:
                    newblock = block[2:6]  # jmin, jmax, kmin, kmax
                    # make blocks shorter if they exceed limit
                    if newblock[1] > limits[1]:
                        newblock[1] = limits[1]
                    if newblock[3] > limits[3]:
                        newblock[3] = limits[3]
                    blocks2d.append(newblock)

    else:
        raise KeyError('Variables must be either "x" or "y" on x-axis, '
                       'and  "y" or "z" on y-axis')
        
    return blocks2d


def magnitude(value):
    if (round(value, 4) == 0):   # setting a threshold for value, everything below is 0.0001 is treated as zero
        return 0
    return int(math.floor(math.log10(abs(value))))


def contouropt(fielddata, cbarmap=plt.cm.get_cmap('RdYlBu_r'), cbartype='extrema'):
    
    # proper rounding for values of all orders of magnitude
    vmax = np.nanmax(fielddata)
    vmin = np.nanmin(fielddata)
    magnitudes = [magnitude(vmin), magnitude(vmax)]
    lowest_magnitude = np.min(magnitudes)
    if lowest_magnitude < 0:
        m = 10**(-lowest_magnitude + 1)
        r = (-lowest_magnitude + 1)
    else:
        m = 10
        r = 2
    cmax = math.ceil(vmax*m)/m
    cmin = math.floor(vmin*m)/m
    tmax = math.floor(vmax*m)/m
    tmin = math.ceil(vmin*m)/m
    
    if cbartype == 'extrema':
        ticks = np.linspace(tmin, tmax, 9)
        ticks = [round(t, r) for t in ticks]
        levels = np.linspace(cmin, cmax, 256)
        contour_opts = {'levels': levels, 'cmap': cbarmap}
    
    elif cbartype == 'equal':
        tex = max(tmax, abs(tmin))
        cex = max(cmax, abs(cmin))
        ticks = np.linspace(-tex, tex, 9)
        ticks = [round(t, r) for t in ticks]
        levels = np.linspace(-cex, cex, 256)
        contour_opts = {'levels': levels, 'cmap': cbarmap}
        
    elif cbartype == 'default':
        ticks = np.linspace(vmin, vmax, 9)
        ticks = [round(t, r) for t in ticks]
        levels = np.linspace(vmin, vmax, 256)
        contour_opts = {'levels': levels, 'cmap': cbarmap}
     
    elif cbartype =='centered':
        if cmin < 0:
            levels1 = np.linspace(cmin, 0, 128)
            levels2 = np.linspace(0, cmax, 128)
            levels = np.concatenate((levels1[:-1], levels2))

            ticks1 = np.linspace(tmin, 0, 5)
            ticks2 = np.linspace(0, tmax, 5)
            ticks = np.concatenate((ticks1[:-1], ticks2))

            midpoint = abs(cmin)/(abs(cmin) + cmax)
        else:
            levels = np.linspace(cmin, cmax, 256)
            ticks = np.linspace(tmin, tmax, 9)
            midpoint = 0.0

        newlevels1 = np.linspace(0.0, midpoint, 128)
        newlevels2 = np.linspace(midpoint, 1.0, 128)
        newlevels = np.concatenate((newlevels1[:-1], newlevels2))

        cdict = {'red': [], 'green': [], 'blue': [],'alpha': []}

        for i, newindex in enumerate(newlevels):
            red, green, blue, alpha = cbarmap(i)

            cdict['red'].append((newindex, red, red))
            cdict['green'].append((newindex, green, green))
            cdict['blue'].append((newindex, blue, blue))
            cdict['alpha'].append((newindex, alpha, alpha))

        shifted_cmap = mpl.colors.LinearSegmentedColormap('shiftedcmap', cdict)
        plt.register_cmap(cmap=shifted_cmap)
                
        ticks = [round(t, r) for t in ticks]
        contour_opts = {'levels': levels, 'cmap': shifted_cmap} 
            
    else:
        raise KeyError('Type must be either "extrema", "equal" or "centered".')
    
    # adjust padding to ticks
    m = magnitude(ticks[0])  # smallest value of ticks
    padding = 20
    # if negative numbers
    if ticks[0] < 0:
        padding = 30
    if m < 0:
        padding += (-m)*10

    return contour_opts, ticks, padding 


# -----------
# 3D block plot support functions

def vertices(blocks):
    vertices = []
    for block in blocks:
        x = block[0:2]
        y = block[2:4]
        z = block[4:6]
        vertices.append(list(itertools.product(x, y, z)))

    return vertices


def edges(blocks):

    edges = []
    for block in blocks:
        ed = []
        xlen = block[1] - block[0]
        ylen = block[3] - block[2]
        zlen = block[5] - block[4]
        verts = vertices([block])
        for s, e in itertools.combinations(np.array(verts[0]), 2):
            if np.sum(np.abs(s - e)) == xlen:
                edge = list(zip(s, e))
                ed.append(edge)
            elif np.sum(np.abs(s - e)) == ylen:
                edge = list(zip(s, e))
                ed.append(edge)
            elif np.sum(np.abs(s - e)) == zlen:
                edge = list(zip(s, e))
                ed.append(edge)
        edges.append(ed)

    return edges


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


def roofs(blocks):

    vertices = []
    
    for block in blocks:
        x = block[0:2]
        y = block[2:4]
        z = block[5]
        # order is important!
        verts = [[x[0], y[0], z], [x[0], y[1], z], [x[1], y[1], z], [x[1], y[0], z]]
        vertices.append(verts)
    
    return np.array(vertices)


def fronts(blocks):

    vertices = []
    
    for block in blocks:
        x = block[0]
        y = block[2:4]
        z = block[4:6]
        # order is important!
        verts = [[x, y[0], z[0]], [x, y[0], z[1]], [x, y[1], z[1]], [x, y[1], z[0]]]
        vertices.append(verts)
    
    return np.array(vertices)


# -----------
# MAIN PLOTTING FUNCTIONS (GENERAL)

def layout(blocks, ax=None, limits=None, **kwargs):
    # block in blocks has form [xmin, xmax, ymin, ymax, ...]
    # limits are plot limits and has form [xmin, xmax, ymin, ymax]
    
    if ax is None:
        ax = plt.axes(projection='3d')

    for block in blocks:
        start = tuple([block[0], block[2]])  # xmin, ymin
        width = block[1] - block[0]  # xmax - xmin
        height = block[3] - block[2]  # ymax - ymin
        p = patches.Rectangle(start, width, height, **kwargs)
        ax.add_patch(p)

    if limits is not None:
        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])
        
    plt.axis('scaled')
    plt.xlabel('x')
    plt.ylabel('y')

    return


def blocks3d(blocks, ax=None, limits=None, **kwargs):
    
    if ax is None:
        ax = plt.axes(projection='3d')
        
    fas = faces(blocks)

    for sides in fas:
        collection = Poly3DCollection(sides, **kwargs)
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
    
    plt.gca().patch.set_facecolor('none')  # sets background

    return


def plot_2dblocks(blocks, **kwargs):
    """Function that adds building addblocks to contour plot"""

    ax = plt.gca()

    for block in blocks:
        start = tuple([block[0], block[2]])  # xmin, ymin
        width = block[1] - block[0]  # xmax - xmin
        height = block[3] - block[2]  # ymax - ymin
        p = patches.Rectangle(start, width, height, **kwargs)
        ax.add_patch(p)

    return


def intersection_plane(blocks, xvar, yvar, position, limits):
    """Function that adds plane around  intersecting blocks"""

    intersecting_blocks = []
    vertices = []
    
    if 'x' in xvar and 'z' in yvar:
        plane = []  # one single intersection plane
        ul_corner = [limits[0], position, limits[5]]
        ll_corner = [limits[0], position, limits[4]]
        lr_corner = [limits[1], position, limits[4]]
        ur_corner = [limits[1], position, limits[5]]
        plane.extend([ul_corner, ll_corner])
        # sort first by x then by y
        sorted_blocks = sorted(blocks, key=operator.itemgetter(0,2))
        for block in sorted_blocks:
            if block[2] <= position <= block[3]:
                intersecting_blocks.append(block)
                x = block[0:2]
                z = block[4:6]
                # order is important!
                verts = [[x[0], position, z[0]], [x[0], position, z[1]], [x[1], position, z[1]], [x[1], position, z[0]]]
                plane.extend(verts)
        plane.extend([lr_corner, ur_corner])
        vertices.append(plane)

    elif 'x' in xvar and 'y' in yvar:
        ul_corner = [limits[0], limits[3], position]
        ll_corner = [limits[0], limits[2], position]
        lr_corner = [limits[1], limits[2], position]
        ur_corner = [limits[1], limits[3], position]
        plane = [ul_corner, ll_corner, lr_corner, ur_corner]
        vertices.append(plane)
        for block in blocks:
            if block[4] <= position <= block[5]:
                intersecting_blocks.append(block)
                x = block[0:2]
                y = block[2:4]
                # order is important!
                verts = [[x[0], y[0], position], [x[0], y[1], position], [x[1], y[1], position], [x[1], y[0], position]]
                vertices.append(verts)

    elif 'y' in xvar and 'z' in yvar:
        ul_corner = [position, limits[2], limits[5]]
        ll_corner = [position, limits[2], limits[4]]
        lr_corner = [position, limits[3], limits[4]]
        ur_corner = [position, limits[3], limits[5]]
        plane = []  # one single intersection plane
        plane.extend([ul_corner, ll_corner])
        # sort first by y then by x
        sorted_blocks = sorted(blocks, key=operator.itemgetter(2,0))
        for block in sorted_blocks:
            if block[0] <= position <= block[1]:
                intersecting_blocks.append(block)
                y = block[2:4]
                z = block[4:6]
                # order is important!
                verts = [[position, y[0], z[0]], [position, y[0], z[1]], [position, y[1], z[1]], [position, y[1], z[0]]]
                plane.extend(verts)
        plane.extend([lr_corner, ur_corner])
        vertices.append(plane)

    else:
        raise KeyError('Variables must be either "x" or "y" on x-axis,'
                       ' and  "y" or "z" on y-axis')
    
    return intersecting_blocks, np.array(vertices)


def plane3d(vertices, ax=None, **kwargs):
    
    if ax is None:
        ax = plt.axes(projection='3d')
    for plane in vertices:
        ax.add_collection3d(Poly3DCollection([plane], **kwargs))
    
    return


# -----------
# MAIN PLOTTING FUNCTIONS (USING FIELD OBJECT)

def minilayout(field_object, ax=None):
    
    if ax is None:
        ax = plt.axes(projection='3d')
       
    # set colours
    try:
        planecolor = mpl.colors.to_rgba(field_object.color, alpha=0.2)
    except Exception:
        planecolor = (0.5, 0.5, 0.5, 0.2)
    # check whether dark or light plotstyle
    axtick_color = ax.xaxis._axinfo['tick']['color']
    # dark plot style
    if axtick_color is 'w':
        block_edgecolor = (1, 1, 1, 0.2)
        block_facecolor = (1, 1, 1, 0.4)
    # light plot style
    else:
        block_edgecolor = (0, 0, 0, 0.1)
        block_facecolor = (0, 0, 0, 0.2) 

    # get intersection plane
    intersecting_blocks, verts = intersection_plane(field_object.miniblocks, 
                                                    field_object.xname, 
                                                    field_object.yname, 
                                                    field_object.minispos, 
                                                    limits=field_object.minilimits)
    # plot plane
    plane = plane3d(verts, ax=ax, facecolors=planecolor, edgecolor=field_object.color)
    # plot all blocks
    allblocks = blocks3d(field_object.miniblocks, ax=ax, limits=field_object.minilimits, 
                         facecolors=(0, 0, 0, 0), edgecolor=block_edgecolor)
    # plot blocks intersecting the plane
    activeblocks = blocks3d(intersecting_blocks, ax=ax, limits=field_object.minilimits, facecolors=block_facecolor, edgecolor=block_edgecolor)

    ax.set_aspect('equal')
    ax.view_init(25, 235)
    ax.dist=12
    ax.grid(False)
    ax.xaxis.pane.set_edgecolor('k')
    ax.yaxis.pane.set_edgecolor('k')
    ax.zaxis.pane.set_edgecolor('k')
    ax.w_xaxis.line.set_color((0, 0, 0, 0.1))
    ax.w_yaxis.line.set_color((0, 0, 0, 0.1))
    ax.w_zaxis.line.set_color((0, 0, 0, 0.1))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.xaxis.labelpad=-15
    ax.yaxis.labelpad=-15
    ax.zaxis.labelpad=-15
    
    try:
        ax.set_title(field_object.reference)
    except Exception:
        pass
    
    return


def contourfield(field_object, ax=None):
    
    if ax is None:
        ax = plt.axes()
        
    # contour options, if defined
    try:
        contour_opts = field_object.contour_opts
    except Exception:
        contour_opts = {}
    # block options, if defined
    try:
        block_opts = field_object.block_opts
    except Exception:
        block_opts = {}

    # plot field
    cs = ax.contourf(field_object.xax, field_object.yax, field_object.data,
                      **contour_opts)
    ax.axis('scaled')
    
    # plot blocks, if defined
    try:
        plot_2dblocks(field_object.blocks2d, **block_opts)
    except Exception:
        pass
    
    # labels, if defined
    try:
        ax.set_xlabel('%s (m)' % (field_object.xname))
        ax.set_ylabel('%s (m)' % (field_object.yname))
        ax.set_title("%s %s" % (field_object.longname, field_object.name))
    except Exception:
        pass

    # colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(cs, cax=cax)
    try:
        cbar.set_ticks(field_object.ticks)
        cbar.set_label(r"$%s$" % field_object.units)
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs, ha='right')
        cbar.ax.yaxis.set_tick_params(pad=field_object.padding)
    except Exception:
        pass

    return


def video_contourfield(field_object, fig=None, ax=None, 
                       show=True, save=False, **save_kwargs):
    
    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = plt.axes()
    
    # create directory to save plots
    if save is True:
        today = datetime.datetime.now()
        timestamp = today.strftime('%Y-%m-%d-%H%M')
        savedir = "./videos"
        if not os.path.exists(savedir):
            os.mkdir(savedir)
        savehere = os.path.join(savedir, timestamp, '')
        os.mkdir(savehere)
        
    # get time range
    timerange = field_object.timepos
    tsteps = field_object.timesteps[timerange]
    
    # contour options, if defined
    try:
        contour_opts = field_object.contour_opts
    except Exception:
        contour_opts = {}
    # block options, if defined
    try:
        block_opts = field_object.block_opts
    except Exception:
        block_opts = {}

    # plot the basic field
    cs = ax.contourf(field_object.xax, field_object.yax, 
                     field_object.data[0], **contour_opts)
    ax.axis('scaled')
    
    # plot blocks, if defined
    try:
        plot_2dblocks(field_object.blocks2d, **block_opts)
    except Exception:
        pass
    
    # labels, if defined
    try:
        ax.set_xlabel('%s (m)' % (field_object.xname))
        ax.set_ylabel('%s (m)' % (field_object.yname))
        ax.set_title("%s %s" % (field_object.longname, field_object.name))
    except Exception:
        pass

    # colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    cbar = plt.colorbar(cs, cax=cax)
    try:
        cbar.set_ticks(field_object.ticks)
        cbar.set_label(r"$%s$" % field_object.units)
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs, ha='right')
        cbar.ax.yaxis.set_tick_params(pad=field_object.padding)
    except Exception:
        pass

    # display video
    for i, field in enumerate(field_object.data):
        ax.collections = []
        cs = ax.contourf(field_object.xax, field_object.yax, 
                         field, **contour_opts)
        ax.axis('scaled')

        # plot blocks, if defined
        try:
            plot_2dblocks(field_object.blocks2d, **block_opts)
        except Exception:
            pass
       
        # labels
        ax.set_xlabel('%s (m)' % (field_object.xname))
        ax.set_ylabel('%s (m)' % (field_object.yname))
        ax.set_title("%s %s at t = %s s" % (field_object.longname, 
                                            field_object.name,
                                            int(tsteps[i]-tsteps[0])))
        # save image
        if save is True:
            plt.savefig("%s/field-%03d" %(savehere, i), 
                        bbox_inches='tight', 
                        **save_kwargs)
        if show is True:
            IPython.display.display(fig)
            IPython.display.clear_output(wait=True)
            time.sleep(0.1)
            
    return fig
        

# -----------
# FIELD OBJECT THAT STORES ALL DATA FOR PLOTS

class Field(object):
    def __init__(self):
        self.info = {}
        
    def __getattr__(self, name):
        if "info" in self.__dict__ and name in self.info:
            return self.info[name]
        else:
            return self.__dict__[name]
    
    def __setattr__(self, name, value):
        if "info" in self.__dict__ and name in self.info:
            self.info[name] = value
        else:
            self.__dict__[name] = value
            
    def add_info(self, name, value):
        self.info[name] = value

    def delete_info(self, name):
        if "info" in self.__dict__ and name in self.info:
            del self.info[name]
            
        
def fieldinfo(field_object, netcdf_variable):
    field_object.add_info('name', netcdf_variable.name)
    field_object.add_info('longname', netcdf_variable.longname)
    field_object.add_info('units', netcdf_variable.units) 
    
    return field_object
    
    
def fieldaxesinfo(field_object, netcdf_variable):
    # number of netcdf coordinates
    ncoor = netcdf_variable.ndim
    
    if ncoor is 4:
        # (time, zaxis, yaxis, xaxis)
        zaxis = netcdf_variable.dimensions[1]
        yaxis = netcdf_variable.dimensions[2]
        xaxis = netcdf_variable.dimensions[3]
        
    elif ncoor is 3:
        # data averaged in y
        # (time, zaxis, xaxis)  
        zaxis = netcdf_variable.dimensions[1]
        xaxis = netcdf_variable.dimensions[2]
        yaxis = ''
        
    elif ncoor is 2:
        # data averaged in x and y
        # (time, zaxis)  
        zaxis = netcdf_variable.dimensions[1]
        yaxis = ''
        xaxis = ''
        
    field_object.add_info('zaxis', zaxis)
    field_object.add_info('yaxis', yaxis)
    field_object.add_info('xaxis', xaxis)
    
    return field_object


def fieldsliceinfo(field_object, fix, spacepos, timepos):
    
    field_object.add_info('spacepos', spacepos)
    field_object.add_info('timepos', timepos)
    field_object.add_info('fix', fix)
    
    if 'zaxis' in fix:
        field_object.add_info('xname', field_object.xaxis)
        field_object.add_info('yname', field_object.yaxis)
                
    elif 'yaxis' in fix:
        field_object.add_info('xname', field_object.xaxis)
        field_object.add_info('yname', field_object.zaxis)
        
    elif 'xaxis' in fix:
        field_object.add_info('xname', field_object.yaxis)
        field_object.add_info('yname', field_object.zaxis)
        
    return field_object


def fieldgetaxes(field_object, exp):

    # netcdf axis
    ncxax = exp.data[field_object.xaxis]
    ncyax = exp.data[field_object.yaxis]
    nczax = exp.data[field_object.zaxis]
    # time axis
    # tricky to know how we get to the right time data...
    if field_object.name.endswith('t'):  # if time averaged variable
        try:
            tsteps = exp.data['tstats']
            field_object.timesteps = tsteps
        except Exception:
            try:
                tsteps = exp.data['time']
                field_object.timesteps = tsteps
            except Exception:
                pass
    else:
        try:
            tsteps = exp.data['tfield']
            field_object.timesteps = tsteps
        except Exception:
            try:
                tsteps = exp.data['tstats']
                field_object.timesteps = tsteps
            except Exception:
                try:
                    tsteps = exp.data['time']
                    field_object.timesteps = tsteps
                except Exception:
                    pass
    
    # get axis data
    xaxis = getaxis(ncxax)
    yaxis = getaxis(ncyax)
    zaxis = getaxis(nczax)
    
    field_object.xgrid = xaxis
    field_object.ygrid = yaxis
    field_object.zgrid = zaxis
    
    return field_object


def fieldgetblocks(field_object, block_indices):
    
    # SOME CORRECTION OF INDICES
    # need to shift indices if on cell centre, i.e. xt, yt, zt
    blocks = copy.deepcopy(block_indices)
    if 't' in field_object.xaxis:
        for block in blocks:
            block[1] -= 1
            
    if 't' in field_object.yaxis:
        for block in blocks:
            block[3] -= 1
            
    if 't' in field_object.zaxis:
        for block in blocks:
            block[5] -= 1

    dimblocks = dap.tools.getdimblocks(blocks, 
                                      field_object.xgrid,
                                      field_object.ygrid,
                                      field_object.zgrid)
    
    field_object.add_info('dimblocks', dimblocks)
    
    return field_object


def fieldsetaxes(field_object, limits=None, sparse=None):
    
    if limits is not None:
        limx = limits[0:2]
        limy = limits[2:4]
    else:
        limx = None
        limy = None
        
    if 'zaxis' in field_object.fix:
        
        fieldfixax = field_object.zgrid.copy()
        fieldxax = getaxis(field_object.xgrid.copy(), limx, sparse)
        fieldyax = getaxis(field_object.ygrid.copy(), limy, sparse)
        
    elif 'yaxis' in field_object.fix:

        fieldfixax = field_object.ygrid.copy()       
        fieldxax = getaxis(field_object.xgrid.copy(), limx, sparse)
        fieldyax = getaxis(field_object.zgrid.copy(), limy, sparse)
        
    elif 'xaxis' in field_object.fix:
        
        fieldfixax = field_object.xgrid.copy()        
        fieldxax = getaxis(field_object.ygrid.copy(), limx, sparse)
        fieldyax = getaxis(field_object.zgrid.copy(), limy, sparse)
        
    # set fixed postion
    field_object.fixax = fieldfixax
    field_object.fixpos = fieldfixax[field_object.spacepos]
    # set plot x and y axes
    field_object.xax = fieldxax
    field_object.yax = fieldyax
    # set plot limits in plot axes
    field_object.limits = [fieldxax[0], fieldxax[-1], 
                           fieldyax[0], fieldyax[-1]]

    return field_object


def fieldgetslice(field_object, netcdf_variable, 
                  limits=None, sparse=None):
    
    indices = dataslice(field_object.xname, field_object.yname, 
                        field_object.spacepos, field_object.timepos, 
                        limits, sparse)
    
    data = netcdf_variable[indices]
    
    field_object.data = data
    
    return field_object
    

def fieldcontours(fields, cbarmap=plt.cm.get_cmap('RdYlBu_r'), cbartype='extrema'):
    
    # if more fields given:
    if isinstance(fields, list):
        # setup contour options equal for all plots
        vmax = 0
        vmin = 0
        for field in fields:
            if np.nanmax(field.data) > vmax:
                vmax = np.nanmax(field.data)       
            if np.nanmin(field.data) < vmin:
                vmin = np.nanmin(field.data)
         
        # get contour plot options
        contour_opts, ticks, padding = contouropt([vmin, vmax], cbarmap, cbartype)
        
        for field in fields:
            field.add_info('ticks', ticks)
            field.add_info('contour_opts',  contour_opts)
            field.add_info('padding', padding)

    else:
        # get contour plot options
        contour_opts, ticks, padding = contouropt(fields.data, cbarmap, cbartype)

        fields.add_info('ticks', ticks)
        fields.add_info('contour_opts',  contour_opts)
        fields.add_info('padding', padding)
    
    return fields


def fieldblocks_2d(field_object, limits=None):
    
    blocks2d = get_2dblocks(field_object.dimblocks, 
                 field_object.xname, 
                 field_object.yname,
                 field_object.fixpos, 
                 field_object.limits)
    
    field_object.add_info('blocks2d', blocks2d)

    return field_object


def fieldblocks_mini3d(field_object, divx=4, divy=2):

    # blocks layout miniplot
    nbl = int(len(field_object.dimblocks)/(divx * divy))
    miniblocks = field_object.dimblocks[:nbl]

    field_object.add_info('miniblocks', miniblocks)

    # domain limits of mini layout
    xlen = round(len(field_object.xgrid)/divx)
    ylen = round(len(field_object.ygrid)/divy)
    
    minilimits=[field_object.xgrid[0], 
                field_object.xgrid[xlen], # maybe + 1?
                field_object.ygrid[0], 
                field_object.ygrid[ylen],
                field_object.zgrid[0], 
                field_object.zgrid[-1]/3]

    # modulo length of one geometry
    # so we can also indicate position corretly 
    # in any repeated layout
    if 'x' in field_object.fix:
        spos = field_object.fixpos%minilimits[1]
    elif 'y' in field_object.fix:
        spos = field_object.fixpos%minilimits[3]
    else:
        spos = field_object.fixpos

    field_object.add_info('minilimits', minilimits)
    field_object.add_info('minispos', spos)

    return field_object


# -----------
# USE DATA FROM FIELD OBJECT TO PREPARE PLOTS

def prepare_fielddata(exp, varname, fix, spos, tpos, limits=None):
    
    # create new field object
    field = Field()
    # get netcdf variable
    ncvariable = exp.data[varname]
    # set color and reference
    field.add_info('color', exp.color)
    field.add_info('reference', exp.reference)
    # get infos from netcdf variable
    field = fieldinfo(field, ncvariable)
    # get axes info from netcdf variable
    field = fieldaxesinfo(field, ncvariable)
    # get info of data slice
    field = fieldsliceinfo(field, fix, spos, tpos)
    # add axes data
    field = fieldgetaxes(field, exp)
    # get the blocks
    field = fieldgetblocks(field, exp.blocks)
    # set plot axes
    field = fieldsetaxes(field, limits)
    # add field data
    field = fieldgetslice(field, ncvariable, limits)
    # add blocks for field in 2D
    field = fieldblocks_2d(field, limits)
        
    return field


def prepare_contourplot(field, cbarmap=plt.cm.get_cmap('RdYlBu_r'), cbartype='centered'):
    
    # add contour plot options
    field = fieldcontours(field, cbarmap, cbartype)
    
    # add block plot options
#     field.block_opts = {'facecolor': 'None', 'edgecolor' : 'k'}
    field.block_opts = {'facecolor': 'w'}

    return field


def prepare_layoutplot(field, divx=4, divy=2):
    
    # get blocks for minilayout
    field = fieldblocks_mini3d(field, divx, divy)

    return field


def prepare_quiverplot(exp, uvar, vvar, wvar, fix,
                       spos, tpos, limits=None, sparse=5):
    
    # first wind speed data for background
        # create new field object
    field = Field()
    # get netcdf variables
    speed1 = exp.data[uvar]
    speed2 = exp.data[vvar]
    # set color and reference
    field.add_info('color', exp.color)
    field.add_info('reference', exp.reference)
    # add info
    field.add_info('name', r'$\sqrt{u^2 + v^2}$')
    field.add_info('longname', 'Wind speed')
    field.add_info('units', 'm/s') 

    # get axes info from u
    field = fieldaxesinfo(field, speed1)
    
    # get info of data slice
    field = fieldsliceinfo(field, fix, spos, tpos)
    # add axes data
    field = fieldgetaxes(field, exp)
    # get the blocks
    field = fieldgetblocks(field, exp.blocks)
    # set plot axes! NO SPARSE FOR FIELD
    field = fieldsetaxes(field, limits, sparse=None)
    # add field data for u
    field = fieldgetslice(field, speed1, limits, sparse=None)
    uwind = field.data.copy()
    # add field data for v
    field = fieldgetslice(field, speed2, limits, sparse=None)
    vwind = field.data.copy()
    speed = [np.sqrt(u**2 + v**2) for u, v in zip(uwind, vwind)]
    field.data = speed
    
    # add blocks for field in 2D
    field = fieldblocks_2d(field, limits)

    # create new field objects for quiver data
    field1 = Field()
    field2 = Field()
    # if axes x, quiver should be u
    # if axes y, quiver should be v
    # if axes z, quiver should be w
    if 'zaxis' in fix:
        # get netcdf variables for quiver
        ncvar1 = exp.data[uvar]  # u
        ncvar2 = exp.data[vvar]  # v
        
    elif 'yaxis' in fix:
        # get netcdf variables for quiver
        ncvar1 = exp.data[uvar]  # u
        ncvar2 = exp.data[wvar]  # w

    elif 'xaxis' in fix:
        # get netcdf variables for quiver
        ncvar1 = exp.data[vvar]  # v
        ncvar2 = exp.data[wvar]  # w
    
    # get axes info from netcdf variable
    field1 = fieldaxesinfo(field1, ncvar1)
    field2 = fieldaxesinfo(field2, ncvar2)
    # get info of data slice
    field1 = fieldsliceinfo(field1, fix, spos, tpos)
    field2 = fieldsliceinfo(field2, fix, spos, tpos)
    
    # add field data! WITH SPARSE
    field1 = fieldgetslice(field1, ncvar1, limits, sparse)
    field2 = fieldgetslice(field2, ncvar2, limits, sparse)
    
    # add axes data
    # can only add axis of one field, let's hope 
    # they are defined on the same axes
    field1 = fieldinfo(field1, ncvar1)
    field1 = fieldgetaxes(field1, exp)
#     field2 = fieldgetaxes(field2, exp)
    # set plot axes
    field1 = fieldsetaxes(field1, limits, sparse)
#     field2 = fieldsetaxes(field2, limits, sparse)
        
    # add quiver axes data to original field object
    field.qxax = field1.xax.copy()
    field.qyax = field1.yax.copy()
    # add quiver field data to original field object
    field.qfield1 = field1.data.copy()
    field.qfield2 = field2.data.copy()
    
    # add quiver plot options
    field.quiver_opts = {'scale' : 0.3, 'scale_units' : 'xy'}
    
    return field


# -----
# PLOT FUNCTIONS WITH FIELD OBJECT

def layout_fieldplot(field, quiver=False):
    
    # plot
    fig = plt.figure(figsize=(8, 6))
    # plot block layout
    ax1 = plt.subplot2grid((1, 3), (0, 0), projection='3d')
    minilayout(field, ax=ax1)
    # plot contour field
    ax2 = plt.subplot2grid((1, 3), (0, 1), colspan=2)
    contourfield(field, ax=ax2)
        
    # quiver plot
    if quiver is True:
        # quiver options, if defined
        try:
            quiver_opts = field_object.quiver_opts
        except Exception:
            quiver_opts = {}
        ax2.quiver(field.qxax, field.qyax, 
                   field.qfield1, field.qfield2,
                   **quiver_opts)
    plt.tight_layout()
    plt.show()
    
    return


def layout_fieldplot_video(field, show=True, save=False, **video_kwargs):
    # plot
    fig = plt.figure(figsize=(8, 6))
    # plot block layout
    ax1 = plt.subplot2grid((1, 3), (0, 0), projection='3d')
    minilayout(field, ax=ax1)
    # plot contour field
    ax2 = plt.subplot2grid((1, 3), (0, 1), colspan=2)
    video_contourfield(field, fig, ax2, show, save, **video_kwargs)
    
    return
