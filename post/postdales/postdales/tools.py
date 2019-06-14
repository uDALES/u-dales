import os
import netCDF4
import f90nml
import glob
import copy
import operator
import numpy as np
import pickle


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


def referencelists():
    referencelist = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)']
    colorlist = ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9', '#DE8F05', '#FBAFE4', '#CA9161', '#8C0900']
    markerlist = ['v', '^', 's', '*', 'h', 'd', 'P',  'p', 'X', 'o']
    
    return referencelist, colorlist, markerlist

    
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


class Exp(object):
    def __init__(self, name):
        self.name = name
        self.data = {}
        self.blockstats = {}
        self.parameters = {}
        
        
    def __getattr__(self, name):
        if "data" in self.__dict__ and name in self.data:
            return self.data[name]
        elif "blockstats" in self.__dict__ and name in self.blockstats:
            return self.blockstats[name]
        elif "parameters" in self.__dict__ and name in self.parameters:
            return self.parameters[name]
        return self.__dict__[name]
    
    
    def __setattr__(self, name, value):
        if "data" in self.__dict__ and name in self.data:
            self.data[name] = value
        elif "blockstats" in self.__dict__ and name in self.blockstats:
            self.blockstats[name] = value
        elif "parameters" in self.__dict__ and name in self.parameters:
            self.parameters[name] = value
        else:
            self.__dict__[name] = value
            
            
    def delete_item(self, name):
        if name in self.__dict__:
            del self.__dict__[name]                    
                
            
    def add_data(self, name, value):
        self.data[name] = value

        
    def delete_data(self, name):
        if "data" in self.__dict__ and name in self.data:
            del self.data[name]
            
            
    def add_blockstats(self, name, value):
        self.blockstats[name] = value
        
        
    def add_parameters(self, dictofdicts):
        for key, val in dictofdicts.items():
        # first dict list, e.g. 'RUN', 'DOMAIN', etc.
            if isinstance(val, dict):
                # second dict list, e.g. 'runtime', 'fieldvars', etc.
                for name, value in val.items():
                    self.parameters[name] = value
        # if normal dict
            else:
                self.parameters[key] = val               

                
    def fulllist(self, only=None):
        if only is None:
            # print everything
            for key, val in vars(self).items():
                if isinstance(val, dict):
                    print(key.upper())
                    for subkey, subval in val.items():
                        if isinstance(subval, (list, tuple, np.ndarray, f90nml.namelist.Namelist, netCDF4._netCDF4.Variable)):
                            print("\t", subkey, ":", type(subval))
                        else:
                            print("\t", subkey, ":", subval)
                else:
                    if isinstance(val, (list, tuple, np.ndarray)):
                        print(key, ":", type(val))
                    else:
                        print(key, ":", val)
        
        elif only in self.__dict__:
            # print only sub dictionary "only"
            for key, val in vars(self).items():
                if key is only:
                    print(key.upper())
                    for subkey, subval in val.items():
                        if isinstance(subval, (list, tuple, np.ndarray, f90nml.namelist.Namelist, netCDF4._netCDF4.Variable)):
                            print("\t", subkey, ":", type(subval))
                        else:
                            print("\t", subkey, ":", subval)

        else:
            print("Key not found. Try only=None to see all data and keys.")
        
        
    def keylist(self, only=None):
        if only is None:
            # return list(vars(self).keys())
            for key, val in vars(self).items():
                if isinstance(val, dict):
                    print(key.upper())
                    for subkey in val.keys():
                            print("\t", subkey)
                else:
                    print(key)
                    
        elif only in self.__dict__:
            # print only sub dictionary "only"
            for key, val in vars(self).items():
                if key is only:
                    print(key.upper())
                    for subkey in val.keys():
                        print("\t", subkey)

        else:
            print("Key not found. Try only=None to see all keys.")

            
def getexp(exps, reference):
    for exp in exps:
        if exp.name == reference:
            return exp
        elif exp.reference == reference:
            return exp

    print("Could not find exp in the list.")
    return
    
    
class ExpStored(object):
    
    def __init__(self):
        pass
    
    @staticmethod
    def create_type(name='DynamicType', dict={}):
        return type(name, (object,), dict)
    
    @staticmethod
    def save(t, fh):
        dict = t.__dict__.copy()
        name = t.name
        for key in dict.keys():
            if key.startswith('__') and key.endswith('__'):
                del dict[key]
        pickle.dump((name, dict), fh)
        
    @classmethod
    def load(cls, fh):
        name, dict = pickle.load(fh)
        return cls.create_type(name, dict)
    
    
def load_saved_exp(expnr, directory):
    
    # check all paths for exp data
    filepaths = "**/*" + expnr + "*"
    searchpath = os.path.join(directory, filepaths)
    filelist = glob.glob(searchpath, recursive=True)
    
    if not filelist:
        # if list of files is empty
        print("Could not find data for", expnr, "Setting up empty exp object.")
        exp = Exp(expnr)
    else:
        # get exp data
        try:
            filename = filelist[0]  # get string only. if several exp files, take the first
            with open(filename, 'rb') as file:
                exp_saved = ExpStored.load(file)
            # copy to new exp object   
            exp = Exp(exp_saved.name)
            for name, value in exp_saved.__dict__.items():
                if not name.startswith("__"):
                    exp.__setattr__(name, value)
        except Exception:
            print("Could not load data for", expnr, "Setting up empty exp object.")
            exp = Exp(expnr)
        
    return exp


def delete_netcdfs(exp):
    netcdfs_add = []
    netcdfs_delete = []    
    for key, val in exp.data.items():
        if isinstance(val, netCDF4._netCDF4.Variable):
            if len(val.shape) > 2:
                netcdfs_delete.append(key)
            else:
                netcdfs_add.append(key)
    [exp.delete_data(d) for d in netcdfs_delete]
    [exp.add_data(a, exp.data[a][:]) for a in netcdfs_add]
    print("Added:", netcdfs_add)
    print("Deleted:", netcdfs_delete)
    
    return exp


def compare_parameters(exps, parameter, position=None, name='reference'):
    """Function outputs different parameters across Exp objects"""
    parameter_dict = {}
    for exp in exps:
        key = exp.__getattr__(name)
        if isinstance(position, int):
            vals = exp.__getattr__(parameter)
            val = vals[position]
        else:
            val = exp.__getattr__(parameter)
        parameter_dict[key] = val
        
    return parameter_dict


def sort_parameters(parameter_dict, reverse=False):
    """Function sorts parameters across Exp objects. Fails if value is not float or int."""
    sorted_by_value = sorted(parameter_dict.items(), key=lambda kv: kv[1], reverse=reverse)
    # print sorted list
    # for key, value in sorted_by_value:
    #    print(key, value)
    return sorted_by_value


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


def generate_canyons(blocks, limits):
    
    testblocks = copy.deepcopy(blocks)
    x0 = limits[0]
    imax = limits[1]
    y0 = limits[2]
    jtot = limits[3]
    canyons = []

    for block in testblocks:
        # 1. compare to other blocks and find all in same row
        samerow = []
        for otherblock in blocks:
            # show which blocks are at same y-heights:
            if block[2] <= otherblock[2] < block[3]:
                samerow.append(otherblock)
            # block[2] truely smaller than otherblock[3] to 
            # not catch the ones that end where my block begins
            elif block[2] < otherblock[3] <= block[3]:
                samerow.append(otherblock)
            elif otherblock[2] <= block[2] < otherblock[3]:
                samerow.append(otherblock)
            elif otherblock[2] < block[3] <= otherblock[3]:
                samerow.append(otherblock)

        # 2. find closest block in y: sort first by y then by x
        samerow = sorted(samerow, key=operator.itemgetter(2,0))
        for otherblock in samerow:
        # show only blocks with larger x start. because we sorted,
        # the first hit should always be the next block in y-line
            if block[1] <= otherblock[0]:
                nextblocky = copy.deepcopy(otherblock)
                # when first block in same row found, break and go to next
                break
            # no block with larger x, we are at end of line
            else:
                # imaginative block with same y coordinates as original block
                # and xcoordinates at domain end
                nextblocky = [imax, imax, block[2], block[3]] 

        # 3. find closest block in x: sort first by x then by y
        samerow = sorted(samerow, key=operator.itemgetter(0,2))
        for otherblock in samerow:
        # show only blocks with larger x start. because we sorted,
        # the first hit should always be the next block in x-line
        # less or equal to also get adjecent blocks
            if block[1] <= otherblock[0]:
                nextblockx = copy.deepcopy(otherblock)
                # when first block in same row found, break and go to next
                break
            # no block with larger x, we are at end of line
            else:
                # imaginative block with same y coordinates as original block
                # and xcoordinates at domain end
                nextblockx = [imax, imax, block[2], block[3]]

        # 4. choose the right next block
        # compare x and y next blocks:
        if nextblockx == nextblocky:
            nextblock = copy.deepcopy(nextblocky)

        # if xblock is same shape as block just shifted then use it
        elif block[2:4] == nextblockx[2:4]:
            nextblock = copy.deepcopy(nextblockx)

        # if x block is adjecent to block or some block overlap
        # choose y block and make it smaller by x block
        elif block[2] < nextblockx[2] and nextblockx[2] < nextblocky[3]:
            nextblock = copy.deepcopy(nextblocky)
            nextblock[3] = nextblockx[2]

        else:
            print("Canyon case that has not been defined!")
            nextx1 = min(nextblockx[0], nextblocky[0])
            nextx2 = min(nextblockx[1], nextblocky[1])
            nexty1 = min(nextblockx[2], nextblocky[2])
            nexty2 = min(nextblockx[3], nextblocky[3])
            nextblock = [nextx1, nextx2, nexty1, nexty2]

        # 5. if ymax of first block is larger define leftover block
        if nextblock[3] < block[3]:
            # adjecent
            if block[1] == nextblockx[0]:
                # and in corner
                if nextblock[3] == nextblockx[2]:
                    leftblock = [block[0], nextblockx[0], nextblockx[3], block[3]]
                else:
                    leftblock = [block[0], nextblockx[0], nextblockx[2], block[3]]
            else:
                leftblock = [block[0], nextblock[0], nextblock[3], block[3]]
                nextblock[3] = block[3]

            # if not just a line:
            if leftblock[3] > leftblock[2]:
                testblocks.append(leftblock)   
        # set nextblock ymax to first block ymax
        elif nextblock[3] > block[3]:
            nextblock[3] = block[3]

        # 6. check if there is a street canyon with no other blocks
        if block[2] < nextblockx[2] and block[2] < nextblocky[2]:
            nextblock = [imax, imax, block[2], nextblock[2]]


        # 7. define canyon in between
        # canyon x coordinates: xmin of first block, xmax of next block
        cxmin = block[1]
        cxmax = nextblock[0]
        # canyon y coordinates: ymin of first block, ymax of next block
        cymin = block[2]
        cymax = nextblock[3]
        canyon = [cxmin, cxmax, cymin, cymax]
        canyons.append(canyon)

    # 8. find canyon onlys
    # now sort first x then y
    allxblocks = sorted(blocks, key=operator.itemgetter(0,2))
    # find all blocks in first x column
    x1col = []
    for block in allxblocks:
        if block[0] == x0 and not all(b == 0 for b in block):
            x1col.append(block)

    cymins = [block[3] for block in x1col]
    ymaxs = [block[2] for block in x1col]
    cymaxs = [*ymaxs[1:], jtot] # discard first ymax and use jtot as final value

    for cymin, cymax in zip(cymins, cymaxs):    
        # check if there are blocks later on in that row
        samerow = []
        for otherblock in allxblocks:
            if cymin <= otherblock[2] < cymax:
                samerow.append(otherblock)
            elif cymin < otherblock[3] <= cymax:
                samerow.append(otherblock)
            elif otherblock[2] <= cymin < otherblock[3]:
                samerow.append(otherblock)
            elif otherblock[2] < cymax <= otherblock[3]:
                samerow.append(otherblock)
        if samerow:
            xmax = min([block[0] for block in samerow])
        else:
            xmax = imax
            
        xmin = x0
        ymin = cymin
        ymax = cymax

        canyon = [xmin, xmax, ymin, ymax]
        canyons.append(canyon)   
        
    return canyons


def getcanyons(blocks, limits):

    if blocks:  # if blocks not empty list
        # move blocks into left corner
        initblocks = sorted(blocks, key=operator.itemgetter(2,0))
        xshift = initblocks[0][0]  # first block xmin
        yshift = initblocks[0][2]  # first block ymin
        allblocks = blockshift(initblocks, limits, xshift, yshift)
        fullcanyons = generate_canyons(allblocks, limits)
        # shift back
        canyons = blockshift(fullcanyons, limits, -xshift, -yshift)
    else:
        canyons = limits
        fullcanyons = limits
        
    return canyons, fullcanyons


def getncpath(datapath, field="*"):
    # specify field if you do not want all files
    
    expnr = datapath[-4:-1]  # expnr is in last 3 characters before '/' of path
    filename = field + "dump." + expnr + ".nc"
    filelist = os.path.join(datapath, filename)
    files = glob.glob(filelist)
    files.sort()
    
    return files


def getncdata(datapaths):
    
    fieldsdata = []
    for file in datapaths:
        fieldsdata.append(netCDF4.Dataset(file))
        
    return fieldsdata


def getncvars(exp, vardict, ncdata=None):
    # exp needs to be exp object and vardict is dict of dicts
    # e.g. selectvars = {'fielddump' : {'tfield' : 'time', 'u' : 'u'}, 'ztdump' : {'tstats' : 'time', 'ubar' : 'uyt'}}
    
    if type(exp) is Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")
        
    if type(vardict) is dict:
        for subdict in vardict.values():
            if type(subdict) is dict:
                pass
            else:
                raise ValueError("Vardict needs to be a dictionary of dictionaries.\nExample of vardict:\n {'fielddump' : {'tfield' : 'time', 'u' : 'u'}, 'ztdump' : {'tstats' : 'time', 'ubar' : 'uyt'}}")
    else:
        raise ValueError("Vardict needs to be a dictionary of dictionaries.\nExample of vardict:\n {'fielddump' : {'tfield' : 'time', 'u' : 'u'}, 'ztdump' : {'tstats' : 'time', 'ubar' : 'uyt'}}")
        
    
    if ncdata is None:
        try:
            ncfields = exp.ncfields
        except Exception:
            print("Ncfields not found. Specify ncdata or exp.ncfields!")
            pass
    else:
        ncfields = ncdata
            
    
    for fieldname, varlist in vardict.items():
        for fdata in ncfields:
            if fdata.title.startswith(fieldname):
                for varname, varkey in varlist.items():
                    try:
                        exp.add_data(varname, fdata.variables[varkey])
                        # print("Variable exp.%s added from variable %s in %s." % (varname, varkey, fieldname))
                    except Exception:
                        # print("Variable %s could not be found in %s." % (varname, fieldname))
                        pass
    
    return exp
      