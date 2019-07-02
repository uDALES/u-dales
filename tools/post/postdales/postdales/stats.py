import math
import numpy as np
from postdales import tools


def first_non_masked_entry(listfloats):
    """Function that returns value of lowest index that is not masked."""
    for item in list(listfloats):
        if np.ma.is_masked(item) == False:
            return item


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


def blockfronts_continuous(blocks, z):
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


def blockstats(dimblocks, a0=None, zm=None, dimcanyons=None):
    blockstatistics = {}
    precision = 3
    
    # if dimblocks defined
    try:
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
        
        # block skewness (??) unweighted
        # in forooghi2017 defined as np.sum and not np.mean but i think this is wrong
        if sigmamean > 0:
            skew_mean = np.mean((blockheights-zmean)**3) / sigmamean**3
        else:
            skew_mean = 0
        blockstatistics['skew_mean'] = round(skew_mean, precision)
    
        # block skewness (??) weighted
        if sigmah > 0:
            skewh = np.sum((blockheights-zh)**3) / sigmah**3
        else:
            skewh = 0
        blockstatistics['skew_h'] = round(skewh, precision)

    except Exception:
        print("Error in dimensionalised block statistics.")
        pass
    
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
    
    # if zm defined
    try:
        # zmax index
        # find index of zm that is closest to zmax:
        zmax_index = np.argmin([abs(zi - zmax) for zi in zm])
        blockstatistics['zmax_index'] = int(zmax_index)
        zmax_mindex = zmax_index  # index of upper boarder of block at zm grid
        zmax_tindex = zmax_index - 1  # last index inside the block at zt grid
        blockstatistics['zmax_mindex'] = int(zmax_mindex)
        blockstatistics['zmax_tindex'] = int(zmax_tindex)
        
        #zh index
        # find index of zm that is closest to zh:
        zh_index = np.argmin([abs(zi - zh) for zi in zm])
        blockstatistics['zh_index'] = int(zh_index)
        zh_mindex = zh_index  # index of upper boarder of block at zm grid
        zh_tindex = zh_index - 1  # last index inside the block at zt grid
        blockstatistics['zh_mindex'] = int(zh_mindex)
        blockstatistics['zh_tindex'] = int(zh_tindex)

        # block function area density covered by blocks per height
        blockfunction_lp = blockareas(dimblocks, zm) / a0
        blockstatistics['blockfunction_lp'] = blockfunction_lp

        mask = blockmask(dimblocks, zm)
        blockstatistics['blockmask'] = mask
 
        blockfunction =  blockfunction_lp
        blockstatistics['blockfunction'] = blockfunction
        
        # correction for ISA data to CSA
        csa = 1 - blockfunction
        blockstatistics['csa'] = csa

        # correction for CSA data to ISA
        isa = 1/(1 - blockfunction)
        blockstatistics['isa'] = isa

        # block function frontal area density covered by blocks per height
        # gives continuous function
        blockfunction_lf = blockfronts_continuous(dimblocks, zm) / a0
        blockstatistics['blockfunction_lf'] = blockfunction_lf

        # block function frontal area similar to lp, with full front
        # gives step function
        blockfunction_lf_discrete = blockfronts_discrete(dimblocks, zm) / a0
        blockstatistics['blockfunction_lf_discrete'] = blockfunction_lf_discrete
        
        # volume of canopy air, normalised
        mask1 = mask.copy()
        # should we correct here for lowest level,
        # which is covered by the street surface?
        # correcting here!
        mask1[0] = 0
        volnorm = np.trapz(mask1*csa, zm)
        blockstatistics['canopyvolume'] = volnorm
        
    except Exception:
        print("No vertical axis z defined.")
        pass
    
    # if dimcanyons defined
    try:
        # exclude continous canyons (main roads)
        ndiff = int(len(dimcanyons) - len(dimblocks))
        # canyon widths
        canyonwidths = [(c[1] - c[0]) for c in dimcanyons]
        if ndiff > 1:
            cws = canyonwidths[:-ndiff]
        else:
            cws = canyonwidths
        blockstatistics['canyonwidths'] = cws
        blockstatistics['canyonwidths_incl_mainroads'] = canyonwidths
        
        # average canyon width excluding main roads
        cw = np.mean(cws)
        blockstatistics['cw'] = round(cw, precision)
        
        # average canyon width including main roads
        cw_alt = np.mean(canyonwidths)
        blockstatistics['cw_incl_mainroads'] = round(cw_alt, precision)

        # block widths
        blockwidths = [(b[1] - b[0]) for b in dimblocks]
        blockstatistics['blockwidths'] = blockwidths
        
        # average block widths
        if blockwidths:
            bw = np.mean(blockwidths)
        else:
            bw = 0
        blockstatistics['bw'] = round(bw, precision)
        
        # average block and canyon width R
        bcw = bw + cw
        blockstatistics['bcw'] = round(bcw, precision)
        
        # average height to canyon width ratio H/W
        # use unweighted average heights here??
        hwratio = zmean/cw
        blockstatistics['hw'] = round(hwratio, precision)
                
        # average canyon width to canyon and block width ratio W/R
        wrratio = cw/bcw
        blockstatistics['wr'] = round(wrratio, precision)

    except Exception:
        print("No dimensionalised canyons defined.")
        pass
        
    return blockstatistics


def height_norm(exp):

    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")
        
    try:
        if exp.zh > 0:
            normh = 1/exp.zh
            norm_zmh = 1/exp.zm[exp.zh_mindex]
            norm_zth = 1/exp.zt[exp.zh_tindex]
            
            normh = round(normh, 3)
            norm_zmh = round(norm_zmh, 3)
            norm_zth = round(norm_zth, 3)
            
        else: 
            # just setting it to 1 is the wrong datatype to multiply it with a netcdf variable later on
            normh = np.array(1)
            norm_zmh = np.array(1)
            norm_zth = np.array(1)
            
        exp.add_data('normh', normh)
        exp.add_data('norm_zmh', norm_zmh)
        exp.add_data('norm_zth', norm_zth)

        if exp.zmax > 0:
            normmax = 1/exp.zmax
            norm_zmmax = 1/exp.zm[exp.zmax_mindex]
            norm_ztmax = 1/exp.zt[exp.zmax_tindex]
            
            normmax = round(normmax, 3)
            norm_zmmax = round(norm_zmmax, 3)
            norm_ztmax = round(norm_ztmax, 3)

        else: 
            normmax = np.array(1)
            norm_zmmax = np.array(1)
            norm_ztmax = np.array(1)
            
        exp.add_data('normmax', normmax)
        exp.add_data('norm_zmmax', norm_zmmax)
        exp.add_data('norm_ztmax', norm_ztmax)
        
    except Exception:
        print("Could not calculate height norm factor.")
        print("The following exp data is required:")
        print("zh, zmax, zh_index, zmax_index, zt, zm")
        pass

    return exp


def momentum_fluxes(exp, varstart=0, varstop=None, method="last", index=None):
    
    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")
        
    try:
        trange = slice(varstart, varstop)
        # add t variable to match var length
        try:
            tsteps = exp.tstats[trange]
            exp.add_data('t', tsteps)
        except Exception:
            print("Could not add stats time range.")
            pass

        try:
            # read in original data
            uavs = exp.ubar[trange, :]
            turbfluxes = -exp.upwp[trange, :]
            dispfluxes = -exp.uw[trange, :]
            # nu = 1e-5
            # dudzs = np.gradient(uavs, axis=1)
            # viscfluxes = nu*dudzs
            # this is actually with nu = nu-SGS
            viscfluxes = exp.usgs[trange, :]
            totfluxes = turbfluxes+dispfluxes+viscfluxes
            
            # add original data
            exp.add_data('uavs_isa', uavs)
            exp.add_data('turbfluxes_isa', turbfluxes)
            exp.add_data('dispfluxes_isa', dispfluxes)
            exp.add_data('viscfluxes_isa', viscfluxes)
            exp.add_data('totfluxes_isa', totfluxes)
            
            # set the standard average procedure
            exp.add_data('uavs', uavs)  # ISA
            # CSA
            exp.add_data('turbfluxes', turbfluxes*exp.csa)
            exp.add_data('dispfluxes', dispfluxes*exp.csa)
            exp.add_data('viscfluxes', viscfluxes*exp.csa)
            exp.add_data('totfluxes', totfluxes*exp.csa)

        except Exception:
            print("Could not read in flux data.")
            print("The following exp data is required:")
            print("ubar, upwp, uw, usgs.")
            pass

        if method in ["last", "final"]:
        # use last time step for fluxes
            uav = exp.uavs[-1]
            turbflux = exp.turbfluxes[-1]
            dispflux = exp.dispfluxes[-1]
            viscflux = exp.viscfluxes[-1]
            totflux = exp.totfluxes[-1]
            
        elif method in ["mean", "average"]:
            uav = np.mean(exp.uavs, axis=0)
            turbflux = np.mean(exp.turbfluxes, axis=0)
            dispflux = np.mean(exp.dispfluxes, axis=0)
            viscflux = np.mean(exp.viscfluxes, axis=0)
            totflux = np.mean(exp.totfluxes, axis=0)
            
        elif method in ["index", "specify"]:
            if isinstance(index, int):
                uav = exp.uavs[index]
                turbflux = exp.turbfluxes[index]
                dispflux = exp.dispfluxes[index]
                viscflux = exp.viscfluxes[index]
                totflux = exp.totfluxes[index]
            else:
                print("Error in calculating fluxes in experiment", exp.reference)
                print("Need to specify index for this method.")
                pass
            
        else:
            print("Error in calculating fluxes in experiment", exp.reference)
            raise ValueError("Method of flux calculation could not be found.\nOptions are: 'last' or 'average'.")
        
        # add data of chosen flux profiles
        # copy data instead of referencing it
        # so that if we patch it later, we do not patch the original
        exp.add_data('uav', uav.copy())
        exp.add_data('turbflux', turbflux.copy())
        exp.add_data('dispflux', dispflux.copy())
        exp.add_data('viscflux', viscflux.copy())
        exp.add_data('totflux', totflux.copy())
        # set switch to show that we (re)set fluxes from original data
        exp.fluxes_patched = False
            
    except Exception:
        print("Could not calculate momentum fluxes.")
        pass

    return exp


def patch_momentum_fluxes(exp):
    
    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")
        
    try:
        # patch only if we have not patched the fluxes already
        if exp.fluxes_patched == False:
            fluxes = [exp.totflux, exp.turbflux, exp.dispflux, exp.viscflux]
            for flux in fluxes:
                # masking of fluxes
                fluxmask = np.ma.getmask(flux)
                # indicies where mask applies
                masked_indices = np.ma.where(fluxmask[fluxmask == True])[0]
                # if there is at least one masked entry
                if len(masked_indices) >= 1:
                    i = masked_indices[-1]
                    flux[i] = 0.
                    exp.fluxes_patched = True

    except Exception:
        print("Could not patch masked entries of momentum fluxes.")
        pass

    return exp        


def estimate_ustar(exp):

    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")
        
    try:
        # use forcing at the bottom of domain
        ustar_sqe = first_non_masked_entry(exp.totstress)
        ustar = math.sqrt(ustar_sqe)
        exp.add_data('ustar', round(ustar, 3))

        ustar_sqe_isa = ustar_sqe * first_non_masked_entry(exp.isa)
        ustar_isa = math.sqrt(ustar_sqe_isa)
        exp.add_data('ustar_isa', round(ustar_isa, 3))
        
    except Exception:
        print("Could not estimate friction velocity ustar.")
        print("The following exp data is required:")
        print("totstress, isa")
        pass

    return exp


def flux_norm(exp):

    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")
        
    try:
        ustar_square = (exp.ustar)**2
        normf = 1/ustar_square
        exp.add_data('normf', round(normf, 3))
        
    except Exception:
        print("Could not calculate flux normalising factor.")
        print("The following exp data is required:")
        print("ustar")
        pass

    return exp


def estimate_inflection_point(exp):
    
    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")

    try:
        # only if we have blocks
        if exp.zmax_index > 0:
            # search for inflection point z_e of uav:

            # find local maximum around z_max (+z_max/2 as margin) as candidate for inflection point
            zmaxx = exp.zmax_index + int(round(exp.zmax_index/2))
            zminn = exp.zh_index - int(round(exp.zh_index/2))
            maxxs = np.argmax(np.gradient(exp.uav[slice(zminn, zmaxx)]))
            ze_index = maxxs + zminn
            # take second derivative and check whether it is close to zero at inflection point
            d2udz2 = np.gradient(np.gradient(exp.uav))
            if d2udz2[ze_index] < 0.1:
                exp.add_data('ze_index', ze_index)
                ze = exp.zt[ze_index]
                exp.add_data('ze', ze)
            else:
                exp.add_data('ze_index', 0)
                exp.add_data('ze', 0)
        else:
            exp.add_data('ze_index', 0)
            exp.add_data('ze', 0)
    except Exception:
        print("Could not estimate inflection point.")
        print("The following exp data is required:")
        print("zmax_index, uav, zt.")
        pass
        
    return exp


def estimate_displacement_length(exp, fitstart=0, fitstop=None, rangestart=0.25, rangestop=45):
    
    # the displacement height is estimated by finding the value zd such that uav(z) and ln(z - zd) has the highest correlation. This yields the value for zd for which ln(z - zd) against uav is closest to a straight line.

    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")

    if exp.zt[fitstop] <= exp.zt[fitstart]:
        print("Error in experiment", exp.reference)
        raise ValueError("Wrong fit range! Fitstart must be smaller index than fitstop.")
            
    try:
        # fitrange defines the range of uav(z) and ln(z - zd) values we consider for calculating the correlation.
        fitrange = slice(fitstart, fitstop)
        exp.add_data("zd_fitstart", fitstart)
        if fitstop < 0:
            fitstopp = len(exp.zt)+fitstop
            exp.add_data("zd_fitstop", fitstopp)
        else:
            exp.add_data("zd_fitstop", fitstop)
        # jrange defines the range of zd values we are testing
        jrange = np.arange(rangestart, rangestop, 0.5)  

        # checks correlations between uav(z) and ln(z - zd)
        corrs = []
        corvar1 = exp.uav[fitrange]
        for j in jrange:
            logzj = np.ma.log(exp.zt - j)
            corvar2 = logzj[fitrange]
            corrmatrix = np.ma.corrcoef(corvar1, corvar2)
            corrs.append(corrmatrix[1, 0])
        maxcorr = np.nanmax(corrs)
        # check that correlation is not nan
        if np.logical_not(np.isnan(maxcorr)):
            maxcorr_index = np.nanargmax(corrs)
            zd = jrange[maxcorr_index]
            exp.add_data("zd_correlation", round(maxcorr, 6))
        else:
            zd = 0
            print("Could not find correlation for zd. Returning zd = 0.")
        exp.add_data("zd", round(zd,3))
        
    except Exception:
        print("Could not estimate displacement length.")
        print("The following exp data is required:")
        print("uav, zt.")
        pass
    
    return exp
    
    
def estimate_roughness_length(exp, fitstart=0, fitstop=None):
    
    # with a value for zd we use the logarithmic equation to estimate u* and z0. the equation is: kappa/ustar * uav(z) = ln(z - zd) - ln(z0).
    # A linear function for ln(z - zd) in uav that obeys this relatino has slope kappa/ustar and y-intercept ln(z0). We use a best fit approximation to find this slope and intercept.
    
    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")

    try:
        if exp.zt[fitstop] <= exp.zt[fitstart]:
            raise ValueError("Wrong fit range! Fitstart must be smaller index than fitstop.")
        fitrange = slice(fitstart, fitstop)
        exp.add_data("z0_fitstart", fitstart)
        exp.add_data("z0_fitstop", fitstop)
        fitvar1 = exp.uav[fitrange]
        fitvar2 = np.ma.log(exp.zt[fitrange] - exp.zd)
        # best fit
        kappa = 0.4
        slope, intercept = np.polyfit(fitvar1, fitvar2, 1)
        z0 = np.exp(intercept)
        exp.add_data('z0', round(z0, 3))
        ustar_loglaw = kappa/slope
        exp.add_data('ustar_loglaw', round(ustar_loglaw, 3))
        
    except Exception:
        print("Could not estimate roughness length.")
        print("The following exp data is required:")
        print("uav, zt, zd.")
        pass
    
    return exp


def estimate_loglaw(exp):
    
    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")

    try:
        kappa = 0.4
        slope=kappa/exp.ustar_loglaw
        
        loglaw = np.exp(slope*exp.uav)*exp.z0 + exp.zd
        exp.add_data("loglaw", loglaw)
        
    except Exception:
        print("Could not estimate log law function.")
        print("The following exp data is required:")
        print("uav, z0, zd, ustar_loglaw.")
        pass
    
    return exp


def estimate_pressure_gradient(exp):
    # forcing (constant pressure gradient in x, <dp/dx>)
    
    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")

    try:
        # use least squares to determine slope of tau
        k, d = np.polyfit(exp.zm[exp.zmax_index+2:], exp.totflux[exp.zmax_index+2:], 1)
        pressure_gradient = - k
        exp.add_data('pressure_gradient', pressure_gradient)
            
    except Exception:
        print("Could not estimate pressure gradient force.")
        print("The following exp data is required:")
        print("totflux, zm, zmax_index.")
        pass
    
    return exp


def estimate_pressure_forcing(exp):
    
    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")

    try:
        # masking of fluxes
        fluxmask = np.ma.getmask(exp.totflux)

        # force as function of height
        forcing_csa = (exp.pressure_gradient)*exp.csa
        # set a mask to where other data (dtau/dz, D) is masked
        forcing_csa = np.ma.masked_where(fluxmask, forcing_csa)

        exp.add_data('forcing', forcing_csa)  # CSA as standard

        # Calculate integrated forcing:
        # differentiating  zgrid: dzm = (zt[i] - zt[i - 1]
        dzs = np.diff(exp.zt)
        dzs = np.append(2*exp.zt[0], dzs)

        lpz = 1 - exp.blockfunction
        lpz[0] = 0  # because blockfunction = 1 there

        # manual integration like np.cumsum, 
        # which does not work for stretched grids
        zcumsum = []
        fi = 0
        for lp, dz in zip(lpz, dzs):
            fi += lp*dz
            zcumsum.append(fi)
        lzstar = fi  # total fluid area
        zcumsum = np.array(zcumsum)  # integrated fluid area function of z

        # stress profile
        stressprofile = lzstar - zcumsum

        # expected total shear stress(CSA)
        ets_csa = (exp.pressure_gradient)*stressprofile
        # set a mask to where other data (dtau/dz, D) is masked
        ets_csa = np.ma.masked_where(fluxmask, ets_csa)

        exp.add_data('totstress', ets_csa)

    except Exception:
        print("Could not estimate total shear stres.")
        print("The following exp data is required:")
        print("totflux, pressure_gradient, blockfunction, zt.")
        pass
    
    return exp


def estimate_flux_gradients(exp):
    
    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")

    try:
        # we need to compute flux derivatives manually because 
        # with np.gradient we get one more NAN entry
        # from central difference scheme.
        # also np.gradient might not work with stretched grid
        # implementing forward difference here:   
        
#         dtau = np.diff(exp.totflux)/np.diff(exp.zm)
#         # extend the last entry to be second to last
#         # to get kmax entries, should be zero anyway
#         dtau = np.ma.append(dtau, dtau[-1])
        
#         dturb = np.diff(exp.turbflux)/np.diff(exp.zm)
#         dturb = np.ma.append(dturb, dturb[-1])

#         ddisp = np.diff(exp.dispflux)/np.diff(exp.zm)
#         ddisp = np.ma.append(ddisp, ddisp[-1])

#         dvisc = np.diff(exp.viscflux)/np.diff(exp.zm)
#         dvisc = np.ma.append(dvisc, dvisc[-1])
        
        # use central difference of np.gradient instead,
        # this gives smoother values around zmax
        dzs = np.gradient(exp.zm)
        
        dtau = np.gradient(exp.totflux)/dzs
        # shift down one NAN value and extend the last entry 
        # to be second to last, should be zero anyway       
        dtau = np.roll(dtau, -1)
        dtau[-1] = dtau[-2]

        dturb = np.gradient(exp.turbflux)/dzs
        dturb = np.roll(dturb, -1)
        dturb[-1] = dturb[-2]
        
        ddisp = np.gradient(exp.dispflux)/dzs
        ddisp = np.roll(ddisp, -1)
        ddisp[-1] = ddisp[-2]
        
        dvisc = np.gradient(exp.viscflux)/dzs
        dvisc = np.roll(dvisc, -1)
        dvisc[-1] = dvisc[-2]
        
        # add data
        exp.add_data('totflux_gradient', dtau)
        exp.add_data('turbflux_gradient', dturb)
        exp.add_data('dispflux_gradient', ddisp)
        exp.add_data('viscflux_gradient', dvisc)

    except Exception:
        print("Could not estimate flux gradients.")
        print("The following exp data is required:")
        print("totflux, turbflux, dispflux, viscflux, zm.")
        pass
    
    return exp


def estimate_drag(exp):
    
    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")

    try:
        # hat D and hat F are defined slightly different. they integrate from top to bottom, tau integrates from bottom to top
        
        # hat D = - tau + hat F
        # drag
        drag_csa = - exp.totflux + exp.totstress
        
        # correct for where drag is nonzero above buildings
        drag_csa[exp.zmax_index + 1:] = 0

        # D = dtau/dz + F
        # we store - D here
        dragforce_csa = - exp.totflux_gradient - exp.forcing

        # add data
        exp.add_data('drag', drag_csa)
        exp.add_data('dragforce', dragforce_csa)
        
    except Exception:
        print("Could not estimate drag force.")
        print("The following exp data is required:")
        print("totflux, totflux_gradient, zmax_index,")
        print("totstress, forcing, isa.")
        pass
    
    return exp


def estimate_mass_centre(exp):
    
    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")

    try:    
        # total drag
        drag0_csa = first_non_masked_entry(exp.drag)
        
        if drag0_csa == 0:
            dm_csa = 0
        else:
            # centre of mass Dm
            dm_csa = np.trapz(exp.drag, exp.zm)/drag0_csa
            dm_csa = round(dm_csa, 2)

        exp.add_data('drag_masscentre', dm_csa)
        exp.add_data('dm', dm_csa)
        
    except Exception:
        print("Could not calculate centre of mass.")
        print("The following exp data is required:")
        print("drag, zm.")
        pass
    
    return exp
                

def pressure_norm(exp):

    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")
        
    try:
        if exp.pressure_gradient is not 0:
            normp = abs(1/(exp.pressure_gradient*exp.zsize))
        else: 
            normp = np.array(1)  # just setting it to 1 is the wrong datatype to multiply it with a netcdf variable later on
        exp.add_data('normp', normp)
        
    except Exception:
        print("Could not calculate pressure gradient norm factor.")
        print("The following exp data is required:")
        print("dpdz, zsize")
        pass

    return exp


def estimate_layer_heights(exp):
    
    if type(exp) is tools.Exp:
        pass
    else:
        raise ValueError("Exp needs to be of type Exp object.")
        
    try:
        # first index where disp flux is greater than threshold:
        # boarder between ISL and RSL
        # approximating from above
        zrsl_index = np.where(abs(exp.dispflux_gradient) > 0.00005)[0][-1]
    except Exception:
        zrsl_index = 0

    exp.add_data('zrsl_index', zrsl_index)
    exp.add_data('zrsl', exp.zm[zrsl_index])
    
    try:
        # first index where drag is below threshold
        # boarder between RSL and UCL
        # approximating from above
        zucl_index = np.where(abs(exp.dragforce) > 0.0001)[0][-1]
    except Exception:
        zucl_index = 0

    exp.add_data('zucl_index', zucl_index)
    exp.add_data('zucl', exp.zm[zucl_index])
        
    return exp


def ustar_norm(exps, expref):
    """Define a reference friction velocity from the experiments"""
    for exp in exps:
        exp.add_data('ustar_smooth', expref.ustar)
        exp.add_data('normu', round(1/(expref.ustar**2), 4))
        
    return exps


def calculate_drag_polynomial_estimation(exps, degree=5):
    """Function that calculates a drag estimation based on the correlation of drag to lambda_f(z). Returns one polynomial for all exps."""
    dragstofit = []
    lfzstofit = []
    for exp in exps:
        # ignore experiments with no blocks
        if exp.lf > 0:
            # test for masked entries
            mask = np.ma.getmask(exp.drag)
            # indicies where mask applies
            masked_indices = np.ma.where(mask[mask == True])[0]
            # if there is at least one masked entry
            if len(masked_indices) >= 1:
                start = masked_indices[-1] + 1
            else:
                start = 0

            # cumulated total drag
            drag0 = exp.drag[start]
            exp.add_data("drag0", drag0)

            drag = exp.drag[start:]/drag0
            blockfunction = exp.blockfunction_lf[start:]/exp.lf
            # add data for fitting the polynomial
            dragstofit.extend(drag)
            lfzstofit.extend(blockfunction)

    # now fit polynomial for all exps
    fit = np.polyfit(lfzstofit, dragstofit, degree)
    p = np.poly1d(fit)

    for exp in exps:
        # ignore experiments with no blocks
        if exp.lf > 0:
            # add general polynomial
            exp.add_data("drag_polynomial", p)
            # add drag estimation from polynomial
            zp = exp.blockfunction_lf/exp.lf
            drag_est = p(zp)*exp.drag0

            # set a mask to where drag is masked
            mask = np.ma.getmask(exp.drag)
            drag_est = np.ma.masked_where(mask, drag_est)
            exp.add_data("drag_estimate", drag_est)

    return exps