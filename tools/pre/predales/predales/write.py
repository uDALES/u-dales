import os
import f90nml


def blocksinp(blocks, filename, parameterlist=None):

    if parameterlist is None:
        z0horiz = 0.01
        z0hhoriz = 0.000067
        Thoriz = 288
        Twest = 288
        Teast = 288
        Tnorth = 288
        Tsouth = 288
        parameterlist = [z0horiz, z0hhoriz, Thoriz, Twest, Teast, Tnorth, Tsouth]
        
    elif len(parameterlist) != 7:
        raise IndexError("Parameter list not the right length.")

    for block in blocks:
        # check that blocks are the right length. For one block only, store as blocks = [[block]].
        if len(block) == 6:
            block.extend(parameterlist)
        else:
            raise IndexError("Blocks not the right length.")

    file = open(filename, 'w')
    file.write("{:12}".format('# Fence location') + '\n')
    file.write("{:100}".format('#  il  iu  jl  ju  kl  ku  z0horiz[m]  z0hhoriz[m]  Thoriz[K]  Twest  Teast  Tnorth  '
                               'Tsouth') + '\n')

    for block in blocks:
        for val in block[0:6]:
            file.write('{:<3.0f} '.format(val))
        for val in block[6:8]:
            file.write('{:<8.6f} '.format(val))
        for val in block[8:12]:
            file.write('{:<4.2f} '.format(val))
        file.write('{:<4.2f}\n'.format(block[12]))

    file.close()
    

def gridinp(xgrid, filename, gridname='x-grid'):

    file = open(filename, 'w')
    file.write("# {:<10}\n".format(gridname))
    file.write("{:<12}\n".format('#'))
    for val in xgrid:
        file.write('{:<20.15f}\n'.format(val))
        
    file.close()
    
    
def lscaleinp(lscales, filename):
        
    if len(lscales[0]) != 10:
        raise IndexError("Lscale list not the right length.")

    file = open(filename, 'w')
    file.write("{:<12}\n".format('# SDBL flow'))
    file.write("{:<60}\n".format('# z uq vq pqx pqy wfls dqtdxls dqtdyls dqtdtls dthlrad'))
    for line in lscales:
        file.write('{:<20.15f} '.format(line[0]))
        for val in line[1:5]:
            file.write('{:<12.6f} '.format(val))
        file.write('{:<15.9f} '.format(line[5]))
        for val in line[6:9]:
            file.write('{:<12.6f} '.format(val))
        file.write('{:<17.12f}\n'.format(line[9]))
        
    file.close()
    
    
def profinp(profs, filename):

    if len(profs[0]) != 6:
        raise IndexError("Prof list not the right length.")

    file = open(filename, 'w')
    file.write("{:<12}\n".format('# SDBL flow'))
    file.write("{:<60}\n".format('# z thl qt u v tke'))
    for line in profs:
        file.write('{:<20.15f} '.format(line[0]))
        for val in line[1:5]:
            file.write('{:<12.6f} '.format(val))
        file.write('{:<12.6f}\n'.format(line[5]))   
        
    file.close() 

    
def scalarinp(scalars, filename):

    if len(scalars[0]) != 5:
        raise IndexError("Scalar list not the right length.")

    file = open(filename, 'w')
    file.write("{:<12}\n".format('# SDBL flow'))
    file.write("{:<60}\n".format('# z sca1 sca2 sca3 sca4'))
    for line in scalars:
        file.write('{:<20.15f} '.format(line[0]))
        for val in line[1:4]:
            file.write('{:<14.10f} '.format(val))
        file.write('{:<14.10f}\n'.format(line[4]))
        
    file.close() 
    
    
def namoptions(parameters, changes, namoptions=None):
    
    # find parameter changes in current parameters
    
    # first dict list, e.g. 'RUN', 'DOMAIN', etc.
    for key, sublist in parameters.items():
    # second dict list, e.g. 'runtime', 'fieldvars', etc.
        for subkey, subval in sublist.items():
    # parameter change list
            for parkey, parval in changes.items():
                if parkey == subkey:
                    parameters[key][subkey] = parval
                    pass

    # change namoptions file
    if namoptions is not None:
        newnamoptions = f90nml.patch(namoptions, parameters, namoptions + "tmp") 
        os.rename(namoptions + "tmp", namoptions)  # rename new namoptions and delete tmp file
    
    return parameters
