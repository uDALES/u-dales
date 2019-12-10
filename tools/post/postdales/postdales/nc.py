import os
import netCDF4
import glob


def getncpath(datadir, field="*"):
    """Get path to netcdf files in datadir. Use
    keyword 'field' to specify only one field."""
    filename = field + "dump.???.nc"  # get only merged files, not from all cores
    filelist = os.path.join(datadir, filename)
    files = glob.glob(filelist)
    files.sort()
    
    return files


def getncdata(datapaths):
    """Get the netcdf data set from files."""    
    fieldsdata = []
    for file in datapaths:
        fieldsdata.append(netCDF4.Dataset(file))
        
    return fieldsdata


def variables_available(list_of_ncfiles):
    """Lists which variables are available from
    the netcdf data sets."""
    print("Available variables:\n")
    for ncfile in list_of_ncfiles:
        variable_names = list(ncfile.variables.keys())
        print(ncfile.title, ":\n", variable_names)
        
    return


def variables_available_1D(list_of_ncfiles):
    """Lists which variables with 1 spatial coordinate
    (vertical profiles) are available from netcdf data sets."""
    print("Available 1D vertical profile variables:\n")
    verticalaxes = ['zt', 'zm']
    for ncfile in list_of_ncfiles:
        print("\n", ncfile.title)
        ncvariables = ncfile.variables.values()
        for variable in ncvariables:
            if any(axis in verticalaxes for axis in variable.dimensions):
                if variable.ndim == 2:  # only time and z dimensions
                    if not variable.name in verticalaxes:  # do not print axis variables themselves
                        print(variable.name)
    return


def variables_available_2D(list_of_ncfiles):
    """Lists which variables with 2 spatial coordinates
    are available from the netcdf data sets."""
    print("Available 2D field variables:")
    fieldaxes = ['xt', 'xm', 'yt', 'ym']
    for ncfile in list_of_ncfiles:
        print("\n", ncfile.title)
        ncvariables = ncfile.variables.values()
        for variable in ncvariables:
            if any(axis in fieldaxes for axis in variable.dimensions):
                if variable.ndim == 3:  # only time and 2 coordinates
                    print(variable.name)
    return


def variables_available_3D(list_of_ncfiles):
    """Lists which variables with 3 spatial coordinates
    are available from the netcdf data sets."""
    print("Available 3D field variables:")
    fieldaxes = ['xt', 'xm', 'yt', 'ym']
    for ncfile in list_of_ncfiles:
        print("\n", ncfile.title)
        ncvariables = ncfile.variables.values()
        for variable in ncvariables:
            if any(axis in fieldaxes for axis in variable.dimensions):
                if variable.ndim == 4:  # time and 3 spatial coordinates
                    print(variable.name)
    return


def select_ncfield(list_of_ncfields, fieldname):
    """Select one netcdf data set from list."""
    for ncfield in list_of_ncfields:
        if ncfield.title.startswith(fieldname):
            return ncfield


def get_coordinates(ncfield):
    """Get coordinates from netcdf data set."""
    coordinates = {}
    for coordinate_name in ncfield.dimensions.keys():
        coordinate = ncfield.variables[coordinate_name]
        coordinates[coordinate_name] = coordinate
    return coordinates

        
def select_variable(variable_name, ncfield):
    """Select a variable from netcdf data set."""
    variable = ncfield.variables[variable_name]
    return variable


def select_coordinates(variable_name, ncfield):
    """Select the coordinates of a variable from 
    the netcdf data set."""
    coords = []
    variable = ncfield.variables[variable_name]
    for dimension in variable.dimensions:
        coordinate = ncfield.variables[dimension]
        coords.append(coordinate)
    return coords
