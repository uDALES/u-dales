import meshio
import numpy as np
import netCDF4 as nc
import sys
import os

file_Sdir = sys.argv[1]
stl_file = sys.argv[2]

mesh = meshio.read(stl_file)
points = mesh.points
cells = mesh.cells

dataset_Sdir = nc.Dataset(file_Sdir)

# netcdf reverses the order!
tSP = np.array(dataset_Sdir['tSP'])
Sdir = np.array(dataset_Sdir['Sdir'])

os.mkdir("./Sdir")

for n in range(len(tSP)):
    print(n)
    filename = "./Sdir/Sdir_" + str(round(tSP[n])) + ".vtk"
    data_dict = {"Sdir": Sdir[n,:]}
    mesh.cell_data = data_dict
    mesh.write(filename, file_format="vtk")
