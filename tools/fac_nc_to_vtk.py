import meshio
import numpy as np
import netCDF4 as nc
import sys
import os

expnr = sys.argv[1]
outdir = os.getenv('DA_WORKDIR')
expdir = outdir + '/' + expnr
stl_file = sys.argv[2]

mesh = meshio.read(stl_file)
points = mesh.points
cells = mesh.cells

file_fac = expdir + "/fac." + expnr + ".nc"

dataset_fac = nc.Dataset(file_fac)

# netcdf reverses the order!
t = np.array(dataset_fac['t'])
tau_x = np.array(dataset_fac['tau_x'])
tau_y = np.array(dataset_fac['tau_y'])
tau_z = np.array(dataset_fac['tau_z'])
pres = np.array(dataset_fac['pres'])
cth = np.array(dataset_fac['cth'])
htc = np.array(dataset_fac['htc'])
os.mkdir(expdir + "/fac")

for n in range(len(t)):
    print(n)
    filename = expdir + "/fac/fac_" + str(round(t[n])) + ".vtk"
    data_dict = {"tau_x": tau_x[n,:], "tau_y": tau_y[n,:], "tau_z": tau_z[n,:], "pres": pres[n,:], "htc": htc[n,:], "cth": cth[n,:]}
    mesh.cell_data = data_dict
    mesh.write(filename, file_format="vtk")

