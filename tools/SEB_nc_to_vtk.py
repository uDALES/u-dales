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

file_emissivity = expdir + "/emissivity.txt"

file_facEB = expdir + "/facEB." + expnr + ".nc"
file_facT = expdir + "/facT." + expnr + ".nc"

emissivity_float = np.genfromtxt(file_emissivity, delimiter=',')

dataset_facEB = nc.Dataset(file_facEB)
dataset_facT = nc.Dataset(file_facT)

# netcdf reverses the order!
tEB = np.array(dataset_facEB['t'])
T = np.array(dataset_facT['T'])
Ts = np.squeeze(T[:,0,:])
H = -np.array(dataset_facEB['hf'])
Knet = np.array(dataset_facEB['netsw'])
Lin = np.array(dataset_facEB['LWin'])
#Lout = 5.67e-8 * np.multiply(emissivity_float, Ts**4)
Lout = np.array(dataset_facEB['LWout'])
Lnet = Lin - Lout

os.mkdir(expdir + "/SEB")

for n in range(len(tEB)):
    print(n)
    filename = expdir + "/SEB/SEB_" + str(round(tEB[n])) + ".vtk"
    data_dict = {"Ts": Ts[n,:], "Knet": Knet[n,:], "H": H[n,:], "Lin": Lin[n,:], "Lout": Lout[n,:]}
    mesh.cell_data = data_dict
    mesh.write(filename, file_format="vtk")
