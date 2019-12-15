# Post-processing

## NetCDF files

The scripts `da_append.sh`, `da_concatenate.sh` in the directory `tools/utils` can be used to merge the output of netCDF files into a single netCDF file.

### Requirements

- nco

### Concatenate NetCDF files

To concatenate the output of serveral cpus from one simulation to a single file, use:
``` sh
da_concatenate.sh <path-to-exp-outputs>
```
The file `da_merge.sh` is used by this script.

### Append NetCDF files

To append the output files of a simulation to the output files of another simulation, use:
``` sh
da_append.sh <path-to-exp-outputs1> <path-to-exp-outputs2>
```

## Data

The python package `postdales` in `tools/postprocessing` provides some post-processing tools for processing, visualising and storing the data.

### Requirements

- Python 3.5 or above
- Python libraries (see requirements.txt)

### Accessing data

The module `nc.py` provides functionality to read in uDALES simulation netCDF data sets and extract variables and coordinate variables from it. Non-netCDF input files such as the Fortran namelists (namoptions) and uDALES block files can be read in using the module `inputs.py`
The jupyter notebook `uDALES_data_handling.ipynb` shows an example of how to read in simulation setups and output data from a uDALES simulation.

### Data processing

Basic data processing is provided by the module `data.py`. Currently this only includes the functionality to convert block indices together with netCDF variables to dimensionalised blocks data. The module `experiments.py` provides an "Exp" data class to simplify the storing and accessing of data associated with a uDALES simulation.
The jupyter notebook `uDALES_data_handling.ipynb` shows how to store simulation data with the "Exp" class object and how to save the data to a file.

### Data visualisation

Simple plotting functions to plot vertical profiles, 2D and 3D contour field plots, and 3D block layouts are available in `plot.py`.
The jupyter notebook `uDALES_output_visualisation.ipynb` shows an example of how to read in output variables and plot vertical profiles and contour fields with them.
