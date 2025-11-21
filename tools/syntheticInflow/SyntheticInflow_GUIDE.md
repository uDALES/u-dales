# Synthetic Inflow Turbulence Generator - Guide

## Overview

This document describes the synthetic turbulence generation workflow for u-DALES and the performance optimization applied to `modSyntheticInflow.f90`. The turbulence is generated at a (y, z) plane inlet based on vertical profiles of mean velocity, Reynolds stresses, and turbulence length/time scales.

---

## Table of Contents

1. [Workflow](#workflow)
2. [Step 1: Generate Input Profiles](#step-1-generate-input-profiles)
3. [Step 2: Generate Synthetic Turbulence](#step-2-generate-synthetic-turbulence)
4. [Optimization Method](#optimization-method)
5. [References](#references)
6. [Changelog](#changelog)

---

## Workflow

The synthetic inflow generation process consists of two main steps:

1. **Generate input profiles**: Create vertical profiles of mean velocity, Reynolds stresses, and turbulence scales
2. **Generate synthetic turbulence**: Run the Fortran code to produce turbulent inflow boundary conditions

---

## Step 1: Generate Input Profiles

### Using write_Reynolds_stress.m

You can generate the required vertical profiles using the MATLAB script `write_Reynolds_stress.m`. There are two options:

1. **Extract from existing simulation**: Use statistics from a previous u-DALES run (mode 1)
2. **Customize profiles**: Define your own target profiles (mode 2)

### Output Files

The MATLAB script generates the following input files:

- `length_time_scales_u.txt` - Length and time scales for u-velocity
- `length_time_scales_v.txt` - Length and time scales for v-velocity  
- `length_time_scales_w.txt` - Length and time scales for w-velocity
- `Reynolds_stress_profiles_velocity.txt` - Mean velocity and Reynolds stresses
- (Optional) `Reynolds_stress_profiles_temp.txt` - Temperature statistics
- (Optional) `Reynolds_stress_profiles_moist.txt` - Moisture statistics

These files define the vertical (z-direction) profiles that characterise the turbulent flow at the inlet.

COPY THE ABOVE FILES TO YOUR /experiments/{expnr}, then go step 2.

---

## Step 2: Generate Synthetic Turbulence

### Method Overview

The `modSyntheticInflow.f90` code implements the **Xie & Castro (2008)** synthetic turbulence generation method. This method generates realistic turbulent fluctuations by:

1. Creating a random field with proper spatial correlations
2. Filtering the random field using prescribed length scales
3. Scaling the result to match target Reynolds stresses
4. Applying mass flux correction to ensure divergence-free inflow

### Running the Generator

**Option 1: Using the shell script**
```bash
bash u-dales/tools/generate_synthetic_inflow.sh experiments/{expnr}
```

**Note**: You must run `write_input` beforehand, as this needs `prof.inp.{expnr}`.

**Option 2: Submit SLURM job (recommended for supercomputers)**

Use the provided `job_synthflow` script. Make sure to configure the settings appropriately for your machine (number of cores, partition, time limit, etc.).

### Output Files

The generator produces driver files for the u-DALES simulation:

- `tdriver_000.{expnr}` - Time stamps
- `udriver_XXX.{expnr}` - U-velocity inflow for each processor
- `vdriver_XXX.{expnr}` - V-velocity inflow for each processor
- `wdriver_XXX.{expnr}` - W-velocity inflow for each processor
- (Optional) `hdriver_XXX.{expnr}` - Temperature inflow
- (Optional) `qdriver_XXX.{expnr}` - Moisture inflow

Where `XXX` is the processor ID (000 to nprocy-1).

---

## Optimization Method
**using FFT to replace the forloop in modSyntheticInflow.f90, e.g., in calc_psi**

### Setting Number of Threads

**Via environment variable** (recommended): `export OMP_NUM_THREADS=128`
**Via SLURM**: Set in the submission script using `OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}`


## References

### Scientific Papers

1. **Xie, Z.-T., & Castro, I. P. (2008)**. "Efficient generation of inflow conditions for large eddy simulation of street-scale flows." *Flow, Turbulence and Combustion*, 81, 449-470.
   - Original synthetic turbulence generation method

2. **Kim, Y., Castro, I. P., & Xie, Z.-T. (2013)**. "Divergence-free turbulence inflow conditions for large-eddy simulations with incompressible flow solvers." *Computers & Fluids*, 84, 56-68.
   - Mass flux correction methodology

### Implementation References

3. **PALM Model System** - Parallelized Large-Eddy Simulation Model
   - Source of the original length/time scale calculation
   - URL: https://palm.muk.uni-hannover.de/

4. **u-DALES** - Urban Dutch Atmospheric Large-Eddy Simulation
   - Main simulation code using this inflow generator
   - URL: https://github.com/uDALES/u-dales

### Algorithm Theory

5. **Separable Filters in Image Processing**
   - Theoretical background on filter decomposition
   - Any standard image processing textbook (e.g., Gonzalez & Woods)

---

## Changelog

### Version 2.0 (November 2025)
- ✅ Implemented separated filtering algorithm
- ✅ Achieved 100-200× speedup for large length scales
- ✅ Increased OpenMP chunk size from 8 to 32
- ✅ Preserved original version as commented reference
- ✅ Added comprehensive documentation

### Version 1.0 (Original)
- Original 4-nested loop implementation
- Based on Xie & Castro (2008) method

---

## License

This code is part of the u-DALES model system and is distributed under the same license as u-DALES.

---

## Contact

For questions or issues related to this optimization:

1. Check the u-DALES GitHub repository issues
2. Consult the u-DALES user manual
3. Contact the u-DALES development team

---

## Acknowledgments

- Original synthetic inflow method: Xie & Castro (2008)
- u-DALES development team
- PALM model developers (for length/time scale methodology)

---

**Last Updated**: November 21, 2025
