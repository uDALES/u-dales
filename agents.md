# uDALES Repository Architecture

## Overview

**uDALES** (urban Dutch Atmospheric Large-Eddy Simulation) is a Large Eddy Simulation (LES) software designed for urban airflow, heat transfer, and pollutant dispersion modeling. The repository contains the complete simulation framework including source code, pre-processing tools, post-processing utilities, documentation, tests, and example cases.

**Current Branch**: `vegetation_3D_SEB`  
**Base**: Fortran (26,450 lines), MATLAB (15,136 lines), Python (9,690 lines), C (10,077 lines)  
**License**: GNU General Public License v3.0  
**Website**: https://udales.github.io/u-dales/

---

## Repository Statistics

### Code Distribution by Language

| Language | Files | Lines | Purpose |
|----------|-------|-------|---------|
| Fortran (.f90, .F) | 49 | 27,533 | Core simulation engine |
| Python (.py) | 33 | 9,690 | Pre/post-processing tools |
| MATLAB (.m) | 83 | 15,136 | Legacy processing tools |
| C (.c, .cpp) | 19 | 10,869 | Support libraries |
| Shell (.sh) | 31 | 1,556 | Build & execution scripts |
| Documentation (.md) | 31 | 3,906 | User guides & tutorials |
| CMake (.cmake) | 9 | 2,301 | Build configuration |

### Input/Data Files
- Boundary condition files (.nbc): 86 files, 100,580 lines
- Boundary index files (.nbi): 86 files, 880 lines
- Sample/example data: 28 files, 1,504 lines
- Various case configurations (.001, .002, .024, .101, .102, .201, .949, .950, .999)

---

## Directory Structure

```
u-dales/
├── src/                       # Fortran core simulation source code (40 modules)
├── tools/                     # Pre/post-processing utilities
│   ├── python/               # Python tools and packages
│   │   ├── udbase.py        # Main post-processing class (2,627 lines)
│   │   ├── udgeom/          # Geometry handling package
│   │   ├── udprep/          # Preprocessing package
│   │   ├── tests/           # Python unit tests
│   │   └── examples/        # Tutorial scripts
│   ├── matlab/              # MATLAB processing scripts
│   ├── IBM/                 # Immersed Boundary Method tools
│   ├── SEB/                 # Surface Energy Balance utilities
│   ├── syntheticInflow/     # Inflow generation
│   ├── View3D/              # 3D visualization
│   └── *.sh                 # Various shell scripts
├── examples/                 # Simulation case studies (001-999)
├── tests/                    # Integration & unit tests
├── docs/                     # Documentation & tutorials
├── build/                    # CMake build artifacts
├── 2decomp-fft/             # Parallelization library (submodule)
├── external/                # External dependencies
└── CMakeLists.txt           # Build configuration
```

---

## Core Components

### 1. Simulation Engine (src/)

The simulation engine consists of **40 Fortran modules** organized by functionality:

#### Main Program & Control
- `program.f90` - Main entry point and simulation loop (271 lines)
- `modglobal.f90` - Global variables and parameters (876 lines)
- `modmpi.f90` - MPI parallelization interface
- `tstep.f90` - Time stepping routines

#### Physics Modules
- **Thermodynamics**: `modthermodynamics.f90`
- **Boundary Conditions**: `modboundary.f90`
- **Subgrid Scale**: `modsubgrid.f90`, `modsubgriddata.f90`
- **Forces**: `modforces.f90` (Coriolis, pressure gradients, etc.)
- **Energy Balance**: `modEB.f90` - Surface energy balance
- **Vegetation**: `vegetation.f90`, `modtrees.f90` - 3D vegetation modeling
- **Chemistry**: `modchem.f90` - Chemical species transport
- **Purifiers**: `modpurifiers.f90` - Air purification systems
- **Heat Pumps**: `heatpump.f90` - Heat pump modeling

#### Numerical Methods
- **Advection**: `advection.f90`, `advec_2nd.f90`, `advec_kappa.f90`, `advec_upw.f90`
- **Poisson Solver**: `modpois.f90`
- **Wall Functions**: `wfmneutral.f90`, `wf_uno.f90`

#### Immersed Boundary Method (IBM)
- `modibm.f90` - IBM implementation
- `modibmdata.f90` - IBM data structures
- `initfac.f90` - Facet initialization

#### Input/Output
- `modstartup.f90` - Initialization and restart
- `readinput.f90` - Input file parsing
- `modsave.f90` - Simulation state saving/loading
- `modfielddump.f90` - Field data output
- `modstatsdump.f90` - Statistics output
- `modstat_nc.f90` - NetCDF statistics output
- `modfields.f90` - Field variables management

#### Driver & Inlet
- `moddriver.f90` - Driver simulations for periodic forcing
- `modinlet.f90`, `modinletdata.f90` - Inflow boundary conditions

#### Testing
- `tests.f90` - Unit tests for specific functionalities
- `modchecksim.f90` - Simulation validation checks

---

### 2. Python Tools (tools/python/)

The Python ecosystem provides modern pre/post-processing capabilities:

#### Core Packages

**UDBase** (`udbase.py` - 2,627 lines)
- Main post-processing class for uDALES simulation data
- Loads and analyzes field data, statistics, facet data
- NetCDF file handling with xarray
- Time series analysis
- Facet and tree data processing
- Compatible with legacy MATLAB workflows

**UDGeom** (`udgeom/` - 1,727 lines total)
- `udgeom.py` (814 lines) - Geometry handling class
- STL file I/O using trimesh
- 3D visualization with matplotlib
- `geometry_generation.py` (392 lines) - Procedural geometry creation
- `split_buildings.py` (184 lines) - Building segmentation
- `calculate_outline.py` (141 lines) - Geometry outline extraction
- `view3d.py` (166 lines) - 3D viewing utilities

**UDPrep** (`udprep/` - 3,302 lines total)
- `udprep.py` (254 lines) - Main preprocessing coordinator
- `udprep_radiation.py` (851 lines) - Radiation calculations
- `directshortwave.py` (1,295 lines) - Direct shortwave radiation (with Numba JIT)
- `solar.py` (1,024 lines) - Solar position and radiation algorithms
- `udprep_grid.py` (251 lines) - Grid generation
- `udprep_init.py` (237 lines) - Initial condition generation
- `udprep_vegetation.py` (226 lines) - Vegetation preprocessing
- `udprep_seb.py` (71 lines) - Surface energy balance setup
- `udprep_ibm.py` (40 lines) - IBM preprocessing
- `udprep_scalars.py` (31 lines) - Scalar field initialization
- `udprep_ic.py` (20 lines) - Initial conditions

#### Python Tests (`tests/` - 441 lines)
- `test_directshortwave.py` (155 lines)
- `test_directshortwave_periodic.py` (83 lines)
- `test_solar.py` (102 lines)
- `test_namelist.py` (72 lines)
- `test_view3d.py` (29 lines)

#### Python Examples (`examples/` - 312 lines)
- `udprep_radiation_tutorial.py` (228 lines)
- `udprep_tutorial.py` (53 lines)
- `udprep_vegetation_tutorial.py` (31 lines)

---

### 3. Tests & Validation (tests/)

Comprehensive testing infrastructure ensures code quality and reproducibility:

#### Integration Tests
- Compare simulation results between branches
- Run on reduced timescales to minimize CPU cost
- Automatic execution in CI/CD pipeline
- Output: boxplots of numerical errors

**Test Cases**:
- Whole cases: `tests/cases/103/`
- Patched examples: `001`, `102`, `201`, `501`, `502`

#### Unit Tests
- **TEST_SPARSE_IJK**: Validates IBM sparse coordinate file reading
  - Tests MPI distribution of solid/fluid boundary points
  - Verifies coordinate remapping between global/local frames
  - Run mode: `runmode = 1004`
  - Location: `tests/integration/sparse_ijk/`

- **TEST_TREES_SPARSE_INPUT**: Tree input file validation
  - Location: `tests/integration/tree_input/`

#### Test Execution
```bash
cd tests/regression
python run_tests.py <branch_a> <branch_b> <build_type>
```

---

### 4. Example Simulations (examples/)

Pre-configured case studies demonstrating various simulation capabilities:

| Case | Description |
|------|-------------|
| 001 | Basic urban geometry |
| 002 | Alternative configuration |
| 024 | Specific test case |
| 101 | IBM demonstration |
| 102 | Extended IBM case |
| 201 | Complex urban setup |
| 949 | Advanced configuration |
| 950 | Specialized scenario |
| 999 | Test/validation case |

Each example contains:
- `namoptions.XXX` - Simulation configuration
- Input files (geometry, boundary conditions, initial conditions)
- README documentation (where applicable)

---

### 5. Documentation (docs/)

Comprehensive user guides and tutorials:

#### Main Documentation
- `index.md` - Landing page
- `udales-installation.md` - Installation guide
- `udales-getting-started.md` - Quick start
- `udales-workflow.md` - Complete workflow guide
- `udales-simulation-setup.md` - Setting up simulations

#### Technical Guides
- `udales-2decomp.md` - Parallelization using 2DECOMP&FFT
- `udales-boundary-conditions.md` - Boundary condition setup
- `udales-namoptions-overview.md` - Configuration file reference
- `udales-driver-simulations.md` - Driver mode documentation

#### Tutorials
- **Geometry**: `udales-geometry-tutorial.md`
- **Fields**: `udales-fields-tutorial.md`
- **Facets**: `udales-facets-tutorial.md`
- **UDBase**: `udales-udbase-tutorial.md`
- **Utilities**: `udales-utility-tutorial.md`

#### Pre/Post-Processing
- `udales-pre-processing.md` - Input preparation
- `udales-post-processing.md` - Output analysis
- `udales-copy-inputs.md` - Input file management

#### Reference
- `udales-how-to-cite.md` - Citation guidelines
- `udales-pub-list.md` - Publications using uDALES
- `references.bib` - Bibliography
- `CONTRIBUTING.md` - Contribution guidelines
- `DEVELOP.md` - Developer notes

#### Media Assets
- `assets/` - Images, GIFs, diagrams
- `*_tutorial_media/` - Tutorial-specific media

---

## Build System

### CMake Configuration (CMakeLists.txt - 127 lines)

**Supported Compilers**:
- GNU Fortran (gfortran) - Primary
- Intel Fortran (ifort) - Secondary

**Build Types**:
- `Release` - Optimized for performance (-O3)
- `Debug` - With debugging symbols, bounds checking, FPE traps

**Dependencies**:
- MPI (OpenMPI or Intel MPI)
- NetCDF & NetCDF-Fortran
- 2DECOMP&FFT library (submodule)
- FFTW (via 2DECOMP)

**Build Options**:
```bash
mkdir -p build/release
cd build/release
cmake ../..
make
```

**Compiler Flags**:
- GNU: `-fdefault-real-8 -ffree-line-length-none -ffpe-trap=invalid,zero,overflow -std=f2008`
- Intel: `-r8 -fpe0 -heap-arrays 10`

---

## Parallelization

### 2DECOMP&FFT Library

uDALES uses the **2DECOMP&FFT library** for domain decomposition:

- **Pencil Decomposition**: 2D decomposition of 3D Cartesian domains
- **MPI parallelization**: Each rank processes one pencil
- **Load balancing**: Decomposition grid must factor mesh exactly

**Communication Patterns**:
1. **Data Transposes**: Z→Y, Y→X, X→Y, Y→Z for FFT-based Poisson solver
2. **Halo Exchanges**: Neighbor data exchange at pencil boundaries

---

## Key Capabilities

### Physical Processes
- ✅ Large Eddy Simulation (LES) of turbulent flows
- ✅ Urban boundary layer dynamics
- ✅ Heat transfer and thermodynamics
- ✅ Pollutant/scalar dispersion
- ✅ Surface Energy Balance (SEB)
- ✅ 3D vegetation modeling with drag and thermal effects
- ✅ Chemical species transport
- ✅ Immersed Boundary Method (IBM) for complex geometries
- ✅ Trees and urban vegetation
- ✅ Air purification systems
- ✅ Heat pump modeling

### Numerical Methods
- ✅ Second-order and higher advection schemes
- ✅ Kappa and upwind advection
- ✅ FFT-based Poisson solver for pressure
- ✅ Runge-Kutta 3rd order time integration
- ✅ Wall functions for boundary layers
- ✅ IBM for buildings and obstacles

### Boundary Conditions
- ✅ Periodic boundaries (x, y, z)
- ✅ Wall functions (neutral and stratified)
- ✅ Prescribed inflow/outflow
- ✅ Synthetic turbulence generation
- ✅ Driver simulations for forcing
- ✅ Time-dependent boundaries

### Input/Output
- ✅ NetCDF field dumps (3D+time)
- ✅ Statistics output (profiles, averages)
- ✅ Facet data for building surfaces
- ✅ Restart files for continuation
- ✅ Tree data output
- ✅ Time series at specified locations

---

## Workflow

### Typical Simulation Workflow

1. **Geometry Preparation**
   - Create or import STL geometry
   - Use `udgeom` to process and validate
   - Generate facet files for IBM

2. **Preprocessing** 
   - Set up computational grid with `udprep`
   - Generate initial conditions
   - Calculate radiation maps (direct/diffuse shortwave)
   - Prepare vegetation maps
   - Configure boundary conditions

3. **Configuration**
   - Edit `namoptions.XXX` file
   - Set domain size, resolution, time steps
   - Enable physics modules (SEB, trees, chemistry, etc.)
   - Configure output options

4. **Execution**
   ```bash
   mpiexec -n <NCPU> <BUILD> namoptions.XXX
   ```

5. **Post-processing**
   - Load results with `UDBase` class
   - Analyze fields, statistics, facets
   - Create visualizations
   - Extract time series and profiles

6. **Validation**
   - Run integration tests if code modified
   - Compare with reference solutions
   - Check conservation properties

---

## Development Branches

Current repository has **4 local branches**:

1. **`master`** - Stable release branch
2. **`enhance-udbase-udgeom-integration`** - Improved Python tool integration
3. **`udbase_py`** - Python port of MATLAB udbase
4. **`vegetation_3D_SEB`** (current) - 3D vegetation with Surface Energy Balance

---

## External Dependencies

### Required
- **Fortran Compiler**: gfortran ≥7.0 or Intel Fortran
- **MPI**: OpenMPI, Intel MPI, or MPICH
- **NetCDF**: NetCDF-C and NetCDF-Fortran
- **CMake**: ≥3.9
- **2DECOMP&FFT**: Included as submodule

### Python Tools (Optional but Recommended)
- **Core**: numpy, xarray, netCDF4
- **Geometry**: trimesh, scipy
- **Visualization**: matplotlib, plotly
- **Performance**: numba (for JIT compilation)
- **Development**: pytest, mkdocs

### MATLAB Tools (Legacy - Optional)
- MATLAB ≥R2016b

---

## Testing Infrastructure

### Continuous Integration (CI)
- GitHub Actions workflow
- Automatic testing on every commit
- Builds Debug and Release configurations
- Runs integration tests comparing branches

### Local Testing
```bash
# Install dependencies
conda env create -f environment.yml
conda activate udales

# Build
mkdir -p build/release
cd build/release
cmake ../.. -DCMAKE_BUILD_TYPE=Release
make

# Run tests
cd ../../tests
python run_tests.py master <your_branch> Release
```

---

## File Formats

### Input Files
- `.stl` - Triangulated surface geometry (ASCII or binary STL)
- `namoptions.XXX` - Fortran namelist configuration
- `.txt` - Sparse IBM coordinate files (solid/fluid boundary points)
- Various ASCII data files for boundary conditions, initial conditions

### Output Files
- `.nc` - NetCDF field data (CF-compliant where possible)
- `*.XXX` - Statistics files (time series, profiles)
- `fac*.XXX` - Facet data (surface temperatures, fluxes, etc.)
- Restart files for simulation continuation

---

## Key Algorithms

### Radiation Calculation
- **Direct Shortwave**: Ray tracing through urban geometry
  - Two methods: `facsec` (default) and `cellsec`
  - Numba JIT compilation for performance
  - Supports periodic boundaries and vegetation
  - Configurable ray density for accuracy/speed tradeoff

- **Solar Position**: NREL SPA algorithm implementation
  - High-accuracy sun position calculations
  - Atmospheric refraction corrections

### IBM (Immersed Boundary Method)
- Facet-based representation of solid boundaries
- Ghost cell approach for boundary conditions
- Wall function integration for turbulent flows
- Sparse storage of solid/fluid boundary points
- MPI-aware domain decomposition

### Vegetation Model
- 3D resolved vegetation effects
- Drag force parameterization
- Thermal effects (latent/sensible heat)
- Leaf Area Density (LAD) distribution

---

## Special Features

### Driver Simulations
- Generate periodic inflow conditions
- Pre-calculate forcing terms
- Accelerate spin-up time for large domains

### Synthetic Inflow
- Turbulence generation for inlet boundaries
- Maintains realistic turbulent structures
- Tools in `tools/syntheticInflow/`

### Surface Energy Balance (SEB)
- Wall, roof, and ground temperature prediction
- Radiant, convective, and conductive heat transfer
- Anthropogenic heat sources
- Facet-based discretization

### Multi-scalar Transport
- Simultaneous transport of multiple scalars
- Passive tracers, pollutants, chemical species
- Source/sink terms

---

## Common Use Cases

1. **Urban Heat Island Studies** - Temperature distributions in cities
2. **Pollutant Dispersion** - Traffic emissions, industrial sources
3. **Wind Comfort Analysis** - Pedestrian-level wind speeds
4. **Building Energy** - Heat fluxes, cooling loads
5. **Green Infrastructure** - Effect of trees and vegetation
6. **Street Canyon Dynamics** - Flow and transport in urban canyons
7. **Urban Planning** - Impact of building layouts

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for:
- How to report bugs
- How to request features  
- Pull request guidelines
- Code style conventions

See [DEVELOP.md](DEVELOP.md) for:
- Development environment setup
- Building from source
- Running tests
- Generating documentation

---

## Publications & Citations

uDALES is an active research code with numerous publications. See:
- [How to Cite](docs/udales-how-to-cite.md)
- [Publication List](docs/udales-pub-list.md)
- [References](docs/references.bib)

---

## License

uDALES is free software licensed under the **GNU General Public License v3.0**.

- Original DALES copyright: Delft University of Technology, Wageningen University, Utrecht University, KNMI (1993-2009)
- uDALES additions: Copyright © 2016-2024 the uDALES Team

See [LICENSE.txt](LICENSE.txt) for full license text.

---

## Contact & Support

- **Website**: https://udales.github.io/u-dales/
- **Repository**: https://github.com/uDALES/u-dales
- **Issues**: https://github.com/uDALES/u-dales/issues

---

## Architecture Summary

uDALES is a mature, well-documented scientific computing framework with:
- **Robust core**: Production-quality Fortran simulation engine
- **Modern tools**: Python ecosystem for workflow automation
- **Active development**: Multiple development branches with new features
- **Comprehensive testing**: Integration and unit tests ensure reliability
- **Rich documentation**: Extensive guides, tutorials, and examples
- **Flexible physics**: Modular design enables urban-specific processes
- **Scalable**: MPI parallelization via 2DECOMP for HPC clusters
- **Open source**: GPL-licensed with active community

The codebase represents a significant investment in urban atmospheric modeling, combining classical LES techniques with modern software engineering practices and specialized urban physics parameterizations.
