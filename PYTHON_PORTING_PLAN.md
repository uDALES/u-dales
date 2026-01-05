# Python Porting Plan for udbase

## Overview
This document outlines the plan to port the MATLAB `udbase` post-processing class and `udgeom` geometry package to Python. The Python implementation will be located in `tools/python/` and provide equivalent functionality for analyzing uDALES simulation outputs.

## File Structure

```
tools/python/
├── udbase.py           # Main post-processing class (equivalent to udbase.m)
├── udgeom.py           # Geometry class (equivalent to +udgeom package)
├── requirements.txt    # Python dependencies
├── README.md          # Usage documentation
└── examples/          # Example scripts
    ├── facets_example.py
    ├── fields_example.py
    └── geometry_example.py
```

## Core Components

### 1. UDBase Class (`udbase.py`)

**Purpose**: Main post-processing class for loading and analyzing uDALES simulation output.

#### Constructor (`__init__`)
- **Inputs**: `expnr` (experiment number), `path` (optional path to experiment directory)
- **Functionality**:
  - Read and parse `namoptions.expnr` file
  - Load grid coordinates from `prof.inp.expnr` or generate uniform grid
  - Load STL geometry if available
  - Load solid point masks (`solid_u/v/w/c.txt`)
  - Load facet data (`facets.inp`, `factypes.inp`, `facetarea.inp`)
  - Load facet section data (`facet_sections_u/v/w/c.txt`, `fluid_boundary_u/v/w/c.txt`)

#### Grid Properties
- `xm`, `xt`: x-coordinates (cell edges and centers)
- `ym`, `yt`: y-coordinates (cell edges and centers)
- `zm`, `zt`: z-coordinates (cell edges and centers)
- `dx`, `dy`: grid spacing in x and y
- `dzm`, `dzt`: variable z-spacing arrays
- `Su`, `Sv`, `Sw`, `Sc`: 3D boolean arrays indicating solid cells

#### Data Loading Methods

**Field Data:**
- `load_field(var=None)` - Load 3D instantaneous fields from `fielddump.expnr.nc`
- `load_stat_xyt(var=None)` - Load slab-averaged statistics from `xytdump.expnr.nc`
- `load_stat_t(var=None)` - Load time-averaged statistics from `tdump.expnr.nc`
- `load_slice(plane, var=None)` - Load 2D slices from `{i,j,k}slicedump.expnr.nc`

**Facet Data:**
- `load_fac_momentum(var=None)` - Load facet momentum data from `fac.expnr.nc`
- `load_fac_eb(var=None)` - Load surface energy balance from `facEB.expnr.nc`
- `load_fac_temperature(var=None)` - Load facet temperatures from `facT.expnr.nc`
- `load_seb()` - Load all surface energy balance terms (returns structure)

#### Facet Analysis Methods
- `assign_prop_to_fac(prop)` - Map facet type properties to individual facets
- `area_average_fac(var, sel=None)` - Area-weighted averaging over facets
- `time_average(var, time, tstart=None, tstop=None)` - Time averaging (static method)
- `area_average_seb(seb)` - Area-average surface energy balance
- `convert_fac_to_field(var, facsec, dz)` - Convert facet variable to 3D density field
- `calculate_frontal_properties()` - Calculate skylines, frontal areas, blockage ratios

#### Visualization Methods
- `plot_fac(var, **kwargs)` - Plot facet variable on 3D mesh
- `plot_fac_type()` - Visualize different surface types

### 2. UDGeom Class (`udgeom.py`)

**Purpose**: Geometry handling and generation for uDALES.

#### UDGeom Class
- `__init__(path=None)` - Constructor, can accept path or triangulation object
- `load(filename)` - Load STL file
- `save(filename)` - Save STL file
- `show(color_buildings=True)` - Visualize geometry with normals

#### Geometry Generation Functions (Static/Module Level)
- `create_flat_surface(xsize, ysize, edgelength)` - Generate flat surface
- `create_canyons(xsize, ysize, B, W, H, shift, edgelength, rotate90=False)` - Generate street canyons
- `create_cubes(xsize, ysize, Hx, Hy, Hz, Cx, Cy, geom_option, edgelength)` - Generate cube arrays
  - Options: 'S' (single), 'AC' (aligned), 'SC' (staggered)
- `create_realistic(building_stl, xsize, ysize, edgelength)` - Add ground to existing buildings

## Implementation Details

### Python Libraries

```python
# Core dependencies
numpy              # Array operations (replaces MATLAB arrays)
xarray             # NetCDF with labeled dimensions
netCDF4            # NetCDF file I/O
scipy              # Sparse matrices, interpolation
trimesh            # STL file handling, mesh operations
matplotlib         # Plotting
pyvista            # Advanced 3D visualization (optional)

# Installation
pip install numpy xarray netCDF4 scipy trimesh matplotlib pyvista
```

### Key Design Decisions

1. **xarray for NetCDF**: Use `xarray.open_dataset()` instead of low-level NetCDF API
   - Provides labeled dimensions
   - Lazy loading for large files
   - Better than MATLAB's NetCDF interface

2. **trimesh for Geometry**: Use `trimesh.Trimesh` objects
   - Direct STL I/O
   - Built-in normals, face centers, areas
   - Boolean operations, contains() method

3. **Properties vs Attributes**: Store all namoptions parameters as instance attributes
   - Use `@property` decorators for computed properties
   - Mimics MATLAB's dynamic properties

4. **Grid Handling**: Store grid arrays as simple numpy arrays
   - Create convenience properties for staggered grids
   - Match MATLAB indexing conventions (adjust for 0-based indexing)

### Facet Data Structures

```python
# Facets structure
self.facs = {
    'typeid': np.array,  # Facet type ID for each facet
    'area': np.array,    # Facet areas
}

# Facet types structure  
self.factypes = {
    'id': np.array,      # Wall type IDs
    'name': list,        # Wall type names
    'lGR': np.array,     # Green facet flag
    'z0': np.array,      # Momentum roughness
    'z0h': np.array,     # Heat roughness
    'al': np.array,      # Albedo
    'em': np.array,      # Emissivity
    'd': np.array,       # Layer thicknesses
    'C': np.array,       # Heat capacities
    'lam': np.array,     # Thermal conductivities
}

# Facet sections structure
self.facsec = {
    'u': {'facid': np.array, 'area': np.array, 'locs': np.array, 'distance': np.array},
    'v': {...},
    'w': {...},
    'c': {...},
}
```

## Implementation Phases

### Phase 1: Core Infrastructure (Priority: High)
**Goal**: Basic simulation loading and grid handling

**Tasks**:
1. Create `udbase.py` skeleton with constructor
2. Implement namoptions parsing
3. Implement grid loading (uniform and stretched)
4. Add file presence checking logic
5. Write unit tests for grid generation

**Deliverable**: Can instantiate UDBase and access grid properties

### Phase 2: Field Data Loading (Priority: High)
**Goal**: Load NetCDF field data

**Tasks**:
1. Implement `load_field()`
2. Implement `load_stat_xyt()`
3. Implement `load_stat_t()`
4. Implement `load_slice()`
5. Add helper method to display NetCDF contents
6. Create example: `examples/fields_example.py`

**Deliverable**: Can load and access field data, matches MATLAB output

### Phase 3: Geometry Handling (Priority: High)
**Goal**: Load and visualize geometry

**Tasks**:
1. Create `udgeom.py` with UDGeom class
2. Implement STL load/save using trimesh
3. Implement `show()` visualization
4. Load geometry in UDBase constructor
5. Load solid point masks
6. Create example: `examples/geometry_example.py`

**Deliverable**: Can load and visualize geometry

### Phase 4: Facet Data (Priority: Medium)
**Goal**: Load and process facet data

**Tasks**:
1. Load facet info files in constructor
2. Implement `load_fac_momentum()`
3. Implement `load_fac_eb()`
4. Implement `load_fac_temperature()`
5. Implement `load_seb()`
6. Implement `assign_prop_to_fac()`
7. Create example: `examples/facets_example.py`

**Deliverable**: Can load facet data

### Phase 5: Facet Analysis (Priority: Medium)
**Goal**: Implement facet analysis methods

**Tasks**:
1. Implement `area_average_fac()`
2. Implement `time_average()` (static method)
3. Implement `area_average_seb()`
4. Implement `convert_fac_to_field()`
5. Implement `calculate_frontal_properties()`
6. Add validation tests against MATLAB outputs

**Deliverable**: Full facet analysis capabilities

### Phase 6: Visualization (Priority: Low)
**Goal**: Plotting and visualization

**Tasks**:
1. Implement `plot_fac()` using matplotlib
2. Implement `plot_fac_type()` 
3. Add optional pyvista support for advanced 3D
4. Create comprehensive visualization examples

**Deliverable**: Can create publication-quality plots

### Phase 7: Geometry Generation (Priority: Low)
**Goal**: Generate synthetic geometries

**Tasks**:
1. Implement `create_flat_surface()`
2. Implement `create_canyons()`
3. Implement `create_cubes()` (single, aligned, staggered)
4. Implement `create_realistic()`
5. Test against MATLAB-generated geometries

**Deliverable**: Can generate geometries for simulations

## Testing Strategy

### Unit Tests
```python
tests/
├── test_udbase.py          # Test UDBase class
├── test_udgeom.py          # Test UDGeom class
├── test_grid.py            # Test grid generation
└── test_facet_analysis.py  # Test facet methods
```

### Regression Tests
- Compare Python outputs with MATLAB outputs for same experiment
- Test with experiments: 065, 110 (from tutorials)
- Verify numerical accuracy (use `np.allclose()`)

### Integration Tests
- Run through complete tutorial workflows in Python
- Verify generated plots match MATLAB versions

## Documentation

### API Documentation
- Docstrings for all public methods (NumPy style)
- Type hints for function signatures
- Example usage in each docstring

### Tutorials
Convert existing MATLAB tutorials to Jupyter notebooks:
- `udbase_tutorial.ipynb`
- `facets_tutorial.ipynb`
- `fields_tutorial.ipynb`
- `geometry_tutorial.ipynb`

### README.md
```markdown
# uDALES Python Tools

## Installation
pip install -r requirements.txt

## Quick Start
from udbase import UDBase

sim = UDBase(expnr=65, path='experiments/065')
u = sim.load_field('u')
print(u.shape)
```

## API Reference

### UDBase Methods Summary

| Method | Purpose | Returns |
|--------|---------|---------|
| `__init__(expnr, path)` | Load simulation | - |
| `load_field(var)` | Load 3D field | xarray.DataArray |
| `load_stat_xyt(var)` | Load slab-averaged stats | xarray.DataArray |
| `load_stat_t(var)` | Load time-averaged stats | xarray.DataArray |
| `load_slice(plane, var)` | Load 2D slice | xarray.DataArray |
| `load_fac_eb(var)` | Load facet energy balance | xarray.DataArray |
| `load_fac_temperature(var)` | Load facet temperatures | xarray.DataArray |
| `load_fac_momentum(var)` | Load facet momentum | xarray.DataArray |
| `load_seb()` | Load full SEB | dict |
| `assign_prop_to_fac(prop)` | Map property to facets | ndarray |
| `area_average_fac(var, sel)` | Area-average facets | ndarray |
| `area_average_seb(seb)` | Area-average SEB | dict |
| `convert_fac_to_field(var, facsec, dz)` | Facet to field | ndarray |
| `calculate_frontal_properties()` | Compute frontal areas | dict |
| `plot_fac(var)` | Plot on mesh | - |
| `plot_fac_type()` | Plot surface types | - |
| `time_average(var, time, t0, t1)` | Time average | ndarray |

### UDGeom Methods Summary

| Method | Purpose | Returns |
|--------|---------|---------|
| `__init__(path)` | Create geometry object | - |
| `load(filename)` | Load STL | - |
| `save(filename)` | Save STL | - |
| `show(color)` | Visualize | - |
| `create_flat_surface(...)` | Generate flat surface | UDGeom |
| `create_canyons(...)` | Generate canyons | UDGeom |
| `create_cubes(...)` | Generate cubes | UDGeom |
| `create_realistic(...)` | Add ground to buildings | UDGeom |

## Compatibility Notes

### MATLAB vs Python Differences

1. **Indexing**: Python uses 0-based indexing, MATLAB uses 1-based
   - Arrays returned match MATLAB shapes but indices differ
   - Solid point files store 1-based indices, subtract 1 when loading

2. **Array Order**: MATLAB uses column-major, NumPy uses row-major by default
   - NetCDF files preserve dimension order
   - Use `.T` for transpose when needed

3. **NetCDF Handling**: xarray is superior to MATLAB
   - Named dimensions instead of numeric indices
   - Lazy loading for memory efficiency
   - Easy subsetting: `data.sel(time=slice(0, 100))`

4. **Visualization**: Matplotlib vs MATLAB plotting
   - Similar API but different defaults
   - Use `from matplotlib import pyplot as plt`
   - Consider pyvista for interactive 3D

## Migration Path

For users transitioning from MATLAB:

```python
# MATLAB
sim = udbase(065, '../experiments/065');
u = sim.load_field('u');

# Python
sim = UDBase(expnr=65, path='../experiments/065')
u = sim.load_field('u')
```

Most method names are identical, just remember:
- Python uses lowercase `True`/`False`
- Optional arguments use `name=value` syntax
- Arrays are 0-indexed

## Success Criteria

The Python implementation is complete when:
1. All methods from `udbase.m` are implemented
2. All methods from `+udgeom` are implemented
3. Unit tests pass with >90% coverage
4. Tutorial examples work in Python
5. Numerical outputs match MATLAB (within floating-point precision)
6. Documentation is complete
7. Can analyze experiments without MATLAB

## Timeline Estimate

- **Phase 1**: 2-3 days
- **Phase 2**: 2-3 days
- **Phase 3**: 2-3 days
- **Phase 4**: 3-4 days
- **Phase 5**: 3-4 days
- **Phase 6**: 2-3 days
- **Phase 7**: 3-4 days
- **Testing & Documentation**: 3-5 days

**Total**: 20-30 days of development time

## Next Steps

1. Review and approve this plan
2. Set up Python development environment
3. Create initial directory structure
4. Begin Phase 1 implementation
5. Iterate with testing and feedback
