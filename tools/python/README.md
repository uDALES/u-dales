# uDALES Python Tools

Python implementation of uDALES post-processing and geometry tools.

## Installation

```bash
cd tools/python
pip install -r requirements.txt
```

## Quick Start

```python
from udbase import UDBase

# Load simulation (automatically loads geometry if available)
sim = UDBase(expnr=65, path='experiments/065')

# Access grid
print(sim.xt, sim.yt, sim.zt)

# Load field data
u = sim.load_field('u')
stats = sim.load_stat_xyt('u')

# Load facet data
seb = sim.load_seb()
seb_avg = sim.area_average_seb(seb)

# Visualize geometry
if sim.geom is not None:
    sim.geom.show()
```

## Features

✅ **Field Data Loading**: Load 3D fields, statistics, and slices from NetCDF files  
✅ **Facet Data Loading**: Load momentum, energy balance, and temperature data  
✅ **Facet Analysis**: Area averaging, time averaging, property assignment  
✅ **xarray Integration**: Labeled dimensions for easy data manipulation  
✅ **Geometry Handling**: Load, create, and visualize STL geometries  
⏳ **Visualization**: Plot facet variables on 3D meshes (coming soon)

## Examples

See `examples/` directory for comprehensive tutorials:
- `fields_example.py` - Working with field data (instantaneous, time-averaged, slices)
- `facets_example.py` - Facet analysis (SEB, momentum, temperatures)
- `geometry_example.py` - Geometry loading, properties, and visualization
- `QUICK_REFERENCE.py` - Quick reference for common operations

Run examples:
```bash
cd examples
python fields_example.py
python facets_example.py
```

## Implementation Status

| Phase | Status | Description |
|-------|--------|-------------|
| Phase 1 | ✅ Complete | Core infrastructure, namoptions, grids |
| Phase 2 | ✅ Complete | Field & facet data loading, analysis methods |
| Phase 3 | ✅ Complete | Geometry handling and visualization |
| Phase 4 | ⏳ In Progress | Facet visualization (plot_fac, plot_fac_type) |
| Phase 5-7 | 📋 Planned | Advanced features |

## Documentation

- [PYTHON_PORTING_PLAN.md](../../PYTHON_PORTING_PLAN.md) - Complete implementation plan
- [PHASE1_COMPLETE.md](PHASE1_COMPLETE.md) - Phase 1 summary
- [PHASE2_COMPLETE.md](PHASE2_COMPLETE.md) - Phase 2 summary
- [PHASE3_COMPLETE.md](PHASE3_COMPLETE.md) - Phase 3 summary
- [QUICK_REFERENCE.py](QUICK_REFERENCE.py) - Quick reference guide
