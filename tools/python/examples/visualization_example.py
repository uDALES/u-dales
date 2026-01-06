"""
Visualization Tutorial - uDALES Python Tools

This example demonstrates how to visualize facet data using the plot_fac
and plot_fac_type methods.

These methods allow you to display scalar data on the 3D geometry surface
and visualize different surface types.
"""

import numpy as np
import sys
from pathlib import Path

# Add parent directory to path to import udbase
sys.path.insert(0, str(Path(__file__).parent.parent))

from udbase import UDBase

# ============================================================================
# Example 1: Plot facet variable - Net Shortwave Radiation
# ============================================================================
print("=" * 70)
print("Example 1: Plot facet variable (Net Shortwave Radiation)")
print("=" * 70)

# Load simulation
sim = UDBase(101, '../../../examples/101')

# Check if geometry and facet data available
if sim.geom is None:
    print("ERROR: Geometry not loaded. Ensure stl_file is in namoptions.")
    sys.exit(1)

print(f"\nLoaded simulation:")
print(f"  Experiment: {sim.expnr}")
print(f"  Geometry: {sim.geom.n_faces} facets")

# Load surface energy balance
try:
    seb = sim.load_seb()
    K = seb['K']  # Net shortwave radiation
    
    print(f"\nLoaded Net Shortwave Radiation:")
    print(f"  Shape: {K.shape}")
    print(f"  Min: {K.min():.2f} W/m²")
    print(f"  Max: {K.max():.2f} W/m²")
    print(f"  Mean: {K.mean():.2f} W/m²")
    
    # Plot at first time step
    print("\nPlotting net shortwave radiation...")
    print("(Close the plot window to continue)")
    
    sim.plot_fac(
        K[:, 0],
        cmap='hot',
        title='Net Shortwave Radiation (W/m²)',
        colorbar=True
    )
    
except Exception as e:
    print(f"Could not load/plot SEB: {e}")

# ============================================================================
# Example 2: Plot facet variable - Temperature
# ============================================================================
print("\n" + "=" * 70)
print("Example 2: Plot facet variable (Surface Temperature)")
print("=" * 70)

try:
    # Load surface temperature
    Ts = sim.load_fac_temperature('Ts')
    
    print(f"\nLoaded Surface Temperature:")
    print(f"  Shape: {Ts.shape}")
    print(f"  Min: {Ts.min():.2f} K")
    print(f"  Max: {Ts.max():.2f} K")
    print(f"  Mean: {Ts.mean():.2f} K")
    
    # Convert to Celsius for plotting
    Ts_C = Ts - 273.15
    
    print("\nPlotting surface temperature...")
    print("(Close the plot window to continue)")
    
    sim.plot_fac(
        Ts_C[:, 0],
        cmap='coolwarm',
        title='Surface Temperature (°C)',
        colorbar=True
    )
    
except Exception as e:
    print(f"Could not load/plot temperature: {e}")

# ============================================================================
# Example 3: Plot facet variable with custom range
# ============================================================================
print("\n" + "=" * 70)
print("Example 3: Plot with custom colormap range")
print("=" * 70)

try:
    # Load momentum flux
    taux = sim.load_fac_momentum('taux')
    
    print(f"\nLoaded momentum flux (taux):")
    print(f"  Shape: {taux.shape}")
    print(f"  Min: {taux.min():.4f} Pa")
    print(f"  Max: {taux.max():.4f} Pa")
    
    print("\nPlotting with symmetric range [-0.5, 0.5]...")
    print("(Close the plot window to continue)")
    
    sim.plot_fac(
        taux[:, 0],
        cmap='RdBu_r',
        vmin=-0.5,
        vmax=0.5,
        title='Momentum Flux τx (Pa)',
        colorbar=True
    )
    
except Exception as e:
    print(f"Could not load/plot momentum: {e}")

# ============================================================================
# Example 4: Plot time-averaged facet variable
# ============================================================================
print("\n" + "=" * 70)
print("Example 4: Plot time-averaged variable")
print("=" * 70)

try:
    # Load net longwave radiation
    L = seb['L']
    time = seb['time']
    
    print(f"\nTime range: {time[0]:.1f} - {time[-1]:.1f} s")
    
    # Time average
    L_avg = sim.time_average(L, time, tstart=time[-1]-3600, tstop=None)
    
    print(f"\nTime-averaged Net Longwave (last hour):")
    print(f"  Min: {L_avg.min():.2f} W/m²")
    print(f"  Max: {L_avg.max():.2f} W/m²")
    print(f"  Mean: {L_avg.mean():.2f} W/m²")
    
    print("\nPlotting time-averaged net longwave...")
    print("(Close the plot window to continue)")
    
    sim.plot_fac(
        L_avg,
        cmap='seismic',
        title='Net Longwave Radiation - Time Averaged (W/m²)',
        colorbar=True
    )
    
except Exception as e:
    print(f"Could not plot time-averaged: {e}")

# ============================================================================
# Example 5: Plot facet types
# ============================================================================
print("\n" + "=" * 70)
print("Example 5: Plot surface types")
print("=" * 70)

try:
    # Check if facet type data available
    if hasattr(sim, 'factypes') and sim.factypes is not None:
        print(f"\nLoaded facet types:")
        print(f"  Number of types: {len(sim.factypes['id'])}")
        print(f"  Type IDs: {sim.factypes['id']}")
        print(f"  Type names: {sim.factypes['name']}")
        
        print("\nPlotting surface types...")
        print("(Close the plot window to continue)")
        
        sim.plot_fac_type(figsize=(12, 10), show_legend=True)
    else:
        print("\nFacet type data not available.")
        print("This requires facets.XXX and factypes.XXX files.")
        
except Exception as e:
    print(f"Could not plot facet types: {e}")

# ============================================================================
# Example 6: Plot computed facet variable
# ============================================================================
print("\n" + "=" * 70)
print("Example 6: Plot computed variable (Net Radiation)")
print("=" * 70)

try:
    # Compute net radiation = K + L
    Q_net = seb['K'] + seb['L']
    
    print(f"\nComputed Net Radiation:")
    print(f"  Min: {Q_net.min():.2f} W/m²")
    print(f"  Max: {Q_net.max():.2f} W/m²")
    print(f"  Mean: {Q_net.mean():.2f} W/m²")
    
    print("\nPlotting net radiation...")
    print("(Close the plot window to continue)")
    
    sim.plot_fac(
        Q_net[:, 0],
        cmap='plasma',
        title='Net Radiation (K + L) (W/m²)',
        colorbar=True
    )
    
except Exception as e:
    print(f"Could not plot computed variable: {e}")

# ============================================================================
# Example 7: Plot with different colormaps
# ============================================================================
print("\n" + "=" * 70)
print("Example 7: Different colormaps for different variables")
print("=" * 70)

print("""
Common colormap choices:

For radiation (0 to positive):
  - 'hot', 'inferno', 'plasma', 'YlOrRd'
  
For temperature:
  - 'coolwarm', 'RdYlBu_r', 'thermal'
  
For diverging data (negative to positive):
  - 'RdBu_r', 'seismic', 'coolwarm', 'PiYG'
  
For general data:
  - 'viridis' (default), 'cividis', 'twilight'

See matplotlib documentation for full list:
https://matplotlib.org/stable/users/explain/colors/colormaps.html
""")

# ============================================================================
# Example 8: Advanced - Area-averaged and plotted
# ============================================================================
print("\n" + "=" * 70)
print("Example 8: Compare area-averaged with visualization")
print("=" * 70)

try:
    # Get sensible heat flux
    H = seb['H']
    
    # Area average
    H_avg_area = sim.area_average_seb(seb)['H']
    
    print(f"\nSensible Heat Flux:")
    print(f"  Instantaneous mean: {H[:, 0].mean():.2f} W/m²")
    print(f"  Area-weighted average: {H_avg_area[0]:.2f} W/m²")
    print(f"  Difference shows importance of area-weighting")
    
    print("\nPlotting sensible heat flux...")
    print("(Close the plot window to continue)")
    
    sim.plot_fac(
        H[:, 0],
        cmap='RdYlBu_r',
        title='Sensible Heat Flux H (W/m²)',
        colorbar=True
    )
    
    print("\nNote: Visualization shows spatial distribution")
    print("      Area average gives single representative value")
    
except Exception as e:
    print(f"Could not complete comparison: {e}")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("SUMMARY: Facet Visualization")
print("=" * 70)

print("""
Key functionality demonstrated:

1. Plot facet variables:
   sim.plot_fac(var, cmap='viridis', title='Variable')

2. Customize colormap range:
   sim.plot_fac(var, vmin=-1, vmax=1, cmap='RdBu_r')

3. Plot surface types:
   sim.plot_fac_type()

4. Time-averaged visualization:
   var_avg = sim.time_average(var, time, tstart, tstop)
   sim.plot_fac(var_avg)

5. Computed variables:
   Q_net = K + L
   sim.plot_fac(Q_net)

Common colormaps:
- Radiation: 'hot', 'plasma', 'inferno'
- Temperature: 'coolwarm', 'RdYlBu_r'
- Diverging: 'RdBu_r', 'seismic'
- General: 'viridis', 'cividis'

Tips:
- Use vmin/vmax to set consistent ranges for comparison
- Time average before plotting to reduce noise
- Choose colormap appropriate for your data type
- Use plot_fac_type() to verify surface type assignments
""")

print("\nVisualization tutorial complete!")
