"""
Advanced Analysis Tutorial - uDALES Python Tools

This example demonstrates advanced facet analysis methods:
- convert_fac_to_field(): Convert facet data to 3D volumetric fields
- calculate_frontal_properties(): Compute skylines, frontal areas, blockage ratios

These methods are useful for urban canopy analysis and flow characterization.
"""

import numpy as np
import sys
from pathlib import Path

# Add parent directory to path to import udbase
sys.path.insert(0, str(Path(__file__).parent.parent))

from udbase import UDBase

# ============================================================================
# Example 1: Calculate frontal properties
# ============================================================================
print("=" * 70)
print("Example 1: Calculate frontal properties")
print("=" * 70)

# Load simulation with geometry and facet data
sim = UDBase(101, '../../../examples/101')

# Check if required data available
if sim.geom is None:
    print("ERROR: Geometry not loaded. Ensure stl_file is in namoptions.")
    sys.exit(1)

if not hasattr(sim, 'facsec') or sim.facsec is None:
    print("ERROR: Facet section data not loaded.")
    print("This requires facet_sections_(u,v,w,c) and fluid_boundary_(u,v,w,c) files.")
    sys.exit(1)

print(f"\nSimulation setup:")
print(f"  Experiment: {sim.expnr}")
print(f"  Domain: {sim.xlen:.1f} x {sim.ylen:.1f} x {sim.zsize:.1f} m")
print(f"  Grid: {sim.itot} x {sim.jtot} x {sim.ktot}")
print(f"  Geometry: {sim.geom.n_faces} facets")

# Calculate frontal properties
print("\nCalculating frontal properties...")
props = sim.calculate_frontal_properties()

print(f"\nResults:")
print(f"  Frontal area (x-direction): {props['Afx']:.1f} m²")
print(f"  Frontal area (y-direction): {props['Afy']:.1f} m²")
print(f"  Blockage ratio (x): {props['brx']:.3f} ({100*props['brx']:.1f}%)")
print(f"  Blockage ratio (y): {props['bry']:.3f} ({100*props['bry']:.1f}%)")

# ============================================================================
# Example 2: Analyze skyline profiles
# ============================================================================
print("\n" + "=" * 70)
print("Example 2: Analyze skyline profiles")
print("=" * 70)

skylinex = props['skylinex']  # Shape: (jtot, ktot)
skyliney = props['skyliney']  # Shape: (itot, ktot)

print(f"\nSkyline shapes:")
print(f"  skylinex: {skylinex.shape} (j, k)")
print(f"  skyliney: {skyliney.shape} (i, k)")

# Calculate blockage at each height
blockage_x_vs_z = np.sum(skylinex, axis=0) / sim.jtot
blockage_y_vs_z = np.sum(skyliney, axis=0) / sim.itot

print(f"\nBlockage fraction at each height:")
print(f"  Bottom (k=0): x={blockage_x_vs_z[0]:.3f}, y={blockage_y_vs_z[0]:.3f}")
print(f"  Mid (k={sim.ktot//2}): x={blockage_x_vs_z[sim.ktot//2]:.3f}, y={blockage_y_vs_z[sim.ktot//2]:.3f}")
print(f"  Top (k={sim.ktot-1}): x={blockage_x_vs_z[-1]:.3f}, y={blockage_y_vs_z[-1]:.3f}")

# ============================================================================
# Example 3: Visualize skylines
# ============================================================================
print("\n" + "=" * 70)
print("Example 3: Visualize skylines")
print("=" * 70)

try:
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Skyline in x-direction (y-z plane)
    im1 = axes[0].imshow(skylinex.T, origin='lower', aspect='auto', cmap='binary')
    axes[0].set_xlabel('j (grid points)')
    axes[0].set_ylabel('k (grid points)')
    axes[0].set_title('Skyline in x-direction\n(View from x, looking at y-z plane)')
    plt.colorbar(im1, ax=axes[0], label='Blocked (1) / Open (0)')
    
    # Skyline in y-direction (x-z plane)
    im2 = axes[1].imshow(skyliney.T, origin='lower', aspect='auto', cmap='binary')
    axes[1].set_xlabel('i (grid points)')
    axes[1].set_ylabel('k (grid points)')
    axes[1].set_title('Skyline in y-direction\n(View from y, looking at x-z plane)')
    plt.colorbar(im2, ax=axes[1], label='Blocked (1) / Open (0)')
    
    plt.tight_layout()
    plt.savefig('skylines.png', dpi=150)
    print("\nSaved skyline visualization to 'skylines.png'")
    print("(Close the plot window to continue)")
    plt.show()
    
except ImportError:
    print("\nmatplotlib not available - skipping visualization")

# ============================================================================
# Example 4: Plot blockage vs height
# ============================================================================
print("\n" + "=" * 70)
print("Example 4: Plot blockage vs height")
print("=" * 70)

try:
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot blockage fraction vs height
    ax.plot(blockage_x_vs_z, sim.zt, 'b-', linewidth=2, label='x-direction')
    ax.plot(blockage_y_vs_z, sim.zt, 'r-', linewidth=2, label='y-direction')
    
    ax.set_xlabel('Blockage Fraction')
    ax.set_ylabel('Height (m)')
    ax.set_title('Vertical Profile of Blockage')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_xlim([0, 1])
    
    plt.tight_layout()
    plt.savefig('blockage_profile.png', dpi=150)
    print("\nSaved blockage profile to 'blockage_profile.png'")
    print("(Close the plot window to continue)")
    plt.show()
    
except ImportError:
    print("\nmatplotlib not available - skipping plot")

# ============================================================================
# Example 5: Convert facet variable to field
# ============================================================================
print("\n" + "=" * 70)
print("Example 5: Convert facet variable to field")
print("=" * 70)

# Create a simple facet variable (e.g., constant = 1 for all facets)
ones = np.ones(sim.geom.n_faces)

print("\nConverting constant facet variable to field...")
field = sim.convert_fac_to_field(ones)

print(f"\nField properties:")
print(f"  Shape: {field.shape}")
print(f"  Min: {field.min():.6f}")
print(f"  Max: {field.max():.6f}")
print(f"  Mean: {field.mean():.6f}")
print(f"  Sum: {field.sum():.2f}")

print("\nThis field represents the density of facet area per unit volume.")
print("High values indicate cells with lots of surface area (e.g., near buildings).")

# ============================================================================
# Example 6: Convert real facet data to field
# ============================================================================
print("\n" + "=" * 70)
print("Example 6: Convert real facet data to field")
print("=" * 70)

try:
    # Load surface temperature
    Ts = sim.load_fac_temperature('Ts')
    
    print(f"\nLoaded surface temperature:")
    print(f"  Shape: {Ts.shape}")
    print(f"  Min: {Ts[:, 0].min():.2f} K")
    print(f"  Max: {Ts[:, 0].max():.2f} K")
    
    # Convert to field (first time step)
    T_field = sim.convert_fac_to_field(Ts[:, 0])
    
    print(f"\nTemperature field:")
    print(f"  Shape: {T_field.shape}")
    print(f"  Non-zero cells: {np.count_nonzero(T_field)}")
    print(f"  Max: {T_field.max():.2f} K")
    
    # This field shows temperature density - where surfaces are and their temps
    
except Exception as e:
    print(f"\nCould not load temperature data: {e}")
    print("This is OK if temperature files don't exist.")

# ============================================================================
# Example 7: Visualize converted field (horizontal slice)
# ============================================================================
print("\n" + "=" * 70)
print("Example 7: Visualize converted field")
print("=" * 70)

try:
    import matplotlib.pyplot as plt
    
    # Take a horizontal slice at mid-height
    k_mid = sim.ktot // 2
    field_slice = field[:, :, k_mid]
    
    print(f"\nHorizontal slice at k={k_mid} (z={sim.zt[k_mid]:.2f} m)")
    print(f"  Non-zero cells: {np.count_nonzero(field_slice)}")
    print(f"  Max value: {field_slice.max():.3f}")
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    im = ax.imshow(field_slice.T, origin='lower', cmap='hot', aspect='equal')
    ax.set_xlabel('i (x-direction)')
    ax.set_ylabel('j (y-direction)')
    ax.set_title(f'Facet Area Density at z={sim.zt[k_mid]:.2f} m')
    plt.colorbar(im, ax=ax, label='Density [1/m]')
    
    plt.tight_layout()
    plt.savefig('field_density.png', dpi=150)
    print("\nSaved field density visualization to 'field_density.png'")
    print("(Close the plot window to continue)")
    plt.show()
    
except ImportError:
    print("\nmatplotlib not available - skipping visualization")

# ============================================================================
# Example 8: Compare different facet sections
# ============================================================================
print("\n" + "=" * 70)
print("Example 8: Compare different facet sections (c vs u)")
print("=" * 70)

# Convert using c-grid (default)
field_c = sim.convert_fac_to_field(ones, facsec=sim.facsec['c'])

# Convert using u-grid
field_u = sim.convert_fac_to_field(ones, facsec=sim.facsec['u'])

print(f"\nField comparison:")
print(f"  c-grid non-zero cells: {np.count_nonzero(field_c)}")
print(f"  u-grid non-zero cells: {np.count_nonzero(field_u)}")
print(f"  c-grid sum: {field_c.sum():.2f}")
print(f"  u-grid sum: {field_u.sum():.2f}")

print("\nDifferent grids capture different aspects of the surface geometry.")

# ============================================================================
# Example 9: Calculate plan area density
# ============================================================================
print("\n" + "=" * 70)
print("Example 9: Calculate plan area density")
print("=" * 70)

# Plan area density λp is the ratio of building plan area to total area
# We can estimate it from the skyline

plan_area_ratio_x = np.sum(skylinex[:, 0]) / (sim.jtot)
plan_area_ratio_y = np.sum(skyliney[:, 0]) / (sim.itot)

print(f"\nPlan area ratio (from ground-level skyline):")
print(f"  From x-direction view: {plan_area_ratio_x:.3f}")
print(f"  From y-direction view: {plan_area_ratio_y:.3f}")
print(f"  Average: {(plan_area_ratio_x + plan_area_ratio_y)/2:.3f}")

# ============================================================================
# Example 10: Urban canopy parameters summary
# ============================================================================
print("\n" + "=" * 70)
print("Example 10: Urban canopy parameters summary")
print("=" * 70)

# Calculate various urban canopy parameters
total_plan_area = sim.xlen * sim.ylen
building_footprint = plan_area_ratio_x * total_plan_area

# Estimate mean building height from frontal area and footprint
if building_footprint > 0:
    mean_height_x = props['Afx'] / (building_footprint ** 0.5 * sim.xlen ** 0.5)
    mean_height_y = props['Afy'] / (building_footprint ** 0.5 * sim.ylen ** 0.5)
else:
    mean_height_x = mean_height_y = 0

print(f"\nUrban canopy parameters:")
print(f"  Domain area: {total_plan_area:.1f} m²")
print(f"  Building footprint (est.): {building_footprint:.1f} m²")
print(f"  Plan area ratio λp: {plan_area_ratio_x:.3f}")
print(f"  Frontal area (x): {props['Afx']:.1f} m²")
print(f"  Frontal area (y): {props['Afy']:.1f} m²")
print(f"  Blockage ratio (x): {props['brx']:.3f}")
print(f"  Blockage ratio (y): {props['bry']:.3f}")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("SUMMARY: Advanced Facet Analysis")
print("=" * 70)

print("""
Key functionality demonstrated:

1. Calculate frontal properties:
   props = sim.calculate_frontal_properties()
   - Frontal areas (Afx, Afy)
   - Blockage ratios (brx, bry)
   - Skyline profiles (skylinex, skyliney)

2. Convert facet to field:
   field = sim.convert_fac_to_field(var)
   - Creates 3D volumetric density field
   - Can specify grid (c, u, v, w)
   - Useful for volume integrals

3. Applications:
   - Urban canopy characterization
   - Flow resistance estimation
   - Building morphology analysis
   - Surface-volume coupling

Frontal Properties:
- Frontal area: Projected area perpendicular to flow
- Blockage ratio: Fraction of domain cross-section blocked
- Skyline: Binary indicator of blocked cells
- Important for drag and flow resistance

Convert to Field:
- Distributes facet data into volume
- Density = area * value / cell_volume
- Non-zero only in cells touching surfaces
- Enables volume-averaged quantities

Physical Interpretation:
- High blockage → Strong flow resistance
- Skyline shows obstruction patterns
- Field density indicates surface concentration
""")

print("\nAdvanced analysis tutorial complete!")
