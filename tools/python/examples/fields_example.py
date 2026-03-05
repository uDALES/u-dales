"""
Field Data Example - Working with uDALES field outputs

This example demonstrates how to:
1. Load instantaneous 3D field data
2. Load time-averaged statistics
3. Load slab-averaged (xy-plane) statistics
4. Load 2D slices
5. Visualize field data

Based on the fields tutorial from the MATLAB version.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from udbase import UDBase


def example_field_loading():
    """Example 1: Load and display field data."""
    
    print("=" * 70)
    print("Example 1: Loading Field Data")
    print("=" * 70)
    
    # Adjust this path to your experiment
    exp_path = Path("../../experiments/110")
    
    if not exp_path.exists():
        print(f"Experiment not found at {exp_path}")
        print("Please adjust the path to a valid experiment.")
        return None
    
    # Load simulation
    sim = UDBase(expnr=110, path=exp_path)
    print(f"\nLoaded simulation: {sim.expnr}")
    print(f"Grid: {sim.itot} x {sim.jtot} x {sim.ktot}")
    print(f"Domain: {sim.xlen} x {sim.ylen} x {sim.zsize} m")
    
    return sim


def example_instantaneous_fields(sim):
    """Example 2: Load instantaneous 3D fields."""
    
    if sim is None:
        return
    
    print("\n" + "=" * 70)
    print("Example 2: Instantaneous 3D Fields")
    print("=" * 70)
    
    # Display available variables
    print("\nAvailable variables in fielddump:")
    sim.load_field()
    
    # Load velocity field
    print("\nLoading u-velocity field...")
    u = sim.load_field('u')
    
    print(f"u shape: {u.shape}")
    print(f"u dimensions: {u.dims}")
    print(f"u range: [{u.min().values:.3f}, {u.max().values:.3f}] m/s")
    
    # Access specific time and level
    if 'time' in u.dims and 'zt' in u.dims:
        u_slice = u.isel(time=0, zt=10)
        print(f"\nHorizontal slice at first time, level k=10:")
        print(f"  Shape: {u_slice.shape}")
        print(f"  Mean: {u_slice.mean().values:.3f} m/s")
    
    return u


def example_time_averaged_stats(sim):
    """Example 3: Load time-averaged 3D statistics."""
    
    if sim is None:
        return
    
    print("\n" + "=" * 70)
    print("Example 3: Time-Averaged Statistics")
    print("=" * 70)
    
    # Display available variables
    print("\nAvailable variables in tdump:")
    sim.load_stat_t()
    
    # Load time-averaged u-velocity
    print("\nLoading time-averaged u-velocity...")
    u_avg = sim.load_stat_t('u')
    
    print(f"u_avg shape: {u_avg.shape}")
    print(f"u_avg dimensions: {u_avg.dims}")
    
    return u_avg


def example_slab_averaged_stats(sim):
    """Example 4: Load slab-averaged (xy-plane) statistics."""
    
    if sim is None:
        return
    
    print("\n" + "=" * 70)
    print("Example 4: Slab-Averaged Statistics")
    print("=" * 70)
    
    # Display available variables
    print("\nAvailable variables in xytdump:")
    sim.load_stat_xyt()
    
    # Load slab-averaged u-velocity
    print("\nLoading slab-averaged u-velocity...")
    u_xyt = sim.load_stat_xyt('u')
    
    print(f"u_xyt shape: {u_xyt.shape}")
    print(f"u_xyt dimensions: {u_xyt.dims}")
    
    return u_xyt


def example_slices(sim):
    """Example 5: Load 2D slices."""
    
    if sim is None:
        return
    
    print("\n" + "=" * 70)
    print("Example 5: 2D Slice Data")
    print("=" * 70)
    
    # Display available variables in horizontal slice
    print("\nAvailable variables in kslicedump (horizontal):")
    sim.load_slice('k')
    
    # Load horizontal slice of u-velocity
    print("\nLoading horizontal slice of u-velocity...")
    u_slice = sim.load_slice('k', 'u')
    
    print(f"u_slice shape: {u_slice.shape}")
    print(f"u_slice dimensions: {u_slice.dims}")
    
    return u_slice


def visualize_vertical_profile(sim):
    """Example 6: Visualize vertical profile."""
    
    if sim is None:
        return
    
    print("\n" + "=" * 70)
    print("Example 6: Vertical Profile Visualization")
    print("=" * 70)
    
    try:
        # Load slab-averaged statistics
        u_xyt = sim.load_stat_xyt('u')
        
        # Extract profile at first time
        if 'time' in u_xyt.dims:
            u_profile = u_xyt.isel(time=0).values
        else:
            u_profile = u_xyt.values
        
        # Plot vertical profile
        plt.figure(figsize=(6, 8))
        plt.plot(u_profile, sim.zt, 'b-', linewidth=2)
        plt.xlabel('u-velocity (m/s)')
        plt.ylabel('Height (m)')
        plt.title('Vertical Profile of u-velocity')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('vertical_profile.png', dpi=150)
        print("\nSaved: vertical_profile.png")
        plt.close()
        
    except Exception as e:
        print(f"Could not create plot: {e}")


def visualize_horizontal_slice(sim):
    """Example 7: Visualize horizontal slice."""
    
    if sim is None:
        return
    
    print("\n" + "=" * 70)
    print("Example 7: Horizontal Slice Visualization")
    print("=" * 70)
    
    try:
        # Load horizontal slice
        u_slice = sim.load_slice('k', 'u')
        
        # Extract first time and level if needed
        if 'time' in u_slice.dims:
            u_2d = u_slice.isel(time=0)
            if 'zt' in u_2d.dims or 'zm' in u_2d.dims:
                u_2d = u_2d.isel(**{list(u_2d.dims)[0]: 0})
        else:
            u_2d = u_slice
        
        # Create meshgrid for plotting
        X, Y = np.meshgrid(sim.xt, sim.yt)
        
        # Plot
        plt.figure(figsize=(10, 8))
        contour = plt.contourf(X, Y, u_2d.values.T, levels=20, cmap='RdBu_r')
        plt.colorbar(contour, label='u-velocity (m/s)')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.title('Horizontal Slice of u-velocity')
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig('horizontal_slice.png', dpi=150)
        print("\nSaved: horizontal_slice.png")
        plt.close()
        
    except Exception as e:
        print(f"Could not create plot: {e}")


def main():
    """Run all examples."""
    
    print("\n" + "=" * 70)
    print("uDALES Field Data Examples")
    print("=" * 70)
    print("\nThese examples show how to work with field data in Python.")
    print("Make sure you have an experiment with field output files.")
    
    # Example 1: Load simulation
    sim = example_field_loading()
    
    if sim is None:
        print("\nExamples require a valid experiment directory.")
        print("Please adjust the path in example_field_loading().")
        return
    
    # Example 2: Instantaneous fields
    u = example_instantaneous_fields(sim)
    
    # Example 3: Time-averaged statistics
    u_avg = example_time_averaged_stats(sim)
    
    # Example 4: Slab-averaged statistics
    u_xyt = example_slab_averaged_stats(sim)
    
    # Example 5: 2D slices
    u_slice = example_slices(sim)
    
    # Example 6: Visualize vertical profile
    visualize_vertical_profile(sim)
    
    # Example 7: Visualize horizontal slice
    visualize_horizontal_slice(sim)
    
    print("\n" + "=" * 70)
    print("✓ All field data examples completed!")
    print("=" * 70)
    print("\nKey takeaways:")
    print("- Use load_field() for instantaneous 3D data")
    print("- Use load_stat_t() for time-averaged 3D data")
    print("- Use load_stat_xyt() for slab-averaged 1D profiles")
    print("- Use load_slice('k') for horizontal slices")
    print("- xarray provides labeled dimensions for easy subsetting")
    print("- Call without var argument to see available variables")


if __name__ == "__main__":
    main()
