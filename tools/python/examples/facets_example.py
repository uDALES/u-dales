"""
Facet Data Example - Working with uDALES surface data

This example demonstrates how to:
1. Load facet momentum data (pressure, shear stress)
2. Load surface energy balance data
3. Load facet temperatures
4. Assign facet type properties
5. Perform area averaging
6. Calculate area-averaged surface energy balance

Based on the facets tutorial from the MATLAB version.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from udbase import UDBase


def example_facet_loading():
    """Example 1: Load simulation with facets."""
    
    print("=" * 70)
    print("Example 1: Loading Simulation with Facets")
    print("=" * 70)
    
    # Adjust this path to your experiment with facets
    exp_path = Path("../../experiments/065")
    
    if not exp_path.exists():
        print(f"Experiment not found at {exp_path}")
        print("Please adjust the path to a valid experiment with facets.")
        return None
    
    # Load simulation
    sim = UDBase(expnr=65, path=exp_path)
    print(f"\nLoaded simulation: {sim.expnr}")
    print(f"Grid: {sim.itot} x {sim.jtot} x {sim.ktot}")
    
    if hasattr(sim, 'nfcts'):
        print(f"Number of facets: {sim.nfcts}")
    
    # Check what facet data is available
    if sim.facs:
        print("\nFacet data loaded:")
        for key in sim.facs:
            print(f"  - {key}")
    
    if sim.factypes:
        print("\nFacet types loaded:")
        print(f"  {len(sim.factypes['id'])} wall types defined")
    
    return sim


def example_facet_properties(sim):
    """Example 2: Work with facet type properties."""
    
    if sim is None or not sim.factypes:
        print("No facet type data available")
        return
    
    print("\n" + "=" * 70)
    print("Example 2: Facet Type Properties")
    print("=" * 70)
    
    # Display facet types
    print("\nWall types in simulation:")
    for i, (id_, name) in enumerate(zip(sim.factypes['id'], sim.factypes['name'])):
        al = sim.factypes['al'][i]
        em = sim.factypes['em'][i]
        print(f"  Type {int(id_):3d}: {name:25s} (albedo={al:.3f}, emissivity={em:.3f})")
    
    # Assign properties to individual facets
    print("\nAssigning albedo to each facet...")
    albedo = sim.assign_prop_to_fac('al')
    print(f"  Albedo array shape: {albedo.shape}")
    print(f"  Albedo range: [{albedo.min():.3f}, {albedo.max():.3f}]")
    
    # Assign emissivity
    emissivity = sim.assign_prop_to_fac('em')
    print(f"\nEmissivity range: [{emissivity.min():.3f}, {emissivity.max():.3f}]")
    
    return albedo, emissivity


def example_facet_momentum(sim):
    """Example 3: Load facet momentum data."""
    
    if sim is None:
        return
    
    print("\n" + "=" * 70)
    print("Example 3: Facet Momentum Data")
    print("=" * 70)
    
    try:
        # Display available variables
        print("\nAvailable variables in fac.nc:")
        sim.load_fac_momentum()
        
        # Load pressure on facets
        print("\nLoading facet pressure...")
        pf = sim.load_fac_momentum('pf')
        print(f"  Shape: {pf.shape}")
        print(f"  Dimensions: {pf.dims}")
        
        return pf
        
    except FileNotFoundError:
        print("Momentum data (fac.nc) not available")
        return None


def example_energy_balance(sim):
    """Example 4: Load surface energy balance data."""
    
    if sim is None:
        return
    
    print("\n" + "=" * 70)
    print("Example 4: Surface Energy Balance")
    print("=" * 70)
    
    try:
        # Display available variables
        print("\nAvailable variables in facEB.nc:")
        sim.load_fac_eb()
        
        # Load individual terms
        print("\nLoading energy balance terms...")
        H = sim.load_fac_eb('hf')  # Sensible heat flux
        K = sim.load_fac_eb('netsw')  # Net shortwave
        
        print(f"  Sensible heat flux shape: {H.shape}")
        print(f"  Net shortwave shape: {K.shape}")
        
        # Load all SEB terms at once
        print("\nLoading complete surface energy balance...")
        seb = sim.load_seb()
        
        print(f"  SEB contains: {list(seb.keys())}")
        print(f"  Time points: {len(seb['t'])}")
        
        return seb
        
    except FileNotFoundError:
        print("Energy balance data (facEB.nc) not available")
        return None


def example_facet_temperature(sim):
    """Example 5: Load facet temperatures."""
    
    if sim is None:
        return
    
    print("\n" + "=" * 70)
    print("Example 5: Facet Temperatures")
    print("=" * 70)
    
    try:
        # Display available variables
        print("\nAvailable variables in facT.nc:")
        sim.load_fac_temperature()
        
        # Load temperature
        print("\nLoading facet temperatures...")
        T = sim.load_fac_temperature('T')
        print(f"  Shape: {T.shape}")
        print(f"  Dimensions: {T.dims}")
        
        if T.ndim >= 2:
            # Surface temperature (first layer)
            T_surf = T[:, 0]
            print(f"\n  Surface temperature range: [{T_surf.min().values:.2f}, {T_surf.max().values:.2f}] K")
        
        return T
        
    except FileNotFoundError:
        print("Temperature data (facT.nc) not available")
        return None


def example_area_averaging(sim, seb):
    """Example 6: Area averaging."""
    
    if sim is None or seb is None:
        return
    
    print("\n" + "=" * 70)
    print("Example 6: Area Averaging")
    print("=" * 70)
    
    # Area-average individual terms
    print("\nArea-averaging net shortwave radiation...")
    K_avg = sim.area_average_fac(seb['Kstar'])
    print(f"  Time-series length: {len(K_avg)}")
    print(f"  Mean value: {K_avg.mean():.2f} W/m²")
    
    # Area-average all SEB terms at once
    print("\nArea-averaging all SEB terms...")
    seb_avg = sim.area_average_seb(seb)
    
    print(f"  Available terms: {list(seb_avg.keys())}")
    
    # Time average
    print("\nTime-averaging (last half of simulation)...")
    t = seb['t']
    t_mid = t[len(t)//2]
    K_time_avg = sim.time_average(seb['Kstar'], t, tstart=t_mid)
    print(f"  Time-averaged net shortwave: {K_time_avg.mean():.2f} W/m²")
    
    return seb_avg


def visualize_seb_timeseries(sim, seb_avg):
    """Example 7: Visualize SEB time series."""
    
    if sim is None or seb_avg is None:
        return
    
    print("\n" + "=" * 70)
    print("Example 7: Surface Energy Balance Time Series")
    print("=" * 70)
    
    try:
        t = seb_avg['t'] / 3600  # Convert to hours
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
        
        # Radiation components
        ax1.plot(t, seb_avg['Kstar'], 'r-', label='K*', linewidth=2)
        ax1.plot(t, seb_avg['Lstar'], 'b-', label='L*', linewidth=2)
        ax1.set_ylabel('Radiation (W/m²)')
        ax1.set_title('Area-Averaged Surface Energy Balance')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Turbulent fluxes and ground heat flux
        ax2.plot(t, seb_avg['H'], 'r-', label='H (sensible)', linewidth=2)
        ax2.plot(t, seb_avg['E'], 'b-', label='E (latent)', linewidth=2)
        ax2.plot(t, seb_avg['G'], 'g-', label='G (ground)', linewidth=2)
        ax2.set_xlabel('Time (hours)')
        ax2.set_ylabel('Flux (W/m²)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('seb_timeseries.png', dpi=150)
        print("\nSaved: seb_timeseries.png")
        plt.close()
        
    except Exception as e:
        print(f"Could not create plot: {e}")


def main():
    """Run all examples."""
    
    print("\n" + "=" * 70)
    print("uDALES Facet Data Examples")
    print("=" * 70)
    print("\nThese examples show how to work with surface/facet data.")
    print("Make sure you have an experiment with facet output files.")
    
    # Example 1: Load simulation
    sim = example_facet_loading()
    
    if sim is None:
        print("\nExamples require a valid experiment directory.")
        print("Please adjust the path in example_facet_loading().")
        return
    
    # Example 2: Facet properties
    albedo, emissivity = example_facet_properties(sim)
    
    # Example 3: Momentum data
    pf = example_facet_momentum(sim)
    
    # Example 4: Energy balance
    seb = example_energy_balance(sim)
    
    # Example 5: Temperature
    T = example_facet_temperature(sim)
    
    # Example 6: Area averaging
    if seb is not None:
        seb_avg = example_area_averaging(sim, seb)
        
        # Example 7: Visualization
        visualize_seb_timeseries(sim, seb_avg)
    
    print("\n" + "=" * 70)
    print("✓ All facet data examples completed!")
    print("=" * 70)
    print("\nKey takeaways:")
    print("- Use load_fac_momentum() for pressure and shear stress")
    print("- Use load_fac_eb() for surface energy balance terms")
    print("- Use load_fac_temperature() for temperatures in facet layers")
    print("- Use load_seb() to load all SEB terms at once")
    print("- Use assign_prop_to_fac() to map wall properties to facets")
    print("- Use area_average_fac() for area-weighted averaging")
    print("- Use time_average() for time averaging")


if __name__ == "__main__":
    main()
