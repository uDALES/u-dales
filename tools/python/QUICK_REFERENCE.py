"""
Quick Start Guide for uDALES Python Tools

A concise reference for common operations.
"""

from udbase import UDBase
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# LOADING SIMULATION
# =============================================================================

# Basic loading
sim = UDBase(expnr=65, path='experiments/065')

# Access properties
print(f"Grid: {sim.itot} x {sim.jtot} x {sim.ktot}")
print(f"Domain: {sim.xlen} x {sim.ylen} x {sim.zsize}")

# Grid arrays
x_centers = sim.xt  # Cell centers in x
x_edges = sim.xm    # Cell edges in x
z_centers = sim.zt  # Cell centers in z
z_edges = sim.zm    # Cell edges in z

# =============================================================================
# FIELD DATA
# =============================================================================

# List available variables
sim.load_field()  # Shows all variables in fielddump

# Load 3D instantaneous field
u = sim.load_field('u')  # Returns xarray.DataArray

# Access specific time/location
u_t0 = u.isel(time=0)           # First time step
u_slice = u.isel(time=0, zt=10)  # Horizontal slice at k=10

# Time-averaged 3D field
u_avg = sim.load_stat_t('u')

# Slab-averaged (vertical profile)
u_profile = sim.load_stat_xyt('u')

# 2D slices
u_horiz = sim.load_slice('k', 'u')  # Horizontal (k = z)
u_vert_x = sim.load_slice('i', 'u')  # Vertical x-normal (i)
u_vert_y = sim.load_slice('j', 'u')  # Vertical y-normal (j)

# =============================================================================
# FACET DATA
# =============================================================================

# Momentum (pressure, shear stress)
pf = sim.load_fac_momentum('pf')     # Pressure
taux = sim.load_fac_momentum('taux')  # Shear stress x

# Surface energy balance
H = sim.load_fac_eb('hf')      # Sensible heat flux
E = sim.load_fac_eb('ef')      # Latent heat flux
K = sim.load_fac_eb('netsw')   # Net shortwave

# Complete SEB
seb = sim.load_seb()  # Dictionary with all terms
# seb['Kstar'], seb['Lstar'], seb['H'], seb['E'], seb['G'], seb['t']

# Facet temperatures
T = sim.load_fac_temperature('T')        # Temperature in layers
dTdz = sim.load_fac_temperature('dTdz')  # Temperature gradient

# =============================================================================
# FACET ANALYSIS
# =============================================================================

# Assign wall properties to facets
albedo = sim.assign_prop_to_fac('al')        # Albedo
emissivity = sim.assign_prop_to_fac('em')    # Emissivity
z0 = sim.assign_prop_to_fac('z0')            # Roughness length
conductivity = sim.assign_prop_to_fac('lam') # Thermal conductivity

# Area-weighted averaging
H_avg = sim.area_average_fac(H)  # Average over all facets

# Area-average specific facets
roof_mask = sim.facs['typeid'] == 1
H_roof = sim.area_average_fac(H, roof_mask)

# Area-average SEB
seb_avg = sim.area_average_seb(seb)

# Time averaging
t = seb['t']
K_time_avg = sim.time_average(seb['Kstar'], t, tstart=3600, tstop=7200)

# =============================================================================
# VISUALIZATION
# =============================================================================

# Vertical profile
plt.figure()
plt.plot(u_profile, sim.zt)
plt.xlabel('u (m/s)')
plt.ylabel('z (m)')
plt.grid(True)
plt.savefig('profile.png')

# Horizontal contour
X, Y = np.meshgrid(sim.xt, sim.yt)
plt.figure()
plt.contourf(X, Y, u_slice.values.T, levels=20)
plt.colorbar(label='u (m/s)')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.savefig('contour.png')

# Time series
plt.figure()
plt.plot(seb_avg['t']/3600, seb_avg['Kstar'], label='K*')
plt.plot(seb_avg['t']/3600, seb_avg['H'], label='H')
plt.xlabel('Time (hours)')
plt.ylabel('Flux (W/m²)')
plt.legend()
plt.savefig('timeseries.png')

# =============================================================================
# WORKING WITH XARRAY
# =============================================================================

# xarray provides labeled dimensions
u = sim.load_field('u')
print(u.dims)    # e.g., ('time', 'zt', 'yt', 'xm')
print(u.shape)   # e.g., (10, 96, 64, 64)

# Select by dimension name
u_t5 = u.sel(time=5.0)           # Select by value
u_z10 = u.isel(zt=10)            # Select by index

# Slicing
u_early = u.sel(time=slice(0, 3600))

# Statistics
u_mean = u.mean(dim='time')      # Time mean
u_max = u.max(dim=['yt', 'xm'])  # Max over horizontal

# Convert to numpy
u_np = u.values  # numpy array

# =============================================================================
# COMMON PATTERNS
# =============================================================================

# Check if facet data is available
if sim.facs and 'area' in sim.facs:
    print(f"Number of facets: {len(sim.facs['area'])}")

# Check if geometry is loaded
if sim.geom is not None:
    print("Geometry available")

# Loop over time in SEB
seb = sim.load_seb()
for i, t in enumerate(seb['t']):
    H_t = seb['H'][:, i]  # Heat flux at time i
    # ... process

# Calculate domain-averaged velocity
u = sim.load_field('u')
u_domain_avg = u.mean(dim=['xm', 'yt', 'zt'])  # Average over space

# =============================================================================
# ERROR HANDLING
# =============================================================================

try:
    u = sim.load_field('u')
except FileNotFoundError:
    print("fielddump not found")

try:
    seb = sim.load_seb()
except FileNotFoundError:
    print("Energy balance files not available")

# =============================================================================
# TIPS
# =============================================================================

# 1. Call load methods without var argument to see available variables
#    sim.load_field()

# 2. xarray is lazy - data not loaded until accessed
#    u = sim.load_field('u')  # Fast, data not loaded
#    u_values = u.values      # Now data is loaded

# 3. Use .isel() for integer indexing, .sel() for coordinate values
#    u.isel(time=0)      # First time step
#    u.sel(time=0.0)     # Time = 0.0 seconds

# 4. Check dimensions before plotting
#    print(u.dims, u.shape)

# 5. For large files, work with subsets
#    u_subset = sim.load_field('u').isel(time=slice(0, 10))
