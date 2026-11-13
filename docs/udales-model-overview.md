# Model overview

uDALES (urban Dutch Atmospheric Large-Eddy Simulation) is an open-source large-eddy simulation model for the urban microclimate. It simulates the flow of air through and above cities at building-resolving resolution, together with the transport of heat, moisture, and pollutants. uDALES originates from the atmospheric LES code DALES and has been extended into a multi-physics urban modelling tool [@Suter2022; @Owens2024].

## Dynamical core

uDALES solves the filtered, incompressible Navier–Stokes equations under the Boussinesq approximation on a staggered Cartesian grid, with grid stretching supported in the vertical. Subgrid-scale turbulence is represented by an eddy-viscosity closure (the Vreman model or a one-equation subgrid TKE model). Several advection schemes are available, including second-order central and flux-limited kappa schemes for scalars; the pressure Poisson equation is solved with an FFT-based method.

## Urban surfaces

Buildings and terrain are represented by an immersed boundary method (IBM). Since uDALES 2.0, the geometry is specified as a triangulated surface (STL), and the IBM is discretely mass-, momentum-, and scalar-conservative [@Owens2024]. Surface fluxes of momentum, heat, and moisture at the resolved surfaces (facets) are computed with urban wall functions.

## Surface energy balance

Each facet can carry a full surface energy balance: shortwave and longwave radiation (with view factors between facets computed during pre-processing), conduction through multi-layer walls, and turbulent sensible and latent heat fluxes [@Suter2022]. This enables coupled simulations of urban airflow and heat transfer, including the effects of shading and radiative trapping in street canyons.

## Vegetation

Trees and other vegetation are represented as porous drag elements that also exchange heat, moisture, and scalars with the air, allowing studies of the effect of urban greening on flow, temperature, and air quality [@Grylls2021a].

## Dispersion and air quality

uDALES transports an arbitrary number of passive or reactive scalars released from point, line, or volume sources. A basic NO~x~–O~3~ chemistry scheme is available for street-scale air quality studies [@Grylls2019], alongside idealised air purifiers and heat-pump sources.

## Inflow and boundary conditions

Simulations can be run with periodic lateral boundaries or with inflow–outflow conditions. Turbulent inflow can be supplied by a precursor (driver) simulation — see [precursor simulations](udales-driver-simulations.md) — or generated synthetically. Details of all boundary condition options are given in the [boundary conditions](udales-boundary-conditions.md) reference.

## Parallelisation

The code is parallelised with MPI using a 2-D pencil domain decomposition based on the 2DECOMP&FFT library, and scales to large domains and core counts on HPC systems — see the [parallelisation notes](udales-2decomp.md).

## Where to read more

The model formulation and validation are described in the uDALES description papers: uDALES 1.0 [@Suter2022] and the conservative IBM of uDALES 2.0 [@Owens2024]. Studies using uDALES are collected in the [publication list](udales-pub-list.md). To get started with the model itself, see the [workflow overview](udales-workflow.md).

## References

\bibliography
