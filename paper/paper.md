---
title: 'uDALES: large-eddy-simulation software for urban flow, dispersion, and microclimate modelling'
tags:
  - fluid dynamics
  - computational fluid dynamics
  - large eddy simulation
  - urban climate
  - fortran
authors:
  - name: Tom Grylls
    orcid: 0000-0002-0948-4932
    affiliation: 1
  - name: Ivo Suter
    orcid: 0000-0002-8612-4033
    affiliation: "1, 2"
  - name: Birgit S. Sützl
    orcid: 0000-0001-9638-5643
    affiliation: 1
  - name: Sam Owens
    orcid: 0000-0003-3253-1896
    affiliation: 1
  - name: David Meyer
    orcid: 0000-0002-7071-7547
    affiliation: "3, 1"
  - name: Maarten van Reeuwijk
    orcid: 0000-0003-4840-5050
    affiliation: 1
affiliations:
 - name: Department of Civil and Environmental Engineering, Imperial College London, London, UK
   index: 1
 - name: Empa, Swiss Federal Laboratories for Materials Science and Technology, Dübendorf, Switzerland
   index: 2
 - name: Department of Meteorology, University of Reading, Reading, UK
   index: 3
date: 10 January 2021
bibliography: paper.bib
---


# Summary

With continuing urbanization, challenges associated with the urban environment such as air quality, heat islands, pedestrian thermal comfort, and wind loads on tall buildings, are increasingly relevant. Our ability to realistically capture processes such as the transport of heat, moisture, momentum and pollutants, and those of radiative transfer in urban environments is key to understanding and facing these challenges [@Oke2017]. The turbulent nature of the urban flow field and the inherent heterogeneity and wide range of scales associated with the urban environment result in a complex modelling problem. Large-eddy simulation (LES) is an approach to turbulence modelling used in computational fluid dynamics to simulate turbulent flows over a wide range of spatial and temporal scales. LES is one of the most promising tools to model the interactions typical of urban areas due to its ability to resolve the urban flow field at resolutions of $O$(1 m, 0.1 s), over spatial domains of $O$(100 m), and time periods of $O$(10 h). Although there are many scalable LES models for atmospheric flows, to our knowledge, only few are capable of explicitly representing buildings and of modelling the full range of urban processes (e.g. PALM-4U @Resler2017; @Maronga2020; or OpenFoam @Weller1998).

uDALES (urban Dutch Atmospheric LES) is an extension of DALES (Dutch Atmospheric LES; @Heus2010, @
2015). It has the additional functionality of modelling buildings within the fluid domain and therefore the capability to model urban environments at the microclimate scale with wet thermodynamics (Table 1). The uDALES framework includes tools to enable users to model a wide variety of idealized and complex urban morphologies [@Sutzl2021]. uDALES uses an Arakawa C-grid and typically uses second-order central-differencing schemes. For scalar quantities, e.g. for pollution concentration, it is possible to use a kappa-scheme for advection to ensure positivity [@Hundsdorfer1995]. A third-order Runge-Kutta time integration scheme is applied. The immersed boundary method, first introduced into DALES by @Pourquie2009 and @Tomas2015, is used to represent buildings (supporting grid-conforming cuboid geometries; @Pourquie2009). Wall functions have been added to calculate the surface scalar and momentum fluxes at the rough immersed boundaries [@Cai2011; @Cai2012; @Uno1995; @Suter2018], based on a local Richardson number. Additional factors are considered to obtain the evapotranspiration fluxes from vegetated surfaces. The code uses fast Fourier transforms to efficiently solve the Poisson equation for pressure and is fully parallelized using MPI (domain decomposition is performed in the spanwise direction, thus producing slabs).

A novel surface energy balance model has been implemented in a two-way-coupled manner [@Suter2018] and includes the effect of turbulent exchange of heat between the surface and the air, as well as radiation and thermal conduction within the surface. uDALES has tools for modelling shortwave and longwave radiative fluxes following a radiosity approach (similar to @Aoyagi2011; @Resler2017), including calculation of direct solar radiation and view factors [@RammohanRao1996], considering shading and multiple reflections. 

uDALES also has the tools necessary to study air quality within cities. Both idealized (point and line) and realistic (street network) sources can be implemented and both passive and reactive scalars can be modelled. The code supports null cycle chemistry (NO-NO$_2$-O$_3$ reactions), and can be easily extended to more sophisticated schemes [@Grylls2019; @Grylls2020]. The high spatial and temporal resolution enables analysis of e.g. real-time pedestrian pollution exposure.

The modelling capabilities of uDALES outlined above, combined with the ability to devise numerous simulation set-ups (e.g. via different lateral boundary conditions: periodic, inflow-outflow and driver) facilitate a plethora of possible studies into the urban environment. Relevant recent publications and their research applications are summarized in Table 1. 

| Research application                                         | Reference                                         |
| ------------------------------------------------------------ | ------------------------------------------------- |
| Urban boundary layers/ boundary-layer meteorology            | @Grylls2019; @Sutzl2021; @Sutzl2021thesis         |
| Urban climate (radiation, green roofs and walls, trees etc.) | @Suter2018; @Suter2021; @Grylls2021               |
| Pollution dispersion/ urban air quality                      | @Grylls2019;  @Grylls2020                         |
| Buoyancy/ convective and stable conditions                   | @Suter2018; @Grylls2020b; @Grylls2020; @Grylls2021 |

Table: uDALES research applications.

Here we present uDALES, a free and open-source large-eddy-simulation software for urban flow, dispersion, and microclimate modelling. In the current platform on GitHub we include: (i) cross-platform support for GNU, Intel, and Cray compilers on Windows and macOS systems, (ii) continuous integration with build and regression tests to check for compilation and simulation errors respectively [@Riechert2019], and (iii) several pre- and post-processing scripts in MATLAB as well as Singularity [@Kurtzer2017] scripts to address issues of scientific reproducibility (e.g. @Meyer2020). Current developments include the role of urban trees [@Grylls2021] and future releases may include the use of PsychroLib [@Meyer2019] to improve the calculation of psychrometric properties of air, and the ability to simulate the diurnal cycle of the urban microclimate in response to solar radiation [@Suter2018]. A detailed description of the model, including a validation study and an example of the surface energy balance is provided in @Suter2021.

# References
