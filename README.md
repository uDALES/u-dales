# uDALES [![CI](https://github.com/uDALES/u-dales/workflows/CI/badge.svg)](https://github.com/uDALES/u-dales/actions)

This is the development repository for the uDALES (urban Dutch Atmospheric Large Eddy Simulation) model. It contains the uDALES software source-code, tests, documentation, and examples with tools to pre- and post-process inputs and outputs.

<div align="center">
<img src="docs/assets/images/fielddump_slice_2D.102.png" alt="uDALES" height="256">
<img src="docs/assets/images/fielddump_slice_3D.102.png" alt="uDALES" height="256"> 
</div>

## Overview

With continuing urbanization, challenges associated with the urban environment such as air quality, heat islands, pedestrian thermal comfort, and wind loads on tall buildings, are increasingly relevant. Our ability to realistically capture processes such as the transport of heat, moisture, momentum and pollutants, and those of radiative transfer in urban environments is key to understanding and facing these challenges.

The turbulent nature of the urban flow field and the inherent heterogeneity and wide range of scales associated with the urban environment result in a complex modelling problem. Large-eddy simulation (LES) is an approach to turbulence modelling used in computational fluid dynamics to simulate turbulent flows over a wide range of spatial and temporal scales. LES is one of the most promising tools to model the interactions typical of urban areas due to its ability to resolve the urban flow field at resolutions of _O_(1 m, 0.1 s), over spatial domains of _O_(100 m), and time periods of _O_(1 h). Although there are many scalable LES models for atmospheric flows, to our knowledge, only few are capable of explicitly representing buildings and of modelling the full range of urban processes.

uDALES is capable of modelling urban environments at the microclimate scale including wet thermodynamics, idealized and complex morphologies, three-dimensional urban surface energy balance and real-time pollution dispersion.

## Validation

The dynamic core, DALES, has been validated extensively and used as part of several atmospheric intercomparison studies over a period of 30 years (Nieuwstadt and Brost, 1986; Heus et al., 2010). The use of the immersed boundary method to model obstacles has been validated against both wind and water tunnel data, for idealised geometries by Tomas et al. 2015 and complex geometries by Grylls, 2020. The wall functions and three-dimensional surface energy balance model was validated by Suter, 2018. The ability to model pollutuon dispersion in cities has been validated by Grylls et al., 2020 and Grylls, 2020. 

## Documentation

**If you are new to uDALES, please follow our [getting started guide](https://udales.github.io/u-dales/udales-getting-started).** For User's guides and general reference documentation, please see the [uDALES website](https://udales.github.io/u-dales/).

If you are a developer, please also refer to the [development notes](DEVELOP.md).

## Contributing

If you are looking to contribute, please read our [Contributors' guide](CONTRIBUTING.md) for details.

## Copyright and License

General DALES copyright applies for any files part of the original DALES distribution and are marked as such at the beginning of each file.

Additional files provided in uDALES are copyrighted "the uDALES Team" and are marked as such at the beginning of each file.

All files are licensed under the GNU General Public License. See [LICENSE.txt](LICENSE.txt).
