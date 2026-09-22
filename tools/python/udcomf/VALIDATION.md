# Thermal-comfort validation

This file records the numerical checks for the `udcomf` thermal-comfort
pipeline. It distinguishes reproducible software verification from physical
validation against measurements.

## Step 12 verification matrix

| Requirement | Automated evidence |
|---|---|
| Wind interpolation at several receptor heights | `test_udcomf_multiheight.py` checks scalar interpolation, saved `ws_local` height selection, `ws_10`, full statistics and k-slice inputs. |
| Nighttime and blocked direct radiation | `test_udcomf_shortwave.py` and `test_udcomf_thermalcomfort.py` check sun visibility, wall shading, direct/non-direct separation and rejection of direct irradiance below the horizon. |
| Open sky and uniform-temperature enclosure | `test_udcomf_shortwave.py`, `test_udcomf_longwave.py` and `test_udcomf_thermalcomfort.py` check analytical open-sky, ground and enclosure limits. |
| Radiation sign and energy consistency | Shortwave tests check albedo conversion and analytical plane irradiance; longwave tests check sky/facet partitioning and uniform emission. Negative source fluxes are rejected. |
| Exact preceding 15-minute means | `test_udcomf_checkpoints.py` checks exact integration of a linear series over `[time-900 s,time)`, missing coverage, cadence gaps, restart checkpoints and multi-height isolation. |
| NetCDF schema and metadata | `test_udcomf_exchange.py` and `test_udcomf_thermalcomfort.py` check `(x,y,time)` order, units, coordinates, masks, time bounds, georeferencing and missing-value handling. |
| Vapour humidity | `test_udcomf_multiheight.py` checks `qv = qt - ql` after height interpolation and rejects negative values without clipping. |
| MRT, PET, UTCI and WBGT references | `test_udcomf_thermalcomfort.py` contains the fixed cases listed below. |
| Model-neutral result | The same test supplies identical exchange data labelled PALM and UrbClim and requires identical index arrays. |
| End-to-end small case | `test_udcomf_end_to_end.py` runs a synthetic one-cell uDALES case through atmospheric extraction, receptor shortwave/longwave calculation, strict exchange export, and all four indices. |

## Fixed numerical references

### MRT

The test uses the documented six-direction standing-person equation with
horizontal weights 0.06, vertical weights 0.22, shortwave absorptivity 0.70,
longwave emissivity 0.97 and a 60 degree solar zenith. Its independently
evaluated weighted irradiances are:

- shortwave: 225.2922678358 W m-2;
- longwave: 395.6 W m-2;
- expected MRT: 314.985841050193 K.

The uniform-radiation end-to-end case also checks the black-body limit. Each
receiving plane's finite angular weights are normalized to the analytical
constant-radiance integral. This removes the 0.1255 K warm MRT bias that the
unnormalized default 8 by 32 quadrature produced in a 25 degree Celsius
uniform enclosure. The test now requires agreement within 0.0002 K.
The radiation checkpoint manifest was advanced to version 2 so checkpoints
created before this correction cannot be mixed with normalized results.

### PET and UTCI

The published `pythermalcomfort` examples are retained as fixed regression
values:

- PET: `tdb=25 C`, `tr=25 C`, `v=0.1 m s-1`, `rh=50%`, `met=1.2`,
  `clo=0.5` gives 24.67 C;
- UTCI: `tdb=25 C`, `tr=25 C`, `v10=1 m s-1`, `rh=50%` gives 24.6 C.

The tests exercise these values through the strict exchange-file pipeline,
not by calling the library directly. The output records the exact
`pythermalcomfort` version.

### WBGT

The ISO 7243 combinations documented by `pythermalcomfort` are reproduced:

- no direct solar load: `Tnwb=25 C`, `Tg=32 C` gives 27.1 C;
- direct solar load: `Tnwb=25 C`, `Tg=32 C`, `Ta=20 C` gives 25.9 C.

A separate directional-radiation case checks the iterative globe and wick
temperatures against roots calculated with an independent bracketed scalar
solver applied to the Liljegren heat balances:

- globe temperature: 311.701170960989 K;
- natural wet-bulb temperature: 296.051700718050 K;
- outdoor WBGT: 26.741424694833 C.

The production fixed-point solver must agree within 0.02 K, its configured
convergence tolerance.

## Scope and remaining evidence

The test suite verifies equations, units, dimensions, numerical convergence,
data flow and limiting cases. The final test uses synthetic solver outputs; it
does not run an LES or establish agreement with a field instrument.

The six-direction urban WBGT radiation treatment deliberately replaces the
flat-ground radiation reconstruction in the published Liljegren program. Its
sensor balances are independently verified, but the directional adaptation
still needs comparison with measured globe and natural-wet-bulb temperatures
before it can be described as measurement-validated. Vegetation radiation is
also rejected explicitly because attenuation and tree longwave emission have
not yet been implemented.

References:

- [Höppe (1999), PET](https://doi.org/10.1007/s004840050118)
- [Bröde et al. (2012), operational UTCI](https://doi.org/10.1007/s00484-011-0454-1)
- [pythermalcomfort model documentation](https://pythermalcomfort.readthedocs.io/en/stable/documentation/models.html)
- [Liljegren WBGT reference source](https://github.com/mdljts/wbgt/blob/master/src/wbgt.c)
- [ISO 7243 WBGT implementation examples](https://pythermalcomfort.readthedocs.io/en/stable/_modules/pythermalcomfort/models/wbgt.html)
