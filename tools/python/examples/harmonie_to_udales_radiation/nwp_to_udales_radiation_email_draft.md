# HARMONIE NWP Radiation to uDALES Inputs

## Summary of Output Mapping

The final mapping is:

```text
HARMONIE strd [J m**-2 accumulated]
  -> time difference and spatial average
  -> downward sky longwave flux [W m**-2]
  -> uDALES timedeplw.inp.<expnr>

HARMONIE ssrd [J m**-2 accumulated]
  -> time difference and spatial average
  -> GHI [W m**-2]
  -> Erbs split into DNI and Dsky [W m**-2]
  -> uDALES facet radiation mapping with shading, SVF, albedo, and reflections
  -> uDALES timedepsw.inp.<expnr> plus Sdir.nc
```

## Variables Used

The relevant HARMONIE variables are accumulated from the start of the forecast
integration and are stored in units of `J m**-2`.

| HARMONIE short name | Full name | Unit in GRIB | Used for | uDALES output |
| --- | --- | --- | --- | --- |
| `ssrd` | Surface solar radiation downwards | `J m**-2`, accumulated | Incoming global horizontal shortwave radiation at the surface | `timedepsw.inp.<expnr>` and companion `Sdir.nc` |
| `strd` | Surface thermal radiation downwards | `J m**-2`, accumulated | Downward sky longwave radiation | `timedeplw.inp.<expnr>` |

We did not use the following available accumulated net radiation fields:

| HARMONIE short name | Full name | Reason not used |
| --- | --- | --- |
| `ssr` | Surface net solar radiation | This already includes the NWP model surface albedo and is not appropriate for uDALES, where net shortwave is recomputed per urban facet using local facet orientation, shading, albedo, sky-view factors, and multiple reflections. |
| `str` | Surface net thermal radiation | uDALES requires the downward sky longwave forcing; the model then applies sky-view factors and facet emissivity internally. |
| `ttr` | Top net thermal radiation | This is a top-of-atmosphere/net column quantity and is not the surface downward longwave forcing required by uDALES. |

## Spatial Averaging Domain

The HARMONIE fields were spatially averaged over the same central Paris domain
used in our NWP-to-uDALES nudging workflow. The default Lambert-93 bounding box
is:

```text
EPSG:2154 easting  = 649375 .. 652600 m
EPSG:2154 northing = 6861175 .. 6864020 m
```

For each forecast time and each radiation variable, we first compute the mean
accumulated radiation over the selected HARMONIE grid cells.

## Accumulated Energy to Flux Conversion

Since `ssrd` and `strd` are accumulated energy per unit area from the start of
the HARMONIE integration, they must be differenced in time to obtain fluxes in
$W\,m^{-2}$.

For a HARMONIE accumulated field $A(x, y, \tau)$ in $J\,m^{-2}$, the spatial
mean over the selected central Paris grid cells $\Omega$ is:

$$
\bar{A}(\tau_k)
= \frac{1}{|\Omega|}
  \sum_{(x,y)\in\Omega} A(x, y, \tau_k)
$$

The instantaneous interval-mean flux is then:

$$
F(\tau_k)
= \max\left(
    0,\,
    \frac{\bar{A}(\tau_k) - \bar{A}(\tau_k - \Delta\tau)}
         {\Delta\tau}
  \right)
$$

where:

$$
\Delta\tau = 900\ \mathrm{s}
$$

for the 15-minute HARMONIE output frequency. The `max(0, ...)` clamp is only a
guard against tiny negative values caused by numerical or GRIB packing noise.
Large negative differences would indicate an inconsistent accumulated field and
are treated as an error.

The resulting native 15-minute flux series is then linearly interpolated to the
uDALES radiation times.

## Longwave Conversion

For longwave radiation, the conversion is direct after differencing:

$$
LW_{\mathrm{sky}}(t)
= \mathrm{interp}_t\left[\frac{d\,strd}{dt}\right]
$$

where $LW_{\mathrm{sky}}(t)$ is in $W\,m^{-2}$.

This is written to:

```text
timedeplw.inp.<expnr>
```

The file contains two header lines followed by:

```text
time_seconds   LWsky_W_m-2
```

For experiment 300, the namelist requires `ntimedeplw = 31`, so we write 31
longwave time points. The values represent the downward sky longwave flux. In
uDALES, this is read as the time-dependent `skyLW` forcing; the surface energy
balance then applies sky-view factors and facet emissivity internally.

## Shortwave Conversion

For shortwave radiation, `ssrd` gives downward global horizontal irradiance
after differencing:

$$
GHI(t)
= \mathrm{interp}_t\left[\frac{d\,ssrd}{dt}\right]
$$

where $GHI(t)$ is in $W\,m^{-2}$.

uDALES needs net shortwave radiation per urban facet, not a single horizontal
grid-cell value. Therefore, we split the HARMONIE-derived `GHI` into direct
normal irradiance and diffuse horizontal sky irradiance, then map these to the
uDALES facets using the existing uDALES radiation preprocessing routines.

### Direct/Diffuse Split

Because the available HARMONIE data contain `ssrd` but not separate direct and
diffuse shortwave components, we use an empirical Erbs clearness-index split.

For solar zenith angle $\theta_z$:

$$
\mu = \cos(\theta_z)
$$

The extraterrestrial horizontal irradiance is:

$$
I_{0h} = S_0 E_0 \mu
$$

$$
E_0 = 1 + 0.033\cos\left(\frac{2\pi\,DOY}{365}\right)
$$

$$
S_0 = 1367\ W\,m^{-2}
$$

The clearness index is:

$$
K_t = \frac{GHI}{I_{0h}}
$$

The diffuse fraction $f_d = D_{\mathrm{sky}} / GHI$ is:

$$
f_d =
\begin{cases}
1 - 0.09K_t,
& K_t \le 0.22 \\
0.9511 - 0.1604K_t + 4.388K_t^2 - 16.638K_t^3 + 12.336K_t^4,
& 0.22 < K_t \le 0.80 \\
0.165,
& K_t > 0.80
\end{cases}
$$

Then:

$$
D_{\mathrm{sky}} = f_d\,GHI
$$

$$
DNI = \frac{GHI - D_{\mathrm{sky}}}{\mu}
$$

At night, or when the sun is below the horizon, both direct and diffuse
shortwave are set to zero. For extremely low sun angles, the positive shortwave
is treated as diffuse to avoid unstable near-horizontal ray tracing.

### Mapping to uDALES Facets

For every uDALES shortwave time step, we compute solar zenith and azimuth from
the case namelist time, latitude, longitude, timezone, elevation, and domain
azimuth. The local solar azimuth used by the facet ray tracer is:

$$
azimuth_{\mathrm{local}}
= azimuth_{\mathrm{solar}} - xazimuth
$$

The direct component is passed to the existing uDALES direct-shortwave
preprocessor, using the configured shortwave method for the case. For experiment
300 this uses the configured scanline/f2py route. This step accounts for facet
orientation and shading by the urban geometry, producing a direct shortwave
field per facet and time:

$$
S_{\mathrm{dir,facet}}(n, t)
$$

The diffuse component is applied using the uDALES sky-view factor and
view-factor/reflection machinery. With surface energy balance enabled, uDALES
computes multiple shortwave reflections using facet albedo and the sparse view
factor matrix. Conceptually, the absorbed/net shortwave per facet is:

$$
netsw_{\mathrm{facet}}(n,t)
= direct_{\mathrm{absorbed}}(n,t)
 + diffuse_{\mathrm{sky,absorbed}}(n,t)
 + reflected_{\mathrm{shortwave,absorbed}}(n,t)
$$

This net facet shortwave matrix is written to:

```text
timedepsw.inp.<expnr>
```

The file format is:

```text
header line
time_1 time_2 ... time_nt
netsw_facet_1_time_1 ... netsw_facet_1_time_nt
netsw_facet_2_time_1 ... netsw_facet_2_time_nt
...
```

For experiment 300, `ntimedepsw = 181` and `dtSP = 600 s`, so the shortwave file
contains 181 time points from 0 s to 108000 s. We also write:

```text
Sdir.nc
```

as a NetCDF companion file containing the time-dependent direct shortwave field
per facet. The uDALES solver uses `timedepsw.inp.<expnr>` at runtime; `Sdir.nc`
is retained for comparison, inspection, and preprocessing consistency.