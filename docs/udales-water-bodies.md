# Water-body facets

uDALES can represent open water — rivers, canals, basins — as a dedicated
facet class within the facet energy-balance framework (`lEB`). A facet whose
wall type lies in the reserved band **(−30, −20]** (by convention `−21`) is
treated as a water surface: its column is held near a prescribed water
temperature and its surface evaporates at the saturated, resistance-free
open-water rate. This page describes the physical model, its configuration,
and its limits.

## Why a water class is needed

Without special treatment, a facet typed as water behaves like thermally thin
dark pavement. With the molecular conductivity of still water
(0.6 W m⁻¹ K⁻¹) the diurnal thermal penetration depth is

$$\delta = \sqrt{\kappa P/\pi} \approx 6\ \mathrm{cm},$$

so only a sliver of the column stores heat and the skin can overheat by tens
of kelvin under summer radiation; and with no latent-heat path the Bowen
ratio is infinite. Such a surface heats the air day and night and humidifies
nothing — the opposite of observed water-body behaviour, in which turbulent
mixing and (for rivers) advection keep the surface near the bulk water
temperature, evaporation dominates the turbulent exchange, and the water body
is a mild warm anomaly only at night (see e.g. Jacobs et al. 2020).

## Physical model

The water facet stays inside the ordinary facet machinery — same radiation,
same wall functions, same implicit layer solver — with three choices that
make it water:

1. **A well-mixed column.** The layer conductivity in `factypes.inp` is set
   to an effective value representing turbulent mixing,

   $$\lambda_{eff} = \rho_w c_w\,\kappa_{eff}, \qquad
     \kappa_{eff} \approx 10^{-3}\ \mathrm{m^2\,s^{-1}}
     \;\Rightarrow\; \lambda_{eff} \approx 4200\ \mathrm{W\,m^{-1}K^{-1}},$$

   the eddy diffusivity shown to reproduce well-mixed urban water columns
   (Jacobs et al. 2020). The steady skin offset under a net surface flux
   imbalance $Q$ across a column of depth $D$ is then $Q D/\lambda_{eff}$ —
   a few tenths of a kelvin at peak summer forcing.

2. **An anchored deep node.** The facet layer solve holds its innermost node
   at a fixed temperature; for water facets this is `waterT`
   (`&ENERGYBALANCE`). Physically the anchor is the advective reservoir of a
   through-flowing river: when the water's residence time in the domain
   ($L/u$, typically an hour for a city-scale river) is short against the
   thermal response time of the mixed column (hours per kelvin), the local
   water temperature is set upstream, not by local fluxes. The diagnosed
   conduction flux into the anchor is the enthalpy exported by the flow, so
   the water surface energy balance closes as
   $R_n = H + LE + G_{advective}$.

3. **Saturated open-water evaporation.** The surface is saturated at the
   skin temperature with no canopy or soil resistance:

   $$LE = \rho L_v\,\frac{q_{sat}(T_s) - q_a}{r_a}, \qquad
     r_a = \frac{1}{C_h |u_{tan}|},$$

   with the transfer coefficient from the standard stability-aware
   Uno/Louis facet wall function. Internally this is the existing ERA40-type
   wet-surface flux evaluated with vegetation fraction 0, relative surface
   humidity 1 and zero canopy/soil resistance; there is no water-content
   depletion. Evaporation only (no dew) is currently allowed, matching the
   vegetated-surface convention. Under periodic boundary conditions the
   moisture source is compensated by the `lperiodicEBcorr` sink, so the
   domain humidity budget does not drift.

Radiatively the facet is opaque, like every other facet: shortwave is
absorbed at the skin with the type's albedo (facet-to-facet reflections come
from preprocessing as usual) and longwave is emitted as
$\varepsilon\sigma T_s^4$. No radiation penetrates the water and no in-water
absorption or reflection is modelled; with the column anchored, absorbed
radiation is exported through the deep node rather than overheating a thin
skin, which is the physically relevant consequence of penetration.

## Configuration

**Namelist** (`&ENERGYBALANCE`):

```fortran
waterT = 297.5   ! water anchor temperature [K]; negative (default) -> use flrT
```

For rivers, take `waterT` from routine river-temperature observations for
the simulated period; large flowing rivers have sub-kelvin diurnal amplitude,
so a constant is usually adequate for day-scale runs. The model response to
the anchor is linear and modest (the skin tracks the anchor one to one; LE
changes by roughly the Clausius–Clapeyron slope over the aerodynamic
resistance, ~10 W m⁻² K⁻¹ at summer temperatures), so an observational
uncertainty of a few tenths of a kelvin is negligible.

**Facet inputs**: type the water facets in the reserved band in
`facets.inp` and give the type a row in `factypes.inp` with, typically,

| property | value | note |
|---|---|---|
| `z0`, `z0h` | 3 mm, 0.03 mm | smooth water surface |
| albedo | 0.06 | mean open-water value |
| emissivity | 0.95–0.97 | 0.97 closer to the open-water literature; the difference is a few W m⁻² in L↑, largely compensated by absorbed sky longwave |
| layer thicknesses | water depth / `nfaclyrs` | |
| heat capacity | 4.18 × 10⁶ J m⁻³ K⁻¹ | liquid water |
| conductivity | ~4200 W m⁻¹ K⁻¹ | λ_eff above; the skin offset scales as Q·D/λ_eff, so values from ~1000 upward differ by well under a kelvin |
| diffusivity columns | λ_eff / (ρc) | keep consistent |

**Format warning**: `factypes.inp` rows must be written in the layout
matching the case's `nfaclyrs` (6 + 4·nfaclyrs + 1 columns). A row written
for a different layer count is read without error but maps thicknesses into
heat capacities, producing a column that silently misbehaves. When editing
facet types, include a molecular-conductivity control case in any test: a
water column that fails to drift under strong forcing with λ = 0.6 W m⁻¹ K⁻¹
indicates a parsing problem, not good behaviour.

## Applicability and limits

- The anchored column is a **flowing-water** model. For stagnant ponds and
  closed basins, whose bulk temperature must drift over multi-day episodes,
  the fixed anchor is not appropriate; a prognostic slab (bulk heat budget
  with optional through-flow relaxation) would be the natural extension.
- No shortwave penetration, in-water radiative transfer, skin (cool-film)
  effect, wave-dependent roughness, spray, or ice.
- Reflected longwave is neglected, consistently with all uDALES facets
  ((1−ε)·L↓, of order 20 W m⁻²).
- Dew/condensation onto the water surface is clipped (evaporation only).
- Water facets do not participate in soil-moisture bookkeeping (`lconstW`,
  green-roof water budgets).

## Verification

The implementation is covered by a self-contained MPI integration test
(`tests/integration/water_facets/`, built on `examples/201`) asserting the
anchored isothermal column at `waterT`, the `flrT` fallback, evaporative
`ef < 0` on water facets and exactly zero latent flux on plain facets.
Longer verification runs reproduce the analytic limits: the skin offset
follows Q·D/λ_eff; under anchor perturbations the surface energy balance
closes identically (dG = −(dLE + dL↑ + dH)) with L↑ responding as 4εσT³;
and the latent flux matches the bulk formula with a realistic open-water
Bowen ratio. A useful post-run acceptance check for any case with water
facets is that the facet-mean SEB residual G = K* + L↓ − L↑ − H − LE stays
within plausible advective-export bounds (order 10² W m⁻²) at all output
times.

## References

- Jacobs, C. et al. (2020): Are urban water bodies really cooling? *Urban
  Climate* 32, 100607. https://doi.org/10.1016/j.uclim.2020.100607
- Mironov, D. et al.: FLake — a two-layer bulk freshwater lake model for NWP.
  https://www.cosmo-model.org/content/model/cosmo/misc/flake/default.htm
- Heiskanen, J. J. et al. (2015): Effects of water clarity on lake
  stratification and lake–atmosphere heat exchange. *JGR-Atmospheres* 120.
  https://doi.org/10.1002/2014JD022938
- Järvi, L., Grimmond, C. S. B., Christen, A. (2011): The Surface Urban
  Energy and Water Balance Scheme (SUEWS). *J. Hydrology* 411.
  https://doi.org/10.1016/j.jhydrol.2011.10.001
- van den Hurk, B. et al. (2000): the ERA-40 surface scheme underlying the
  uDALES vegetated-facet evaporation.
