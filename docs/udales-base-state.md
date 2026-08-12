# Hydrostatic base state and bulk flow quantities

This page documents the mathematical formulation of the hydrostatic base state
(`src/modbasestate.f90`), the anchoring of the pressure column when terrain covers part or
all of the domain bottom, the treatment of horizontal slab averages at levels without fluid
cells, and the computation of bulk (volume- and plane-integrated) flow quantities. These
pieces interact in ways that are easy to get wrong; the conventions below were introduced
while resolving issues [#299](https://github.com/uDALES/u-dales/issues/299) and
[#302](https://github.com/uDALES/u-dales/issues/302).

## Role of the base state

uDALES solves the Boussinesq equations, so the dynamics respond only to *differences* of
virtual potential temperature from a horizontal mean — an additive shift of the reference
column does not accelerate the flow. A reference (base) state is nevertheless required for:

- converting liquid-water potential temperature $\theta_l$ to absolute temperature and
  pressure (the Exner function), and hence for the **saturation adjustment** in moist runs,
  where the *absolute* pressure matters;
- anchoring the hydrostatic integration of the slab-mean pressure column each time step;
- providing well-defined values at vertical levels that contain **no fluid cells** (slabs
  buried inside terrain);
- the constant reference $\theta_v$ in subgrid buoyancy terms and the TKE-budget statistics.

The base state is a fixed 1-D column $\big(\theta_{l,b},\, q_{t,b},\, \theta_{v,b},\, p_b,\,
\pi_b\big)(k)$, built **once per start** (cold, warm, or stratified start) from the initial
profiles in `prof.inp` together with the namelist surface pressure `ps`. It is *not* updated
during the run; the evolving slab-mean column (recomputed every RK3 substep, see below)
carries the actual thermodynamic state of the simulation.

## Discrete hydrostatics in Exner form

With the ideal-gas law and the definition of virtual potential temperature,
$T_v = \theta_v \left(p/p_{00}\right)^{\kappa}$ with $\kappa = R_d/c_p$ and
$p_{00}$ = `pref0`, hydrostatic balance

$$ \frac{\mathrm{d}p}{\mathrm{d}z} = -\rho g = -\frac{p\, g}{R_d T_v} $$

becomes, after multiplying by $\kappa\, p^{\kappa-1}$, **linear in $p^{\kappa}$**:

$$ \frac{\mathrm{d}\,p^{\kappa}}{\mathrm{d}z} = -\frac{g\, p_{00}^{\kappa}}{c_p\, \theta_v}. $$

For a layer over which $\theta_v$ is constant this integrates *exactly*:

$$ p^{\kappa}(z + \Delta z) = p^{\kappa}(z) - \frac{g\, p_{00}^{\kappa}\, \Delta z}{c_p\, \theta_v}. $$

The discrete column uses this relation on the staggered grid, marching half-level pressures
$p_h$ with the full-level $\theta_v(k)$ (constant over the cell of thickness
$\Delta z_f(k)$ = `dzf(k)`), and full-level pressures $p_f$ with the interface value

$$ \theta_{v,h}(k) \;=\; \frac{\theta_v(k)\,\Delta z_f(k\!-\!1) + \theta_v(k\!-\!1)\,\Delta z_f(k)}{2\,\Delta z_h(k)}, $$

the distance-weighted interpolation to the half level ($\Delta z_h$ = `dzh`). The Exner
function follows as $\pi = (p/p_{00})^{\kappa}$. The same recurrence is used in
`initbasestate` (once, on the initial profiles) and in `fromztop` (every substep, on the
evolving slab means), so the two columns are identical whenever the inputs are identical —
this is what the buried-slab substitution below relies on.

## The anchor: `ps` at the lowest fluid level

The integration constant is supplied by the namelist parameter `ps`, imposed at the **lowest
vertical slab that contains at least one fluid cell**:

$$ k_{ps} \;=\; \min\{\,k \;:\; N_f(k) > 0\,\}, \qquad p_h(k_{ps}) = p_s \;\text{ at } z = z_h(k_{ps}), $$

where $N_f(k)$ = `IIcs(k)` is the global count of fluid cells in slab $k$. For every case
whose bottom slab contains fluid — all standard cases with streets or open ground at $z=0$ —
this reduces to $k_{ps} = k_b$ and $z_h(k_b) = 0$: the anchor is the domain bottom, exactly
as before. Only when solid geometry fully covers the domain bottom does the anchor move up
to the first level the atmosphere actually occupies.

The reason for this choice is best seen in $p^{\kappa}$ space. If the anchor sat at $z = 0$
below fully-buried slabs, the pressure at the first fluid level would be

$$ p^{\kappa}\!\left(z_h(k_{ps})\right) \;=\; p_s^{\kappa} \;-\; \sum_{k=k_b}^{k_{ps}-1} \frac{g\, p_{00}^{\kappa}\, \Delta z_f(k)}{c_p\, \theta_{v,b}(k)}, $$

i.e. it would depend on the $\theta_v$ assigned to a column of levels that contain **no
atmosphere** — whatever `prof.inp` happens to say inside the terrain. For dry Boussinesq
dynamics this offset is a pure gauge (buoyancy involves only differences from the slab
mean), but in moist runs the absolute pressure enters the Clausius–Clapeyron saturation
curve, so the invented buried column would change $q_l$ at real fluid levels. Anchoring at
$z_h(k_{ps})$ makes the pressure the flow actually feels a direct user input, and renders
the sub-surface part of `prof.inp` dynamically inert.

The startup log reports the anchor, e.g.
`Base state: thv_b = 288.000 K, ps = 101325.0 Pa at z = 20.0 m (k = 3, lowest fluid level)`.
When specifying `ps` for a terrain case, use the pressure at that height (for example, a
station pressure reduced to the elevation of the lowest open level), not sea-level pressure.

## The buried-column continuation

Below $k_{ps}$ the base-state pressures are still defined — several arrays span the full
column and the statistics files dump all levels — by marching the *inverse* recurrences
downward from the anchor:

$$ p_h^{\kappa}(k) \;=\; p_h^{\kappa}(k+1) + \frac{g\, p_{00}^{\kappa}\, \Delta z_f(k)}{c_p\, \theta_{v,b}(k)}, \qquad p_f^{\kappa}(k) \;=\; p_f^{\kappa}(k+1) + \frac{g\, p_{00}^{\kappa}\, \Delta z_h(k+1)}{c_p\, \theta_{v,h,b}(k+1)}. $$

Because the upward and downward steps are exact algebraic inverses in $p^{\kappa}$, the
telescoped identity

$$ p_h^{\kappa}(k_b) \;=\; p_s^{\kappa} \;+\; \sum_{k=k_b}^{k_{ps}-1} \frac{g\, p_{00}^{\kappa}\, \Delta z_f(k)}{c_p\, \theta_{v,b}(k)} $$

holds to machine precision, and re-marching upward from $p_h(k_b)$ recovers `ps` at the
anchor exactly. These sub-anchor pressures are **derived diagnostics** — a smooth reference
continuation, not user input — and nothing dynamical depends on them. This is verified by
the unit-test runmodes 1006 (base state, including an elevated-anchor case) and 1007
(buried-slab pressure continuation).

## The runtime column and fully-solid slabs

Every RK3 substep, `diagfld` recomputes the slab-mean state with fluid-masked averages,

$$ \langle \varphi \rangle_f(k) \;=\; \frac{1}{N_f(k)} \sum_{i,j} \varphi(i,j,k)\, \mathbb{I}(i,j,k), $$

and re-integrates the pressure column from `ps` at $k_{ps}$ using the recurrences above.
The average is **undefined** where $N_f(k) = 0$; the averaging routine (`avexy_ibm`) writes
the missing-data marker `nodata = -999.` there. That marker must never participate in
arithmetic, so every consumer substitutes a well-defined value at fully-solid slabs:

| Quantity | Empty-slab value | Rationale |
| --- | --- | --- |
| `thl0av`, `qt0av` | $\theta_{l,b}(k)$, $q_{t,b}(k)$ | base-state continuation; keeps the hydrostatic march exact |
| `ql0av` | $0$ | no cloud water in the base state |
| `thvf`, `thvh` | $\theta_{v,b}(k)$ | reference for SGS buoyancy; avoids $\sqrt{g/\theta_v}$ and divisions hitting the marker |
| `u0av`, `v0av`, `sv0av` | $0$ | the fluid-volume mean of velocity/scalars in a slab with no fluid vanishes; keeps finite differences in `lstend`, nudging and outflow sums clean |

With these substitutions the marker survives only in the statistics output (`xytdump`),
where the NetCDF `_FillValue` attribute declares it as missing data at exactly the
fully-solid levels — it is a *label*, never an operand.

## Bulk velocities: integrals, not sums of slab means

Several forcings and boundary conditions need a bulk velocity. Composing it from slab means,

$$ U_{\mathrm{wrong}} \;=\; \frac{1}{H} \sum_k \langle u \rangle_f(k)\, \Delta z_f(k), $$

is subtly wrong with an immersed boundary: writing $\phi(k) = N_f(k)/(N_x N_y)$ for the
fluid (open-area) fraction of slab $k$, the true fluid-volume mean is

$$ U \;=\; \frac{\displaystyle\sum_k \langle u \rangle_f(k)\, \phi(k)\, \Delta z_f(k)}{\displaystyle\sum_k \phi(k)\, \Delta z_f(k)}, $$

and the two agree **only when $\phi$ is constant with height**. In an urban canopy
($\phi < 1$ near the ground, $\phi = 1$ aloft) the naive sum overweights the blocked slabs —
and where $\phi(k) = 0$ it adds the `nodata` marker into the sum outright. uDALES therefore
computes bulk quantities as masked integrals that never pass through slab averages:

- **Convective outflow velocity** (`uouttot`, `vouttot`): the fluid-masked volume flux
  through the outlet plane divided by its free (fluid) area,

$$ u_{\mathrm{out}} \;=\; \frac{1}{A_f} \int_{\Omega_f} u \,\mathrm{d}A, \qquad A_f = \int_{\Omega_f} \mathrm{d}A, $$

  where $\Omega_f$ is the fluid part of the outlet plane. A fully-solid slab contributes
  zero to a sum (rather than an undefined value to an average), so no marker can appear and
  each slab is weighted by its actual open area.

- **Volume-flow-rate forcing**: the fluid-volume mean above, evaluated from per-slab sums as
  $\sum_k S(k)\, \Delta z_f(k) \,/\, \sum_k N_f(k)\, \Delta z_f(k)$ with
  $S(k) = \langle u \rangle_f(k)\, N_f(k)$ the global slab sum (the horizontal cell area
  cancels between numerator and denominator on the equidistant grid).

The same principle applies generally: **averages need an "undefined" state; sums do not.**
Whenever a bulk quantity can be written as (masked sum)/(masked sum), that form is both
exact under variable fluid fraction and immune to missing-data markers by construction.
