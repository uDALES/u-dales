# One-way nesting by velocity imposition and relaxation zones

**Status:** implemented on branch `nesting` (PR #376); validation campaign in progress, see §10.5.
**Code inspected:** uDALES `master` @ `1f8ff3e8`; DALES `v4.4_openBC` @ `aadd296d`; PALM `v23.04`.
**Scope of v1:** aligned Cartesian grids, offline parent data, one-way, velocity only; scalars keep the existing inlet and outflow treatment. Refinement
ratios $r>1$ are supported by the writer and covered by its unit tests, but are **not yet validated
end to end** — every system test so far runs at $r=1$. See the note under §10.4.

---

## 0. What "correct" means here

The zone is not required to be accurate. Errors near an imposed boundary are unavoidable — the
parent cannot know what the child is doing — and the design accepts them.

**v1 limitation: only the velocity is nested.** Scalars (temperature, humidity, passive scalars)
are not imposed from the parent: under `BCxm_nesting`/`BCym_nesting` they keep the profile inlet
plus convective outflow treatment they had before, with the single domain-mean outflow velocity
`uouttot`, so a face where the parent flow reverses extrapolates scalars rather than importing them.

> **Revised after V1 (both runs, §10.5).** This section originally named a single **decay length**
> — how far into the domain one must go before the statistics are indistinguishable from an
> unnested reference. V1 shows that conflates two error behaviours that point in opposite
> directions, and that no single number can serve both. The criterion is now split.

**Criterion A — interior fidelity, for the mean flow.** The mean-flow error does *not* decay
inward, because the boundary is where the solution is imposed and the interior is the only place
the child is free to differ. Measured, it is smallest at the wall and *grows* inward to a plateau:
$0.003 \to 0.028$ normalised across the V1 child. The right criterion is therefore a bound on the
free interior, not a distance:

$$\max_{\text{interior}} \frac{|\langle q\rangle_{\rm child} - \langle q\rangle_{\rm parent}|}{q_{\rm ref}} \le 0.05 .$$

V1 meets this ($0.028$ plateau, RMS $0.008\,u_\star$ against a $0.120\,u_\star$ sampling floor).

**Criterion B — recovery fetch, for the turbulence.** Second-order statistics *do* behave the way
the original criterion assumed: the resolved-TKE error peaks at the inner edge of the zone and
decays monotonically inward. Here a distance is the right measure — the fetch $L_{\rm rec}$ at
which the error falls to the reference's own sampling noise and stays there.

V1 measures the peak at $0.37$ just inside the zone, crossing the noise floor by $\approx2.5h$, and
still falling at $0.087$ where the far zone truncates it 13$h$ later. **The fetch is longer than the
V1 child domain**, so $L_{\rm rec}$ is bounded below but not yet measured.

**What the child is recovering from — settled by C0 (§10.5).** The V1 child carried a $\approx10\,\%$
deficit of resolved turbulence above the canopy, concentrated in the 8–64 m band. It was first read
as a limitation of the fetch alone. It is not: **the loss is at the boundary and the fetch is only
the recovery.** Boundary data stored every $\Delta t_P$ and interpolated in time carries nothing
below the wavelength $2U\Delta t_P$ (Taylor's hypothesis), and attenuates the octave above it. At
the V1 cadence of 3 s that is 21 m at $z/h=2$ — the depleted band exactly — and re-running the same
child at cadences from 0.5 s to 9 s moves the deficit monotonically from $-2\,\%$ to $-24\,\%$. The
governing number is a **dump Courant number**

$$C_{\rm dump} = \frac{U\,\Delta t_P}{\Delta x_P},$$

evaluated at the largest wind the zone sees: nothing the parent resolved is lost when
$C_{\rm dump}\le2$ (i.e. $2U\Delta t_P\le4\Delta x_P$). V1 ran at 5.4 at $z/h=2$ and $\approx7.5$
at the top. The same rule covers refinement: an $r=4$ parent supplies nothing below $4\Delta x_P$
either, and the child regenerates neither band within $13h$ (V0, §10.5). So the statement that
survives all four campaigns is: **a band that is missing from the boundary data, whatever removed
it, comes back over a fetch that grows with its wavelength; the zone width has nothing to do with
it (V2), and the child's interior extent is the only lever after the cadence itself.**
$L_{\rm rec}(\lambda)$ is the quantity a user needs, and §10.5 gives its first points.

Urban geometry helps, and helps specifically. Buildings are a strong, distributed, always-on
generator of turbulence, so structural errors injected at the boundary — reflections, interpolation
noise, spurious gradients, a near-surface shear out of equilibrium — are shredded and mixed out
over a short fetch. This is a real advantage over the flat dry convective case the DALES open-BC
paper studied, where the only turbulence source is buoyancy and the measured development fetches
are correspondingly long. **Urban nesting should be easier than that benchmark, not harder.**

The distinction matters because it sorts the concerns of §7 into two piles:

| Turbulence **does** disperse these | Turbulence does **not** disperse these |
|---|---|
| reflections off the imposed faces | net mass flux error ($\Phi\ne0$) — a global constraint |
| interpolation noise, sub-parent-scale divergence | bulk momentum error — a wrong imposed bulk flux is simply advected through |
| the *fluctuating* part of the zone's pressure signature | the *mean* part of that signature |
| a near-surface profile out of equilibrium with the child's surface | a mean scalar or temperature bias |

The errors that dilute are exactly the ones that were hardest to bound analytically; the ones that
survive are exactly the ones that can be enforced exactly and cheaply — $\Phi=0$ offline, and a
conservative flux interpolation. That is a comfortable place to be, and it is why C2 remains the
top concern while C1 and C6 relax.

---

## 1. The scheme

### 1.1 Forcing

For each velocity component $q\in\{u,v,w\}$ at its own staggered location $\mathbf{x}_q$, let
$\tilde q(\mathbf{x}_q,t)$ be the parent field interpolated onto that location and to time $t$,
and $W_q(\mathbf{x}_q)\in[0,1]$ a static weight. The continuous forcing is

$$
S_q \;=\; \frac{W_q(\mathbf{x}_q)}{\tau}\,\bigl[\tilde q(\mathbf{x}_q,t) - q(\mathbf{x}_q,t)\bigr],
\qquad
\frac{\partial q}{\partial t} \;=\; \mathcal{N}_q \;+\; S_q ,
$$

with $\mathcal{N}_q$ all existing tendencies (advection, subgrid, IBM, buoyancy, …).

$W_q$ is built per face $f$ from the inward normal distance $s_f$ **measured at $q$'s own
location**, with a strong-imposition width $L_{\rm imp}$ and a relaxation width $L_{\rm rel}$:

$$
W_f(s)=
\begin{cases}
1, & s\le L_{\rm imp},\\[2pt]
\tfrac12\bigl[1+\cos(\pi\xi)\bigr], \quad \xi=\dfrac{s-L_{\rm imp}}{L_{\rm rel}},
  & L_{\rm imp}<s<L_{\rm imp}+L_{\rm rel},\\[6pt]
0, & s\ge L_{\rm imp}+L_{\rm rel},
\end{cases}
\qquad
W_q=1-\prod_f\bigl(1-W_f\bigr).
$$

The raised cosine is $C^1$ with zero slope at both ends, so neither the guard strip nor the
trusted interior sees a kink; a quintic $1-(6\xi^5-15\xi^4+10\xi^3)$ is the $C^2$ alternative.
The product form is the bounded union: exactly $1$ if any face demands $1$, never $>1$, and
$\ge\max_f W_f$ at corners.

Finally $W_q \leftarrow W_q\cdot\mathbb{1}[\text{fluid}]$, using $II_u/II_v/II_w$ at the matching
stagger, with an additional erosion of $n_{\rm wall}$ cells away from any solid point (§5).

### 1.2 Discrete update

uDALES integrates with Wicker–Skamarock RK3 from the start-of-step state $u^n$ (`um`), with
$\Delta t_s = \Delta t/(4-s)$ for substep $s=1,2,3$ ([modtstep.f90:191](../src/modtstep.f90#L191)).
Writing $q^\ast = q^n + \Delta t_s\,\mathcal{N}_q$ for the predicted value carried in
`qm + rk3coef*qp`, the forcing is applied **exactly over the substep** rather than as an explicit
tendency:

$$
q^{\rm new} \;=\; \tilde q \;+\; \bigl(q^\ast-\tilde q\bigr)\,
   \exp\!\left(-\frac{W_q\,\Delta t_s}{\tau}\right),
\qquad
q_p \;\leftarrow\; \frac{q^{\rm new}-q^n}{\Delta t_s}.
$$

Properties:

* unconditionally stable for any $W_q,\tau,\Delta t$ — no explicit relaxation stability limit;
* $\tau\to0$ (or $W_q\Delta t_s/\tau\to\infty$) gives $q^{\rm new}=\tilde q$ exactly: the guard
  strip is the same formula, not a special case. Note the corollary — since the rate is $W_q/\tau$,
  $\tau=0$ imposes Dirichlet **wherever $W_q>0$**, not only where $W_q=1$, so it defeats the ramp
  entirely. $\tau=0$ is legal only for a pure guard strip with $L_{\rm rel}=0$, and
  `checkinitvalues` enforces that;
* $W_q\Delta t_s/\tau\ll1$ recovers the linear form $q_p \mathrel{+}= (W_q/\tau)(\tilde q-q^\ast)$;
* $W_q=0$ leaves $q_p$ **bit-identical**, so the trusted interior is untouched;
* within the zone the update is first order in time, not RK3. This is deliberate, and is the same
  compromise `bcpup` already makes for `BCxm_profile`/`BCxm_driver`.

Because the update *overwrites* $q_p$ rather than adding to it, it must be the last modification
of $q_p$ before `poisson`. See §4.

**When the target is evaluated.** $\tilde q$ is the parent interpolated at `timee`, and uDALES
advances `timee` to $t^{n+1}$ at the *first* substep (`tstep_update`), so all three substeps of
step $n$ relax towards the same target $\tilde q(t^{n+1})$ — not towards targets at the substep
times $t^n+\Delta t_s$ the formula above might suggest. This is the convention `timedep` and the
driver inflow already use, it keeps the three substeps consistent with the single projection target
`nesting_bcpup` imposes, and its error is $O(\Delta t)$ inside a zone whose update is first order in
time anyway (see the last property above). Test I6 (restart parity) pins the convention.

### 1.3 Interpolating the parent

Two choices here have more consequence than they look.

**In space — interpolate fluxes, not point values.** For aligned grids with an integer refinement
ratio, distribute each parent *face flux* over the child faces it contains rather than sampling
the parent velocity pointwise:

$$
\tilde u_{\rm child}\,A_{\rm child}\;=\;\frac{A_{\rm child}}{A_{\rm parent}}\;u_P A_{\rm parent}.
$$

for the tangential distribution, combined with **linear interpolation in each component's own
normal direction**. That combination — piecewise-constant tangentially, linear normally — is the
standard divergence-preserving prolongation, and is stronger than merely conservative: within a
parent cell each of $\partial u/\partial x$, $\partial v/\partial y$, $\partial w/\partial z$ is
constant and equal to the parent's, so

$$
(\nabla\!\cdot\mathbf{u})_{\rm child}\big|_{\rm cell}
 = (\nabla\!\cdot\mathbf{u})_{\rm parent}\big|_{\rm containing\ cell}.
$$

**A discretely solenoidal parent therefore gives an exactly solenoidal child target**, to
round-off — not merely one conservative in the parent-cell mean. Verified numerically: a parent
field built from a vector potential (div $\sim2\times10^{-17}$) interpolates to a child field with
div $\sim8\times10^{-17}$, relative $2\times10^{-16}$. Pointwise sampling gives neither property.

The consequence is worth stating plainly: **the interpolation contributes nothing for the pressure
solver to clean up.** All divergence the projection sees comes from the *blend* between parent and
child where $0<W<1$ — not from the interpolation. That materially reduces C1.

> **This guarantee describes `prolongation="constant"`, which is the writer's default — restored to
> it on 2026-09-07 after the Codex review.**
> Everything above holds only while the tangential reconstruction is piecewise constant. V0 found
> that scheme leaves a staircase in the child's mean profile — the reconstruction is flat across the
> $r$ child faces inside one parent cell where the true profile is sloped, so the error alternates
> about the parent-cell mean with period $r$ — and that the staircase is large enough to fail
> criterion A on its own ($0.126\,u_\star$ at $r=2$, $0.275\,u_\star$ at $r=4$; §10.5). W8 therefore
> made the tangential reconstruction **linear** and shipped it as the default.
> Linear tangential slopes are not divergence-preserving: the $x$-, $y$- and $z$-derivative terms no
> longer cancel cell by cell, so the identity above fails and the child's projection has to remove
> the remainder.
>
> **Measured, on the V0 coarse arm's own 4 m LES field** (a 128 m window, $|u|_{\max}=7.2$ m/s), the
> two reconstructions differ by seven orders of magnitude in the child target's divergence:
>
> | | child divmax, $r=2$ | as a fraction of $u/\Delta x_{\rm child}$ |
> |---|---|---|
> | parent itself | $1.27\times10^{-7}$ | — |
> | `constant` | $1.27\times10^{-7}$ (identical) | $3.5\times10^{-8}$ |
> | `linear` (shipped) | $1.74\times10^{-1}$ | $4.8\times10^{-2}$ |
>
> ($r=4$: $2.6\times10^{-1}$, $3.6\times10^{-2}$. The parent's own $1.27\times10^{-7}$ is the
> single-precision floor of the dumps, not a physical divergence.) So the linear scheme asks the
> child's pressure solve to absorb a divergence source of order 5 % of $u/\Delta x$ inside the zone,
> every timestep. Parent-face conservation, $\Phi$ and the post-projection `divmax` are all still
> clean, which is exactly why none of the *flux* assertions notice — but
> `test_v0_tiny.py::test_the_prolongation_reproduces_the_parent_divergence` and
> `test_a_solenoidal_field_stays_solenoidal_through_the_round_trip` do, and both fail with the linear
> default. They were not run when it was introduced.
>
> **The default is therefore back to `constant`** while the trade is settled: a default must not
> violate the contract this document states and the suite asserts, and `constant` is additionally
> the only one of the two validated at production size. `linear` remains available and unit-tested.
>
> **Neither scheme is yet validated end to end as shipped:** V0's production numbers were produced
> with `constant` at the old 3 s linear-in-time cadence, and `linear` has never been run at
> production size at all. The two are on opposite sides of a real trade — one exact in the
> divergence and wrong in the mean, the other the reverse — and the choice is being settled by
> measurement (plan §7, R2), not by argument. Until then this section should be read as describing
> the *option*, not the default. The proper resolution, if exact local solenoidality is judged
> necessary, is a constrained (Balsara-style) linear reconstruction that is divergence-preserving by
> construction; that is scoped as separate work rather than rushed.

**In time — the cadence is a low-pass filter, and it is the dominant error of the whole scheme.**
This subsection originally discussed the time interpolant only as a smoothness question (linear is
$C^0$, cubic is $C^1$, the jump in $\partial\tilde q/\partial t$ at each crossing gives a pressure
transient at $1/\Delta t_P$). That is true and secondary. The primary effect is that data stored every
$\Delta t_P$ and interpolated between levels cannot carry frequencies above $1/(2\Delta t_P)$, and
attenuates those below it by the interpolant's transfer function — $\mathrm{sinc}^4(f\Delta t_P)$ for
linear interpolation. Under Taylor's hypothesis, $f = U/\lambda$, so at a face where the wind is $U$
every wavelength below $2U\Delta t_P$ is absent from the target and the octave above it is damped.
Computed from the V1 parent's own spectra at 3 s: at $z/h=2$ the target keeps 0 % of the 8–16 m
variance, 21 % of 16–32 m, 66 % of 32–64 m and 90 % of 64–128 m, and the child's interior ratios
0.82, 0.85, 0.95, 1.04 order band by band with that. **The requirement is $C_{\rm dump}=U\Delta t_P/\Delta x_P\le2$**
(§0), which for uDALES-to-uDALES nesting at $\Delta x_P=2$ m and $U=5$ m/s means $\Delta t_P\le0.8$ s —
every one or two parent steps, not the tens of seconds the DALES open-BC guidance ("$\le30\times$
temporal refinement") suggested. That guidance came from a convective boundary layer whose large
eddies are slow; it does not transfer to a sheared urban layer.

**The interpolant matters more than expected, and Catmull–Rom is now the default.** Measured on the
same child at the same 3 s cadence (C0, §10.5): switching `nest_timeinterp` from linear to the
unlimited cubic halves the interior deficit ($-9.9\to-5.7\,\%$) and lifts the 8–16 m ratio from 0.83
to 0.89 — a band that is *entirely* above the target's Nyquist frequency and so cannot be helped by a
flatter passband directly. What the cubic supplies is the 16–64 m band, which the cascade then uses
to rebuild 8–16 m faster. The cubic costs one extra buffer slot and nothing in I/O.

**The Hermite slopes must be unlimited, and this is a constraint, not a preference.** §3.1(2)
needs the interpolant to be *linear in the data*: if
$\tilde q_f(t)=\sum_k c_k(t)\,q_{f,k}$ with coefficients $c_k$ that depend only on $t$ and the
level spacings — the same for every boundary face $f$ — then

$$\Phi(t)=\sum_f A_f\,\hat n\!\cdot\!\tilde q_f(t)=\sum_k c_k(t)\,\Phi(t_k)=0,$$

because the offline correction already zeroed $\Phi$ at every stored level. Linear interpolation
and unlimited Catmull–Rom both have this form. A *monotone* limiter does not: the Fritsch–Carlson
slope is a weighted harmonic mean of the one-sided slopes, plus a sign-switch branch, so $c_k$
becomes a function of the local data and each face gets its own weighting. The cancellation then
collapses. Measured, with every stored level corrected to $|\Phi|\approx10^{-14}$:

| interpolant | linear in the data? | $\max|\Phi|$ across the interval |
|---|---|---|
| linear | yes | $1.1\times10^{-14}$ |
| Catmull–Rom (unlimited Hermite) | yes | $1.4\times10^{-14}$ |
| Fritsch–Carlson (monotone) | **no** | $1.0\times10^{0}$ |

In the solver the monotone version produced a runtime $\Phi$ of $5.2\times10^{-5}$ and a `divtot`
of $1.7$, and aborted at the first substep under `nest_lfluxassert`. Monotonicity is the wrong
property to ask for here anyway — velocity is sign-unconstrained, and a small overshoot in an
imposed boundary value is harmless, whereas losing compatibility is not. Test **U44** pins the
linearity directly.

### 1.4 How wide does the zone have to be?

Three criteria; take the maximum. Throughout, $N$ denotes a width in child cells,
$\tau=n_\tau\Delta t$, and $C=U\Delta t/\Delta x$ is the local advective Courant number
(uDALES's `courant` is a *sum-of-directions* limit of 1.1–1.5, so a representative $C$ for the mean
flow near the boundary is ≈0.5).

**(a) Imposition strip — numerics only.** Momentum advection is second-order central and needs one
halo ([modadvection.f90:158](../src/modadvection.f90#L158)); the subgrid terms need one. Two cells
suffice; **use $N_{\rm imp}=3$** — one for the stencil, one of margin, one to cover $ihc=2$ if
scalars are added later. More buys nothing numerically and costs domain.

**(b) Absorption — easily satisfied, and *not* the binding constraint.** A disturbance advected out
through the zone is attenuated by $\exp(-\mathcal{D})$ with optical depth

$$
\mathcal{D}=\frac{1}{U\tau}\int W\,ds=\frac{L_{\rm imp}+\tfrac12 L_{\rm rel}}{U\tau}
\quad\Longrightarrow\quad
\mathcal{D}=\frac{N_{\rm imp}+\tfrac12 N_{\rm rel}}{n_\tau C},
$$

(the raised cosine has mean $\tfrac12$). A reflected signal makes the round trip, so its amplitude
at the inner edge is $\sim\exp(-2\mathcal{D})$; 99% suppression needs only $\mathcal{D}\gtrsim2.3$.
With $n_\tau=4$, $C=0.5$ that is $N_{\rm imp}+\tfrac12 N_{\rm rel}\gtrsim5$, i.e. $N_{\rm rel}\gtrsim4$.
**Absorption is cheap.** What actually sets the width is (c).

**(c) Not scattering off the zone's own inner edge.** A ramp of length $L$ scatters a wave of
wavelength $\lambda$ only weakly when $L\gtrsim\lambda$. Two scales matter: the shortest resolved
motions ($\lambda\sim4\Delta x$, giving $N_{\rm rel}\gtrsim8$) and the energetic eddies near the
boundary ($\lambda\sim h$, the building height, or the shear-layer depth). So

$$
L_{\rm rel}\;\gtrsim\;\max\bigl(8\Delta x,\;\tfrac12\ \text{to}\ 1\times h,\;2\Delta x_{\rm parent}\bigr),
$$

the last term ensuring the transition is resolved in the *parent's* own terms (with refinement
$r\le4$ that is $\le8$ cells, so rarely binding).

**Recommended starting point:** $N_{\rm imp}=3$, $N_{\rm rel}=8$–$12$, i.e. a **total of about
11–15 cells**, with 12 as the default. Reassuringly, that is the same order as decades of
mesoscale practice — WRF's `spec_bdy_width` is 5, COSMO-style Davies zones are ~10–15 — for
essentially the same job.

**The real knob is $\tau$, not the width.** Because the update of §1.2 is unconditionally stable,
$\tau$ can be made arbitrarily small, so width and stiffness are freely interchangeable: a narrow
stiff zone is numerically fine. The cost of stiffness is *physical* — more turbulence damping and
a larger misfit at the inner edge — so **choose $\tau$ from how much damping you can tolerate,
then read the width off (b) and (c).**

**Domain cost, given a building-free zone.** The margin is
$N_{\rm imp}+N_{\rm rel}+n_{\rm wall}\approx13$–$16$ cells per side, i.e. ~26–32 cells consumed in
each direction: ~5–6% of a 512-cell domain, ~10–12% of 256, but 20–25% of a 128-cell domain — at
which point nesting has stopped paying for itself. Rule of thumb: **the child should be at least
~20× the zone width across, i.e. $\gtrsim250$ cells.**

**Specify it in metres, not cells.** The zone width is set by the eddy scale, so halving $\Delta x$
leaves it unchanged in metres and doubles the cell count. `nesting_init` should report the
equivalent cell count and warn if it falls below 6 cells (under-resolved transition) or exceeds
15% of the domain.

**On buildings in the zone.** The design rule is **no buildings anywhere $W>0$**, plus the
$n_{\rm wall}$ erosion — so the first solid cell may appear at $N_{\rm imp}+N_{\rm rel}+1$. The
underlying requirement is really that *the parent must be consistent with whatever geometry sits
in the zone*; for a mesoscale parent, which resolves no buildings at all, that is equivalent to
the rule above. (Strictly, self-nesting from a uDALES parent resolving the same geometry would
permit buildings in the zone — but the general rule is kept, since it is the only one that also
holds for a coarse parent and it keeps the Big Brother validation clean.)

---

## 2. uDALES's discrete operators

Let $V_k = \Delta x\,\Delta y\,\Delta z_f(k)$ be the cell volume. Define the discrete divergence
$\mathcal{D}$, gradient $\mathcal{G}$ and the solver's Laplacian $\mathcal{L}$ as they are
actually coded:

$$
(\mathcal{D}\mathbf{u})_{ijk}=\frac{u_{i+1}-u_{i}}{\Delta x}
 +\frac{v_{j+1}-v_{j}}{\Delta y}
 +\frac{w_{k+1}-w_{k}}{\Delta z_f(k)}
\tag{fillps, modpois.f90:966--973}
$$

$$
(\mathcal{L}p)_{ijk}=\rho_f(k)\!\left[\frac{\delta^2_x p}{\Delta x^2}+\frac{\delta^2_y p}{\Delta y^2}\right]
+\frac{1}{\Delta z_f(k)}\!\left[\rho_h(k{+}1)\frac{p_{k+1}-p_k}{\Delta z_h(k{+}1)}
-\rho_h(k)\frac{p_k-p_{k-1}}{\Delta z_h(k)}\right]
\tag{modpois.f90:154--157, 202--211}
$$

$$
(\mathcal{G}p)_x = \frac{p_i-p_{i-1}}{\Delta x},\quad\text{etc.}
\tag{tderive, modpois.f90:1044--1055}
$$

The substep solves, with $\mathbf{u}^\ast$ the predicted velocity after all forcing,

$$
\mathcal{L}p=\frac{1}{\Delta t_s}\mathcal{D}\mathbf{u}^\ast,
\qquad
\mathbf{u}^{n+1}=\mathbf{u}^\ast-\Delta t_s\,\mathcal{G}p .
$$

With $\rho\equiv1$ (see F1 below) $\mathcal{L}=\mathcal{D}\mathcal{G}$ exactly, so
$\mathcal{D}\mathbf{u}^{n+1}=0$ wherever the solve is exact.

**Four facts from the source that the rest of this document depends on.**

**F1.** $\rho_f=\rho_h\equiv1$: `rhobf`/`rhobh` are allocated to `1.`
([modfields.f90:571](../src/modfields.f90#L571)) and **never assigned anywhere else in `src/`**.
uDALES is Boussinesq in practice. Note also that $\mathcal{L}$ carries $\rho$ but $\mathcal{D}$
does not — inert today, silently wrong the day anelasticity is switched on. Keep the $\rho$
factors in all new code; raise the `fillps` inconsistency as a separate issue against `master`.

**F2.** *Imposed boundary-normal velocity survives the projection exactly.* `bcp` sets
$p_{0,j,k}=p_{1,j,k}$ and $p_{itot+1}=p_{itot}$ for non-periodic $x$
([modboundary.f90:1377](../src/modboundary.f90#L1377)), so $(\mathcal{G}p)_x=0$ on the boundary
face; and `tderive`'s loop runs $i=ib..ie$, never touching $u_{itot+1}$. Homogeneous Neumann
pressure $\Leftrightarrow$ the boundary flux is *locked*.

**F3.** *The singular system is silently regularised.* For non-periodic laterals the horizontal
eigenvalues are Neumann (DCT, `FFTW_REDFT10`, [modpois.f90:111](../src/modpois.f90#L111)) and the
vertical operator is Neumann–Neumann, so $\mathcal{L}$ is singular with null space the constants.
uDALES removes the null space by replacing the top row of the **horizontal-mean mode only** with
a Dirichlet-across-the-top-cell coefficient ([modpois.f90:205](../src/modpois.f90#L205)):

$$
b^{\rm N}_{ktot}=-a_{ktot}, \qquad b^{\rm D}_{ktot}=-a_{ktot}-2c_{ktot},
\qquad \tilde b_{ktot}=b^{\rm D}_{ktot}\ \ \text{iff}\ \ \lambda_{ij}=0 .
$$

The modified system is non-singular, so **an incompatible right-hand side produces no error and
no warning.**

**F4.** `masscorr` is already a no-op under any non-periodic BC (every branch guarded by
`.not.linoutflow`, [modforces.f90:352](../src/modforces.f90#L352)), and `linoutflow` is set by
`checkinitvalues` for non-periodic $BC_{xm}/BC_{ym}$. No conflict with the zone forcing.

---

## 3. Mass compatibility, and what the lid does

### 3.1 The compatibility condition

Scale row $(i,j,k)$ of $\mathcal{L}$ by $V_k$. Every term becomes a difference of face fluxes with
a shared face coefficient, so $A=V\mathcal{L}$ is symmetric and, with Neumann on all six faces,
$A\mathbf{1}=0$. Hence $\operatorname{range}(A)=\mathbf{1}^\perp$ and the left null vector of the
unscaled operator is $V_k$. Solvability requires

$$
\sum_{ijk} V_k\,(\mathcal{D}\mathbf{u}^\ast)_{ijk}=0 .
$$

The divergence telescopes, leaving only domain-boundary faces:

$$
\Phi \;\equiv\;
\Delta y\!\sum_{j,k}\!\Delta z_f\bigl[u^\ast_{itot+1}-u^\ast_{1}\bigr]
+\Delta x\!\sum_{i,k}\!\Delta z_f\bigl[v^\ast_{jtot+1}-v^\ast_{1}\bigr]
+\Delta x\Delta y\!\sum_{i,j}\bigl[w^\ast_{ktot+1}-w^\ast_{1}\bigr]
\;=\;0 .
$$

**$\Phi$ is the net volume flux of the predicted velocity through the domain boundary, and
nothing else.** Two consequences that shape the whole architecture:

1. **The interior relaxation zone cannot break compatibility, however strong it is.** Only the
   values on the boundary faces can. So the constraint belongs in `bcpup`, where those values are
   set explicitly and exactly.
2. $\Phi$ is a **linear** functional of the boundary data. Therefore if every stored parent time
   level satisfies $\Phi=0$, so does every linearly time-interpolated target, and **no
   per-substep correction is required.** Correct the input once, at initialisation. (This is why
   DALES needs its per-substep `radcorrection` and we do not: their outflow faces are free
   — radiation — whereas ours are imposed.)

With immersed boundaries, solid faces carry $u\equiv0$ (`ibmnorm`,
[modibm.f90:716](../src/modibm.f90#L716)), so both $\Phi$ and any correction must be summed over
**fluid faces only** and normalised by the **fluid** face area. PALM does exactly this — its flux
sum, its correction and its `face_area` all carry the same `topo_flags` mask
(`pmc_interface_mod.f90:6238`, `:6379`, `:2444`). Masking the correction but dividing by the
geometric area leaves a residual.

### 3.2 What happens when $\Phi\ne0$ — the lid analysis

This is the decision that needs settling, and it turns entirely on whether $w_{ktot+1}$ is a free
variable.

Consider the horizontal-mean mode $\bar p_k$, which is the only mode with $\lambda_{ij}=0$ and
hence the only one whose operator F3 modifies. Its vertical system in flux form is

$$
\frac{1}{\Delta z_f(k)}\bigl[F_{k+1/2}-F_{k-1/2}\bigr]=\bar R_k,
\qquad F_{k+1/2}=\rho_h(k{+}1)\frac{\bar p_{k+1}-\bar p_k}{\Delta z_h(k{+}1)} ,
$$

with $F_{1/2}=0$ at the bottom ($a_1=0$, Neumann). Summing over $k$ weighted by $V_k$:

$$
\sum_k V_k \bar R_k \;=\; \Delta x\,\Delta y\; F_{ktot+1/2}.
$$

Under the **Neumann** top row $b^{\rm N}$ we would have $F_{ktot+1/2}=0$, hence the solvability
condition $\Phi=0$ — and no solution at all when $\Phi\ne0$. Under the **Dirichlet** top row
$b^{\rm D}$ actually used, the ghost value is $\bar p_{ktot+1}=-\bar p_{ktot}$, i.e. $\bar p=0$ on
the top *face*, giving

$$
F_{ktot+1/2}=-2c_{ktot}\,\Delta z_f(ktot)\,\bar p_{ktot}
\quad\Longrightarrow\quad
\bar p_{ktot} \;\propto\; -\,\Phi .
$$

So the solver always returns an answer, and the incompatibility is expressed as a **non-zero mean
pressure at the lid**.

Formally: with $\tilde{\mathcal{L}}$ the modified operator, $\mathcal{L}=\mathcal{DG}$ under the
BCs actually applied to the velocity, and $\tilde{\mathcal{L}}p=\mathcal{D}\mathbf{u}^\ast/\Delta t_s$,
the post-projection divergence is
$\mathcal{D}\mathbf{u}^{n+1}=\Delta t_s(\tilde{\mathcal{L}}-\mathcal{L})p$, which is non-zero **only** in
the horizontal-mean mode of the top cell layer. (The same pin is active in ordinary periodic runs,
where $\Phi\equiv0$, so $\bar p_{ktot}=0$ and nothing happens — it is not a nesting-specific artefact.)

Whether that pressure moves any mass depends on the top BC:

**Case A — `BCtopm_freeslip` (rigid lid).** `boundary` forces $w_{ktot+1}\equiv0$
([modboundary.f90:174](../src/modboundary.f90#L174)), `bcpup` sets $w^\ast_{ktot+1}=0$, and
`tderive` updates $w_p$ only for $k=2..ktot$. The flux $F_{ktot+1/2}$ that the pressure solution
implies is therefore **never realised as a velocity**. The projection removes the divergence
everywhere except the top cell layer, which retains a horizontally uniform residual

$$
\bigl(\mathcal{D}\mathbf{u}^{n+1}\bigr)_{ij,ktot}
=\frac{\Phi}{\Delta x\,\Delta y\,\Delta z_f(ktot)\;itot\;jtot} \;\ne\;0 .
$$

There is **no breathing at all**: the imbalance becomes a spurious volume source in the top layer,
regenerated every substep, silently (F3). The scheme is exact if and only if $\Phi=0$ is enforced.

**Case B — `BCtopm_pressure` (leaky lid).** `bcpup` sets
$w^\ast_{ktot+1}=w^n_{ktot+1}/\Delta t_s+2\langle p^{\rm acc}\rangle_{ktot}/\Delta z_h(ktot{+}1)$
and `tderive` adds the matching increment to $w_p$
([modboundary.f90:1236](../src/modboundary.f90#L1236),
[modpois.f90:1057](../src/modpois.f90#L1057)), where $p^{\rm acc}$ is the *accumulated* pressure
`pres0`. So

$$
w_{ktot+1}\;\propto\;\langle p^{\rm acc}\rangle_{ktot},
\qquad
p^{\rm acc}\;=\;\textstyle\sum_{\rm substeps} p ,
$$

i.e. mass excess $\to$ mean pressure increment $\to$ accumulated pressure $\to$ outflow at the
lid. This is genuine breathing, and it is the mechanism `checkinitvalues` relies on when it
force-switches every non-periodic run to `BCtopm_pressure`
([modstartup.f90:845](../src/modstartup.f90#L845)). But it is an *integrator driving a flux*: a
second-order feedback whose timescale and damping are set by $\Delta z_h(ktot{+}1)$ and the domain,
not by us, and which can ring. And physically the child then exchanges mass with a fictitious
reservoir rather than with its parent.

**Case C — nest the lid.** Impose the parent's $w$ at $k=ktot{+}1$ as a fifth Dirichlet face.
Then $\Phi$ includes a known top term and the correction is distributed over all five faces —
PALM's choice (`pmc_interface_mod.f90:6340`–`6375`). Physically the most consistent for a nesting
scheme; requires the parent's $w$ at the child lid to be trustworthy, and interacts with `grwdamp`
(§7 C8).

**Recommendation, for discussion.** A for v1, because the point of nesting is that the parent sets
the flux, and A makes the constraint exact, cheap and enforceable offline; with a hard runtime
assertion on $\Phi$ (F3 guarantees the solver will not complain on its own). B stays available as
a namelist option for a mesoscale parent whose column budget is not closed. C is the natural v2
once the lateral scheme is trusted. This choice is not load-bearing for anything else in the
design — A, B and C differ only in which faces enter $\Phi$ and which are corrected.

### 3.3 The three options in the brief

1. **Adjust the imposed normal velocity so $\Phi=0$** — yes, but *once at initialisation, per
   input time*, distributing a uniform normal-velocity increment over the fluid lateral faces in
   proportion to fluid face area. This is DALES's `openboundary_divcorr`
   (`modopenboundary.f90:444`–`576`). By §3.1(2) that suffices for all times.
2. **Subtract the global mean of the RHS** — no. The operator is *not* singular (F3), so this
   would not be a projection onto the range; it would add a spurious uniform source. Use
   $\sum V_k R$ only as a diagnostic.
3. So: physical flux correction offline, numerical null-space treatment left alone, both
   monitored at runtime.

**Do we need to project the target field on the child grid?** No — and with the
divergence-preserving interpolation of §1.3 the reason is stronger than first thought: the target
is *already* discretely solenoidal wherever the parent was, so there is nothing to project. What
remains is the divergence generated by blending parent and child across $0<W<1$. Make the *magnitude* of
that cleanup a first-class diagnostic (§6). If it becomes comparable to the physical pressure
signal — most likely with a strongly coarsened parent — pre-project the target in the Python
preprocessing tool (one Poisson solve per input time, offline), rather than adding machinery to
the solver.

---

## 4. Where it hooks into the timestep

```
tstep_update                                                     program.f90:134
timedep                                                                     :136
+ nesting_update_target        advance the two-slab buffer; interpolate in time
advection  subgrid  bottom  coriolis  forces  lstend  nudge             :142-164
ibmwallfun  periodicEBcorr  masscorr(no-op, F4)  ibmnorm               :166-171
EB  vegetation  heatpump  scalsource  fixuinf2  fixuinf1               :173-186
grwdamp                                                                     :191
+ nesting_apply                <== zone relaxation of up,vp,wp (§1.2), masked
poisson                                                                     :193
  fillps -> bcpup  <== + nesting_bcpup: Dirichlet on the outer faces (F2)
tstep_integrate  halos  checksim  fielddump  statsdump                 :197-205
boundary           <== + case(BCxm_nesting): ghost planes from the parent    :207
```

**Why there.** After all physical tendencies, so imposition is authoritative. After `grwdamp`,
which relaxes toward *slab means* ([modboundary.f90:1447](../src/modboundary.f90#L1447)) and would
otherwise fight the zone. After `ibmnorm`, so `um` is already zero in solids — but mask
explicitly regardless. Before `poisson`, so the projection cleans the divergence the forcing
introduces.

**Units and bookkeeping.** At that point $q_p$ is a tendency in m s$^{-2}$ accumulated since the
last `tstep_integrate`, and $q^\ast=q^n+\Delta t_s q_p$. Since §1.2 *overwrites* $q_p$, nothing may
be inserted between `nesting_apply` and `poisson`. Comment this at both call sites.

**Halos.** Momentum advection is second-order central
([modadvection.f90:158](../src/modadvection.f90#L158)) and needs **one** ghost plane of
$u_0,v_0,w_0$ — $i_b{-}1$, $i_e{+}1$, $j_b{-}1$, $j_e{+}1$ — filled by the nesting branch of
`boundary`, exactly as `xmi_driver` does ([modboundary.f90:720](../src/modboundary.f90#L720)).
Before the projection, `fillps` additionally needs the far faces $pup_{ie+1}$, $pvp_{je+1}$;
interior halos are exchanged inside `bcpup` ([modboundary.f90:1218](../src/modboundary.f90#L1218)),
the domain-boundary values are not and must be set by the new branch.

**Reusing the existing BC machinery.** Add `BCxm_nesting = 4`, `BCym_nesting = 3` and one branch
each in `boundary` (ghosts), the outflow section (impose rather than convect), and `bcpup`
(faces). No new mechanism is required — `BCxm_profile` ([modboundary.f90:1256](../src/modboundary.f90#L1256))
and `BCxm_driver` ([modboundary.f90:1283](../src/modboundary.f90#L1283)) already do exactly this
job for other data sources. One trap: `checkinitvalues` force-switches
`BCtopm := BCtopm_pressure` for those two modes ([modstartup.f90:845](../src/modstartup.f90#L845));
the nesting mode must not inherit that if Case A is chosen.

**Solver restriction.** Only `ipoiss = POISS_FFT2D` supports non-periodic laterals:
`POISS_FFT2D_2DECOMP` uses an r2c transform in $x$ and a periodic eigenvalue set in $y$ (its own
comment says the non-periodic branch "is not correct", [modpois.f90:246](../src/modpois.f90#L246));
`POISS_FFT3D` is triply periodic; `POISS_CYC` is commented out. Add a hard `checkinitvalues`
guard, not a silent fallback.

---

## 5. Immersed boundaries

There are **no conservative fluid fractions in uDALES**. $II_c/II_u/II_v/II_w$ are integer 0/1
masks built from solid-point lists ([modibm.f90:2122](../src/modibm.f90#L2122)); `mask_u/v/w/c`
are their real-typed twins; `secareas` are areas of *facet sections*, not cell face fractions;
`lconservativeibm` selects a scalar-advection correction, not a volume fraction. The brief's
premise does not hold, so:

1. mask $W_q$ by $II_u/II_v/II_w$ at the matching stagger;
2. additionally zero $W_q$ within $n_{\rm wall}$ cells of any solid point (default 1) — a fluid
   face adjacent to a wall carries that wall's shear signature, and forcing it toward a parent
   value that knows nothing about the wall injects spurious momentum exactly where `ibmwallfun`
   is setting the stress;
3. prefer a **building-free zone**: count solid points inside the zone at init and warn loudly
   (error under a strict flag). For an urban child this is a modelling constraint as much as a
   numerical one — the region where the parent is trusted should not contain geometry the parent
   cannot represent;
4. run $\Phi$ and its correction over fluid faces with fluid face areas (§3.1).

`ibmnorm` already zeroes both `um` and `up` at solid points, so an error in (1) would be partly
and silently masked; (2) is not covered by it at all.

---

## 6. Parent data, I/O, and diagnostics

### 6.1 What has to be in the file (container-independent)

The *contract* matters; the container does not. Required content, whatever the format:

* the parent velocity on **zone slabs only** — four lateral slabs of the zone thickness, plus an
  optional top slab (Case C) — each component sampled at **its own stagger**, plus an optional
  full-3D block used once for cold-start initialisation;
* a time axis, and enough of the child grid (`xh, xf, yh, yf, zf, zh`, `itot, jtot, ktot`,
  `xlen, ylen`) to be **validated against the run at init and abort on mismatch**;
* $\rho_f,\rho_h$ — the density the fluxes were built with;
* an explicit stagger tag per variable, so the file cannot be silently misread;
* provenance: `parent_model`, `parent_dx`, `parent_dt`, `child_origin_x/y`, `rotation_deg`
  (constrained to 0 in v1, marking the extension point);
* `divergence_corrected` and the residual $\Phi(t)$ left by the offline correction (§3.3).

The validate-against-the-run requirement is the important one. The existing driver path has no
such check and instead prints warnings telling the user to make `ylen`, `jtot` and `nprocy` match
by hand ([moddriver.f90:765](../src/moddriver.f90#L765)); that is the failure mode to design out.

### 6.2 Volume and read frequency

Per parent time level, with $n_z$ the zone thickness in cells and 8-byte reals,

$$
V \;\approx\; 3\,\bigl[\,2 n_z\,j_{tot}\,k_{tot} + 2 n_z\,i_{tot}\,k_{tot}\bigr]\times 8\ \text{B}
 \;=\; 48\,n_z\,k_{tot}\,(i_{tot}+j_{tot})\ \text{B}.
$$

For $i_{tot}=j_{tot}=512$, $k_{tot}=128$, $n_z=16$: **≈100 MB per time level**, so ≈6 GB per hour
at $\Delta t_P=60$ s, ≈36 GB at 10 s, ≈360 GB at 1 s. **The parent time interval, not the spatial
extent, is what sets the cost — and the interval is not free to choose.** §1.3 and §0 set
$C_{\rm dump}\le2$, i.e. $\Delta t_P\le2\Delta x_P/U_{\max}$: 0.8 s for this case at 5 m/s, so
≈450 GB per simulated hour of boundary data, and the design's original 60 s sizing (143 GB for a
production window) was under-sampled by a factor of 75. The $\le30\times$ temporal-refinement
guidance imported from the DALES open-BC paper is withdrawn for the reason given in §1.3. Storing
single precision halves all of this at no meaningful accuracy cost (the data are a boundary
condition and are about to be interpolated); promote to double on read. uDALES already has a
`SINGLE_PRECISION_OUTPUT` build option. Three consequences:

- **Full-domain field dumps are not a production path.** V1's 170 GB of 3-D dumps covered 3 h at
  3 s; at 0.5 s the same window is 1 TB. The parent must write the zone slabs and the initial block
  only (a `lnestdump` switch in the parent, item D1 in the plan), which is $n_z(i_{tot}+j_{tot})$
  cells per level against $i_{tot}j_{tot}$ — a factor of $\approx15$ for the production case.
- **A coarser cadence is a choice with a measured price**, not a free parameter: the deficit curve
  of §10.5 (C0) and the recovery fetch of §10.5 (V2) together say how much interior a child needs
  for a given $C_{\rm dump}$.
- **A mesoscale parent at 10–60 s output cannot drive this scheme** at LES scales; see §9.4.

Crucially, a read happens **once per parent interval, not per timestep**: with $\Delta t\approx0.3$ s
and $\Delta t_P=60$ s that is one read per ≈200 steps, ≈1.6 MB per rank on 64 ranks. The
arithmetic is comfortable. The design target is therefore *"do nothing pathological"*, not
*"maximise bandwidth"* — and the pathologies are specific and avoidable:

1. **Non-contiguous per-rank reads.** Order the bulk data so that each rank's slab is a **single
   contiguous run** — i.e. the decomposed index ($j$ for the west/east slabs, $i$ for
   south/north) must be the slowest-varying spatial index, giving a layout
   `[face][component][time][j][k][n_z]`. This is the single most important performance property
   and it is independent of the container.
2. **A metadata storm.** $N$ ranks each opening and seeking in one file on a shared Lustre/GPFS
   filesystem can cost far more than the bytes do. `AGENTS.md` already warns about this for large
   NetCDF/HDF5 work on this cluster.
3. **Reading on demand.** Trigger the read for level $n{+}2$ as soon as the model *enters*
   interval $n\!\to\!n{+}1$, into a third buffer, so the latency is hidden behind a whole parent
   interval rather than stalling one substep. Double buffering already gives the slots; only the
   trigger point has to be early.

**As implemented (D1, `src/modnestdump.f90`).** The parent writes the band with `&NESTDUMP`
(`lnestdump`, `tnestdump`, the child box `nestdump_x0/y0/xsize/ysize` in parent coordinates on
parent faces, `nestdump_nzone` parent cells inside each lateral face, `nestdump_linit`): each rank
writes the intersection of its subdomain with each of the four strips of the band, `u0/v0/w0` in
single precision with the upper staggered faces included, to its own `nestdump.<ipx>.<ipy>.<expnr>.nc`,
plus its part of the whole box once to `nestdump_init.<ipx>.<ipy>.<expnr>.nc` at the first dump
(spec section 9). The band is the child's guard + ramp on the parent grid, rounded up, **plus one
cell**, which is what the linear tangential reconstruction of §1.3 needs at the outermost zone cell
of a refined child. The data are raw: the flux correction of §3 and the initial-condition projection
stay in `udprep.nesting`, which `make_child_case --source nestdump` feeds exactly as it feeds the
full dumps -- on the `tiny` preset the two sources give bit-identical nesting files, the tiny
geometry's band files are 4.0x smaller than its field dumps (5.1x by cell count; the ~15x above is
the production geometry), and the parent's own accounting printed 35 MB in 0.06 s of write calls
for 40 dumps on 4 ranks (`tests/validation/nesting/test_nestdump_tiny.py`).

### 6.3 Container options

| Option | Speed | Introspectable | Notes |
|---|---|---|---|
| Rank-0 read + `MPI_Scatterv` | Good | yes | One file handle, no metadata storm, container-agnostic, easiest to get right. Costs one copy and a collective per parent interval — irrelevant at this frequency. |
| Raw stream binary + small self-describing sidecar | Best | via sidecar + Python reader | One `MPI_File_read_at_all` per rank per level; trivially contiguous. Not readable with `ncdump`. |
| NetCDF-4, chunk shape matched to the decomposition | Good if chunked right, poor if not | yes | Parallel collective read only if the build has parallel NetCDF; otherwise independent reads and pathology 2. |

**Recommendation at design time** was rank-0 read + scatter behind a container-agnostic
`nesting_read_slab` interface. **What is implemented** is the third row: NetCDF-4, every rank
opening the file read-only and reading its own contiguous hyperslab of each slab with
`nf90_get_var` (`modnestingio`), with no scatter and no parallel NetCDF. Four levels are buffered
for the Hermite stencil; at a parent-interval crossing the slots are rolled and the one new level is
read *synchronously, at the crossing, in the substep that crosses* — there is no read-ahead or
prefetch, contrary to what §6.2 item 3 asks for. This was measured rather than argued (§10.7
item 3): the per-crossing read is small against a parent interval at 64 ranks, so the "metadata
storm" of §6.2 did not materialise at that scale and the simpler design stands. Stage the file to
node-local/ephemeral storage before the run when it fits, as `AGENTS.md` already recommends.

**Measure it.** `nesting_stats` reports the read wall time per parent interval and its fraction of
total runtime, with a warning above 1%. That turns "does the I/O slow us down?" from a design
argument into a number.

### 6.4 Diagnostics

Nesting fails quietly (F3), so instrument from the first commit:
$\Phi$ (must be at round-off — hard assertion); `divtot`/`divmax`, which `chkdiv` already computes
([modchecksim.f90:161](../src/modchecksim.f90#L161)) and whose drift is the Case-A symptom;
zone energy injection $\sum_{\rm zone}\mathbf{u}\cdot\mathbf{S}\,dV$ split guard/relaxation;
$\|\mathcal{G}p\|$ in the zone relative to the interior (the answer to §3.3);
the misfit $\mathrm{rms}(\tilde q-q)$ at the inner zone edge; and resolved TKE on planes normal to
the inflow, to measure the turbulence development fetch.

The **headline acceptance metrics are the two of §0**, and which one applies depends on the
statistic. First-order quantities (mean profiles, facet stresses) take criterion A, an interior
fidelity bound — their error grows inward from the imposed boundary, so a distance is meaningless
for them. Second-order quantities (resolved TKE, spectra) take criterion B, the recovery fetch, in
metres and in building heights. Report both per validation case. Everything else in this table is a
means of diagnosing *why* those numbers are what they are.

Reporting a single decay length across both, as this document originally specified, produces a
number that looks excellent and means nothing: in V1 the metric returned "one cell past the ramp"
on seven of eight curves purely because the error was below the sampling floor everywhere, while
the underlying shapes were opposite.

---

## 7. Main concerns

Ordered by how much they worry me, not by how easy they are to write down.

**C1 — the zone is not a firewall; the pressure response is global.** The forcing is local, but
the projection that cleans up after it is elliptic. Every substep the zone injects a divergent
velocity field and the Poisson solve distributes the correction over the *whole* domain, trusted
interior included. This is intrinsic to relaxation-zone nesting — PALM and WRF have it too — but
it should be stated plainly rather than assumed away: the interior is insulated from the zone's
*momentum* forcing, not from its *pressure* signature. Mitigation (i) is now settled: the §1.3
interpolation makes the target *exactly* solenoidal, so the interpolation contributes no divergence
at all and the only source is the parent/child blend across $0<W<1$. That leaves (ii), treating $\|\mathcal{G}p\|_{\rm zone}/\|\mathcal{G}p\|_{\rm interior}$ as a headline
diagnostic, not a footnote. Per §0 this concern splits: the *fluctuating* part of the zone's
pressure signature is mixed away by canopy turbulence within the recovery fetch (§0, criterion B),
so it is bounded by measurement rather than by construction; the *mean* part is not dispersed and
must be controlled by the conservative interpolation and $\Phi=0$. Watch the ratio, but judge it on
the fetch, not on its value in the zone. V1 measured the ratio at 0.16 with buildings present —
the projection works *less* hard in the zone than in the interior.

**C2 — silent failure (F3).** The top concern, and §0 is why: a net mass error is a global
constraint violation that turbulence cannot disperse — it accumulates. This is also a property of
the *code*, not of the scheme, which is what makes it dangerous: an incompatible right-hand side produces no error, no warning, and no crash —
just a slow drift and a residual divergence buried in the top cell layer. Any bug in the flux
correction, the IBM face masking or the input file presents this way. The runtime assertion on
$\Phi$ must be **on by default**, not opt-in, and `divtot` must be watched over long runs rather
than spot-checked.

**C3 — warm restarts.** Not yet designed, and a real gap. On restart the buffer must be
repositioned to the correct parent interval, the interpolation weights recomputed for the restart
time, and the initial state must be consistent with the target in the zone — otherwise the first
substeps see a large misfit and the imposition kicks. This is exactly why the DALES prototype has
its "read the full 3-D field" hack (`modnudgeboundary.f90:274-296`). Design it properly: store the
parent interval index and the target time in the restart file, and validate on read.

**C4 — the zone costs domain.** A building-free zone of $L_{\rm imp}+L_{\rm rel}$ on each side is
$2(L_{\rm imp}+L_{\rm rel})$ of non-urban buffer in each direction. At 20 cells out of 512 that is
~8% per side and unremarkable; with a coarse parent needing a wide zone it can become a large
fraction of an expensive domain — and C1 says a wider zone also means more pressure contamination.
The zone width is therefore a genuine three-way trade (reflection vs turbulence damping vs domain
economics) that only V1/V5 can settle.

**C5 — turbulence damping and the development fetch.** Small for a matched LES parent; grows once
the parent is coarsened. The DALES paper's fetch measurements are, however, an upper bound for our
case: they are for a flat dry convective boundary layer whose only turbulence source is buoyancy,
whereas an urban child regenerates turbulence continuously in building wakes (§0). Expect a
substantially shorter fetch, and measure it (**V5**) rather than inheriting their numbers. The
scale-selective forcing of §8.4(ii) remains the real fix for a strongly coarsened parent and
belongs on the v2 roadmap.

**C6 — reflections at the imposed faces.** Fully Dirichlet lateral faces cannot let disturbances
out, so outgoing structures reflect. **We are deliberately not using a convective or radiation
outflow to fix this**: with buoyancy in play, a radiation condition lets fluid flow back in through
a face that has flipped from outflow to inflow, and degrades the solution in a way that is worse
and harder to diagnose than reflection. (DALES's own top radiation condition has to add an
explicit buoyancy damping term for precisely this reason, and their in/outflow switch is per-cell
and per-substep.) The remedy within this design is the zone itself: disturbances are damped
*before* they reach the face, and canopy turbulence disperses whatever survives (§0). If reflection
appears the answer is a wider or stronger zone, not a different boundary condition. Reflections are
a *structural* error, the kind that dilutes fastest, so this concern is milder than it looks — but
it is still the one that decides the zone width, which is what **V2** and **V5** measure.

**C7 — buildings near the zone.** §5. Warning plus wall erosion is the mechanism; a building-free
zone is the actual requirement.

**C8 — `grwdamp` versus a nested lid.** The sponge relaxes toward slab means over the top quarter
by default ([modboundary.f90:44](../src/modboundary.f90#L44)). Require `igrw_damp = 0` if Case C is
ever enabled.

**C9 — the $\rho$ inconsistency in `fillps`** (F1). Latent, not active. Separate issue against
`master`; do not fix it in this branch.

**C10 — RK order in the zone.** Accepted, and better-behaved than it first looks: because
`tstep_integrate` always integrates from $u^n$ and the third substep has $\Delta t_3=\Delta t$, the
**net relaxation over a full timestep is exactly $\exp(-W\Delta t/\tau)$, independent of the
substepping**. The first two substeps only relax the intermediate states used to evaluate the
other tendencies. First-order in the zone, RK3 everywhere else.

**Still open:** zone width and $\tau$ defaults (empirical, **V1**/**V2**); lid Case A/B/C (§3.2); whether $w$
is relaxed pointwise or only imposed in the guard strip — the DALES prototype hard-wires
$\tilde w=0$ (`modnudgeboundary.f90:410`), wrong for any convective or sloped parent; scalars and
temperature out of scope for v1 but the schema should leave room.

## 8. How PALM and DALES actually do it

### 8.1 PALM — online self-nesting (PMC)

Parent and child run **concurrently** as separate MPI process groups joined by inter-communicators
with pre-computed index maps. Per coupled step:

1. `pmci_datatrans` / `pmci_child_datatrans` — the parent sends its prognostic fields over the
   child's boundary region; the child receives into buffer arrays.
2. `pmci_interpolation` → `pmci_interp_lr`, `_sn`, `_t` — for each child face, interpolate the
   parent values onto the child's boundary *and* ghost cells **at each variable's own stagger**,
   using precomputed maps (`jfl/jfu`, `kfl/kfu`: the child index range covered by each parent
   cell). Boundary values are set, not relaxed: PALM has no relaxation zone here.
3. `pmci_boundary_conds` — re-apply PALM's *own* boundary conditions on top of the interpolated
   values, so the nested data are made consistent with the model's BC framework rather than
   overriding it.
4. `pmci_ensure_nest_mass_conservation` — integrate the boundary-normal volume flux over all five
   faces with the topography mask, `MPI_ALLREDUCE`, then add a **uniform** velocity correction
   distributed over the faces in proportion to **fluid** face area so the net flux is zero. Called
   immediately before `pres`.
5. The pressure solver runs and projects the child interior.
6. *Two-way only:* `pmci_anterpolation` averages the child solution back onto the parent grid in
   the interior, explicitly **excluding a buffer of cells near the child boundary**
   (`anterpolation_buffer_width`) because the near-boundary child solution is not trustworthy.

PALM's **offline** nesting (`nesting_offl_mod`, Kadasch et al. 2021) is the closer analogue to
what we want: read a coarse `_dynamic` driver file, buffer two time levels, interpolate linearly
in time, impose Dirichlet boundary values, apply the same mass-flux correction, and optionally add
synthetic turbulence. It has **no relaxation zone**; Kadasch et al. note one may be added.

*Transferable:* the topography-masked flux/area treatment (§3.1); interpolation at each variable's
native stagger; correction immediately before the pressure solve; and the anterpolation buffer,
whose rationale is exactly the rationale for our guard strip. *Not transferable:* the PMC
infrastructure itself, which exists to support concurrent two-way coupling we do not want.

### 8.2 DALES — open boundary conditions (`modopenboundary`)

At initialisation: read the *whole* `openboundaries.inp.XXX.nc` (all times, per-rank hyperslabs);
`openboundary_divcorr` loops over every input time, computes $\oint\rho\mathbf{u}_B\!\cdot\!\hat n\,dS$
over the five faces and distributes a uniform correction over the lateral faces to null it; then
copies the $t_0$ boundary values into `u0/um` for a cold start.

Per RK substep, at the **top** of the loop (`program.f90:226-228`):

1. `openboundary_turb` — synthetic turbulence (substep 1 only).
2. `openboundary_ghost` — halo exchange on non-open sides, then `applyboundaryf` fills ghost cells
   for cell-centred variables and *tangential* velocities: homogeneous Neumann where outflow, a
   **Robin** condition (paper Eq. 17) where inflow, blending Dirichlet→Neumann as $u_n\to0$ via a
   timescale $\tau_0[1+(u_s/u_n)^p]$ built from the subgrid velocity scale.
3. `openboundary_tend` → `applyboundaryh` — *sets* (does not add to) the tendency of the
   **boundary-normal** velocity on each face: Sommerfeld **radiation** with a Hedley–Yau phase
   speed where $u_B\!\cdot\!\hat n\ge0$, **nudging** toward the input over $\tau=\max(\tau_m,\Delta t)$
   where $u_B\!\cdot\!\hat n<0$. Then `radcorrection` per face: for each patch $S_{\rm int}$, add a
   uniform $\epsilon$ forcing the patch-integrated normal tendency to equal the patch-integrated
   parent tendency $\partial u_B/\partial t$ — so mass is conserved per *parent cell*, while the
   child stays free to generate turbulence on smaller scales.
4. Advection then runs with its **start index shifted** (`sx=3`, `modadvection.f90:53-56`) so that
   advection does not overwrite the tendency just prescribed on the boundary face.
5. `poisson` — `fillps` exchanges `up/vp/um/vm` across open boundaries and builds the
   $\rho$-weighted divergence; `tderive` sets homogeneous Neumann $p$ at open faces
   ($p_1=p_2$ etc.), so the projection leaves the prescribed boundary-normal velocity untouched.
6. `openboundary_phasevelocity` — **after** the projection, estimate the phase speed for the next
   substep from the current tendency and the velocity gradient one cell inside, patch-averaged.
7. `tstep_integrate`, then `openboundary_ghost` again.

Note the structural contrast in step 4: uDALES needs no equivalent change, because `bcpup`
overwrites the face value regardless. The uDALES idiom is that `boundary` sets the face *value*
in `u0`/`um` while `bcpup` sets it in the *predicted* field and zeroes the tendency so the RK
update cannot move it — "u(ib) only evolves according to pressure correction"
([modboundary.f90:1288](../src/modboundary.f90#L1288)). Our nesting branch must follow that idiom
exactly.

### 8.3 DALES — the relaxation prototype (`modnudgeboundary`)

Conceptually the closest to this proposal, and instructive mostly as a list of things to do
differently: Gaussian weights that are **not monotone** (peak at the offset, so the outermost
cells are forced *less* than cells further in, `:104-107`); one cell-centred weight applied to
`u`, `v`, `w`, `thl`, `qt` alike (`:415-427`), a half-cell error on every normal component;
`min(1, Σ)` corner saturation needing a special-cased radial patch; a hard-wired $\tilde w=0$
(`:410`); per-rank binary stream input keyed on `myidx/myidy` (`:276-317`); a "read the full 3-D
initial field" hack (`:274-296`); a requirement that the zone fit inside one MPI tile
(`:160-166`); explicit relaxation with $\tau_i=1/\Delta t$ as the default; and no divergence
correction of the target at all. It is applied on top of *periodic* BCs — it is a body force, not
a boundary scheme.

### 8.4 What to take, and what not to

uDALES is a DALES fork: the two still share the staggered layout, the RK3 tendency-plus-projection
structure, the `fillps`/`tderive` split, and — per Liqui Lung et al. (2024) §2.1.3, for exactly
the same reason — **homogeneous Neumann pressure with a cosine-transform solver**. Anything
written in that idiom ports.

The fragile parts of `modopenboundary` are all on the outermost face: the phase-speed estimator
with `sign(max(abs(·),1e-10))` guards, CFL clamps and patch averaging; the per-cell in/outflow
switch; and a Robin inflow needing a subgrid velocity from the SFS-TKE scheme that uDALES does not
always run (`loneeqn`). Their only test case is a flat dry convective boundary layer — none of it
has met an immersed boundary, and in our domains that face is exactly where the buildings are. It
also imports concepts uDALES lacks (`commrow`/`commcol` sub-communicators, a second patch
decomposition, `field_r`/`pois_r` kinds).

Worth noting that the relaxation-zone-plus-prescribed-BC approach proposed here is the
**mainstream** one — Moeng et al. (2007), Zhu et al. (2010), Heinze et al. (2017), as the paper's
own introduction states — and it is *simpler* than DALES's: no phase speed, no in/outflow branch,
no per-cell estimator.

| Take | Where |
|---|---|
| The input contract: self-describing, decomposition-free, hyperslab-read | `modopenboundary.f90:268-442` |
| Init-time divergence correction of the input, all times | `modopenboundary.f90:444-576` |
| Homogeneous-Neumann-$p$ + velocity correction (not an inhomogeneous Neumann $p$ BC) | paper §2.1.3 |
| Refinement limits $\le4$ spatial, $\le30$ temporal | paper, Conclusions |
| Topography-masked flux and **fluid** face area | PALM `pmc_interface_mod.f90:6207-6430`, `:2426-2520` |
| The near-boundary buffer rationale | PALM `anterpolation_buffer_width` |
| `modsynturb` — later, for a mesoscale parent | DALES `modsynturb.f90` |

| Do **not** take | Why |
|---|---|
| Radiation outflow + phase-speed estimator | Fragile in a canopy; the zone replaces its function |
| Robin inflow condition | Needs a subgrid velocity scale we may not have |
| Per-cell in/outflow switching | Constant switching under boundary recirculation |
| Patch-based `radcorrection` | Unnecessary once all faces are imposed (§3.1(2)) |
| The PMC coupler | Built for concurrent two-way coupling we do not want |
| Anything from `modnudgeboundary` verbatim | §8.3 |

**LES-to-LES vs mesoscale.** These are different modes, not one compromise. For a matched or
$\le4\times$ parent, impose the full pointwise field with a narrow zone and short $\tau$. For a
mesoscale/RANS parent there is no resolved turbulence to impose, and pointwise forcing actively
destroys the child's; then (i) widen the zone and scale $\tau$ to the parent's time resolution,
(ii) make the forcing **scale-selective**,

$$
S_q=\frac{W_q}{\tau}\Bigl[\,\mathcal{I}_{P\to C}\,q_P-\langle q_C\rangle_P\Bigr],
$$

with $\langle\cdot\rangle_P$ the box average over one parent cell, so that only the scales the
parent actually resolves are constrained and the child keeps its own sub-parent-scale
fluctuations; and (iii) inject synthetic turbulence at the guard strip. (ii) is the
highest-value v2 feature.

### 8.5 Side by side

| | DALES `openBC` | PALM (online PMC / offline) | This plan |
|---|---|---|---|
| Where the parent acts | outermost face + ghosts | outermost face + ghosts | **face (Dirichlet) + interior relaxation zone** |
| Outflow treatment | Sommerfeld radiation, Hedley–Yau phase speed | radiation (online: interpolated parent) | **none — imposed; the zone absorbs outgoing structures** |
| In/outflow switching | per cell, per substep | per boundary | **none** |
| Tangential + scalars at inflow | Robin (blends Dirichlet→Neumann) | Dirichlet | Dirichlet (velocity only in v1) |
| Relaxation zone | separate prototype on periodic BCs | none in production; noted as future work | **the core mechanism** |
| Weight function | Gaussian, cell-centred, non-monotone | n/a | raised cosine, **per-stagger**, bounded union |
| Mass compatibility | offline `divcorr` **+** per-substep patch `radcorrection` | uniform correction over all fluid faces every coupled step | **offline correction only** (linearity of $\Phi$, §3.1(2)) |
| Interpolation | pointwise | pointwise, native stagger | **conservative face fluxes**, native stagger (§1.3) |
| Time interpolation | linear, 2 levels | linear, 2 levels | linear or $C^1$ Hermite, 2–3 levels |
| Pressure BC | homogeneous Neumann + DCT | solver's own | homogeneous Neumann + DCT — **unchanged** |
| Advection at the boundary | start index shifted to `sx=3` | — | **unchanged** (`bcpup` overwrites the face) |
| Immersed boundaries | untested | `topo_flags` masks throughout | `II_u/II_v/II_w` + wall erosion + building-free zone |
| Input | NetCDF, all times read at init | NetCDF `_dynamic` driver | slabs, 2-level rolling, container-agnostic |
| Two-way | no | yes (anterpolation + buffer) | no |

In one sentence: **PALM's mass treatment, DALES's input contract, the relaxation zone neither of
them runs in production, and none of the estimator machinery.** The parts that carry risk in their
codes — phase speeds, in/outflow switches, Robin timescales, patch corrections — are all
consequences of letting the boundary be free. Imposing it and putting the compliance inside the
domain removes them all, and moves the entire risk budget onto one question: is the zone wide
enough (§7 C6)?

## 9. Implementation plan

### 9.1 Files touched

| File | Change |
|---|---|
| `src/modnesting.f90` | **new** — the whole scheme (§9.2) |
| `src/modglobal.f90` | `BCxm_nesting = 4`, `BCym_nesting = 3`; six `TEST_NESTING_*` runmode constants |
| `src/modstartup.f90` | `&NESTING` namelist + broadcasts; `checkinitvalues` guards; `call nesting_init`; `&NESTDUMP` namelist + broadcasts (D1) |
| `src/program.f90` | three call sites (§4); runmode dispatch for the new tests; `initnestdump`/`nestdump`/`exitnestdump` next to the fielddump calls (D1) |
| `src/modboundary.f90` | `case(BCxm_nesting)`/`case(BCym_nesting)` in `boundary`, in the outflow block, and in `bcpup` — each a thin delegation to `modnesting` |
| `src/tests.f90` | five new in-solver test entry points (§10.1) |
| `src/modsave.f90` | *(unchanged — the restart state is reconstructed, not stored; §9.5)* |
| `tools/python/udprep/nesting.py` | **new** — writer, conservative interpolation, divergence correction, validation |
| `tools/python/namelists.json` | `&NESTING` and `&NESTDUMP` metadata |
| `src/modnestdump.f90` | **new** (D1) -- the parent-side zone dump of §6.2: `&NESTDUMP`, per-rank band and initial-block files (spec section 9) |
| `tests/validation/nesting/caselib.py`, `make_parent_case.py`, `make_child_case.py`, `config.py` | D1: `NestDump` reader, `Preset.parent_output` / `nestdump_sections`, the `nestdump` driving source; `test_nestdump_tiny.py` closes it |
| `tests/test_suites.yml` | new `nesting` group |
| `docs/udales-boundary-conditions.md` | document `BCxm = 4`, `BCym = 3` |

`CMakeLists.txt` needs no change (`file(GLOB_RECURSE ... CONFIGURE_DEPENDS "src/*.f90")`,
[CMakeLists.txt:103](../CMakeLists.txt#L103)).

### 9.2 `modnesting` — structure

```fortran
module modnesting
  implicit none;  save;  private
  public :: nesting_init, nesting_update_target, nesting_apply, &
            nesting_boundary, nesting_bcpup, nesting_stats, nesting_finalize
  ! test hooks, public so tests.f90 can reach them without duplicating logic
  public :: nest_shape_fn, nest_union, nest_stagger_coord, nest_flux_residual

  !------------------------------------------------------------------ namelist
  logical :: lnesting        = .false.
  character(len=256) :: nestfile = ''        ! default nesting.inp.<expnr>
  real    :: nest_guardwidth = 0.            ! [m]  L_imp
  real    :: nest_zonewidth  = 0.            ! [m]  L_rel
  real    :: nest_tau        = 0.            ! [s]  0 => Dirichlet in the strip
  integer :: nest_shape      = 1             ! 1 raised cosine, 2 quintic
  logical :: nest_lateral(4) = .true.        ! W, E, S, N
  logical :: nest_top        = .false.       ! Case C
  integer :: nest_timeinterp = 2             ! 1 linear, 2 cubic Hermite (Catmull-Rom, unlimited)
  integer :: nest_nwall      = 1             ! wall erosion, cells
  logical :: nest_lparentgeom= .false.       ! parent resolves the child geometry
  real    :: nest_fluxtol    = 1.e-10        ! abort threshold on Phi
  logical :: nest_lfluxassert= .true.        ! assertion ON by default (C2)
  logical :: nest_lfluxcheckall = .false.    ! recompute Phi from the slabs at init (10.6.3)
  logical :: nest_linitfromparent = .false.  ! cold start from the full-3D block (10.6.4)

  !------------------------------------------------------------- zone geometry
  type zone_type                             ! one per staggered component
    integer              :: npts = 0
    integer, allocatable :: ijk(:,:)         ! (npts,3) local indices
    real,    allocatable :: w(:)             ! weight, IBM-masked and eroded
    integer, allocatable :: src(:)           ! offset into the slab buffer
  end type
  type(zone_type) :: zone_u, zone_v, zone_w

  !------------------------------------------------------------- parent buffer
  real,    allocatable :: bufu(:,:), bufv(:,:), bufw(:,:)   ! (npts, nslot)
  real,    allocatable :: tgtu(:), tgtv(:), tgtw(:)         ! (npts) this substep
  real,    allocatable :: facu_b(:,:,:), facu_e(:,:,:)      ! imposed face values
  real,    allocatable :: facv_b(:,:,:), facv_e(:,:,:)      !   (.,.,nslot)
  real    :: tnest(3)                                       ! buffered times
  integer :: it_lo = 0, nslot = 3
  real    :: phi_last = 0., tread_total = 0.
end module
```

Storing the zone as a **point list with weights**, rather than three more 3-D arrays, follows the
established `solid_info_type` idiom ([modibm.f90:59](../src/modibm.f90#L59)), keeps memory
proportional to the zone rather than the domain, and makes the inner loop a single flat sweep.

### 9.3 Routine contracts

| Routine | Called from | Contract |
|---|---|---|
| `nesting_init` | `program.f90`, after `createmasks`/`calcfluidvolumes`, before `readinitfiles` | Read and **validate** the file header against the run (abort on mismatch); build $W_q$ and the three zone lists; verify $\Phi=0$ for every stored time and abort if not; load the first `nslot` time levels; optionally initialise `u0/um` from the parent. |
| `nesting_update_target` | `program.f90`, after `timedep` | If `timee` has crossed into a new parent interval, roll the buffer and trigger the read for the level **one interval ahead** (§6.2). Evaluate `tgtu/tgtv/tgtw` and the imposed face values at the current time. |
| `nesting_apply` | `program.f90`, between `grwdamp` and `poisson` | The §1.2 update over the three zone lists. **Overwrites** `up/vp/wp` at zone points; nothing may be inserted between this and `poisson`. |
| `nesting_boundary` | `modboundary::boundary` | Fill the ghost planes of `u0/um`, `v0/vm`, `w0/wm` from the parent, in the `xmi_driver` pattern ([modboundary.f90:720](../src/modboundary.f90#L720)). |
| `nesting_bcpup(pup,pvp,pwp,rk3coef)` | `modboundary::bcpup` | Set `pup(ib)`, `pup(ie+1)`, `pvp(jb)`, `pvp(je+1)` to the imposed values$/\Delta t_s$; zero the matching tendencies so the RK update cannot move them; recompute $\Phi$ and assert. |
| `nesting_stats` | `program.f90`, with `statsdump` | $\Phi$, `divtot`, zone energy injection, $\|\mathcal{G}p\|$ ratio, inner-edge misfit, read time fraction. |
| `nesting_finalize` | `program.f90`, at the end of the run | Close the parent file and free the buffers. There is no restart record: `nesting_init` reconstructs the buffer state exactly from `timee` (§9.5). |
| `nest_shape_fn`, `nest_union`, `nest_stagger_coord`, `nest_flux_residual` | tests | Pure/near-pure entry points, public **so the unit tests exercise the production code rather than a copy**. |

### 9.4 When parent and child geometry differ

The parent may resolve different buildings from the child's, and V4 shows the child's interior
canopy then follows its *own* geometry (§10.5). A parent that resolves **no** buildings is
**out of scope**, decided on 2026-09-07 on the strength of V3. It is mechanically supported and
numerically clean, but V3 measures it as unusable — the child's canopy wind is wrong by a factor of
two at the first building row and by $38$–$57\,\%$ at the last, and no child equilibrated within
$26h$ of fetch, while the same child driven by a parent that *does* resolve buildings is inside the
sampling spread by $5h$ (§10.5). The underlying reason is not specific to nesting: a building-free
domain has far more difficulty generating and sustaining canopy turbulence in the first place, so
there is nothing in the imposed flow for the child to inherit. **The parent must resolve the
canopy. This is a requirement of the scheme, not a quality of implementation, and configurations
that violate it are not supported and will not be tuned for.** With that settled, three
consequences follow for the geometry that remains. It is **not** the same as saying that a mesoscale parent can drive the child:
as built this is an LES-to-LES tool. A parent whose output carries no resolved turbulence at the
scales the child needs — any parent at $C_{\rm dump}\gg2$, which a mesoscale model's 10–60 s output
always is — would need a turbulence-generating inflow the scheme does not have (§6.2), and the
child would spend its whole interior regenerating what the boundary never supplied. The earlier
wording here ("the expected case for a mesoscale parent") is withdrawn.

1. **The zone must be building-free** (§1.4), which makes the lateral boundary faces entirely
   fluid. The masked and unmasked forms of $\Phi$ then coincide — but keep the masked form, so the
   scheme stays correct if the rule is ever relaxed and so the floor is handled uniformly.

   **This is a constraint on the child alone, and that is worth stating explicitly because it is
   easy to over-read.** The *parent* may have buildings anywhere, the zone included; only the
   child's own solid mask has to be clear there, which is exactly what `nest_lparentgeom = .false.`
   asserts. The parent's buildings are still imprinted on the flow that arrives at the boundary —
   their wakes are in the imposed velocity field — so the child's zone receives physically
   meaningful forcing while containing no solid cells of its own. **But something is lost, and V2
   measured it** (§10.5): the wakes those cubes would have shed *inside* the zone are absent, so
   the first building rows of the interior receive an inflow with too little canopy-layer deficit,
   and the interior mean wind runs $0.05$–$0.07\,u_\star$ low up to $z/h\approx1.5$ — enough to fail
   criterion A for a child that clears 12–20 cubes, where a child whose zone was building-free in
   the parent too passes. Clearing the child is correct at the boundary; the interior then starts
   its own canopy adjustment at the inner edge, which is point 2 below seen from the other side.

   The practical consequence is that a building-free zone never requires modifying the parent.
   V1 carved a plaza out of the parent geometry to achieve it, which worked but was unnecessary,
   and which constrained the child sizes and zone widths a later sweep could reach — a parent
   modified for one child fits only that child. Clear the child instead.

   **A matched-geometry control is not always constructible, and that is a geometric fact rather
   than an oversight.** Clearing the child requires dropping every cube within the zone depth of a
   face, so a dropped cube may not occupy more than the building period allows: at V1's spacing
   that is 26 m of clearance plus an 8 m half-width at each end, against a 32 m period. A regular
   array satisfies it; a **staggered** array does not, at any child size, because its two column
   families sit exactly half a period apart, so one family always lands in a blocked window and
   clearing it would remove buildings from the analysis interior rather than only from the zone.
   `Preset.validate()` refuses such a configuration rather than silently producing one. The
   consequence for V4 is that its staggered-parent case has no matched-geometry control of its own
   and must be judged against the V1 child.
2. **The imposed near-surface profile will not be in equilibrium with the child's surface.** The
   parent's wind near the ground reflects its own roughness and (absent) canopy; the child's zone
   has only ground roughness, and the interior has buildings. An internal boundary layer must
   therefore develop between the inner edge of the zone and the first building row. This is
   physics, not a bug: it is the same adjustment fetch that any nested LES pays. It must be
   *measured* (validation **V3**) and budgeted for in domain layout — but note that
   **the buildings should start immediately at the inner edge of the zone, not after a standoff.**
   A building-free standoff is actively counterproductive: the flow there adjusts only to the
   ground roughness through weak shear-driven mixing, and then has to adjust a *second* time on
   reaching the canopy. Starting the buildings at $W=0$ gives one adjustment instead of two, driven
   by vigorous wake turbulence. The first few building rows will carry visibly wrong facet
   stresses; that is accepted (§0).
   Do **not** try to fix it by tapering $W$ with height: that would make the imposed lateral flux
   height-dependent in a way not reflected in the offline correction, breaking $\Phi=0$.

   **The standoff argument above was refuted by V3 (§10.5), for a configuration now out of
   scope.** With a parent resolving no buildings, every standoff of 5, 15 or 40 cells beat a zero
   standoff by more than the sampling spread, monotonically at $20h$. The "two adjustments are
   worse than one" reasoning missed that the first adjustment sheds some of the excess near-surface
   momentum such a parent imposes. Since that parent is no longer supported, the *advice* stands
   unchanged — start the buildings at the zone edge — but it now rests on the cleared-parent-cubes
   result, where a child at zero standoff off a proper parent is inside the sampling spread by
   $5h$, rather than on the argument given above, which is withdrawn.
3. **`nest_lparentgeom` records which case applies.** When `.false.` (parent has no matching
   geometry), `nesting_init` *requires* a building-free zone and errors otherwise; when `.true.`
   (self-nesting on identical geometry) it downgrades to a warning.

### 9.5 Restart (C3)

The zone state is entirely reconstructible from the input file plus the time — unlike the DALES
prototype's hack (`modnudgeboundary.f90:274-296`), which re-reads full fields. **In the
implementation this turns out to be enough on its own, so the restart record is not written.**
`nesting_init` is called *after* `readinitfiles`, by which point `timee` holds the restart time
(cold start: `modstartup.f90:1203`; warm start: `readrestartfiles`); `set_interval(timee)` is
deterministic in `timee`, so the reconstructed `it_lo` and target are identical to the ones a
continuous run would have held. Test **I6** verifies this: 100 steps versus 50 + restart + 50 is
bitwise identical, for both a mid-interval restart and one exactly on a parent level.

There is no nesting restart record. An earlier revision carried `nesting_restart_write`/
`nesting_restart_read` (and a unit test, U20) against the day the buffer state might be pinned in
the `initd` file rather than reconstructed; they had no production call site and were removed in
the 2026-09 review pass, since reconstruction from `timee` is exact and I6 pins it. `nesting_finalize`
is called from `program.f90` at the end of the run and closes the parent file.

**The call order is load-bearing.** `nesting_init` reads `timee`, and `readinitfiles` is what
assigns it — `real :: timee` at `modglobal.f90:441` has no initialiser. Calling `nesting_init`
first, as the first implementation did, leaves `set_interval`/`eval_target` operating on an
undefined value: on a Debug build (`-init=snan -fpe0`, `CMakeLists.txt:62-66`) every nested case
trapped immediately in `eval_target`, and on a Release build a warm start silently read `timee` as
0 from static storage and positioned the buffer at $t=0$, breaking restart parity at
$\sim10^{-7}$.

### 9.6 Staging

| Stage | Content | Gate |
|---|---|---|
| **M0** | `modnesting` skeleton: namelist, `nest_shape_fn`, `nest_union`, `nest_stagger_coord`, zone lists, diagnostics. `nesting_apply` a no-op. BC constants, `checkinitvalues` guards. Runmodes `TEST_NESTING_WEIGHTS`, `TEST_NESTING_GEOMETRY`. | **U1–U14**, **I1** |
| **M1** | Python writer: conservative interpolation, divergence correction, schema, validation. | **P1–P11** |
| **M2** | Reader, time buffer, prefetch, restart state. Runmode `TEST_NESTING_IO`. | **U15–U22** |
| **M3** | `nesting_bcpup`, `nesting_boundary`, $\Phi$ assertion. Runmode `TEST_NESTING_FLUX`. | **U23–U28**, **I2** |
| **M4** | `nesting_apply` (§1.2). Runmode `TEST_NESTING_UPDATE`. | **U29–U34**, **I3–I4** |
| **M5** | Full path; decomposition, restart and IBM parity. | **I5–I8** |
| **M6** | Validation campaign. | **V1–V7** |

**Status: M0–M5 complete, and the §10.6 follow-ups with them.** 43 unit
assertions across six runmodes on four decompositions, 59 Python writer tests,
and the I1–I10 matrix all pass against Release; the matrix also passes against
Debug (`-init=snan -fpe0`), with the two I1-vs-baseline comparisons skipped
since they need a Release binary. M5 found one real defect, the `nesting_init`
call order (§9.5), which is fixed. M6 has not been started.

---

## 10. Test suite

The design principle is that **every mechanism is testable in isolation, by a test that fails only
if that mechanism is wrong.** The integration tests then check composition, not correctness of the
parts. Unit tests call the *production* routines through the public test hooks of §9.3 — never a
reimplementation.

### 10.1 Unit — in-solver runmodes

Following the existing `TEST_SPARSE_IJK` pattern ([tests.f90](../src/tests.f90),
[program.f90:238](../src/program.f90#L238)): a `logical function` per test, dispatched by
`runmode`, exiting 0/1. Cheap, MPI-aware, CI-able. `tests/unit/` currently holds only a README —
these are its first real occupants.

**`TEST_NESTING_WEIGHTS = 1006`** — pure algebra, no grid, no I/O.

| ID | Isolates | Pass criterion |
|---|---|---|
| U1 | shape function value | $W(0)=1$; $W=1$ for $s\le L_{\rm imp}$; $W=0$ for $s\ge L_{\rm imp}+L_{\rm rel}$ |
| U2 | range and monotonicity | $0\le W\le1$ and $W$ non-increasing, over 10⁴ samples |
| U3 | $C^1$ continuity | numerical $dW/ds$ continuous, and **zero**, at both $s=L_{\rm imp}$ and $s=L_{\rm imp}+L_{\rm rel}$, to $O(\Delta s^2)$ |
| U4 | quintic variant | as U1–U3, plus $d^2W/ds^2$ continuous |
| U5 | corner union bounds | $W_\cup\le1$; $W_\cup\ge\max_f W_f$; $W_\cup=1$ iff some $W_f=1$; symmetric under face permutation |
| U6 | union degeneracy | one active face $\Rightarrow W_\cup=W_f$ to within one ulp (the union evaluates $1-(1-a)$, which is not bitwise $a$; measured $5.6\times10^{-17}$) |
| U7 | integral identity | $\int W\,ds = L_{\rm imp}+\tfrac12 L_{\rm rel}$ to quadrature error — the quantity §1.4(b) depends on |

**`TEST_NESTING_GEOMETRY = 1007`** — stagger, indexing, masking. Run on 1×1, 2×1, 1×2, 2×2.

| ID | Isolates | Pass criterion |
|---|---|---|
| U8 | stagger coordinates | for every local point, `nest_stagger_coord` returns `xh(ig)`/`yf(jg)`/`zf(k)` for $u$, `xf`/`yh`/`zf` for $v$, `xf`/`yf`/`zh` for $w$, matching an independent global formula to round-off. **This is the half-cell-error test** — the defect present in the DALES prototype. |
| U9 | zone membership | the set of global points with $W>0$ is identical across all four decompositions |
| U10 | weight invariance | $\sum W$ (MPI-reduced) identical across decompositions to round-off |
| U11 | IBM masking | no solid point has $W>0$, for a synthetic solid list covering faces, edges and corners |
| U12 | wall erosion | no point within `nest_nwall` of a solid has $W>0$; count matches an independently computed dilation |
| U13 | building-free enforcement | with `nest_lparentgeom=.false.` and a solid inside the zone, init **errors**; with `.true.`, warns and continues |
| U14 | width reporting | the cell-equivalent width reported from a width in metres is correct, and the <6-cell / >15%-of-domain warnings fire |

**`TEST_NESTING_IO = 1008`** — reader and time buffer. The file is generated by the Python writer
from an analytic field, so this test also pins the writer/reader contract.

| ID | Isolates | Pass criterion |
|---|---|---|
| U15 | spatial read | every buffered value equals $f(x_q,y_q,z_q)$ for a non-separable analytic $f$ (so a transposed or off-by-one index cannot pass), to round-off |
| U16 | decomposition independence | buffers gathered from 1×1, 2×1, 1×2, 2×2 are bitwise identical |
| U17 | time interpolation exactness | for a field **linear in $t$**, the target at 20 arbitrary times is exact to round-off (linear and Hermite modes both) |
| U18 | Hermite $C^1$ | for a smooth-in-$t$ field, $\partial\tilde q/\partial t$ is continuous across an interval crossing; linear mode shows the expected jump — so the test documents the difference rather than hiding it |
| U19 | buffer roll | stepping across several interval boundaries gives a target identical to a run that read every level eagerly |
| U20 | *(deleted)* | tested the restart record's repositioning; removed with `nesting_restart_write/read`, since §9.5 reconstructs the buffer from `timee` and there is no state to save (I6 pins the restart bitwise) |
| U21 | prefetch does not change results | prefetch on/off give bitwise-identical targets |
| U22 | header validation | a file with mismatched `itot`, `xlen`, `zf` or stagger tag **aborts with a specific message** — one subtest per field |

**`TEST_NESTING_FLUX = 1009`** — the compatibility machinery.

| ID | Isolates | Pass criterion |
|---|---|---|
| U23 | $\Phi$ correctness | for an imposed face field with analytically known net flux, `nest_flux_residual` matches to round-off |
| U24 | $\Phi$ decomposition invariance | identical across 1×1, 2×1, 1×2, 2×2 |
| U25 | $\Phi$ with a masked face | with synthetic solids on a boundary face, $\Phi$ uses fluid faces only (compare against a hand-summed value) |
| U26 | corrected input | $\Phi=0$ to round-off for every time level of a writer-corrected file |
| U27 | assertion fires | a deliberately uncorrected file trips `nest_lfluxassert` and aborts — **the guard against C2 must itself be tested** |
| U28 | linearity of $\Phi$ in time | $\Phi$ at an interpolated time equals the interpolation of the endpoint values, confirming §3.1(2) in the code and not just on paper |

**`TEST_NESTING_UPDATE = 1010`** — the relaxation update, on synthetic single points.

| ID | Isolates | Pass criterion |
|---|---|---|
| U29 | no-op off-zone | $W=0\Rightarrow$ `up` **bitwise** unchanged |
| U30 | Dirichlet limit | $\tau\to0$ or $W\to\infty\Rightarrow q^{\rm new}=\tilde q$ exactly |
| U31 | linear limit | $W\Delta t_s/\tau\le10^{-4}\Rightarrow$ agrees with $q_p{+}{=}(W/\tau)(\tilde q-q^\ast)$ to $O(\epsilon^2)$ |
| U32 | full-step composition | driving one point through all three substeps with $\mathcal{N}=0$ gives exactly $\tilde q+(q^n-\tilde q)e^{-W\Delta t/\tau}$ — the §7 C10 property |
| U33 | stability and monotonicity | for $\Delta t/\tau\in[10^{-3},10^{6}]$, no overshoot, no sign change, no NaN |
| U34 | solid points untouched | `up` at solid points is unchanged by `nesting_apply` |

**`TEST_NESTING_FLUX = 1009`, continued** — §10.6 items 3 and 5.

| ID | Isolates | Pass criterion |
|---|---|---|
| U35 | the lid split | `nest_flux_split` gives $\Phi$, $\Phi_{\rm lid}$ and $\Phi-\Phi_{\rm lid}$ each equal to an independently summed global reference over the six / top / five faces |
| U36 | a closed lid | with $w^\ast=0$ at the top, $\Phi_{\rm lid}$ is **exactly** zero and the asserted quantity is the old six-face $\Phi$ — so case A is unchanged |
| U37 | the stored residual is real | `flux_residual` equals the residual recomputed from the four boundary-normal slabs, read straight through `modnestingio`, for every stored level; and `fluid_lateral_area` is this grid's |
| U38 | schema 1 still works | a schema 1 file loads, warns that it predates schema 2, recomputes, and gives $\Phi=0$ at run time |
| U39 | what the cheap check costs | a file whose stored residual lies is accepted at init and caught by the per-substep assertion instead — the trade-off is pinned rather than implied. Two abort cases close it: `nest_lfluxcheckall` catches the same file, and an `fluid_lateral_area` mismatch forces the recompute on its own |

**`TEST_NESTING_INIT = 1011`** — cold-start initialisation from the parent (§10.6 item 4). Run on
1×1, 2×1, 1×2, 2×2.

| ID | Isolates | Pass criterion |
|---|---|---|
| U40 | the block is read correctly | every point of `u0/um`, `v0/vm`, `w0/wm` that the reader must fill equals the stored block — a non-separable analytic 3-D field — to round-off, including the far faces `u(ie+1)`, `v(je+1)` and `w(ke+1)` on the ranks that own them. **This is the $t=0$ statement**, and it is made here rather than in an integration test because a restart file can only be written after a step |
| U41 | the switch is a switch | with `nest_linitfromparent = .false.` and the same file, the fields are **bitwise** unchanged |
| U42 | warm starts are not touched | with the switch on and `lwarmstart`, the fields are bitwise unchanged and the run says so |
| U43 | the failure modes abort | switch on with no block; a block at the wrong shape; a block at the wrong stagger — one subtest per case, each its own process |

### 10.2 Unit — Python (`tools/python/tests/`)

| ID | Isolates | Pass criterion |
|---|---|---|
| P1 | conservative interpolation, constants | a uniform parent field maps to the identical uniform child field |
| P2 | conservative interpolation, flux identity | $\sum$ child face fluxes over each parent face $=$ the parent face flux, to round-off — **the defining property** |
| P3 | divergence preservation | a discretely solenoidal parent field gives a child target whose divergence, integrated over each parent cell, is zero to round-off |
| P4 | exact divergence preservation | child-cell divergence equals that of the containing parent cell to round-off; a solenoidal parent gives a solenoidal child target (§1.3) |
| P5 | divergence correction | $\Phi=0$ after correction, for random parent data, on every time level |
| P6 | correction is minimal and shape-preserving | it adds a constant per unit fluid area — differences between faces are unchanged |
| P7 | idempotence | correcting an already-corrected file changes nothing to round-off |
| P8 | masked correction | with a partially solid boundary, correction uses fluid area only and still gives $\Phi=0$ |
| P9 | schema round-trip | write→read reproduces every field and attribute; a missing required attribute is rejected |
| P10 | refinement guard | spatial ratio >4 or temporal >30 is refused without an explicit override |
| P11 | container equivalence | the raw-binary and NetCDF back-ends produce byte-identical buffers on read |
| P12 | the stored residual | `flux_residual` is the residual *after* correction, not before; an uncorrected file stores its real one; `fluid_lateral_area` is written and is the geometric (not $\rho$-weighted) masked area |
| P13 | the initial-condition block | round-trips bitwise on both back-ends; carries the contract's dimensions and stagger; a wrong stagger, an undeclared block, a declared-but-missing block and a schema 1 file asked to carry one are each rejected |
| P14 | the projection works | a random closed-box field with $O(1)$ divergence comes back at $<10^{-12}$, on a uniform and on a stretched vertical |
| P15 | the projection is minimal | every boundary-normal velocity is bitwise unchanged, while the interior demonstrably moves |
| P16 | idempotence | projecting an already-solenoidal field changes nothing to round-off; an incompatible field is **refused**, not absorbed |
| P17 | schema 1 compatibility | a schema 1 file writes, validates and reads back with the same slab bits, carries none of the schema 2 items, and an unknown schema is still rejected |

### 10.3 Integration (`tests/integration/nesting/`)

| ID | Isolates | Method | Pass criterion |
|---|---|---|---|
| I1 | **no-op guarantee** | existing case, `lnesting=.false.`; the pre-branch binary is built by the driver from `origin/master` with the compiler and build type of the build under test (`tests/integration/nesting/_baseline.py`), so the small case runs in CI | bitwise identical to the pre-branch binary on a small periodic case. Larger cases cannot be held to bitwise: every FFT is planned with `FFTW_MEASURE` (`modpois.f90:110-191`, `2decomp-fft/src/fft_fftw3.f90:26`), which selects the algorithm by run-time timing, so *the same binary* differs from itself at $\sim5\times10^{-12}$ relative. There, judge against the measured baseline-vs-baseline self-noise. |
| I2 | **face imposition survives projection** (F2) | impose a known $u$ on the faces, one substep | $u$ at the boundary faces after `poisson`+`tstep_integrate` equals the imposed value to round-off |
| I3 | **uniform flow** | parent $=(U,0,0)$, full path | field preserved exactly; $\Phi$, `divtot`, `divmax` and $p$ all at round-off |
| I4 | **manufactured solenoidal field** | two fields, $W\equiv1$. (a) $u=\sin\frac{2\pi x}{L}\cos\frac{2\pi y}{L}$, $v=-\cos\frac{2\pi x}{L}\sin\frac{2\pi y}{L}$, $w=0$; (b) $\psi=\sin\frac{2\pi x}{L}\sin\frac{4\pi y}{L}$, $u=\partial_y\psi$, $v=-\partial_x\psi$ | (a) is *exactly* discretely solenoidal, not $O(h^2)$ — see the note below — so it must sit at round-off on every grid and cannot measure a rate. (b) has genuinely $O(h^2)$ discrete divergence and carries the convergence check: $\|\mathcal{G}p\|$ at second order over three grids at fixed $\Delta t$ |
| I5 | **decomposition parity** | I3 and I4 on 1×1, 2×1, 1×2, 2×2; and, in CI, `tests/cases/064` (a cube) with a 3 m guard, a 20 m ramp and $\tau=4$ s on 1×1 and 2×2 (`TestI5CubeParity2x2`) | fields agree to $10^{-9}$, mirroring `tests/integration/processor_boundaries/`; measured exactly 0 on both |
| I6 | **restart parity** | 100 steps vs 50 + restart + 50 on one rank; 8 vs 3 + 5 and 5 + 3 on 2×2 ranks with the cube of case 064 (`TestI6CubeRestartParity2x2`) | bitwise identical, including mid-parent-interval restarts and restarts exactly on an interval boundary; on 2×2 all 40 per-rank restart records |
| I7 | **zone isolation** | two runs differing only in the interior, identical in the zone | the difference stays confined as expected; quantifies C1's global pressure response rather than assuming it away |
| I8 | **IBM interaction** | buildings adjacent to the zone edge | facet stresses on the building match a no-nesting reference to a stated tolerance: `tau_y`, `tau_z`, `pres` on the windward face, `tau_x`, `tau_z`, `pres` on the side faces (`tau_x` is identically zero on a face whose normal is $x$) |
| I9 | **cold start from the parent** | ZONED case, `prof.inp` carrying `u = 0` against a parent carrying `u = U`, run with and without `nest_linitfromparent` | the run is **bitwise identical** to one whose `prof.inp` carries the same field — two runs of this solver cannot agree bit for bit unless they started from the same bits; the first `divmax` is at round-off; the first substep's $\|\mathcal{G}p\|$ in the interior is $4\times10^{-6}$ of the `prof.inp` control's |
| I10 | **leaky lid** (case B) | `BCtopm_pressure`, warm-started from a restart whose `pres0` carries a uniform offset, with `nest_lfluxassert = .true.` | the run completes; $\Phi_{\rm lid}$ is eight orders of magnitude above `nest_fluxtol`, so the pre-fix assertion **would** have fired; $\Phi$ over the closed faces stays at round-off; `divmax` stays at round-off |

**Why the I4(a) field is exact, not second order.** On the staggered grid $u$
sits at $(x_h,y_f)$ and $v$ at $(x_f,y_h)$, and $x_h(i)+h/2=x_f(i)$. With
$k=2\pi/L$,

$$\frac{u_{i+1,j}-u_{i,j}}{h}=\frac{2\sin\frac{kh}{2}}{h}\cos k x_f\,\cos k y_f,
\qquad
\frac{v_{i,j+1}-v_{i,j}}{h}=-\frac{2\sin\frac{kh}{2}}{h}\cos k x_f\,\cos k y_f,$$

which cancel identically for **any** $h$: the two terms carry the same
$\sin(kh/2)$ factor because $u$ and $v$ use the same wavenumber. Measured
$\max|\mathcal{D}\mathbf{u}|$ is $4.2\times10^{-16}$, $8.3\times10^{-16}$,
$1.8\times10^{-15}$ at $h=2,1,0.5$ m — round-off, and *growing* with
resolution. Asserting second order on this field would be asserting a property
of the round-off. Giving $u$ and $v$ different wavenumbers, as in (b), breaks
the cancellation and leaves the expected
$\mathcal{D}\mathbf{u}=\frac{h^2}{24}k_xk_y(k_x^2-k_y^2)\cos k_xx_f\cos k_yy_f+O(h^4)$.

### 10.4 System / validation (`heavy`)

| ID | Question | Method | Deliverable |
|---|---|---|---|
| V1 | Does matched LES-to-LES nesting reproduce the parent? | Big Brother: periodic parent, writer dumps zone slabs, sub-domain child at matched resolution | **DONE — §10.5.** Mean flow yes ($0.008\,u_\star$); canopy turbulence yes (1–2 %); above the canopy a real $\approx10\,\%$ resolved-TKE deficit from insufficient fetch |
| V2 | Does the zone width behave as §1.4 predicts? | V1 repeated over $N_{\rm rel}\in\{4,9,12,16\}$ at fixed child size, plus a child-size arm (interior $5h$, $9h$, $13h$) at fixed zone width, every child clearing its own zone (§9.4) | **DONE — §10.5.** Zone width moves the deficit by 0.16 % over the whole range (0.03 of a sampling spread); interior extent moves it from $-13.5$ to $-9.9\,\%$. The deficit is a recovery over fetch, and §10.5's note on *what* is being recovered from applies |
| V3 | **Parent without buildings** (**closed — out of scope**) | parent resolves no geometry; child has buildings starting **at** the inner zone edge, compared against 0/5/15/40-cell standoffs; plus a **cleared-parent-cubes** arm at standoff 0 — the identical child, but its parent resolves V1's own aligned canopy where the child's zone sits and the child clears it (`nest_lparentgeom = .false.`) | **DONE — §10.5.** The adjustment length is **not measurable**: with a building-free parent no child equilibrated within $26h$, and the canopy wind is $38$–$57\,\%$ wrong at the last row. A parent that resolves buildings puts the same child inside the sampling spread by $5h$. §9.4's prediction that a standoff lengthens adjustment is **REFUTED**, for a configuration now out of scope. **Decision 2026-09-07: building-free parents are not supported and this row is closed** — a building-free domain struggles to generate and sustain canopy turbulence at all, so there is nothing for the child to inherit. The cleared-parent-cubes arm isolates what a genuinely absent boundary condition costs from what "parent had cubes, child removed them" costs (§0, "New from V2") |
| V4 | **Different parent geometry** | parent with a different building layout, child identical to V1's | **DONE — §10.5.** Turbulence: the child's canopy is its own ($+46\,\%$ against its parent, 3–4$\sigma$) and aloft it matches its parent to $1.1\,\%$. Mean flow: **inherited, not re-established** — the child sits at $0.80$ of the V1 child, tracking the $0.83$ bulk ratio of the two parents. The V1-child reference is valid for turbulence only |
| **V0** | **Does a child at higher resolution than its parent reproduce it?** | genuinely coarse parent grid ($r = 2, 4$), child at $\Delta x$; the writer's conservative interpolation carries the refinement | **DONE, both arms — §10.5**, and restated after the TKE-estimator correction. It reproduces its *parent*: with genuine 4 m / 8 m parents the child carries the parent's own $0.9$–$1.5\,u_\star$ mean-flow bias through the whole interior. Above the canopy it recovers at most a few points on its parent within $13h$; in the canopy it *overshoots* the truth by $+18$–$28\,\%$, because the too-fast imposed momentum drives too much production. Mean flow is not readable at all until the prolongation question (R2) is settled and the case re-run |
| V5 | How far can the parent be coarsened? | parent smoothed at 2/4/8 in space, 10/30/60 in time | **Superseded**: the time axis is C0 (0.5–9 s, §10.5) and the space axis is V0; the ≤4/≤30 guidance is withdrawn (§1.3). C6 is a go, on the terms of §6.2 |
| V6 | Does mass drift over long runs? | 10⁵-step run | `divtot` bounded, not drifting |
| V7 | Does the I/O cost anything? | production-sized case | read time <1% of runtime; if not, switch container (§6.3) |

**Refinement ratio 1 is a deliberate simplification, and a temporary one.** V1 and V2 both run
parent and child on the *same* grid. That is the right way to start — at $r=1$ the prolongation is
the identity, so anything the experiment finds is attributable to the scheme rather than to the
interpolation, and the interpolation has its own unit coverage in P1–P17. But the point of nesting
is to run the child at *higher* resolution than the parent, and **no end-to-end result exists yet
for that**: the V0 harness is built and validated at tiny scale, but its production run has not
completed. V1 and V2 must not be read as validating refinement.

V0 above is that test, and it is listed first because it is the most important one outstanding, not
because it is next in sequence. Two of its results are not predictable from what has been measured
so far: the deficit V1 found is in the 8–64 m band, which at $r>1$ straddles the parent's filter
scale, and the child would there be asked to *generate* structure the parent never resolved rather
than to reproduce structure the parent had. That could go either way — the child has more capacity
to build its own inertial range, or it has a larger gap to bridge.


### 10.5 V1 results

Two runs of the Big Brother experiment, `production` (job 3990472) and `converged` (job 3991175),
identical in every physical parameter and differing only in schedule. Full write-up and plots in
the published report; the numbers that change this document are below.

Parent $256\times256\times64$ at $\Delta x = 2$ m with 228 cubes of $h=16$ m; child
$128\times128\times64$ centred, leaving $104\times104$ cells = $13h\times13h$ of interior outside
the zone. Zone 3 imposed + 9 relaxed cells per side, $\tau = 1$ s, linear time interpolation, child
cold-started from the parent's 3-D block. 64 cores; 53 min and 4 h 08 respectively.

| | `production` | `converged` |
|---|---|---|
| statistics window | 1491 s (497 samples) | 10 191 s (3397 samples) |
| $\|\Phi\|/A$ | $2.85\times10^{-14}$ | $2.98\times10^{-14}$ |
| $\langle u\rangle$ RMS difference | $0.0129\,u_\star$ | $0.0083\,u_\star$ |
| $\langle u\rangle$ sampling floor | $0.165\,u_\star$ | $0.120\,u_\star$ |
| TKE RMS difference | $0.114\,u_\star^2$ | $0.121\,u_\star^2$ |
| TKE sampling floor | $0.158\,u_\star^2$ | $0.036\,u_\star^2$ |

**The deficit is real.** The sampling floor fell by $4.4\times$ between runs and the difference did
not move; a sampling artefact would have fallen with the floor. Above $z/h=2$ the deficit is
$-9.9\,\%$ in the mean, reaching $4$–$8\sigma$ at $z/h = 3$–$5$ against per-height spreads of
1.1–1.5 %. Below $z/h=1$ it is $+1.5\,\%$ against a 0.7 % spread — the canopy layer is clean, as
§1.4 assumed when the zone width was chosen.

**The deficit recovers over fetch — but see the correction below on where it comes from.** The deficit is scale-selective: at
$z/h=2.06$ the child/parent spectral ratio is $1.03$ for $\lambda > L/4$, $0.875$ for $16$–$64$ m,
$0.827$ for $8$–$16$ m, and $1.05$ below $4\Delta x$. The imposed large scales come through
*strong*; what is missing is the energy-containing middle, which the child must regenerate itself.
Those two intermediate bands reproduce to four parts in a thousand across runs sharing no samples,
which is what makes the result deterministic rather than statistical. In physical space the same
thing appears as a TKE error decaying monotonically from $0.37$ at the zone edge to $0.087$ at the
far side of the interior — *still falling* where the far zone truncates it.

**One correction to an earlier reading.** The `production` TKE trace appeared to climb steadily and
was read as an unfinished spin-up. Over the longer window the domain is seen to carry a slow
oscillation of a few thousand seconds with $\pm10\,\%$ excursions and no net trend: the short run
had landed on a rising limb. The lesson is about window length, not spin-up — a window must span
several of these cycles. The child follows the parent's slow modes with a lag of several hundred
seconds.

**C0 — the cadence discriminator (jobs 3993705 and 3993706, 2026-09-07).** Where the deficit comes
from. The paragraph above names the recovery and misses the loss: the boundary data is sampled
every 3 s and interpolated in time, which removes every wavelength below $2U\Delta t_P$ (§1.3).
Two arms. **C0a** re-ran the V1 child from the same 3 s dumps subsampled to 6 s and 9 s, and once
more at 3 s with the Catmull–Rom interpolant, over the full 10 191 s window. **C0b** warm-started
the converged parent from its end-of-spin-up restart ($t=10\,800$ s) for 2 400 s dumping every
0.5 s (4 800 levels, 240 GB), and sliced those dumps to six cadences, so one parent realisation
drives the whole ladder and the comparison is paired. Pre-registered predictions are in
`tests/validation/nesting/README.md`. Deficit is resolved TKE above $z/h=2$; band ratios at
$z/h=2.06$.

| cadence | interpolant | $C_{\rm dump}$ at $u_0$ | deficit | 8–16 m | 16–64 m | samples |
|---|---|---|---|---|---|---|
| 0.5 s | linear | 0.75 | $-2.1\,\%$ | 0.935 | 0.971 | 600 |
| 1 s | linear | 1.5 | $-3.3\,\%$ | 0.920 | 0.960 | 599 |
| 1.5 s | linear | 2.25 | $-4.9\,\%$ | 0.914 | 0.950 | 599 |
| 3 s | linear | 4.5 | $-11.3\,\%$ | 0.831 | 0.873 | 597 |
| 3 s (V1) | linear | 4.5 | $-9.9\,\%$ | 0.833 | 0.869 | 3 397 |
| **3 s (C0a)** | **Catmull–Rom** | 4.5 | $\mathbf{-5.7\,\%}$ | **0.892** | **0.930** | 3 397 |
| 6 s | linear | 9 | $-20.5\,\%$ | 0.713 | 0.725 | 594 |
| 6 s (C0a) | linear | 9 | $-18.1\,\%$ | 0.713 | 0.726 | 3 394 |
| 9 s | linear | 13.5 | $-24.5\,\%$ | 0.663 | 0.650 | 592 |
| 9 s (C0a) | linear | 13.5 | $-21.5\,\%$ | 0.664 | 0.654 | 3 391 |

**The cause is the cadence.** The band ratios fall monotonically with $\Delta t_P$ at every height,
and the 3 s, 6 s and 9 s points agree between the two arms — different parent realisations, windows
differing by $5.7\times$ — to 0.002 in the 8–16 m band. Against the pre-registered table the
"cadence causes it" row is met on 16–64 m ($0.971\ge0.97$) and just short on 8–16 m (0.935), where
the linear interpolant's own passband still attenuates an 8 m eddy at 3.6 m/s by
$\mathrm{sinc}^4(0.225)\approx0.84$ even at 0.5 s: the residual is the interpolant, not the scheme,
and no relaxation-time arm is needed. The criterion-A flags on the 600-sample rows are the
mean-flow sampling floor at 1 800 s (spread 8.2 %), not a result; the same child passes at
3 397 samples.

**One prediction was wrong, and usefully.** Catmull–Rom at 3 s was predicted to move 16–64 m a
little and 8–16 m not at all, since the latter lies wholly above the target's Nyquist frequency.
It moved 8–16 m from 0.833 to 0.892 and halved the deficit. The cubic's flatter passband supplies
more of 16–64 m, and the cascade rebuilds 8–16 m from it. `nest_timeinterp = 2` is therefore the
default (§1.3), and its own cadence ladder (C0c, jobs listed in the plan) gives the operating curve
for the interpolant the design recommends.

**The deficit curve is the deliverable, and C0c gives it for the recommended interpolant** (job
3994598, the same six cadences off the same 960 dumps with `nest_timeinterp = 2`):

| cadence | 0.5 s | 1 s | 1.5 s | 3 s | 6 s | 9 s |
|---|---|---|---|---|---|---|
| deficit, linear | $-2.1$ | $-3.3$ | $-4.9$ | $-11.3$ | $-20.5$ | $-24.5\,\%$ |
| deficit, Catmull–Rom | $-1.6$ | $-1.8$ | $-2.4$ | $-6.4$ | $-13.3$ | $-15.6\,\%$ |
| 8–16 m, Catmull–Rom | 0.944 | 0.937 | 0.936 | 0.897 | 0.814 | 0.768 |
| 16–64 m, Catmull–Rom | 0.974 | 0.971 | 0.975 | 0.940 | 0.825 | 0.750 |

The cubic halves the cost at every cadence from 3 s up, and at 1.5 s it matches what the linear
interpolant needs 0.5 s for — three times the storage. Below 1.5 s both flatten at $-2\,\%$, the
floor of this child's own fetch. One caveat the ladder exposes: at $C_{\rm dump}\ge9$ the cubic
*overshoots*, putting $8$–$14\,\%$ excess energy into the scales above $L/4$ and $+3$–$6\,\%$ into
canopy TKE, so a coarsely sampled parent should not be read as merely "smoothed" under the
cubic. The operating rule that follows: **Catmull–Rom, $C_{\rm dump}\le2$ for a deficit under
3 %, $C_{\rm dump}\le4.5$ for one under 7 %.** With the V2 size arm this gives the first points of
$L_{\rm rec}(\lambda)$: at $C_{\rm dump}=5.4$ (linear) the 8–16 m band recovers from 0.62 at $5h$
to 0.83 at $13h$ and is still rising.

**V2 results (job 3992816, 3 h 47).** Six children of the converged parent, five run and the V1
child reused, all clearing their own zone (`nest_lparentgeom = .false.`), 3397 samples each.
Deficit is resolved TKE above $z/h=2$; spectral ratios at $z/h=2.06$.

| arm | point | interior | $N_{\rm rel}$ | deficit | 16–64 m | 8–16 m | $>L/4$ | criterion A |
|---|---|---|---|---|---|---|---|---|
| zone | `nrel4` | $14.25h$ | 4 | $-9.85\,\%$ | 0.880 | 0.841 | 1.025 | $0.055\,u_\star$ |
| both | `ref` (V1) | $13h$ | 9 | $-9.91\,\%$ | 0.869 | 0.833 | 1.033 | $0.039\,u_\star$ |
| zone | `nrel12` | $12.25h$ | 12 | $-9.95\,\%$ | 0.863 | 0.826 | 1.016 | $0.033\,u_\star$ |
| zone | `nrel16` | $11.25h$ | 16 | $-10.01\,\%$ | 0.867 | 0.815 | 1.005 | $0.030\,u_\star$ |
| size | `size96` | $9h$ | 9 | $-11.49\,\%$ | 0.806 | 0.750 | 0.972 | $0.143\,u_\star$ |
| size | `size64` | $5h$ | 9 | $-13.52\,\%$ | 0.706 | 0.621 | 0.811 | $0.135\,u_\star$ |

P1 holds: the zone arm spans 0.16 % against per-point spreads of 4–5 %, and the band ratios
drift by 0.02 in the direction of a *wider* zone costing slightly more, i.e. the ramp is not
where the energy comes back. P2 holds: the size arm moves 3.6 % over $8h$ of interior, and the
8–16 m ratio climbs from 0.62 to 0.83, still rising at $13h$. Read together with the cadence
correction above, this is the recovery curve of the band the 3 s boundary removed, at
$C_{\rm dump}=5.4$; C0 supplies the other axis.

**An unplanned finding: clearing the child's zone has a mean-flow cost.** The two smaller
children fail criterion A at $0.14\,u_\star$ where every full-size child passes. Their profiles
show why: `size64` and `size96` remove 12 and 20 of the parent's cubes from their own zone
(§9.4), and their canopy-layer mean wind runs $0.05$–$0.07\,u_\star$ below the parent's up to
$z/h\approx1.5$, decaying to $<0.01\,u_\star$ by $z/h=3$. The parent's cubes still imprint on the
imposed velocity, but the *wakes* those cubes shed inside the zone are absent, so the first
building rows of the interior see an inflow with too little canopy-layer deficit. This is the
adjustment length V3 was designed to measure, seen from the other side; it does not appear in
the full-size children because their zone was building-free in the parent too. It bears on the
choice in §9.4: clearing the child's zone is cheap and correct for the *boundary*, but the
interior then starts its own canopy adjustment at the inner edge.

**V0 filtered arm (job 3993078, 2 h 25).** The child of V1 at $\Delta x = 2$ m, driven by the
converged truth box-filtered onto 4 m ($r=2$) and 8 m ($r=4$) grids at the same 3 s cadence, so
the driving parent is a *perfect* coarse LES. $|\Phi|$ and `divmax` at $10^{-15}$ in both.

| $r$ | TKE above $2h$ | filtered parent's own | canopy TKE | above parent filter | straddling it | below it | criterion A |
|---|---|---|---|---|---|---|---|
| 2 | $-13.1\,\%$ | $-6.9\,\%$ | $-0.1\,\%$ | 1.00 | 0.73 (8–16 m) | 0.87 ($<8$ m) | $0.126\,u_\star$ |
| 4 | $-19.9\,\%$ | $-20.2\,\%$ | $-9.1\,\%$ | 0.98 | 0.65 (16–32 m) | 0.63 ($<16$ m) | $0.275\,u_\star$ |

*The criterion-A failures are a prolongation staircase.* The mean-wind error alternates in sign
between adjacent child levels inside one parent cell — period 2 at $r=2$, period 4 at $r=4$ —
which is the signature of piecewise-constant conservative interpolation of a curved profile:
with $\partial u/\partial z\approx0.07\ \mathrm{s^{-1}}$ at $z/h=2$ a 4 m parent cell puts
$\pm0.07$ m/s $=\pm0.18\,u_\star$ on the boundary, and the interior keeps about half. The writer's
prolongation is being replaced by a conservative piecewise-linear reconstruction (unlimited, so it
stays linear in the data and preserves every parent-face integral; W8 in the plan). V0's
mean-flow row is not readable until that is in.

*Above the canopy the child recovers little of what an $r=4$ parent never resolved, over this
fetch.* Its deficit is $-19.9\,\%$ against its parent's own $-22.8\,\%$: about three points
better than the parent, so a refined child **can** improve on its parent, but only marginally
within $13h$. At $r=2$ it regenerates the sub-filter band to 0.87 yet ends worse than its parent
($-13.1$ against $-9.5\,\%$), because the 8–16 m band is lost twice, to the filter and to the 3 s
cadence (§10.5 correction above). The canopy layer recovers fully at $r=2$ and not at $r=4$,
where 8 m cells sit on 16 m streets. Read with V2: **a band that is missing from the boundary
data, whatever removed it, comes back over a fetch that grows with its wavelength**, and $13h$
is not enough for 16–32 m. The recovery length $L_{\rm rec}(\lambda)$ is the quantity that
governs cadence and refinement alike, and it is what C0 and the V2 size arm together measure.

**An estimator correction that changed the V0 canopy conclusions (2026-09-07).** The
driving-parent profile and the reference profile were not the same statistic. `make_child_case`
averaged each velocity component over horizontal planes *first* and then formed
$\tfrac12\sum(\langle q^2\rangle-\langle q\rangle^2)$, so its variance was taken about the combined
space-and-time mean and therefore included the **dispersive** part — the spatial variance of the
time-mean field. `analyse.Bundle` instead subtracts each cell's own time mean and only then averages
horizontally. `analyse_v0` compared the two directly. Neither definition is wrong: plane averaging is
standard where the flow is horizontally homogeneous, and around buildings the dispersive term is a
real quantity. Comparing one against the other is what was wrong.

The manifest now reports all three — `tke` (temporal, consistent with the reference), `tke_dispersive`
and `tke_total` (the old number) — and a band-only source, which has no interior and one stored
block, reports the TKE comparison as unavailable rather than emitting a figure the analysis would
accept. Measured on V0's own parents, the dispersive term is **3 % of the temporal one above
$z/h=2$ and 85–145 % of it inside the canopy**, which is exactly where buildings lock a spatial
pattern in place. So the correction barely moves the results aloft (deficits go 2–3 points more
negative) and **inverts the sign of every canopy comparison**. The rewritten paragraphs below are
the corrected ones; the child-versus-truth numbers were computed by `analyse.py` on both sides
throughout and never moved.

**V0 coarse arm (job 3994025).** The same children driven by *genuine* 4 m and 8 m LES runs of
the same domain, same forcing, same 3 s cadence — the real use case, unpaired. The result is
dominated by something the filtered arm could not show: **the coarse parents are wrong in the
mean, and the child transmits that error intact.** Against the truth, the 4 m parent's wind is
8–20 % high (about $0.9\,u_\star$ at the child's heights); the 8 m parent, with two cells per cube,
is 10–70 % high. The child's mean wind then runs $0.9$–$1.1\,u_\star$ (r = 2) and
$1.4$–$1.6\,u_\star$ (r = 4) above the truth through the whole interior, so criterion A reads
$1.34$ and $1.80\,u_\star$ — not a nesting error but the parent's, advected through exactly as §0's
table said a bulk-momentum error would be.

The turbulence needs the correction of §10.5's estimator note to read at all. Above the canopy the
child's deficit is close to its parent's own ($-22.0\,\%$ against $-14.5\,\%$ at r = 2; $-18.1$
against $-18.7$ at r = 4, i.e. equal within a point). **In the canopy the earlier reading was
inverted by the statistic, and is now the opposite.** The coarse parents were reported as carrying
$+60$–$100\,\%$ *excess* canopy turbulence, and the child as pulling that back to $+18$ and
$+28\,\%$; measured consistently, the parents carry a canopy *deficit* ($-2.9\,\%$ at 4 m,
$-14.1\,\%$ at 8 m), and the child **overshoots** the truth by those same $+18$ and $+28\,\%$. The
coherent reading is the momentum one: the coarse parent feeds the child too much near-surface
momentum, the child's own resolved buildings turn that into canopy shear production, and the canopy
ends too energetic. Refinement does not rescue the canopy here; it produces a canopy that is too
active because what it is fed is too fast. The paired "cost of a real parent beyond a perfect
filtered one" is $-8.9\,\%$ of TKE at r = 2 and $+1.8\,\%$ at r = 4, against a criterion-A
difference of $1.2$–$1.5\,u_\star$. **The conclusion for the intended use:
a child at higher resolution than its parent reproduces its parent, not the truth; the parent's
mean flow at the boundary is the limiting factor, and an LES with fewer than four cells per
building is not a usable parent for this geometry.**

**V4 — a parent whose geometry differs from the child's (job 3996514, 5 h 05).** Parent: a
*staggered* 16 m array. Child: the V1 *aligned* array, driven through the zone dump at 0.5 s with
the cubic interpolant ($C_{\rm dump}=0.75$). Two references, and they answer different questions.

*Against its own parent* — the perfect-model test — the child is excellent above the canopy and
deliberately different inside it:

| | mean relative difference | spread | significance |
|---|---|---|---|
| resolved TKE, $z/h>2$ | $-1.1\,\%$ | $6.2\,\%$ | $0.18$ |
| resolved TKE, canopy | $+45.7\,\%$ | $15.0\,\%$ | $3.1$–$3.9$ |

The $-1.1\,\%$ aloft is an **independent confirmation of the C0 operating rule on a different
geometry**: at $C_{\rm dump}\le2$ with the cubic, the boundary delivers the parent's turbulence and
the child keeps it, where V1 at $C_{\rm dump}=5.4$ lost $10\,\%$. The $+46\,\%$ in the canopy is the
answer V4 was built to get: **the child's canopy turbulence is generated by the child's own
buildings, not inherited from the parent's**, and at 3–4$\sigma$ it is unambiguous.

*Against the V1 child* — the reference §10.4 nominated — the mean flow appears to fail badly, RMS
$1.88\,u_\star$, criterion A$'$ $=2.85$. It is not a nesting error. Both runs are driven by the same
fixed $\mathrm{d}p/\mathrm{d}x$, and a staggered array has more drag than an aligned one, so the two
parents equilibrate at different bulk winds: $3.46$ m/s staggered against $4.13$–$4.17$ m/s for the
aligned parents. The V4 child's interior wind sits at $0.80$ of the V1 child's, against a parent
bulk ratio of $0.83$. **The child inherits its parent's bulk momentum and does not re-establish the
bulk flow its own geometry would produce** — the same conclusion V0 reached from the other
direction, and the sharpest statement of the scheme's one-way nature.

So the two halves separate cleanly: **turbulence structure adjusts to the child's geometry within
the fetch; bulk momentum does not adjust at all.** §10.4's V4 row asked for "interior statistics
against the V1 child", which conflates them — the mean must be compared normalised by bulk velocity
(or against the child's own parent), and only the turbulence against V1. Two confounds are recorded
rather than corrected: the V1 child ran at 3 s with linear interpolation while this child ran at
0.5 s with the cubic, which C0 prices at $+8$–$9\,\%$ of TKE aloft and accounts for essentially all
of the $+9.3\,\%$ "mismatch" difference the summary table reports; and the staggered layout leaves
$-8$ m of building clearance in the child's zone (§9.4's geometric limit), so the cleared cubes
encroach on the interior by one cell — to be checked before this run is quoted on canopy numbers.

**V3 — a parent that resolves no buildings (job 3996513, 5 h 40).** Five children off two
parents, at 0.5 s with the cubic interpolant. Four have a **building-free parent** and differ only
in the standoff between the inner zone edge and the first building row (0, 5, 15, 40 cells); the
fifth, `cleared-parent-cubes`, has a parent that **does** resolve buildings, including where the
child's zone lies, and the child clears them from its own mask (§9.4). The measure is the residual
canopy-velocity error against the periodic equilibrium, at fixed stations from the inner zone edge:

| child | parent | 5$h$ | 10$h$ | 20$h$ | error at first row | at last row |
|---|---|---|---|---|---|---|
| `standoff0` | no buildings | $91.8\,\%$ | $49.2\,\%$ | $37.6\,\%$ | $193\,\%$ | $92\,\%$ |
| `standoff5` | no buildings | $63.9\,\%$ | $28.2\,\%$ | $23.3\,\%$ | $166\,\%$ | $108\,\%$ |
| `standoff15` | no buildings | $62.7\,\%$ | $17.3\,\%$ | $11.1\,\%$ | $136\,\%$ | $79\,\%$ |
| `standoff40` | no buildings | — | $33.1\,\%$ | $7.9\,\%$ | $115\,\%$ | $57\,\%$ |
| **`cleared-parent-cubes`** | **buildings** | $\mathbf{2.2\,\%}$ | $\mathbf{3.7\,\%}$ | $\mathbf{3.6\,\%}$ | $9.3\,\%$ | $7.1\,\%$ |

Sampling spread $4$–$6\,\%$.

**The supported configuration passes, and passes quickly.** `cleared-parent-cubes` — a parent that
resolves buildings everywhere, including under the child's zone, with the child clearing them from
its own mask — puts the child's canopy inside the sampling spread by $5h$ of fetch and holds it
there ($2.2$, $3.7$, $3.6\,\%$). That is V3's answer for the case the scheme is for, and it is the
row to quote.

**A parent that resolves no buildings cannot drive an urban child at this domain size.** No
building-free-parent child reached equilibrium anywhere in $26h$ of available fetch, so the
adjustment length §9.4 asked for is not measurable — it is longer than the domain. The canopy wind
is wrong by a factor of two at the first building row and still wrong by $38$–$57\,\%$ at the last.
The one child whose parent resolved buildings is inside the sampling spread by $5h$. §9.4 called a
building-free parent "explicitly supported"; mechanically it is, and the flux and divergence
diagnostics are as clean here as anywhere, but **it is not usable**, and the design should say so:
the parent must resolve the canopy, and it is the *parent's* buildings, not the child's own, that
put the canopy deficit into the imposed flow. This is the same one-way property V0 and V4 found —
the child inherits its parent's bulk momentum — seen at its most extreme, because a building-free
parent has no canopy deficit to inherit.

**§9.4's standoff argument is refuted.** It reasoned that a building-free standoff is
counterproductive: the flow adjusts once to ground roughness and then again on reaching the canopy,
where starting the buildings at the zone edge would give one adjustment instead of two. Measured,
every standoff beats zero by more than the sampling spread, monotonically at $20h$
($37.6\to23.3\to11.1\to7.9\,\%$ for 0, 5, 15, 40 cells), and the error at the first building row
falls the same way ($193\to115\,\%$). The second adjustment is not wasted: with a building-free
parent the imposed flow arrives with far too much near-surface momentum, and the building-free
approach sheds some of it before the canopy. Note the scope of the refutation — it is established
only for a parent that resolves no buildings, which is the configuration this row was built around
and, per the paragraph above, the one nobody should use. Whether a standoff helps when the parent
*does* resolve buildings is untested; the `cleared-parent-cubes` arm runs at standoff 0 only.

**Was the building-free parent turbulent at all? Yes, but weakly, and that is the mechanism.**
These runs are neutral — `ltempeq`, `lmoist` and `lbuoyancy` are all `.false.` — so there is no
convection, and with no buildings the only source of turbulence is shear over the ground roughness.
Measured over the canopy depth, against the two parents that do resolve buildings:

| parent | canopy $\langle u\rangle$ | canopy resolved TKE | resolved $\overline{u'w'}$ at roof height |
|---|---|---|---|
| `921`, no buildings | $2.84$ m/s | $0.143$ m$^2$/s$^2$ | $0.034$ m$^2$/s$^2$ |
| `926`, buildings | $0.87$ m/s | $0.343$ m$^2$/s$^2$ | $0.117$ m$^2$/s$^2$ |
| `920`, reference | $0.86$ m/s | $0.358$ m$^2$/s$^2$ | $0.120$ m$^2$/s$^2$ |

So the building-free parent is genuinely turbulent in the ordinary wall-turbulence sense — it is not
laminar — but it carries **2.4$\times$ less resolved turbulent energy and 3.5$\times$ less resolved
momentum flux**, and its near-surface wind is **3.3$\times$ faster** because nothing is extracting
momentum there. That is the whole of the V3 result restated as a property of the *parent*: what the
boundary imposes is not a weakly-wrong version of the child's flow, it is a different flow, and the
child is asked to remove some 70 % of the near-surface momentum before its canopy can be right.
Over $26h$ it does not finish. Nothing about the nesting scheme would fix this, which is why §9.4
now treats a canopy-resolving parent as a requirement.

One asymmetry to note when reading the table: `921` is driven at a fixed volume flow rate, set to
the reference's bulk velocity (4.134 m/s) so that the two are compared at matched bulk flow, while
`926` and `920` are driven at fixed $\mathrm{d}p/\mathrm{d}x$. That is deliberate — it is the only
way to make a building-free domain and a building-resolving one comparable at all — but it means
the two are not at matched *stress*, and the stress difference is exactly the point.

*An analysis defect found here and fixed.* The standoff comparison was keyed by standoff cells, and
`cleared-parent-cubes` also has standoff 0, so it silently overwrote `standoff0` in the reported
table and entered the "beats zero" list as though it were a standoff — while being driven by a
different parent. The sweep is now split by parent and keyed by child (`run_geometry.summarise_v3`).
The verdict is unchanged in direction, but the numbers it reported were one child's under another's
name.

### 10.6 Wiring

Two groups in `tests/test_suites.yml`. `nesting-unit` (included by `supported` and
`supported-macos`, so it runs on all four CI legs): the unit runmodes on one rank, the abort and
message cases, and the one multi-rank case in CI, I5 on 2×2 with a ramp, $\tau>0$ and a cube. The I1
small case is in `supported` directly (Linux legs), against a baseline the driver builds from
`origin/master`. Every suite in the gate carries `UDALES_REQUIRE_LAUNCHER=1`: an unusable MPI
launcher there is a failure, not a skip. `nesting` holds the rest of the integration matrix as
`class: experimental`, `platform: hpc` -- verified on the cluster, not in CI; `run_tests.py --list`
shows which is which. System/validation is `heavy`, never in GitHub Actions. Fixtures under
`tests/cases/`; the analytic-field generator lives in `tools/python/udprep/nesting.py` so the tests
and the production writer cannot drift apart.

Ordering note: **U15–U22 and P1–P11 must pass before I3 is attempted.** The uniform-flow test is
seductive because it exercises everything at once — which also means that when it fails it tells
you nothing. The unit layer exists so that it never has to be the first thing that fails.

---

## 10.7 Open items after the first implementation pass

Recorded here rather than in a tracker so they travel with the design.

1. ~~**`stagger` attribute is written but never checked.**~~ **DONE.** `nestio_validate` now
   compares the `stagger` attribute of every slab variable present against the layout the reader
   assumes, and aborts on a mismatch or a missing attribute. Spec §6 updated; U22 promoted from an
   informational line to an abort case.
2. ~~**The C1 diagnostics are not implemented.**~~ **DONE.** `nesting_stats` now takes `p` as an
   argument — passed from `program.f90`, which already holds it, so the cycle
   `modnesting → modpois → modboundary → modnesting` is avoided — and reports
   $\|\mathcal{G}p\|_{\rm zone}$, $\|\mathcal{G}p\|_{\rm interior}$ and their ratio, plus the
   energy injection split between guard strip and relaxation ramp (accumulated in `nesting_apply`).
3. ~~**Init recomputes $\Phi$ from the boundary slabs for every stored time.**~~ **DONE, measured.**
   Schema 2 stores `flux_residual(time)` — the residual of the data *as stored* — next to the
   `fluid_lateral_area` it was summed over, and `check_stored_flux` validates that instead of
   reading anything. The full recompute is kept behind `nest_lfluxcheckall` (default `.false.`) and
   the reader falls back to it, with a named warning, for a schema 1 file or when the file's fluid
   lateral area is not this run's — so the cheap path can never silently paper over a writer/solver
   mask mismatch. Note that §6.1 always asked for "the residual $\Phi(t)$ **left by** the offline
   correction"; the v1 file stored the residual *before* it instead, so this closes a gap between
   the contract and the implementation rather than extending the contract.

   **Measured** on CX3, Release build, `nesting_init` timed around `check_stored_flux`, on two
   13 GB files built by the production writer:

   | Case | full recompute, **cold** | full recompute, cached | stored residual |
   |---|---|---|---|
   | $256^2\times128$, $n_z=16$, 256 levels; 17.8 MB read per level, 4.6 GB total, 1024 reads | **35.4 s** | 1.4–1.6 s (0.62 s on 2×2) | **9.6 µs** |
   | $128^2\times64$, $n_z=8$, 2048 levels; 295 kB read per level, 604 MB total, 8192 reads | **38.7 s** | 1.5–1.7 s (0.70–0.81 s on 2×2) | **12 µs** |

   "Cold" is the first read after the file has left the filesystem client's cache; "cached" is any
   read after that. **The distinction is the whole story, and it is easy to measure the wrong one**
   — the first numbers taken here were the cached ones, straight after writing the fixture, and they
   made the check look cheap. Cold, the two cases cost 130 MB/s and 16 MB/s of *effective*
   bandwidth respectively: 34 ms per 4.5 MB read in the first, 4.7 ms per 74 kB read in the second.
   The second is not bandwidth at all, it is ~5 ms of latency per `nf90_get_var` — design §6.2's
   pathology 2, seen in the wild.

   A production file is worse on both axes. $512^2\times128$, $n_z=16$ with 4000 levels reads
   35.7 MB per level, 143 GB over 16000 calls; at the cold rates above that is **≈20 minutes** of
   initialisation before the first timestep, on every restart of every run. The stored residual
   removes it: ~10 µs, no slab read, cold or warm, and independent of `ntime` in practice.
4. ~~**Cold-start initialisation of `u0`/`um` from the parent is not implemented.**~~ **DONE.**
   Schema 2 carries an optional full-domain block `u_init`/`v_init`/`w_init` at the first stored
   time, and `nest_linitfromparent` (default `.false.`) fills `u0/um`, `v0/vm`, `w0/wm` from it on a
   **cold start only** — a warm start already holds a consistent state and overwriting it would
   break restart parity (I6), so the switch is ignored there with a message. The writer makes the
   block usable rather than merely present: it takes the boundary-normal velocities from the
   *corrected* slabs at that time, closes the floor and the lid ($w=0$, case A), and then projects
   the whole 3-D field onto the discretely solenoidal subspace with the solver's own operators —
   a DCT-II in $x$ and $y$ and a tridiagonal sweep in the stretched vertical, homogeneous Neumann
   pressure on all six faces, so every boundary-normal velocity survives untouched (F2) and an
   already-solenoidal parent comes back unchanged (bit for bit when its discrete divergence is
   exactly zero, to $10^{-15}$ when it is merely at round-off). Errors, not warnings, on: the switch set
   with no block; a block at the wrong shape; a block at the wrong stagger.

   One constraint fell out of the derivation and is now enforced: the projection uses the solver's
   **density-free** divergence (F1: `fillps` carries no $\rho$) while the slab flux correction is
   $\rho$-weighted, so a file carrying an initial condition must have `rhobf == rhobh == 1`. That is
   always true in uDALES; the writer refuses anything else rather than storing a block whose
   boundary flux does not close.
5. ~~**Case B (`BCtopm_pressure`) trips the flux assertion.**~~ **DONE**, and the design's own
   analysis of §3.2 settles it. Under a leaky lid `bcpup` sets $w^\ast_{ktot+1}$ from the accumulated
   pressure and `tderive` adds the matching increment $2\langle p\rangle_{ktot}/\Delta z_h$ — which is
   *exactly* the Dirichlet-in-the-mean-mode row the solver pins (F3). So $\tilde{\mathcal{L}} =
   \mathcal{DG}$ under the BCs actually applied, the projection is complete for **any** $\Phi$, and
   $\Phi=0$ is not a solvability requirement in case B at all: the lid flux is the child breathing
   against its reservoir, not an error. Asserting on the six-face $\Phi$ there was asserting a
   property the scheme is not required to have.

   `nest_flux_split` therefore reports $\Phi$ and the part the lid carries, and `nesting_bcpup`
   asserts on the *closed* faces — $\Phi$ under a rigid lid, $\Phi-\Phi_{\rm lid}$ under a leaky one.
   Under a rigid lid `bcpup` forces $w^\ast=0$ at the top, $\Phi_{\rm lid}$ is identically zero and
   nothing changes. `nest_lfluxassert` is now on by default in both cases.

   Worth recording, because it made the test hard to write: with a flux-balanced parent and a cold
   start, case B's lid **never** moves. The lid velocity is driven by $\langle p^{\rm acc}\rangle_{ktot}$,
   the pin makes that proportional to $-\Phi$, and $\Phi$ is in turn the lid flux — an autonomous
   linear feedback started from rest, so it stays at round-off for ever ($\Phi_{\rm lid}\sim8\times10^{-18}$
   measured). The assertion only ever fires when the child *arrives* with a column-pressure excess,
   which is what I10 constructs by offsetting `pres0` in a restart file: there $\Phi_{\rm lid}=1.0\times10^{-2}$,
   eight orders of magnitude above `nest_fluxtol`, while $\Phi_{\rm closed}=3\times10^{-17}$ and
   `divmax` stays at $8.9\times10^{-16}$ — case B's substantive claim, confirmed.
6. ~~**`nesting_init` ran before `readinitfiles`**~~ **DONE.** It reads `timee`, which
   `readinitfiles` is what assigns; the two calls are now in that order (`program.f90:103-108`).
   Found by I6; see §9.5 for the failure modes it caused. Nothing between the old and new call
   sites needs nesting — `nesting_boundary`/`_apply`/`_bcpup`/`_stats` all return while `linit` is
   false, and every `call boundary` inside `readinitfiles` is commented out — and `nesting_init`'s
   own prerequisites `createmasks`/`calcfluidvolumes` still precede it.
7. **I1 cannot be a bitwise test on a non-trivial case.** `FFTW_MEASURE` makes the solver
   irreproducible against itself at $\sim5\times10^{-12}$; see the I1 row in §10.3. If a genuinely
   bitwise regression gate is ever wanted, `FFTW_ESTIMATE` behind a build flag would give it, at a
   performance cost that has not been measured. A second, distinct noise source showed up on the
   small `lnesting = .false.` case (`TestI1NoOpSmallCase`, ubuntu-latest Release): the baseline and
   branch binaries are two *different* compiles (`origin/master` vs. this branch), and even on the
   shared, unchanged no-op code path, cross-build `-O3` codegen drift (vectorisation, instruction
   scheduling, FMA contraction) moved results by $\sim10^{-15}$ to $2.3\times10^{-12}$ relative — a
   same-binary self-noise reference cannot see this, since it is one binary run twice. The test's
   bound was floored at this file's own round-off constant, $1\times10^{-9}$, instead.

---

## 11. Code reference index

**uDALES** `master` @ `1f8ff3e8`

| Topic | Location |
|---|---|
| time loop / hook points | `src/program.f90:132-222` |
| RK3 coefficient, integration, tendency reset | `src/modtstep.f90:191`, `:206-231`, `:322-338` |
| $\mathcal{D}$, predicted velocity | `src/modpois.f90:911-998` (esp. `:942`, `:963`, `:966-973`) |
| $\mathcal{G}$, projection, top-BC increment | `src/modpois.f90:1001-1105` (esp. `:1044-1055`, `:1057-1067`) |
| $\mathcal{L}$: Neumann eigenvalues, tridiagonal, top pin | `src/modpois.f90:111-146`, `:148-176`, `:205-215` |
| backend capability limits | `src/modpois.f90:226-253`, `:298-352` |
| `bcpup` face imposition (profile / driver) | `src/modboundary.f90:1191-1341`, `:1256`, `:1283` |
| `bcp` Neumann pressure | `src/modboundary.f90:1344-1430` |
| ghost filling, inflow profile / driver | `src/modboundary.f90:688-763` |
| convective outflow (present, deliberately **not** used — §7 C6) | `src/modboundary.f90:908-926` |
| sponge layer | `src/modboundary.f90:44-52`, `:1447-1491` |
| `masscorr` and its `linoutflow` guard | `src/modforces.f90:328-497` |
| local→global index idiom | `src/modforces.f90:953-980`, `src/modpois.f90:200-207` |
| IBM solid handling and masks | `src/modibm.f90:716-764`, `:767-845`, `:2122-2253` |
| $\rho_f,\rho_h$ (identically 1) | `src/modfields.f90:380-381`, `:571-572` |
| grid, BC constants, index ranges | `src/modglobal.f90:97-176`, `:611-661`, `:718-780` |
| namelists, BC validation, `BCtopm` auto-switch | `src/modstartup.f90:105-176`, `:694-941`, `:845` |
| decomposition setup | `src/modstartup.f90:652-691` |
| divergence diagnostics | `src/modchecksim.f90:161-203` |
| per-rank binary driver I/O (anti-pattern) | `src/moddriver.f90:750-931` |
| NetCDF read precedent, linkage | `src/initfac.f90:265-271`; `CMakeLists.txt:100-101`, `:133` |
| decomposition-parity test pattern | `tests/integration/processor_boundaries/` |

**DALES** `franslql/dales` @ `v4.4_openBC`

| Topic | Location |
|---|---|
| NetCDF boundary input, per-rank hyperslab | `src/modopenboundary.f90:268-442` |
| init-time divergence correction | `src/modopenboundary.f90:444-576` |
| ghost filling | `src/modopenboundary.f90:578-621` |
| normal-velocity tendencies (radiation / nudging) | `src/modopenboundary.f90:884-1019` |
| per-substep patch mass correction | `src/modopenboundary.f90:1021-1176` |
| phase-speed estimator | `src/modopenboundary.f90:648-743` |
| anelastic Poisson RHS, Neumann $p$ | `src/modpois.f90:161-243`, `:276-309` |
| call ordering | `src/program.f90:213-295` |
| relaxation prototype | `src/modnudgeboundary.f90` (weights `:57-112`, per-rank I/O `:276-317`, forcing `:359-440`) |

**PALM** `v23.04`: `pmc_interface_mod.f90:6207-6430` (masked mass conservation),
`:2426-2520` (fluid face areas), `:6041-6197` (post-interpolation BCs).

**Papers.** Liqui Lung, Jakob, Siebesma & Jansson (2024), GMD **17**, 4053–4076,
doi:10.5194/gmd-17-4053-2024 — §2.1.3 mass conservation and homogeneous Neumann $p$, §2.2.2 Robin
inflow, Conclusions for the $\le4$/$\le30$ limits. Hellsten et al. (2021), GMD **14**, 3185–3214,
doi:10.5194/gmd-14-3185-2021. Craske & van Reeuwijk (2013), *Comput. Fluids* **86**, 284–297 — the
basis for Neumann-on-everything-but-the-normal-velocity.
