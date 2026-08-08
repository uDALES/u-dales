# Code architecture

This page is a tour of the uDALES source code for developers who are about to modify the solver: where things live, how the Fortran modules are organised, and what happens — in order — during one timestep.

## Repository layout

- `src/` — the Fortran MPI solver (one flat directory of `mod*.f90` files plus `program.f90`).
- `tools/python/` — the current Python tooling: `udbase` (postprocessing), `udprep` (preprocessing), `udgeom` (geometry), plus tests and example notebooks.
- `tools/matlab/` — legacy MATLAB postprocessing and geometry helpers.
- `tools/IBM/` — immersed-boundary preprocessing utilities.
- `tools/SEB/` — surface energy balance preprocessing utilities.
- `tools/syntheticInflow/` — synthetic-inflow generation tooling.
- `tools/View3D/` — external view-factor calculator (git submodule) used for radiation preprocessing.
- `tools/*.sh` — build and run wrappers (`build_executable.sh`, `hpc_execute.sh`, ...).
- `examples/` — runnable example cases (`001`, `002`, ..., `999`) with reference inputs.
- `tests/` — unit, system, integration and regression tests; see `tests/README.md`.
- `docs/` — this documentation (MkDocs).
- `2decomp-fft/` — the 2DECOMP&FFT parallel decomposition/FFT library (git submodule).

## Module organisation

Everything in `src/` compiles into a single executable whose entry point is `program.f90`. The modules fall into a few natural groups.

### State and configuration

- `modglobal.f90` — global constants, grid dimensions and indices (`ib:ie`, `jb:je`, `kb:ke`), halo widths, switches, and run-mode constants; `initglobal` derives grid/timestep settings.
- `modfields.f90` — declares and allocates the 3-D prognostic and diagnostic field arrays (`initfields`).
- `modmpi.f90` — MPI layer: communicators, rank/neighbour info, reduction helpers.

### Startup and input

- `modstartup.f90` — reads the `namoptions` namelists (`readnamelists`), initialises 2DECOMP (`init2decomp`), sanity-checks settings (`checkinitvalues`), and reads initial/restart fields and profiles (`readinitfiles`).
- `readinput.f90` — generic readers for the sparse (i,j,k)-indexed input file formats.
- `initfac.f90` — reads the facet input files (facet geometry, types, properties).

### Dynamics and numerics

- `modadvection.f90` — advection schemes for momentum and scalars.
- `modsubgrid.f90` / `modsubgriddata.f90` — subgrid-scale diffusion (Vreman, Smagorinsky, one-equation TKE models) and its variables.
- `modforces.f90` — remaining momentum sources: body forces, Coriolis, large-scale forcing, nudging, mass-flow correction.
- `modpois.f90` — the FFT-based Poisson solver for the pressure correction.
- `modtstep.f90` — adaptive timestepping and the third-order Runge–Kutta integration (Wicker & Skamarock, 2002).

### Physics

- `modthermodynamics.f90` — diagnostic thermodynamics (liquid water, buoyancy-related quantities).
- `modEB.f90` — the facet surface energy balance.
- `vegetation.f90` — sparse vegetation: momentum drag, canopy energy balance, scalar deposition.
- `modchem.f90` — NO/NO2/O3 null-cycle chemistry, applied once per full timestep.
- `modscalsource.f90` — point/line/planar scalar sources.
- `modpurifiers.f90`, `heatpump.f90` — air purifiers and heat-pump exhaust.
- `modsurfdata.f90` — surface-model variables.

### Boundaries and immersed boundary method (IBM)

- `modboundary.f90` — all domain boundary conditions (except immersed boundaries), top-of-domain damping, and halo updates.
- `modibm.f90` / `modibmdata.f90` — immersed-boundary forcing: wall shear (`ibmwallfun`), zero normal velocities (`ibmnorm`), floor treatment (`bottom`).
- `modwallfunctions.f90` — wall-function fluxes (`wfuno`, `wfmneutral`) used by the IBM and facet code.
- `modinlet.f90` / `modinletdata.f90` — turbulent inflow generation by rescale–recycle (Lund et al.).
- `moddriver.f90` — driver simulations: writing and replaying precursor inlet planes.
- `modtimedep.f90` — time-dependent boundary values and forcings.

### Output and monitoring

- `modsave.f90` — restart files (`writerestartfiles`).
- `modfielddump.f90` — instantaneous 3-D field output (`fielddump*.nc`).
- `modstatsdump.f90` / `modstatistics.f90` — time- and space-averaged statistics output.
- `modstat_nc.f90` — low-level NetCDF writing routines used by the dump modules.
- `modchecksim.f90` — monitors Courant/Péclet numbers and divergence during the run.
- `tests.f90` — in-solver test routines dispatched via special `runmode` values.

## Life of a timestep

`program.f90` is short and readable; it is the best map of the solver. It has three phases.

**Startup.** In order: `initmpi`, `readnamelists`, `init2decomp`, `checkinitvalues`, `initglobal`, then the `init*` routines of the other modules (`initfields`, `initboundary`, `initthermodynamics`, `initsubgrid`, `initdriver`, `initpois`), the IBM/facet setup (`readfacetfiles`, `initibm`, `createmasks`, `calcfluidvolumes`), the initial fields (`readinitfiles`, `createscals`), and finally the output and add-on modules (`initstatsdump`, `initEB`, `inittimedep`, `initfielddump`, `boundary`, `init_vegetation`, `createpurifiers`, `init_heatpump`).

**Time loop.** The loop is `do while ((timeleft>0) .or. (rk3step < 3))`, and — important — **one pass through the loop is one RK3 substep**, not one full timestep. Each substep, every process adds its contribution to the tendency arrays (`up`, `vp`, `wp`, `thlp`, ...), and `tstep_integrate` applies them at the end. The ordered calls, with the module implementing each:

1. `tstep_update` (`modtstep`) — advance `rk3step` (1→2→3→1...); at the start of substep 1 choose a new `dt` from the CFL and diffusion-number limits (if `ladaptive`).
2. `timedep` (`modtimedep`) — update time-dependent forcings.
3. `advection` (`modadvection`) — advection tendencies (includes the predicted pressure-gradient term).
4. `shiftedPBCs` (`modforces`) — shifted periodic boundary conditions.
5. `subgrid` (`modsubgrid`) — subgrid-scale diffusion tendencies.
6. `bottom` (`modibm`) — surface layer at the domain floor.
7. `coriolis`, `forces`, `lstend`, `nudge` (`modforces`) — Coriolis, remaining body forces, large-scale forcings, top-cell nudging.
8. `ibmwallfun` (`modibm`) — immersed-boundary wall shear (via `modwallfunctions`), then `periodicEBcorr` (`modforces`).
9. `masscorr` (`modforces`) — correct the predicted velocity for the prescribed mass flow.
10. `ibmnorm` (`modibm`) — set velocity components normal to immersed boundaries to zero.
11. `EB` (`modEB`), `vegetation_forcing` (`vegetation`), `heatpump` (`modheatpump`), `scalsource` (`modscalsource`) — facet energy balance, vegetation, heat pump, and scalar-source terms.
12. `fixuinf2`, `fixuinf1` (`modforces`) — free-stream velocity controllers.
13. `grwdamp` (`modboundary`) — gravity-wave damping near the top of the domain.
14. `poisson` (`modpois`) — solve the Poisson equation for the pressure correction and project the velocity field to be divergence-free.
15. `purifiers` (`modpurifiers`).
16. `tstep_integrate` (`modtstep`) — the RK3 update: `u0 = um + rk3coef * up` with `rk3coef = dt/(4 - rk3step)`, then reset the tendencies to zero; on substep 3, copy `u0` into `um` (chemistry via `modchem` is also applied here, once per full step).
17. `halos` (`modboundary`) — exchange halo data between neighbouring pencils.
18. `checksim` (`modchecksim`), `fielddump` (`modfielddump`), `statsdump` (`modstatsdump`) — monitoring and output.
19. `boundary` (`modboundary`) — apply domain boundary conditions to the updated fields.
20. `thermodynamics` (`modthermodynamics`) — diagnostic thermodynamic fields.
21. `writerestartfiles` (`modsave`) — restart output when due.

**Finalisation.** `exitfielddump`, `exitstatsdump`, `exit_heatpump`, `exitmpi`.

When adding a new physical process, the usual pattern is: a module with an `init*` routine called from the startup phase, and a forcing routine called from the loop that adds to the tendency arrays *before* `poisson`/`tstep_integrate`.

## Field arrays and grid

uDALES uses a staggered grid: scalars (`thl`, `qt`, `sv`, pressure) live at cell centres, and each velocity component on the cell face in its own direction (see the [fields tutorial](udales-fields-tutorial.md) for the exact layout).

Each prognostic variable exists at three time levels (declared in `modfields.f90`):

- `um`, `vm`, `wm`, `thlm`, ... — the field at the previous full timestep (the RK3 base state);
- `u0`, `v0`, `w0`, `thl0`, ... — the current state, updated after every RK3 substep;
- `up`, `vp`, `wp`, `thlp`, ... — the accumulated tendencies for the current substep.

Physics and dynamics routines only ever *add* to the `*p` arrays; `tstep_integrate` combines them as `u0 = um + rk3coef*up` and zeroes them again.

Arrays are allocated per MPI rank over the local pencil `(ib:ie, jb:je, kb:ke)` (all starting at 1) extended by halo (ghost) cells of width `ih`, `jh`, `kh` — e.g. `up(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh)`. The halo widths are set in `modglobal.f90` based on the widest advection stencil in use (up to 3 cells in x and y). The domain is split among MPI ranks with a 2-D pencil decomposition, and halo cells are filled by 2DECOMP&FFT exchange routines — see [parallelisation](udales-2decomp.md) for the decomposition and its constraints on grid loops.
