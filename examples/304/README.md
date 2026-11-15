# 304: dry convective boundary layer over a floor of fixed temperature

Validation case for the Monin-Obukhov facets with a facet temperature (`iwallmom = 4`, `iwalltemp = 2`), the
uDALES analogue of the DALES surface with a prescribed surface temperature (`isurf = 2`). Not a published case:
the setup is chosen to give a clear, quickly developing convective boundary layer driven by a surface flux that
follows from Monin-Obukhov similarity (not prescribed), small enough for a few PBS hours per seed ensemble.
Run on CX3 with `tools/hpc_execute.sh examples/304` (16 ranks, 4 x 4); compare with
`tools/python/examples/dales_cases/compare_dales.py 304 ...` (see the README there). Results: [RESULTS.md](RESULTS.md).
Outputs of interest: `bldump_t.304.nc` (`ustar`, `obl`, `wtheta`, `zi`), `bldump_z.304.nc`.

Inputs are written by `python tools/python/examples/dales_cases/make_cbl304.py examples/304
tools/python/examples/dales_cases/dales_inputs/cbl304` (profiles, forcing, the flat-floor IBM files with
`Tfacinit.inp.304`, and the DALES `prof.inp.001`, `lscale.inp.001`, `baseprof.inp.001`); the namoptions files of
both codes are written by hand.

## Floor

One horizontal facet covering the domain (`libm = .true.`, `lbottom = .false.`, `nfcts = 1`; `udprep.flat_floor`
with `--z1 12.5 --z0 0.1 --z0h 0.01 --no-lgr --tfacinit 302`). The facet temperature is 302 K (`Tfacinit.inp.304`)
and stays there: `lEB = .false.`, and without the energy balance nothing changes `facT`. The facet type has
`lGR = 0` (dry case, no moisture exchange). `iwallmom = 4` with `iwalltemp = 2`: per facet section, L from the bulk
Richardson number between the first-level air and the facet, `u*^2 = Cm U^2`, heat flux `Cs U (Tfacet - thl)`.

## Settings and sources

| setting | value | source |
|---|---|---|
| grid | 96 x 96 x 80, dx = dy = 50 m, dz = 25 m (z up to 2000 m) | case choice: zi grows from 500 to about 800 m, well below the sponge (top quarter) |
| duration, step | 14400 s (4 h), fixed `dtmax = 1.0` s (`ladaptive = .false.`) | case choice; DALES the same fixed step |
| initial thl | 300 K up to 500 m, +6 K/km above | case choice |
| moisture | dry (`lmoist = .false.`) | DALES: `lmoist = .true.` with `qt = 0` and a dry surface `rsisurf2 = 1e10 s/m`, because DALES `isurf = 2` sets `qskin` only in `qtsurf` (with `lmoist`): dry DALES would use an uninitialised `qskin` in the Richardson number |
| wind | u = ug = 5 m/s, v = vg = 0, `lcoriol`, `xlat = 45` | case choice |
| surface | `thls` = facet temperature 302 K, z0 = 0.1 m, z0h = 0.01 m | DALES `isurf = 2`, `thls = 302`, `z0mav`, `z0hav`; `lmostlocal` default (`.true.`) |
| `ps` | 1e5 Pa | case choice; DALES `baseprof.inp.001` (constant density, `ibas_prf = 5`) from the uDALES hydrostatic integral |
| SGS | `loneeqn`, `ldelta = .false.` | DALES default one-equation model, `ldelta = .false.` |
| advection | cd2 for all variables | DALES `iadv_* = 2` |
| sponge | `igrw_damp = 2`, `ksp` default | DALES default |
| initial e12 | 0.3 m/s at all levels | both codes read `e12` (column 6) |
| noise | `randthl = 0.1` K on all 80 levels, no velocity or qt noise, seeds 43-45 | DALES `randthl`, `krandumax = 0`, `randqt = 0` |
| Galilean shift | none | DALES `cu = cv = 0` |
| decomposition | `nprocx = nprocy = 4` | DALES the same |
| bldump | every 60 s | DALES `ltimestat`/`lstat`, 60 s samples, 600 s profile blocks |
