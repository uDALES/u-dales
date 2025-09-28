# uDALES Variable Index Report

This report shows the status of variables across different aspects of the uDALES code:

- **Namelist**: Variable is defined in a Fortran namelist (✓/✗)
- **JSON Read**: Variable is read from JSON configuration (count or ✗)
- **Broadcast**: Variable is broadcast via MPI (count or ✗)
- **Schema**: Variable is present in the JSON schema (✓/✗)

## BC

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `bcbotm` | ✓ | ✓ | ✓ | ✓ |
| `bcbotq` | ✓ | ✓ | ✓ | ✓ |
| `bcbots` | ✓ | ✓ | ✓ | ✓ |
| `bcbott` | ✓ | ✓ | ✓ | ✓ |
| `bcqfxm` | ✓ | ✓ | ✓ | ✓ |
| `bcqfxp` | ✓ | ✓ | ✓ | ✓ |
| `bcqfym` | ✓ | ✓ | ✓ | ✓ |
| `bcqfyp` | ✓ | ✓ | ✓ | ✓ |
| `bcqfz` | ✓ | ✓ | ✓ | ✓ |
| `bctfxm` | ✓ | ✓ | ✓ | ✓ |
| `bctfxp` | ✓ | ✓ | ✓ | ✓ |
| `bctfym` | ✓ | ✓ | ✓ | ✓ |
| `bctfyp` | ✓ | ✓ | ✓ | ✓ |
| `bctfz` | ✓ | ✓ | ✓ | ✓ |
| `bctopm` | ✓ | ✓ | ✓ | ✓ |
| `bctopq` | ✓ | ✓ | ✓ | ✓ |
| `bctops` | ✓ | ✓ | ✓ | ✓ |
| `bctopt` | ✓ | ✓ | ✓ | ✓ |
| `bcxm` | ✓ | ✓ | ✓ | ✓ |
| `bcxq` | ✓ | ✓ | ✓ | ✓ |
| `bcxs` | ✓ | ✓ | ✓ | ✓ |
| `bcxt` | ✓ | ✓ | ✓ | ✓ |
| `bcym` | ✓ | ✓ | ✓ | ✓ |
| `bcyq` | ✓ | ✓ | ✓ | ✓ |
| `bcys` | ✓ | ✓ | ✓ | ✓ |
| `bcyt` | ✓ | ✓ | ✓ | ✓ |
| `bczp` | ✓ | ✓ | ✓ | ✓ |
| `ds` | ✓ | ✓ | ✓ | ✓ |
| `qt_top` | ✓ | ✓ | ✓ | ✓ |
| `qts` | ✓ | ✓ | ✓ | ✓ |
| `thl_top` | ✓ | ✓ | ✓ | ✓ |
| `thls` | ✓ | ✓ | ✓ | ✓ |
| `wqsurf` | ✓ | ✓ | ✓ | ✓ |
| `wsvsurfdum` | ✓ | ✓ | ✓ | ✓ |
| `wsvtopdum` | ✓ | ✓ | ✓ | ✓ |
| `wtsurf` | ✓ | ✓ | ✓ | ✓ |
| `wttop` | ✓ | ✓ | ✓ | ✓ |
| `z0` | ✓ | ✓ | ✓ | ✓ |
| `z0h` | ✓ | ✓ | ✓ | ✓ |

## CHEMISTRY

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `jno2` | ✓ | ✓ | ✓ | ✓ |
| `k1` | ✓ | ✓ | ✓ | ✓ |
| `lchem` | ✓ | ✓ | ✓ | ✓ |

## DOMAIN

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `itot` | ✓ | ✓ | ✓ | ✓ |
| `jtot` | ✓ | ✓ | ✓ | ✓ |
| `ksp` | ✓ | ✓ | ✓ | ✓ |
| `ktot` | ✓ | ✓ | ✓ | ✓ |
| `xday` | ✓ | ✓ | ✓ | ✓ |
| `xlat` | ✓ | ✓ | ✓ | ✓ |
| `xlen` | ✓ | ✓ | ✓ | ✓ |
| `xlon` | ✓ | ✓ | ✓ | ✓ |
| `xtime` | ✓ | ✓ | ✓ | ✓ |
| `ylen` | ✓ | ✓ | ✓ | ✓ |

## DRIVER

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `chunkread_size` | ✓ | ✓ | ✓ | ✓ |
| `driverjobnr` | ✓ | ✓ | ✓ | ✓ |
| `driverstore` | ✓ | ✓ | ✓ | ✓ |
| `dtdriver` | ✓ | ✓ | ✓ | ✓ |
| `iangledeg` | ✓ | ✓ | ✓ | ✓ |
| `idriver` | ✓ | ✓ | ✓ | ✓ |
| `iplane` | ✓ | ✓ | ✓ | ✓ |
| `lchunkread` | ✓ | ✓ | ✓ | ✓ |
| `tdriverstart` | ✓ | ✓ | ✓ | ✓ |

## DYNAMICS

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `iadv_mom` | ✓ | ✓ | ✓ | ✓ |
| `iadv_qt` | ✓ | ✓ | ✓ | ✓ |
| `iadv_sv` | ✓ | ✓ | ✓ | ✓ |
| `iadv_thl` | ✓ | ✓ | ✓ | ✓ |
| `iadv_tke` | ✓ | ✓ | ✓ | ✓ |
| `ipoiss` | ✓ | ✓ | ✓ | ✓ |
| `lqlnr` | ✓ | ✓ | ✓ | ✓ |

## ENERGYBALANCE

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `bldt` | ✓ | ✓ | ✓ | ✓ |
| `dteb` | ✓ | ✓ | ✓ | ✓ |
| `flrt` | ✓ | ✓ | ✓ | ✓ |
| `fraction` | ✓ | ✓ | ✓ | ✓ |
| `grlai` | ✓ | ✓ | ✓ | ✓ |
| `lconstw` | ✓ | ✓ | ✓ | ✓ |
| `leb` | ✓ | ✓ | ✓ | ✓ |
| `lfactlyrs` | ✓ | ✓ | ✓ | ✓ |
| `lperiodicebcorr` | ✓ | ✓ | ✓ | ✓ |
| `lvfsparse` | ✓ | ✓ | ✓ | ✓ |
| `lwriteebfiles` | ✓ | ✓ | ✓ | ✓ |
| `nfaclyrs` | ✓ | ✓ | ✓ | ✓ |
| `nnz` | ✓ | ✓ | ✓ | ✓ |
| `rsmin` | ✓ | ✓ | ✓ | ✓ |
| `sinkbase` | ✓ | ✓ | ✓ | ✓ |
| `skylw` | ✓ | ✓ | ✓ | ✓ |
| `wfc` | ✓ | ✓ | ✓ | ✓ |
| `wgrmax` | ✓ | ✓ | ✓ | ✓ |
| `wsoil` | ✓ | ✓ | ✓ | ✓ |
| `wwilt` | ✓ | ✓ | ✓ | ✓ |

## HEATPUMP

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `lfan_hp` | ✓ | ✓ | ✓ | ✓ |
| `lheatpump` | ✓ | ✓ | ✓ | ✓ |
| `nhppoints` | ✓ | ✓ | ✓ | ✓ |
| `q_dot_hp` | ✓ | ✓ | ✓ | ✓ |
| `qh_dot_hp` | ✓ | ✓ | ✓ | ✓ |

## INFO

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `dtin` | ✓ | ✓ | ✓ | ✓ |
| `jgtotinl` | ✓ | ✓ | ✓ | ✓ |
| `kmaxin` | ✓ | ✓ | ✓ | ✓ |
| `nprocsinl` | ✓ | ✓ | ✓ | ✓ |
| `totalreadu` | ✓ | ✓ | ✓ | ✓ |
| `wtop` | ✓ | ✓ | ✓ | ✓ |

## INLET

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `di` | ✓ | ✓ | ✓ | ✓ |
| `dti` | ✓ | ✓ | ✓ | ✓ |
| `inletav` | ✓ | ✓ | ✓ | ✓ |
| `lfixinlet` | ✓ | ✓ | ✓ | ✓ |
| `lfixutauin` | ✓ | ✓ | ✓ | ✓ |
| `linletra` | ✓ | ✓ | ✓ | ✓ |
| `lreadminl` | ✓ | ✓ | ✓ | ✓ |
| `lstoreplane` | ✓ | ✓ | ✓ | ✓ |
| `lwallfunc` | ✓ | ✓ | ✓ | ✓ |
| `uinf` | ✓ | ✓ | ✓ | ✓ |
| `vinf` | ✓ | ✓ | ✓ | ✓ |

## NAMSUBGRID

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `c_vreman` | ✓ | ✓ | ✓ | ✓ |
| `cf` | ✓ | ✓ | ✓ | ✓ |
| `cn` | ✓ | ✓ | ✓ | ✓ |
| `cs` | ✓ | ✓ | ✓ | ✓ |
| `lbuoycorr` | ✓ | ✓ | ✓ | ✓ |
| `ldelta` | ✓ | ✓ | ✓ | ✓ |
| `lmason` | ✓ | ✓ | ✓ | ✓ |
| `loneeqn` | ✓ | ✓ | ✓ | ✓ |
| `lsmagorinsky` | ✓ | ✓ | ✓ | ✓ |
| `lvreman` | ✓ | ✓ | ✓ | ✓ |
| `nmason` | ✓ | ✓ | ✓ | ✓ |
| `prandtl` | ✓ | ✓ | ✓ | ✓ |
| `rigc` | ✓ | ✓ | ✓ | ✓ |

## OUTPUT

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `fieldvars` | ✓ | ✓ | ✓ | ✓ |
| `islice` | ✓ | ✓ | ✓ | ✓ |
| `jslice` | ✓ | ✓ | ✓ | ✓ |
| `kslice` | ✓ | ✓ | ✓ | ✓ |
| `lfielddump` | ✓ | ✓ | ✓ | ✓ |
| `lislicedump` | ✓ | ✓ | ✓ | ✓ |
| `ljslicedump` | ✓ | ✓ | ✓ | ✓ |
| `lkslicedump` | ✓ | ✓ | ✓ | ✓ |
| `lmintdump` | ✓ | ✓ | ✓ | ✓ |
| `ltdump` | ✓ | ✓ | ✓ | ✓ |
| `ltkedump` | ✓ | ✓ | ✓ | ✓ |
| `lxydump` | ✓ | ✓ | ✓ | ✓ |
| `lxytdump` | ✓ | ✓ | ✓ | ✓ |
| `lydump` | ✓ | ✓ | ✓ | ✓ |
| `lytdump` | ✓ | ✓ | ✓ | ✓ |
| `slicevars` | ✓ | ✓ | ✓ | ✓ |
| `tcheck` | ✓ | ✓ | ✓ | ✓ |
| `tfielddump` | ✓ | ✓ | ✓ | ✓ |
| `tsample` | ✓ | ✓ | ✓ | ✓ |
| `tstatsdump` | ✓ | ✓ | ✓ | ✓ |
| `tstatstart` | ✓ | ✓ | ✓ | ✓ |

## PHYSICS

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `dpdx` | ✓ | ✓ | ✓ | ✓ |
| `ifixuinf` | ✓ | ✓ | ✓ | ✓ |
| `igrw_damp` | ✓ | ✓ | ✓ | ✓ |
| `lbuoyancy` | ✓ | ✓ | ✓ | ✓ |
| `lcoriol` | ✓ | ✓ | ✓ | ✓ |
| `lmoist` | ✓ | ✓ | ✓ | ✓ |
| `lnudge` | ✓ | ✓ | ✓ | ✓ |
| `lnudgevel` | ✓ | ✓ | ✓ | ✓ |
| `lprofforc` | ✓ | ✓ | ✓ | ✓ |
| `ltempeq` | ✓ | ✓ | ✓ | ✓ |
| `ltimedeplw` | ✓ | ✓ | ✓ | ✓ |
| `ltimedepnudge` | ✓ | ✓ | ✓ | ✓ |
| `ltimedepsurf` | ✓ | ✓ | ✓ | ✓ |
| `ltimedepsw` | ✓ | ✓ | ✓ | ✓ |
| `luoutflowr` | ✓ | ✓ | ✓ | ✓ |
| `luvolflowr` | ✓ | ✓ | ✓ | ✓ |
| `lvinf` | ✓ | ✓ | ✓ | ✓ |
| `lvoutflowr` | ✓ | ✓ | ✓ | ✓ |
| `lvvolflowr` | ✓ | ✓ | ✓ | ✓ |
| `nnudge` | ✓ | ✓ | ✓ | ✓ |
| `ntimedeplw` | ✓ | ✓ | ✓ | ✓ |
| `ntimedepnudge` | ✓ | ✓ | ✓ | ✓ |
| `ntimedepsurf` | ✓ | ✓ | ✓ | ✓ |
| `ntimedepsw` | ✓ | ✓ | ✓ | ✓ |
| `ps` | ✓ | ✓ | ✓ | ✓ |
| `tnudge` | ✓ | ✓ | ✓ | ✓ |
| `tscale` | ✓ | ✓ | ✓ | ✓ |
| `uflowrate` | ✓ | ✓ | ✓ | ✓ |
| `vflowrate` | ✓ | ✓ | ✓ | ✓ |

## PURIFS

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `epu` | ✓ | ✓ | ✓ | ✓ |
| `lpurif` | ✓ | ✓ | ✓ | ✓ |
| `npurif` | ✓ | ✓ | ✓ | ✓ |
| `qpu` | ✓ | ✓ | ✓ | ✓ |

## RUN

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `author` | ✓ | ✓ | ✓ | ✓ |
| `courant` | ✓ | ✓ | ✓ | ✓ |
| `diffnr` | ✓ | ✓ | ✓ | ✓ |
| `dtmax` | ✓ | ✓ | ✓ | ✓ |
| `iexpnr` | ✓ | ✓ | ✓ | ✓ |
| `irandom` | ✓ | ✓ | ✓ | ✓ |
| `krand` | ✓ | ✓ | ✓ | ✓ |
| `ladaptive` | ✓ | ✓ | ✓ | ✓ |
| `libm` | ✓ | ✓ | ✓ | ✓ |
| `ljson_input` | ✓ | ✓ | ✓ | ✓ |
| `lles` | ✓ | ✓ | ✓ | ✓ |
| `lper2inout` | ✓ | ✓ | ✓ | ✓ |
| `lrandomize` | ✓ | ✓ | ✓ | ✓ |
| `lreadmean` | ✓ | ✓ | ✓ | ✓ |
| `lstratstart` | ✓ | ✓ | ✓ | ✓ |
| `lwalldist` | ✓ | ✓ | ✓ | ✓ |
| `lwarmstart` | ✓ | ✓ | ✓ | ✓ |
| `nprocx` | ✓ | ✓ | ✓ | ✓ |
| `nprocy` | ✓ | ✓ | ✓ | ✓ |
| `randqt` | ✓ | ✓ | ✓ | ✓ |
| `randthl` | ✓ | ✓ | ✓ | ✓ |
| `randu` | ✓ | ✓ | ✓ | ✓ |
| `runmode` | ✓ | ✓ | ✓ | ✓ |
| `runtime` | ✓ | ✓ | ✓ | ✓ |
| `startfile` | ✓ | ✓ | ✓ | ✓ |
| `trestart` | ✓ | ✓ | ✓ | ✓ |

## SCALARS

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `lreadscal` | ✓ | ✓ | ✓ | ✓ |
| `lscasrc` | ✓ | ✓ | ✓ | ✓ |
| `lscasrcl` | ✓ | ✓ | ✓ | ✓ |
| `lscasrcr` | ✓ | ✓ | ✓ | ✓ |
| `nscasrc` | ✓ | ✓ | ✓ | ✓ |
| `nscasrcl` | ✓ | ✓ | ✓ | ✓ |
| `nsv` | ✓ | ✓ | ✓ | ✓ |

## TREES

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `cd` | ✓ | ✓ | ✓ | ✓ |
| `dec` | ✓ | ✓ | ✓ | ✓ |
| `dqdt` | ✓ | ✓ | ✓ | ✓ |
| `lad` | ✓ | ✓ | ✓ | ✓ |
| `lsize` | ✓ | ✓ | ✓ | ✓ |
| `ltreedump` | ✓ | ✓ | ✓ | ✓ |
| `ltrees` | ✓ | ✓ | ✓ | ✓ |
| `ntrees` | ✓ | ✓ | ✓ | ✓ |
| `qstar` | ✓ | ✓ | ✓ | ✓ |
| `r_s` | ✓ | ✓ | ✓ | ✓ |
| `ud` | ✓ | ✓ | ✓ | ✓ |

## WALLS

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `dtfac` | ✓ | ✓ | ✓ | ✓ |
| `fkar` | ✓ | ✓ | ✓ | ✓ |
| `iwallmoist` | ✓ | ✓ | ✓ | ✓ |
| `iwallmom` | ✓ | ✓ | ✓ | ✓ |
| `iwallscal` | ✓ | ✓ | ✓ | ✓ |
| `iwalltemp` | ✓ | ✓ | ✓ | ✓ |
| `lbottom` | ✓ | ✓ | ✓ | ✓ |
| `lnorec` | ✓ | ✓ | ✓ | ✓ |
| `lwritefac` | ✓ | ✓ | ✓ | ✓ |
| `nblocks` | ✓ | ✓ | ✓ | ✓ |
| `nbndpts_c` | ✓ | ✓ | ✓ | ✓ |
| `nbndpts_u` | ✓ | ✓ | ✓ | ✓ |
| `nbndpts_v` | ✓ | ✓ | ✓ | ✓ |
| `nbndpts_w` | ✓ | ✓ | ✓ | ✓ |
| `nfcts` | ✓ | ✓ | ✓ | ✓ |
| `nfctsecs_c` | ✓ | ✓ | ✓ | ✓ |
| `nfctsecs_u` | ✓ | ✓ | ✓ | ✓ |
| `nfctsecs_v` | ✓ | ✓ | ✓ | ✓ |
| `nfctsecs_w` | ✓ | ✓ | ✓ | ✓ |
| `nsolpts_c` | ✓ | ✓ | ✓ | ✓ |
| `nsolpts_u` | ✓ | ✓ | ✓ | ✓ |
| `nsolpts_v` | ✓ | ✓ | ✓ | ✓ |
| `nsolpts_w` | ✓ | ✓ | ✓ | ✓ |
| `prandtlturb` | ✓ | ✓ | ✓ | ✓ |

## Summary Statistics

- **Total namelist variables**: 245
- **Total JSON-readable variables**: 245
- **Total broadcast variables**: 245
- **Total schema variables**: 250
- **Total namelists found**: 17

## Potential Issues

### Variables in schema but not in namelists:
- `fact`
- `qt0`
- `thl0`
- `u0`
- `v0`
