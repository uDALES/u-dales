# uDALES Variable Index Report

This report shows the status of variables across different aspects of the uDALES code:

- **Namelist**: Variable is defined in a Fortran namelist (✓/✗)
- **JSON Read**: Variable is read from JSON configuration (count or ✗)
- **Broadcast**: Variable is broadcast via MPI (count or ✗)
- **Schema**: Variable is present in the JSON schema (✓/✗)

## BC Namelist *(defined in modstartup.f90)*

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

## CHEMISTRY Namelist *(defined in modstartup.f90)*

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `jno2` | ✓ | ✓ | ✓ | ✗ |
| `k1` | ✓ | ✓ | ✓ | ✗ |
| `lchem` | ✓ | ✓ | ✓ | ✗ |

## DOMAIN Namelist *(defined in modstartup.f90)*

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

## DRIVER Namelist *(defined in modstartup.f90)*

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `chunkread_size` | ✓ | ✓ | ✓ | ✓ |
| `driverjobnr` | ✓ | ✓ | ✓ | ✓ |
| `driverstore` | ✓ | ✓ | ✓ | ✓ |
| `dtdriver` | ✓ | ✓ | ✓ | ✓ |
| `iangledeg` | ✓ | ✓ | ✓ | ✓ |
| `idriver` | ✓ | ✓ | ✓ (2x) | ✓ |
| `iplane` | ✓ | ✓ | ✓ | ✓ |
| `lchunkread` | ✓ | ✓ | ✓ | ✓ |
| `tdriverstart` | ✓ | ✓ | ✓ | ✓ |

## DYNAMICS Namelist *(defined in modstartup.f90)*

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `iadv_mom` | ✓ | ✓ | ✓ | ✓ |
| `iadv_qt` | ✓ | ✓ | ✓ | ✓ |
| `iadv_sv` | ✓ | ✓ | ✓ (2x) | ✓ |
| `iadv_thl` | ✓ | ✓ | ✓ | ✓ |
| `iadv_tke` | ✓ | ✓ | ✓ | ✓ |
| `ipoiss` | ✓ | ✓ | ✓ | ✓ |
| `lqlnr` | ✓ | ✓ | ✓ | ✓ |

## ENERGYBALANCE Namelist *(defined in modstartup.f90)*

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

## HEATPUMP Namelist *(defined in modstartup.f90)*

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `lfan_hp` | ✓ | ✓ | ✓ | ✓ |
| `lheatpump` | ✓ | ✓ | ✓ | ✓ |
| `nhppoints` | ✓ | ✓ | ✓ | ✓ |
| `q_dot_hp` | ✓ | ✓ | ✓ | ✓ |
| `qh_dot_hp` | ✓ | ✓ | ✓ | ✓ |

## INFO Namelist *(defined in modinlet.f90)*

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `dtin` | ✓ | ✓ | ✓ (2x) | ✗ |
| `jgtotinl` | ✓ | ✓ | ✓ (2x) | ✗ |
| `kmaxin` | ✓ | ✓ | ✓ | ✗ |
| `nprocsinl` | ✓ | ✓ | ✓ (2x) | ✗ |
| `totalreadu` | ✓ | ✓ | ✓ (2x) | ✗ |
| `wtop` | ✓ | ✓ | ✓ (2x) | ✗ |

## INLET Namelist *(defined in modstartup.f90)*

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

## NAMCHECKSIM Namelist *(defined in modchecksim.f90)*

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `tcheck` | ✓ | ✓ | ✓ | ✗ |

## NAMSTATSDUMP Namelist *(defined in modstatsdump.f90)*

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `anymore` | ✓ | ✗ | ✗ | ✗ |
| `is` | ✓ | ✗ | ✗ | ✗ |
| `khigh` | ✓ | ✓ | ✓ (3x) | ✗ |
| `klow` | ✓ | ✓ | ✓ (3x) | ✗ |
| `lmintdump` | ✓ | ✓ | ✓ (2x) | ✗ |
| `ltdump` | ✓ | ✓ | ✓ (2x) | ✗ |
| `ltkedump` | ✓ | ✓ | ✓ | ✗ |
| `ltreedump` | ✓ | ✓ | ✓ (2x) | ✗ |
| `lxydump` | ✓ | ✓ | ✓ | ✗ |
| `lxytdump` | ✓ | ✓ | ✓ | ✗ |
| `lydump` | ✓ | ✓ | ✓ | ✗ |
| `lytdump` | ✓ | ✓ | ✓ | ✗ |
| `maybe` | ✓ | ✗ | ✗ | ✗ |
| `namstatsdump` | ✓ | ✗ | ✗ | ✗ |
| `removed` | ✓ | ✗ | ✗ | ✗ |
| `tsample` | ✓ | ✓ | ✓ | ✗ |
| `tstatsdump` | ✓ | ✓ | ✓ | ✗ |

## NAMSUBGRID Namelist *(defined in modsubgrid.f90)*

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `c_vreman` | ✓ | ✓ | ✓ (2x) | ✗ |
| `cf` | ✓ | ✓ | ✓ (2x) | ✗ |
| `cn` | ✓ | ✓ | ✓ (2x) | ✗ |
| `cs` | ✓ | ✓ | ✓ (2x) | ✗ |
| `lbuoycorr` | ✓ | ✓ | ✓ (2x) | ✗ |
| `ldelta` | ✓ | ✓ | ✓ (2x) | ✗ |
| `lmason` | ✓ | ✓ | ✓ (2x) | ✗ |
| `loneeqn` | ✓ | ✓ | ✓ (2x) | ✗ |
| `lsmagorinsky` | ✓ | ✓ | ✓ (2x) | ✗ |
| `lvreman` | ✓ | ✓ | ✓ (2x) | ✗ |
| `nmason` | ✓ | ✓ | ✓ (2x) | ✗ |
| `prandtl` | ✓ | ✓ | ✓ (2x) | ✗ |
| `rigc` | ✓ | ✓ | ✓ (2x) | ✗ |

## OUTPUT Namelist *(defined in modstartup.f90)*

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `fieldvars` | ✓ | ✓ | ✓ | ✓ |
| `islice` | ✓ | ✓ | ✓ | ✓ |
| `jslice` | ✓ | ✓ | ✓ | ✓ |
| `kslice` | ✓ | ✓ | ✓ | ✓ |
| `lfielddump` | ✓ | ✓ | ✓ (2x) | ✓ |
| `lislicedump` | ✓ | ✓ | ✓ | ✓ |
| `ljslicedump` | ✓ | ✓ | ✓ | ✓ |
| `lkslicedump` | ✓ | ✓ | ✓ | ✓ |
| `lmintdump` | ✓ | ✓ | ✓ (2x) | ✓ |
| `ltdump` | ✓ | ✓ | ✓ (2x) | ✓ |
| `ltkedump` | ✓ | ✓ | ✓ | ✓ |
| `lxydump` | ✓ | ✓ | ✓ | ✓ |
| `lxytdump` | ✓ | ✓ | ✓ | ✓ |
| `lydump` | ✓ | ✓ | ✓ | ✓ |
| `lytdump` | ✓ | ✓ | ✓ | ✓ |
| `slicevars` | ✓ | ✗ | ✗ | ✗ |
| `tfielddump` | ✓ | ✓ | ✓ | ✓ |
| `tsample` | ✓ | ✓ | ✓ | ✓ |
| `tstatsdump` | ✓ | ✓ | ✓ | ✓ |
| `tstatstart` | ✓ | ✓ | ✓ | ✓ |

## PHYSICS Namelist *(defined in modstartup.f90)*

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
| `ltimedeplw` | ✓ | ✓ | ✓ (2x) | ✓ |
| `ltimedepnudge` | ✓ | ✓ | ✓ (2x) | ✓ |
| `ltimedepsurf` | ✓ | ✓ | ✓ (2x) | ✓ |
| `ltimedepsw` | ✓ | ✓ | ✓ (2x) | ✓ |
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

## PURIFS Namelist *(defined in modstartup.f90)*

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `epu` | ✓ | ✓ | ✓ | ✓ |
| `lpurif` | ✓ | ✓ | ✓ | ✓ |
| `npurif` | ✓ | ✓ | ✓ (2x) | ✓ |
| `qpu` | ✓ | ✓ | ✓ | ✓ |

## RUN Namelist *(defined in modstartup.f90)*

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
| `runmode` | ✓ | ✓ | ✓ | ✗ |
| `runtime` | ✓ | ✓ | ✓ | ✓ |
| `startfile` | ✓ | ✓ | ✓ | ✓ |
| `trestart` | ✓ | ✓ | ✓ | ✓ |

## SCALARS Namelist *(defined in modstartup.f90)*

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `lreadscal` | ✓ | ✓ | ✓ | ✓ |
| `lscasrc` | ✓ | ✓ | ✓ | ✓ |
| `lscasrcl` | ✓ | ✓ | ✓ | ✓ |
| `lscasrcr` | ✓ | ✓ | ✓ | ✓ |
| `nscasrc` | ✓ | ✓ | ✓ | ✓ |
| `nscasrcl` | ✓ | ✓ | ✓ | ✓ |
| `nsv` | ✓ | ✓ | ✓ | ✓ |

## TREES Namelist *(defined in modstartup.f90)*

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `cd` | ✓ | ✓ | ✓ | ✓ |
| `dec` | ✓ | ✓ | ✓ | ✓ |
| `dqdt` | ✓ | ✓ | ✓ | ✓ |
| `lad` | ✓ | ✓ | ✓ | ✓ |
| `lsize` | ✓ | ✓ | ✓ | ✓ |
| `ltreedump` | ✓ | ✓ | ✓ (2x) | ✓ |
| `ltrees` | ✓ | ✓ | ✓ | ✓ |
| `ntrees` | ✓ | ✓ | ✓ (2x) | ✓ |
| `qstar` | ✓ | ✓ | ✓ | ✓ |
| `r_s` | ✓ | ✓ | ✓ | ✓ |
| `ud` | ✓ | ✓ | ✓ | ✓ |

## WALLS Namelist *(defined in modstartup.f90)*

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

- **Total namelist variables**: 261
- **Total JSON-readable variables**: 245
- **Total broadcast variables**: 352
- **Total schema variables**: 224
- **Total namelists found**: 19

## Potential Issues

### Variables in schema but not in namelists:
- `fact`
- `qt0`
- `thl0`
- `u0`
- `v0`

### Variables in namelists but not in schema:
- `anymore`
- `c_vreman`
- `cf`
- `cn`
- `cs`
- `dtin`
- `is`
- `jgtotinl`
- `jno2`
- `k1`
- `khigh`
- `klow`
- `kmaxin`
- `lbuoycorr`
- `lchem`
- `ldelta`
- `lmason`
- `lmintdump`
- `loneeqn`
- `lsmagorinsky`
- `ltdump`
- `ltkedump`
- `ltreedump`
- `lvreman`
- `lxydump`
- `lxytdump`
- `lydump`
- `lytdump`
- `maybe`
- `namstatsdump`
- `nmason`
- `nprocsinl`
- `prandtl`
- `removed`
- `rigc`
- `runmode`
- `slicevars`
- `tcheck`
- `totalreadu`
- `tsample`
- `tstatsdump`
- `wtop`
