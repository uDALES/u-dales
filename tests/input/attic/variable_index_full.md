# uDALES Variable Index Report

This report shows the status of variables across different aspects of the uDALES code:

- **Namelist**: Variable is defined in a Fortran namelist
- **JSON Read**: Variable is read from JSON configuration
- **Broadcast**: Variable is broadcast via MPI
- **Schema**: Variable is present in the JSON schema

## BC Namelist

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
| `wsvsurfdum` | ✓ | ✗ | ✓ | ✓ |
| `wsvtopdum` | ✓ | ✗ | ✓ | ✓ |
| `wtsurf` | ✓ | ✓ | ✓ | ✓ |
| `wttop` | ✓ | ✓ | ✓ | ✓ |
| `z0` | ✓ | ✓ | ✓ | ✓ |
| `z0h` | ✓ | ✓ | ✓ | ✓ |

## CHEMISTRY Namelist

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `jno2` | ✓ | ✗ | ✓ | ✗ |
| `k1` | ✓ | ✗ | ✓ | ✗ |
| `lchem` | ✓ | ✗ | ✓ | ✗ |

## DOMAIN Namelist

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

## DRIVER Namelist

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

## DYNAMICS Namelist

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `iadv_mom` | ✓ | ✓ | ✓ | ✓ |
| `iadv_qt` | ✓ | ✓ | ✓ | ✓ |
| `iadv_sv` | ✓ | ✓ | ✓ | ✓ |
| `iadv_thl` | ✓ | ✓ | ✓ | ✓ |
| `iadv_tke` | ✓ | ✓ | ✓ | ✓ |
| `ipoiss` | ✓ | ✓ | ✓ | ✓ |
| `lqlnr` | ✓ | ✓ | ✓ | ✓ |

## ENERGYBALANCE Namelist

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

## HEATPUMP Namelist

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `lfan_hp` | ✓ | ✗ | ✓ | ✓ |
| `lheatpump` | ✓ | ✗ | ✓ | ✓ |
| `nhppoints` | ✓ | ✗ | ✓ | ✓ |
| `q_dot_hp` | ✓ | ✗ | ✓ | ✓ |
| `qh_dot_hp` | ✓ | ✗ | ✓ | ✓ |

## INFO Namelist

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `dtin` | ✓ | ✗ | ✓ | ✗ |
| `inl` | ✓ | ✗ | ✗ | ✗ |
| `jgtotinl` | ✓ | ✗ | ✓ | ✗ |
| `kmaxin` | ✓ | ✗ | ✗ | ✗ |
| `namezinlet` | ✓ | ✗ | ✗ | ✗ |
| `nprocsinl` | ✓ | ✗ | ✓ | ✗ |
| `totalreadu` | ✓ | ✗ | ✓ | ✗ |
| `wtop` | ✓ | ✗ | ✓ | ✗ |
| `zgrid` | ✓ | ✗ | ✗ | ✗ |

## INLET Namelist

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

## INPS Namelist

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `fact` | ✗ | ✓ | ✓ | ✓ |
| `qt0` | ✗ | ✓ | ✗ | ✓ |
| `thl0` | ✗ | ✓ | ✗ | ✓ |
| `u0` | ✗ | ✓ | ✗ | ✓ |
| `v0` | ✗ | ✓ | ✗ | ✓ |

## NAMCHECKSIM Namelist

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `tcheck` | ✓ | ✗ | ✓ | ✗ |

## NAMSTATSDUMP Namelist

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `khigh` | ✓ | ✗ | ✓ | ✗ |
| `klow` | ✓ | ✗ | ✓ | ✗ |
| `lmintdump` | ✓ | ✓ | ✓ | ✗ |
| `ltdump` | ✓ | ✓ | ✓ | ✗ |
| `ltkedump` | ✓ | ✓ | ✓ | ✗ |
| `ltreedump` | ✓ | ✓ | ✓ | ✗ |
| `lxydump` | ✓ | ✓ | ✓ | ✗ |
| `lxytdump` | ✓ | ✓ | ✓ | ✗ |
| `lydump` | ✓ | ✓ | ✓ | ✗ |
| `lytdump` | ✓ | ✓ | ✓ | ✗ |
| `tsample` | ✓ | ✓ | ✓ | ✗ |
| `tstatsdump` | ✓ | ✓ | ✓ | ✗ |

## NAMSUBGRID Namelist

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `c_vreman` | ✓ | ✗ | ✓ | ✗ |
| `cf` | ✓ | ✗ | ✓ | ✗ |
| `cn` | ✓ | ✓ | ✓ | ✗ |
| `cs` | ✓ | ✓ | ✓ | ✗ |
| `lbuoycorr` | ✓ | ✗ | ✓ | ✗ |
| `ldelta` | ✓ | ✗ | ✓ | ✗ |
| `lmason` | ✓ | ✓ | ✓ | ✗ |
| `loneeqn` | ✓ | ✗ | ✓ | ✗ |
| `lsmagorinsky` | ✓ | ✓ | ✓ | ✗ |
| `lvreman` | ✓ | ✗ | ✓ | ✗ |
| `nmason` | ✓ | ✗ | ✓ | ✗ |
| `prandtl` | ✓ | ✓ | ✓ | ✗ |
| `rigc` | ✓ | ✓ | ✓ | ✗ |

## OUTPUT Namelist

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
| `slicevars` | ✓ | ✗ | ✓ | ✗ |
| `tfielddump` | ✓ | ✓ | ✓ | ✓ |
| `tsample` | ✓ | ✓ | ✓ | ✓ |
| `tstatsdump` | ✓ | ✓ | ✓ | ✓ |
| `tstatstart` | ✓ | ✓ | ✓ | ✓ |

## PHYSICS Namelist

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

## PURIFS Namelist

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `epu` | ✓ | ✓ | ✓ | ✓ |
| `lpurif` | ✓ | ✓ | ✓ | ✓ |
| `npurif` | ✓ | ✓ | ✓ | ✓ |
| `qpu` | ✓ | ✓ | ✓ | ✓ |

## RUN Namelist

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
| `runmode` | ✓ | ✓ | ✓ | ✓ |
| `runtime` | ✓ | ✓ | ✓ | ✓ |
| `startfile` | ✓ | ✓ | ✓ | ✓ |
| `trestart` | ✓ | ✓ | ✓ | ✓ |

## SCALARS Namelist

| Variable | Namelist | JSON Read | Broadcast | Schema |
|----------|----------|-----------|-----------|--------|
| `lreadscal` | ✓ | ✓ | ✓ | ✓ |
| `lscasrc` | ✓ | ✓ | ✓ | ✓ |
| `lscasrcl` | ✓ | ✓ | ✓ | ✓ |
| `lscasrcr` | ✓ | ✓ | ✓ | ✓ |
| `nscasrc` | ✓ | ✓ | ✓ | ✓ |
| `nscasrcl` | ✓ | ✓ | ✓ | ✓ |
| `nsv` | ✓ | ✓ | ✓ | ✓ |

## TREES Namelist

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

## WALLS Namelist

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

- **Total namelist variables**: 259
- **Total JSON-readable variables**: 306
- **Total broadcast variables**: 351
- **Total schema variables**: 225
- **Total namelists found**: 20

## Potential Issues

### Variables in schema but not in namelists:
- `fact`
- `qt0`
- `thl0`
- `u0`
- `v0`

### Variables in namelists but not in schema:
- `c_vreman`
- `cf`
- `cn`
- `cs`
- `dtin`
- `inl`
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
- `namezinlet`
- `nmason`
- `nprocsinl`
- `prandtl`
- `rigc`
- `slicevars`
- `tcheck`
- `totalreadu`
- `tsample`
- `tstatsdump`
- `wtop`
- `zgrid`
