# Namelist variable report

Sections found: bc, chemistry, domain, driver, dynamics, energybalance, heatpump, info, inlet, namchecksim, namstatsdump, namsubgrid, output, physics, purifs, run, scalars, trees, walls

## bc 
Read from JSON: bcbotm, bcbotq, bcbots, bcbott, bcqfxm, bcqfxp, bcqfym, bcqfyp, bcqfz, bctfxm, bctfxp, bctfym, bctfyp, bctfz, bctopm, bctopq, bctops, bctopt, bcxm, bcxq, bcxs, bcxt, bcym, bcyq, bcys, bcyt, bczp, ds, qt_top, qts, thl_top, thls, wqsurf, wsvsurfdum, wsvtopdum, wtsurf, wttop, z0, z0h

Broadcast via MPI_BCAST: bcbotm, bcbotq, bcbots, bcbott, bcqfxm, bcqfxp, bcqfym, bcqfyp, bcqfz, bctfxm, bctfxp, bctfym, bctfyp, bctfz, bctopm, bctopq, bctops, bctopt, bcxm, bcxq, bcxs, bcxt, bcym, bcyq, bcys, bcyt, bczp, ds, qt_top, qts, thl_top, thls, wqsurf, wsvsurfdum, wsvtopdum, wtsurf, wttop, z0, z0h
## chemistry 
Read from JSON: jno2, k1, lchem

Broadcast via MPI_BCAST: jno2, k1, lchem
## domain 
Read from JSON: itot, jtot, ksp, ktot, xday, xlat, xlen, xlon, xtime, ylen

Broadcast via MPI_BCAST: itot, jtot, ksp, ktot, xday, xlat, xlen, xlon, xtime, ylen
## driver 
Read from JSON: chunkread_size, driverjobnr, driverstore, dtdriver, iangledeg, idriver, iplane, lchunkread, tdriverstart

Broadcast via MPI_BCAST: chunkread_size, driverjobnr, driverstore, dtdriver, iangledeg, idriver, iplane, lchunkread, tdriverstart
## dynamics 
Read from JSON: iadv_mom, iadv_qt, iadv_sv, iadv_thl, iadv_tke, ipoiss, lqlnr

Broadcast via MPI_BCAST: iadv_mom, iadv_qt, iadv_sv, iadv_thl, iadv_tke, ipoiss, lqlnr
## energybalance 
Read from JSON: bldt, dteb, flrt, fraction, grlai, lconstw, leb, lfactlyrs, lperiodicebcorr, lvfsparse, lwriteebfiles, nfaclyrs, nnz, rsmin, sinkbase, skylw, wfc, wgrmax, wsoil, wwilt

Broadcast via MPI_BCAST: bldt, dteb, flrt, fraction, grlai, lconstw, leb, lfactlyrs, lperiodicebcorr, lvfsparse, lwriteebfiles, nfaclyrs, nnz, rsmin, sinkbase, skylw, wfc, wgrmax, wsoil, wwilt
## heatpump 
Read from JSON: lfan_hp, lheatpump, nhppoints, q_dot_hp, qh_dot_hp

Broadcast via MPI_BCAST: lfan_hp, lheatpump, nhppoints, q_dot_hp, qh_dot_hp
## info 
Read from JSON: dtin, jgtotinl, kmaxin, nprocsinl, totalreadu, wtop

Broadcast via MPI_BCAST: dtin, jgtotinl, kmaxin, nprocsinl, totalreadu, wtop
## inlet 
Read from JSON: di, dti, inletav, lfixinlet, lfixutauin, linletra, lreadminl, lstoreplane, lwallfunc, uinf, vinf

Broadcast via MPI_BCAST: di, dti, inletav, lfixinlet, lfixutauin, linletra, lreadminl, lstoreplane, lwallfunc, uinf, vinf
## namchecksim 
Read from JSON: tcheck

Broadcast via MPI_BCAST: tcheck
## namstatsdump 
Read from JSON: khigh, klow

Broadcast via MPI_BCAST: khigh, klow
## namsubgrid 
Read from JSON: c_vreman, cf, cn, cs, lbuoycorr, ldelta, lmason, loneeqn, lsmagorinsky, lvreman, nmason, prandtl, rigc

Broadcast via MPI_BCAST: c_vreman, cf, cn, cs, lbuoycorr, ldelta, lmason, loneeqn, lsmagorinsky, lvreman, nmason, prandtl, rigc
## output 
Read from JSON: fieldvars, islice, jslice, kslice, lfielddump, lislicedump, ljslicedump, lkslicedump, lmintdump, ltdump, ltkedump, lxydump, lxytdump, lydump, lytdump, tfielddump, tsample, tstatsdump, tstatstart

Broadcast via MPI_BCAST: fieldvars, islice, jslice, kslice, lfielddump, lislicedump, ljslicedump, lkslicedump, lmintdump, ltdump, ltkedump, lxydump, lxytdump, lydump, lytdump, tfielddump, tsample, tstatsdump, tstatstart
## physics 
Read from JSON: dpdx, ifixuinf, igrw_damp, lbuoyancy, lcoriol, lmoist, lnudge, lnudgevel, lprofforc, ltempeq, ltimedeplw, ltimedepnudge, ltimedepsurf, ltimedepsw, luoutflowr, luvolflowr, lvinf, lvoutflowr, lvvolflowr, nnudge, ntimedeplw, ntimedepnudge, ntimedepsurf, ntimedepsw, ps, tnudge, tscale, uflowrate, vflowrate

Broadcast via MPI_BCAST: dpdx, ifixuinf, igrw_damp, lbuoyancy, lcoriol, lmoist, lnudge, lnudgevel, lprofforc, ltempeq, ltimedeplw, ltimedepnudge, ltimedepsurf, ltimedepsw, luoutflowr, luvolflowr, lvinf, lvoutflowr, lvvolflowr, nnudge, ntimedeplw, ntimedepnudge, ntimedepsurf, ntimedepsw, ps, tnudge, tscale, uflowrate, vflowrate
## purifs 
Read from JSON: epu, lpurif, npurif, qpu

Broadcast via MPI_BCAST: epu, lpurif, npurif, qpu
## run 
Read from JSON: author, courant, diffnr, dtmax, iexpnr, irandom, krand, ladaptive, libm, lles, lper2inout, lrandomize, lreadmean, lstratstart, lwalldist, lwarmstart, nprocx, nprocy, randqt, randthl, randu, runmode, runtime, startfile, trestart

Broadcast via MPI_BCAST: author, courant, diffnr, dtmax, iexpnr, irandom, krand, ladaptive, libm, lles, lper2inout, lrandomize, lreadmean, lstratstart, lwalldist, lwarmstart, nprocx, nprocy, randqt, randthl, randu, runmode, runtime, startfile, trestart
## scalars 
Read from JSON: lreadscal, lscasrc, lscasrcl, lscasrcr, nscasrc, nscasrcl, nsv

Broadcast via MPI_BCAST: lreadscal, lscasrc, lscasrcl, lscasrcr, nscasrc, nscasrcl, nsv
## trees 
Read from JSON: cd, dec, dqdt, lad, lsize, ltreedump, ltrees, ntrees, qstar, r_s, ud

Broadcast via MPI_BCAST: cd, dec, dqdt, lad, lsize, ltreedump, ltrees, ntrees, qstar, r_s, ud
## walls 
Read from JSON: dtfac, fkar, iwallmoist, iwallmom, iwallscal, iwalltemp, lbottom, lnorec, lwritefac, nblocks, nbndpts_c, nbndpts_u, nbndpts_v, nbndpts_w, nfcts, nfctsecs_c, nfctsecs_u, nfctsecs_v, nfctsecs_w, nsolpts_c, nsolpts_u, nsolpts_v, nsolpts_w, prandtlturb

Broadcast via MPI_BCAST: dtfac, fkar, iwallmoist, iwallmom, iwallscal, iwalltemp, lbottom, lnorec, lwritefac, nblocks, nbndpts_c, nbndpts_u, nbndpts_v, nbndpts_w, nfcts, nfctsecs_c, nfctsecs_u, nfctsecs_v, nfctsecs_w, nsolpts_c, nsolpts_u, nsolpts_v, nsolpts_w, prandtlturb

# Comparison with canonical namelist index

## bc

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `bcbotm` | ✓ | ✓ | ✓ |
| `bcbotq` | ✓ | ✓ | ✓ |
| `bcbots` | ✓ | ✓ | ✓ |
| `bcbott` | ✓ | ✓ | ✓ |
| `bcqfxm` | ✓ | ✓ | ✓ |
| `bcqfxp` | ✓ | ✓ | ✓ |
| `bcqfym` | ✓ | ✓ | ✓ |
| `bcqfyp` | ✓ | ✓ | ✓ |
| `bcqfz` | ✓ | ✓ | ✓ |
| `bctfxm` | ✓ | ✓ | ✓ |
| `bctfxp` | ✓ | ✓ | ✓ |
| `bctfym` | ✓ | ✓ | ✓ |
| `bctfyp` | ✓ | ✓ | ✓ |
| `bctfz` | ✓ | ✓ | ✓ |
| `bctopm` | ✓ | ✓ | ✓ |
| `bctopq` | ✓ | ✓ | ✓ |
| `bctops` | ✓ | ✓ | ✓ |
| `bctopt` | ✓ | ✓ | ✓ |
| `bcxm` | ✓ | ✓ | ✓ |
| `bcxq` | ✓ | ✓ | ✓ |
| `bcxs` | ✓ | ✓ | ✓ |
| `bcxt` | ✓ | ✓ | ✓ |
| `bcym` | ✓ | ✓ | ✓ |
| `bcyq` | ✓ | ✓ | ✓ |
| `bcys` | ✓ | ✓ | ✓ |
| `bcyt` | ✓ | ✓ | ✓ |
| `bczp` | ✓ | ✓ | ✓ |
| `ds` | ✓ | ✓ | ✓ |
| `qt_top` | ✓ | ✓ | ✓ |
| `qts` | ✓ | ✓ | ✓ |
| `thl_top` | ✓ | ✓ | ✓ |
| `thls` | ✓ | ✓ | ✓ |
| `wqsurf` | ✓ | ✓ | ✓ |
| `wsvsurfdum` | ✓ | ✓ | ✓ |
| `wsvtopdum` | ✓ | ✓ | ✓ |
| `wtsurf` | ✓ | ✓ | ✓ |
| `wttop` | ✓ | ✓ | ✓ |
| `z0` | ✓ | ✓ | ✓ |
| `z0h` | ✓ | ✓ | ✓ |

- Canonical variables in index: 39
- JSON reads detected: 39 (including 39 in index)
- MPI_BCAST detected: 39 (including 39 in index)

## chemistry

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `jno2` | ✓ | ✓ | ✓ |
| `k1` | ✓ | ✓ | ✓ |
| `lchem` | ✓ | ✓ | ✓ |

- Canonical variables in index: 3
- JSON reads detected: 3 (including 3 in index)
- MPI_BCAST detected: 3 (including 3 in index)

## domain

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `itot` | ✓ | ✓ | ✓ |
| `jtot` | ✓ | ✓ | ✓ |
| `ksp` | ✓ | ✓ | ✓ |
| `ktot` | ✓ | ✓ | ✓ |
| `xday` | ✓ | ✓ | ✓ |
| `xlat` | ✓ | ✓ | ✓ |
| `xlen` | ✓ | ✓ | ✓ |
| `xlon` | ✓ | ✓ | ✓ |
| `xtime` | ✓ | ✓ | ✓ |
| `ylen` | ✓ | ✓ | ✓ |

- Canonical variables in index: 10
- JSON reads detected: 10 (including 10 in index)
- MPI_BCAST detected: 10 (including 10 in index)

## driver

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `chunkread_size` | ✓ | ✓ | ✓ |
| `driverjobnr` | ✓ | ✓ | ✓ |
| `driverstore` | ✓ | ✓ | ✓ |
| `dtdriver` | ✓ | ✓ | ✓ |
| `iangledeg` | ✓ | ✓ | ✓ |
| `idriver` | ✓ | ✓ | ✓ |
| `iplane` | ✓ | ✓ | ✓ |
| `lchunkread` | ✓ | ✓ | ✓ |
| `tdriverstart` | ✓ | ✓ | ✓ |

- Canonical variables in index: 9
- JSON reads detected: 9 (including 9 in index)
- MPI_BCAST detected: 9 (including 9 in index)

## dynamics

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `iadv_mom` | ✓ | ✓ | ✓ |
| `iadv_qt` | ✓ | ✓ | ✓ |
| `iadv_sv` | ✓ | ✓ | ✓ |
| `iadv_thl` | ✓ | ✓ | ✓ |
| `iadv_tke` | ✓ | ✓ | ✓ |
| `ipoiss` | ✓ | ✓ | ✓ |
| `lqlnr` | ✓ | ✓ | ✓ |

- Canonical variables in index: 7
- JSON reads detected: 7 (including 7 in index)
- MPI_BCAST detected: 7 (including 7 in index)

## energybalance

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `bldt` | ✓ | ✓ | ✓ |
| `dteb` | ✓ | ✓ | ✓ |
| `flrt` | ✓ | ✓ | ✓ |
| `fraction` | ✓ | ✓ | ✓ |
| `grlai` | ✓ | ✓ | ✓ |
| `lconstw` | ✓ | ✓ | ✓ |
| `leb` | ✓ | ✓ | ✓ |
| `lfactlyrs` | ✓ | ✓ | ✓ |
| `lperiodicebcorr` | ✓ | ✓ | ✓ |
| `lvfsparse` | ✓ | ✓ | ✓ |
| `lwriteebfiles` | ✓ | ✓ | ✓ |
| `nfaclyrs` | ✓ | ✓ | ✓ |
| `nnz` | ✓ | ✓ | ✓ |
| `rsmin` | ✓ | ✓ | ✓ |
| `sinkbase` | ✓ | ✓ | ✓ |
| `skylw` | ✓ | ✓ | ✓ |
| `wfc` | ✓ | ✓ | ✓ |
| `wgrmax` | ✓ | ✓ | ✓ |
| `wsoil` | ✓ | ✓ | ✓ |
| `wwilt` | ✓ | ✓ | ✓ |

- Canonical variables in index: 20
- JSON reads detected: 20 (including 20 in index)
- MPI_BCAST detected: 20 (including 20 in index)

## heatpump

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `lfan_hp` | ✓ | ✓ | ✓ |
| `lheatpump` | ✓ | ✓ | ✓ |
| `nhppoints` | ✓ | ✓ | ✓ |
| `q_dot_hp` | ✓ | ✓ | ✓ |
| `qh_dot_hp` | ✓ | ✓ | ✓ |

- Canonical variables in index: 5
- JSON reads detected: 5 (including 5 in index)
- MPI_BCAST detected: 5 (including 5 in index)

## info

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `dtin` | ✓ | ✓ | ✓ |
| `jgtotinl` | ✓ | ✓ | ✓ |
| `kmaxin` | ✓ | ✓ | ✓ |
| `nprocsinl` | ✓ | ✓ | ✓ |
| `totalreadu` | ✓ | ✓ | ✓ |
| `wtop` | ✓ | ✓ | ✓ |

- Canonical variables in index: 6
- JSON reads detected: 6 (including 6 in index)
- MPI_BCAST detected: 6 (including 6 in index)

## inlet

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `di` | ✓ | ✓ | ✓ |
| `dti` | ✓ | ✓ | ✓ |
| `inletav` | ✓ | ✓ | ✓ |
| `lfixinlet` | ✓ | ✓ | ✓ |
| `lfixutauin` | ✓ | ✓ | ✓ |
| `linletra` | ✓ | ✓ | ✓ |
| `lreadminl` | ✓ | ✓ | ✓ |
| `lstoreplane` | ✓ | ✓ | ✓ |
| `lwallfunc` | ✓ | ✓ | ✓ |
| `uinf` | ✓ | ✓ | ✓ |
| `vinf` | ✓ | ✓ | ✓ |

- Canonical variables in index: 11
- JSON reads detected: 11 (including 11 in index)
- MPI_BCAST detected: 11 (including 11 in index)

## namchecksim

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `tcheck` | ✓ | ✓ | ✓ |

- Canonical variables in index: 1
- JSON reads detected: 1 (including 1 in index)
- MPI_BCAST detected: 1 (including 1 in index)

## namstatsdump

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `anymore` | ✓ |  |  |
| `is` | ✓ |  |  |
| `khigh` | ✓ | ✓ | ✓ |
| `klow` | ✓ | ✓ | ✓ |
| `lmintdump` | ✓ |  |  |
| `ltdump` | ✓ |  |  |
| `ltkedump` | ✓ |  |  |
| `ltreedump` | ✓ |  |  |
| `lxydump` | ✓ |  |  |
| `lxytdump` | ✓ |  |  |
| `lydump` | ✓ |  |  |
| `lytdump` | ✓ |  |  |
| `maybe` | ✓ |  |  |
| `namstatsdump` | ✓ |  |  |
| `removed` | ✓ |  |  |
| `tsample` | ✓ |  |  |
| `tstatsdump` | ✓ |  |  |

- Canonical variables in index: 17
- JSON reads detected: 2 (including 2 in index)
- MPI_BCAST detected: 2 (including 2 in index)

## namsubgrid

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `c_vreman` | ✓ | ✓ | ✓ |
| `cf` | ✓ | ✓ | ✓ |
| `cn` | ✓ | ✓ | ✓ |
| `cs` | ✓ | ✓ | ✓ |
| `lbuoycorr` | ✓ | ✓ | ✓ |
| `ldelta` | ✓ | ✓ | ✓ |
| `lmason` | ✓ | ✓ | ✓ |
| `loneeqn` | ✓ | ✓ | ✓ |
| `lsmagorinsky` | ✓ | ✓ | ✓ |
| `lvreman` | ✓ | ✓ | ✓ |
| `nmason` | ✓ | ✓ | ✓ |
| `prandtl` | ✓ | ✓ | ✓ |
| `rigc` | ✓ | ✓ | ✓ |

- Canonical variables in index: 13
- JSON reads detected: 13 (including 13 in index)
- MPI_BCAST detected: 13 (including 13 in index)

## output

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `fieldvars` | ✓ | ✓ | ✓ |
| `islice` | ✓ | ✓ | ✓ |
| `jslice` | ✓ | ✓ | ✓ |
| `kslice` | ✓ | ✓ | ✓ |
| `lfielddump` | ✓ | ✓ | ✓ |
| `lislicedump` | ✓ | ✓ | ✓ |
| `ljslicedump` | ✓ | ✓ | ✓ |
| `lkslicedump` | ✓ | ✓ | ✓ |
| `lmintdump` | ✓ | ✓ | ✓ |
| `ltdump` | ✓ | ✓ | ✓ |
| `ltkedump` | ✓ | ✓ | ✓ |
| `lxydump` | ✓ | ✓ | ✓ |
| `lxytdump` | ✓ | ✓ | ✓ |
| `lydump` | ✓ | ✓ | ✓ |
| `lytdump` | ✓ | ✓ | ✓ |
| `slicevars` | ✓ |  |  |
| `tfielddump` | ✓ | ✓ | ✓ |
| `tsample` | ✓ | ✓ | ✓ |
| `tstatsdump` | ✓ | ✓ | ✓ |
| `tstatstart` | ✓ | ✓ | ✓ |

- Canonical variables in index: 20
- JSON reads detected: 19 (including 19 in index)
- MPI_BCAST detected: 19 (including 19 in index)

## physics

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `dpdx` | ✓ | ✓ | ✓ |
| `ifixuinf` | ✓ | ✓ | ✓ |
| `igrw_damp` | ✓ | ✓ | ✓ |
| `lbuoyancy` | ✓ | ✓ | ✓ |
| `lcoriol` | ✓ | ✓ | ✓ |
| `lmoist` | ✓ | ✓ | ✓ |
| `lnudge` | ✓ | ✓ | ✓ |
| `lnudgevel` | ✓ | ✓ | ✓ |
| `lprofforc` | ✓ | ✓ | ✓ |
| `ltempeq` | ✓ | ✓ | ✓ |
| `ltimedeplw` | ✓ | ✓ | ✓ |
| `ltimedepnudge` | ✓ | ✓ | ✓ |
| `ltimedepsurf` | ✓ | ✓ | ✓ |
| `ltimedepsw` | ✓ | ✓ | ✓ |
| `luoutflowr` | ✓ | ✓ | ✓ |
| `luvolflowr` | ✓ | ✓ | ✓ |
| `lvinf` | ✓ | ✓ | ✓ |
| `lvoutflowr` | ✓ | ✓ | ✓ |
| `lvvolflowr` | ✓ | ✓ | ✓ |
| `nnudge` | ✓ | ✓ | ✓ |
| `ntimedeplw` | ✓ | ✓ | ✓ |
| `ntimedepnudge` | ✓ | ✓ | ✓ |
| `ntimedepsurf` | ✓ | ✓ | ✓ |
| `ntimedepsw` | ✓ | ✓ | ✓ |
| `ps` | ✓ | ✓ | ✓ |
| `tnudge` | ✓ | ✓ | ✓ |
| `tscale` | ✓ | ✓ | ✓ |
| `uflowrate` | ✓ | ✓ | ✓ |
| `vflowrate` | ✓ | ✓ | ✓ |

- Canonical variables in index: 29
- JSON reads detected: 29 (including 29 in index)
- MPI_BCAST detected: 29 (including 29 in index)

## purifs

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `epu` | ✓ | ✓ | ✓ |
| `lpurif` | ✓ | ✓ | ✓ |
| `npurif` | ✓ | ✓ | ✓ |
| `qpu` | ✓ | ✓ | ✓ |

- Canonical variables in index: 4
- JSON reads detected: 4 (including 4 in index)
- MPI_BCAST detected: 4 (including 4 in index)

## run

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `author` | ✓ | ✓ | ✓ |
| `courant` | ✓ | ✓ | ✓ |
| `diffnr` | ✓ | ✓ | ✓ |
| `dtmax` | ✓ | ✓ | ✓ |
| `iexpnr` | ✓ | ✓ | ✓ |
| `irandom` | ✓ | ✓ | ✓ |
| `krand` | ✓ | ✓ | ✓ |
| `ladaptive` | ✓ | ✓ | ✓ |
| `libm` | ✓ | ✓ | ✓ |
| `lles` | ✓ | ✓ | ✓ |
| `lper2inout` | ✓ | ✓ | ✓ |
| `lrandomize` | ✓ | ✓ | ✓ |
| `lreadmean` | ✓ | ✓ | ✓ |
| `lstratstart` | ✓ | ✓ | ✓ |
| `lwalldist` | ✓ | ✓ | ✓ |
| `lwarmstart` | ✓ | ✓ | ✓ |
| `nprocx` | ✓ | ✓ | ✓ |
| `nprocy` | ✓ | ✓ | ✓ |
| `randqt` | ✓ | ✓ | ✓ |
| `randthl` | ✓ | ✓ | ✓ |
| `randu` | ✓ | ✓ | ✓ |
| `runmode` | ✓ | ✓ | ✓ |
| `runtime` | ✓ | ✓ | ✓ |
| `startfile` | ✓ | ✓ | ✓ |
| `trestart` | ✓ | ✓ | ✓ |

- Canonical variables in index: 25
- JSON reads detected: 25 (including 25 in index)
- MPI_BCAST detected: 25 (including 25 in index)

## scalars

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `lreadscal` | ✓ | ✓ | ✓ |
| `lscasrc` | ✓ | ✓ | ✓ |
| `lscasrcl` | ✓ | ✓ | ✓ |
| `lscasrcr` | ✓ | ✓ | ✓ |
| `nscasrc` | ✓ | ✓ | ✓ |
| `nscasrcl` | ✓ | ✓ | ✓ |
| `nsv` | ✓ | ✓ | ✓ |

- Canonical variables in index: 7
- JSON reads detected: 7 (including 7 in index)
- MPI_BCAST detected: 7 (including 7 in index)

## trees

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `cd` | ✓ | ✓ | ✓ |
| `dec` | ✓ | ✓ | ✓ |
| `dqdt` | ✓ | ✓ | ✓ |
| `lad` | ✓ | ✓ | ✓ |
| `lsize` | ✓ | ✓ | ✓ |
| `ltreedump` | ✓ | ✓ | ✓ |
| `ltrees` | ✓ | ✓ | ✓ |
| `ntrees` | ✓ | ✓ | ✓ |
| `qstar` | ✓ | ✓ | ✓ |
| `r_s` | ✓ | ✓ | ✓ |
| `ud` | ✓ | ✓ | ✓ |

- Canonical variables in index: 11
- JSON reads detected: 11 (including 11 in index)
- MPI_BCAST detected: 11 (including 11 in index)

## walls

| Variable | In index | JSON read | MPI_BCAST |
|---|---:|:---:|:---:|
| `dtfac` | ✓ | ✓ | ✓ |
| `fkar` | ✓ | ✓ | ✓ |
| `iwallmoist` | ✓ | ✓ | ✓ |
| `iwallmom` | ✓ | ✓ | ✓ |
| `iwallscal` | ✓ | ✓ | ✓ |
| `iwalltemp` | ✓ | ✓ | ✓ |
| `lbottom` | ✓ | ✓ | ✓ |
| `lnorec` | ✓ | ✓ | ✓ |
| `lwritefac` | ✓ | ✓ | ✓ |
| `nblocks` | ✓ | ✓ | ✓ |
| `nbndpts_c` | ✓ | ✓ | ✓ |
| `nbndpts_u` | ✓ | ✓ | ✓ |
| `nbndpts_v` | ✓ | ✓ | ✓ |
| `nbndpts_w` | ✓ | ✓ | ✓ |
| `nfcts` | ✓ | ✓ | ✓ |
| `nfctsecs_c` | ✓ | ✓ | ✓ |
| `nfctsecs_u` | ✓ | ✓ | ✓ |
| `nfctsecs_v` | ✓ | ✓ | ✓ |
| `nfctsecs_w` | ✓ | ✓ | ✓ |
| `nsolpts_c` | ✓ | ✓ | ✓ |
| `nsolpts_u` | ✓ | ✓ | ✓ |
| `nsolpts_v` | ✓ | ✓ | ✓ |
| `nsolpts_w` | ✓ | ✓ | ✓ |
| `prandtlturb` | ✓ | ✓ | ✓ |

- Canonical variables in index: 24
- JSON reads detected: 24 (including 24 in index)
- MPI_BCAST detected: 24 (including 24 in index)

