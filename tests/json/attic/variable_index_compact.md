# uDALES Variable Index - Compact Summary

## Overview

- **Namelists found**: 20
- **Total namelist variables**: 259
- **Variables with JSON support**: 306
- **Variables broadcast via MPI**: 351
- **Variables in JSON schema**: 225

## Namelist Summary

| Namelist | Total Vars | JSON Read | Broadcast | In Schema |
|----------|------------|-----------|-----------|-----------|
| BC | 39 | 37 | 39 | 39 |
| CHEMISTRY | 3 | 0 | 3 | 0 |
| DOMAIN | 10 | 10 | 10 | 10 |
| DRIVER | 9 | 9 | 9 | 9 |
| DYNAMICS | 7 | 7 | 7 | 7 |
| ENERGYBALANCE | 20 | 20 | 20 | 20 |
| HEATPUMP | 5 | 0 | 5 | 5 |
| INFO | 9 | 0 | 5 | 0 |
| INLET | 11 | 11 | 11 | 11 |
| INPS | 5 | 5 | 1 | 5 |
| NAMCHECKSIM | 1 | 0 | 1 | 0 |
| NAMSTATSDUMP | 12 | 10 | 12 | 0 |
| NAMSUBGRID | 13 | 6 | 13 | 0 |
| OUTPUT | 20 | 19 | 20 | 19 |
| PHYSICS | 29 | 29 | 29 | 29 |
| PURIFS | 4 | 4 | 4 | 4 |
| RUN | 25 | 25 | 25 | 25 |
| SCALARS | 7 | 7 | 7 | 7 |
| TREES | 11 | 11 | 11 | 11 |
| WALLS | 24 | 24 | 24 | 24 |

## Variables with Full Support
*(Namelist + JSON + Broadcast + Schema)*

- `BC.bcbotm`
- `BC.bcbotq`
- `BC.bcbots`
- `BC.bcbott`
- `BC.bcqfxm`
- `BC.bcqfxp`
- `BC.bcqfym`
- `BC.bcqfyp`
- `BC.bcqfz`
- `BC.bctfxm`
- `BC.bctfxp`
- `BC.bctfym`
- `BC.bctfyp`
- `BC.bctfz`
- `BC.bctopm`
- `BC.bctopq`
- `BC.bctops`
- `BC.bctopt`
- `BC.bcxm`
- `BC.bcxq`
- `BC.bcxs`
- `BC.bcxt`
- `BC.bcym`
- `BC.bcyq`
- `BC.bcys`
- `BC.bcyt`
- `BC.bczp`
- `BC.ds`
- `BC.qt_top`
- `BC.qts`
- `BC.thl_top`
- `BC.thls`
- `BC.wqsurf`
- `BC.wtsurf`
- `BC.wttop`
- `BC.z0`
- `BC.z0h`
- `DOMAIN.itot`
- `DOMAIN.jtot`
- `DOMAIN.ksp`
- `DOMAIN.ktot`
- `DOMAIN.xday`
- `DOMAIN.xlat`
- `DOMAIN.xlen`
- `DOMAIN.xlon`
- `DOMAIN.xtime`
- `DOMAIN.ylen`
- `DRIVER.chunkread_size`
- `DRIVER.driverjobnr`
- `DRIVER.driverstore`
- `DRIVER.dtdriver`
- `DRIVER.iangledeg`
- `DRIVER.idriver`
- `DRIVER.iplane`
- `DRIVER.lchunkread`
- `DRIVER.tdriverstart`
- `DYNAMICS.iadv_mom`
- `DYNAMICS.iadv_qt`
- `DYNAMICS.iadv_sv`
- `DYNAMICS.iadv_thl`
- `DYNAMICS.iadv_tke`
- `DYNAMICS.ipoiss`
- `DYNAMICS.lqlnr`
- `ENERGYBALANCE.bldt`
- `ENERGYBALANCE.dteb`
- `ENERGYBALANCE.flrt`
- `ENERGYBALANCE.fraction`
- `ENERGYBALANCE.grlai`
- `ENERGYBALANCE.lconstw`
- `ENERGYBALANCE.leb`
- `ENERGYBALANCE.lfactlyrs`
- `ENERGYBALANCE.lperiodicebcorr`
- `ENERGYBALANCE.lvfsparse`
- `ENERGYBALANCE.lwriteebfiles`
- `ENERGYBALANCE.nfaclyrs`
- `ENERGYBALANCE.nnz`
- `ENERGYBALANCE.rsmin`
- `ENERGYBALANCE.sinkbase`
- `ENERGYBALANCE.skylw`
- `ENERGYBALANCE.wfc`
- `ENERGYBALANCE.wgrmax`
- `ENERGYBALANCE.wsoil`
- `ENERGYBALANCE.wwilt`
- `INLET.di`
- `INLET.dti`
- `INLET.inletav`
- `INLET.lfixinlet`
- `INLET.lfixutauin`
- `INLET.linletra`
- `INLET.lreadminl`
- `INLET.lstoreplane`
- `INLET.lwallfunc`
- `INLET.uinf`
- `INLET.vinf`
- `OUTPUT.fieldvars`
- `OUTPUT.islice`
- `OUTPUT.jslice`
- `OUTPUT.kslice`
- `OUTPUT.lfielddump`
- `OUTPUT.lislicedump`
- `OUTPUT.ljslicedump`
- `OUTPUT.lkslicedump`
- `OUTPUT.lmintdump`
- `OUTPUT.ltdump`
- `OUTPUT.ltkedump`
- `OUTPUT.lxydump`
- `OUTPUT.lxytdump`
- `OUTPUT.lydump`
- `OUTPUT.lytdump`
- `OUTPUT.tfielddump`
- `OUTPUT.tsample`
- `OUTPUT.tstatsdump`
- `OUTPUT.tstatstart`
- `PHYSICS.dpdx`
- `PHYSICS.ifixuinf`
- `PHYSICS.igrw_damp`
- `PHYSICS.lbuoyancy`
- `PHYSICS.lcoriol`
- `PHYSICS.lmoist`
- `PHYSICS.lnudge`
- `PHYSICS.lnudgevel`
- `PHYSICS.lprofforc`
- `PHYSICS.ltempeq`
- `PHYSICS.ltimedeplw`
- `PHYSICS.ltimedepnudge`
- `PHYSICS.ltimedepsurf`
- `PHYSICS.ltimedepsw`
- `PHYSICS.luoutflowr`
- `PHYSICS.luvolflowr`
- `PHYSICS.lvinf`
- `PHYSICS.lvoutflowr`
- `PHYSICS.lvvolflowr`
- `PHYSICS.nnudge`
- `PHYSICS.ntimedeplw`
- `PHYSICS.ntimedepnudge`
- `PHYSICS.ntimedepsurf`
- `PHYSICS.ntimedepsw`
- `PHYSICS.ps`
- `PHYSICS.tnudge`
- `PHYSICS.tscale`
- `PHYSICS.uflowrate`
- `PHYSICS.vflowrate`
- `PURIFS.epu`
- `PURIFS.lpurif`
- `PURIFS.npurif`
- `PURIFS.qpu`
- `RUN.author`
- `RUN.courant`
- `RUN.diffnr`
- `RUN.dtmax`
- `RUN.iexpnr`
- `RUN.irandom`
- `RUN.krand`
- `RUN.ladaptive`
- `RUN.libm`
- `RUN.lles`
- `RUN.lper2inout`
- `RUN.lrandomize`
- `RUN.lreadmean`
- `RUN.lstratstart`
- `RUN.lwalldist`
- `RUN.lwarmstart`
- `RUN.nprocx`
- `RUN.nprocy`
- `RUN.randqt`
- `RUN.randthl`
- `RUN.randu`
- `RUN.runmode`
- `RUN.runtime`
- `RUN.startfile`
- `RUN.trestart`
- `SCALARS.lreadscal`
- `SCALARS.lscasrc`
- `SCALARS.lscasrcl`
- `SCALARS.lscasrcr`
- `SCALARS.nscasrc`
- `SCALARS.nscasrcl`
- `SCALARS.nsv`
- `TREES.cd`
- `TREES.dec`
- `TREES.dqdt`
- `TREES.lad`
- `TREES.lsize`
- `TREES.ltreedump`
- `TREES.ltrees`
- `TREES.ntrees`
- `TREES.qstar`
- `TREES.r_s`
- `TREES.ud`
- `WALLS.dtfac`
- `WALLS.fkar`
- `WALLS.iwallmoist`
- `WALLS.iwallmom`
- `WALLS.iwallscal`
- `WALLS.iwalltemp`
- `WALLS.lbottom`
- `WALLS.lnorec`
- `WALLS.lwritefac`
- `WALLS.nblocks`
- `WALLS.nbndpts_c`
- `WALLS.nbndpts_u`
- `WALLS.nbndpts_v`
- `WALLS.nbndpts_w`
- `WALLS.nfcts`
- `WALLS.nfctsecs_c`
- `WALLS.nfctsecs_u`
- `WALLS.nfctsecs_v`
- `WALLS.nfctsecs_w`
- `WALLS.nsolpts_c`
- `WALLS.nsolpts_u`
- `WALLS.nsolpts_v`
- `WALLS.nsolpts_w`
- `WALLS.prandtlturb`

## Variables Missing from JSON Schema
*(In namelists but not in schema)*

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

## Schema Variables Missing from Namelists
*(In schema but not in namelists)*

- `fact`
- `qt0`
- `thl0`
- `u0`
- `v0`