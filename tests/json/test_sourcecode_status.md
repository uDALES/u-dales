# uDALES Variable Status Report

This report shows the status of variables with color coding:

- 🟢 GREEN: Full support (Namelist + JSON + Broadcast + Schema)
- 🟠 ORANGE: Warnings (Schema mismatches, duplicates)
- 🔴 RED: Errors (JSON/Broadcast issues)

## Executive Summary

- **Total Namelist Variables**: 244
- **JSON Support**: 244 variables
- **Namelist Broadcasts**: 243 variables
- **Non-namelist Broadcasts**: 1 variables (internal/computed)
- **Schema Coverage**: 224 variables

**Status Distribution:**
- 🟢 **219 variables** with full support
- 🟠 **25 variables** with warnings
- 🔴 **0 variables** with errors

## Variables by Namelist
*Only showing variables that are defined in Fortran namelists*

### BC

🟢 **Full Support**: `bcbotm`, `bcbotq`, `bcbots`, `bcbott`, `bcqfxm`, `bcqfxp`, `bcqfym`, `bcqfyp`, `bcqfz`, `bctfxm`, `bctfxp`, `bctfym`, `bctfyp`, `bctfz`, `bctopm`, `bctopq`, `bctops`, `bctopt`, `bcxm`, `bcxq`, `bcxs`, `bcxt`, `bcym`, `bcyq`, `bcys`, `bcyt`, `bczp`, `ds`, `qt_top`, `qts`, `thl_top`, `thls`, `wqsurf`, `wsvsurfdum`, `wsvtopdum`, `wtsurf`, `wttop`, `z0`, `z0h`

### CHEMISTRY

🟠 **Warnings**:
- `jno2` (missing from schema)
- `k1` (missing from schema)
- `lchem` (missing from schema)

### DOMAIN

🟢 **Full Support**: `itot`, `jtot`, `ksp`, `ktot`, `xday`, `xlat`, `xlen`, `xlon`, `xtime`, `ylen`

### DRIVER

🟢 **Full Support**: `chunkread_size`, `driverjobnr`, `driverstore`, `dtdriver`, `iangledeg`, `idriver`, `iplane`, `lchunkread`, `tdriverstart`

### DYNAMICS

🟢 **Full Support**: `iadv_mom`, `iadv_qt`, `iadv_sv`, `iadv_thl`, `iadv_tke`, `ipoiss`, `lqlnr`

### ENERGYBALANCE

🟢 **Full Support**: `bldt`, `dteb`, `flrt`, `fraction`, `grlai`, `lconstw`, `leb`, `lfactlyrs`, `lperiodicebcorr`, `lvfsparse`, `lwriteebfiles`, `nfaclyrs`, `nnz`, `rsmin`, `sinkbase`, `skylw`, `wfc`, `wgrmax`, `wsoil`, `wwilt`

### HEATPUMP

🟢 **Full Support**: `lfan_hp`, `lheatpump`, `nhppoints`, `q_dot_hp`, `qh_dot_hp`

### INFO

🟠 **Warnings**:
- `dtin` (missing from schema)
- `jgtotinl` (missing from schema)
- `kmaxin` (missing from schema)
- `nprocsinl` (missing from schema)
- `totalreadu` (missing from schema)
- `wtop` (missing from schema)

### INLET

🟢 **Full Support**: `di`, `dti`, `inletav`, `lfixinlet`, `lfixutauin`, `linletra`, `lreadminl`, `lstoreplane`, `lwallfunc`, `uinf`, `vinf`

### NAMSUBGRID

🟠 **Warnings**:
- `c_vreman` (missing from schema)
- `cf` (missing from schema)
- `cn` (missing from schema)
- `lbuoycorr` (missing from schema)
- `ldelta` (missing from schema)
- `lmason` (missing from schema)
- `loneeqn` (missing from schema)
- `lsmagorinsky` (missing from schema)
- `lvreman` (missing from schema)
- `nmason` (missing from schema)
- `prandtl` (missing from schema)
- `rigc` (missing from schema)
- `sg_cs` (missing from schema)

### OUTPUT

🟢 **Full Support**: `fieldvars`, `islice`, `jslice`, `kslice`, `lfielddump`, `lislicedump`, `ljslicedump`, `lkslicedump`, `lmintdump`, `ltdump`, `ltkedump`, `lxydump`, `lxytdump`, `lydump`, `lytdump`, `tfielddump`, `tsample`, `tstatsdump`, `tstatstart`

🟠 **Warnings**:
- `slicevars` (missing from schema)
- `tcheck` (missing from schema)

### PHYSICS

🟢 **Full Support**: `dpdx`, `ifixuinf`, `igrw_damp`, `lbuoyancy`, `lcoriol`, `lmoist`, `lnudge`, `lnudgevel`, `lprofforc`, `ltempeq`, `ltimedeplw`, `ltimedepnudge`, `ltimedepsurf`, `ltimedepsw`, `luoutflowr`, `luvolflowr`, `lvinf`, `lvoutflowr`, `lvvolflowr`, `nnudge`, `ntimedeplw`, `ntimedepnudge`, `ntimedepsurf`, `ntimedepsw`, `ps`, `tnudge`, `tscale`, `uflowrate`, `vflowrate`

### PURIFS

🟢 **Full Support**: `epu`, `lpurif`, `npurif`, `qpu`

### RUN

🟢 **Full Support**: `author`, `courant`, `diffnr`, `dtmax`, `iexpnr`, `irandom`, `krand`, `ladaptive`, `libm`, `lles`, `lper2inout`, `lrandomize`, `lreadmean`, `lstratstart`, `lwalldist`, `lwarmstart`, `nprocx`, `nprocy`, `randqt`, `randthl`, `randu`, `runtime`, `startfile`, `trestart`

🟠 **Warnings**:
- `runmode` (missing from schema)

### SCALARS

🟢 **Full Support**: `lreadscal`, `lscasrc`, `lscasrcl`, `lscasrcr`, `nscasrc`, `nscasrcl`, `nsv`

### TREES

🟢 **Full Support**: `cd`, `dec`, `dqdt`, `lad`, `lsize`, `ltreedump`, `ltrees`, `ntrees`, `qstar`, `r_s`, `ud`

### WALLS

🟢 **Full Support**: `dtfac`, `fkar`, `iwallmoist`, `iwallmom`, `iwallscal`, `iwalltemp`, `lbottom`, `lnorec`, `lwritefac`, `nblocks`, `nbndpts_c`, `nbndpts_u`, `nbndpts_v`, `nbndpts_w`, `nfcts`, `nfctsecs_c`, `nfctsecs_u`, `nfctsecs_v`, `nfctsecs_w`, `nsolpts_c`, `nsolpts_u`, `nsolpts_v`, `nsolpts_w`, `prandtlturb`

## Recommendations

### Medium Priority
3. **Update JSON schema** to include missing namelist variables
4. **Review schema variables** that don't correspond to namelists

### General
5. **Improve JSON coverage**: Currently 244/244 namelist variables support JSON
6. **Focus on namelist variables**: 1 namelist variables lack JSON support despite being broadcast