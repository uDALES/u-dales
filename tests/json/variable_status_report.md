# uDALES Variable Status Report

This report shows the status of variables with color coding:

- 🟢 **GREEN**: Full support (Namelist + JSON + Broadcast + Schema)
- 🟠 **ORANGE**: Warnings (Schema mismatches, duplicates)
- 🔴 **RED**: Errors (JSON/Broadcast issues)

## Executive Summary

- **Total Namelist Variables**: 251
- **JSON Support**: 245 variables
- **Namelist Broadcasts**: 245 variables
- **Non-namelist Broadcasts**: 107 variables (internal/computed)
- **Schema Coverage**: 224 variables

**Status Distribution:**
- 🟢 **207 variables** with full support
- 🟠 **54 variables** with warnings
- 🔴 **0 variables** with errors

## Variables by Namelist
*Only showing variables that are defined in Fortran namelists*

### BC Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `bcbotm`, `bcbotq`, `bcbots`, `bcbott`, `bcqfxm`, `bcqfxp`, `bcqfym`, `bcqfyp`, `bcqfz`, `bctfxm`, `bctfxp`, `bctfym`, `bctfyp`, `bctfz`, `bctopm`, `bctopq`, `bctops`, `bctopt`, `bcxm`, `bcxq`, `bcxs`, `bcxt`, `bcym`, `bcyq`, `bcys`, `bcyt`, `bczp`, `ds`, `qt_top`, `qts`, `thl_top`, `thls`, `wqsurf`, `wsvsurfdum`, `wsvtopdum`, `wtsurf`, `wttop`, `z0`, `z0h`

### CHEMISTRY Namelist *(defined in modstartup.f90)*

🟠 **Warnings**:
- `jno2` (missing from schema)
- `k1` (missing from schema)
- `lchem` (missing from schema)

### DOMAIN Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `itot`, `jtot`, `ksp`, `ktot`, `xday`, `xlat`, `xlen`, `xlon`, `xtime`, `ylen`

### DRIVER Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `chunkread_size`, `driverjobnr`, `driverstore`, `dtdriver`, `iangledeg`, `iplane`, `lchunkread`, `tdriverstart`

🟠 **Warnings**:
- `idriver` (duplicate operations: Broadcast: 2x)

### DYNAMICS Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `iadv_mom`, `iadv_qt`, `iadv_thl`, `iadv_tke`, `ipoiss`, `lqlnr`

🟠 **Warnings**:
- `iadv_sv` (duplicate operations: Broadcast: 2x)

### ENERGYBALANCE Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `bldt`, `dteb`, `flrt`, `fraction`, `grlai`, `lconstw`, `leb`, `lfactlyrs`, `lperiodicebcorr`, `lvfsparse`, `lwriteebfiles`, `nfaclyrs`, `nnz`, `rsmin`, `sinkbase`, `skylw`, `wfc`, `wgrmax`, `wsoil`, `wwilt`

### HEATPUMP Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `lfan_hp`, `lheatpump`, `nhppoints`, `q_dot_hp`, `qh_dot_hp`

### INFO Namelist *(defined in modinlet.f90)*

🟠 **Warnings**:
- `dtin` (missing from schema)
- `jgtotinl` (missing from schema)
- `kmaxin` (missing from schema)
- `nprocsinl` (missing from schema)
- `totalreadu` (missing from schema)
- `wtop` (missing from schema)

### INLET Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `di`, `dti`, `inletav`, `lfixinlet`, `lfixutauin`, `linletra`, `lreadminl`, `lstoreplane`, `lwallfunc`, `uinf`, `vinf`

### NAMCHECKSIM Namelist *(defined in modchecksim.f90)*

🟠 **Warnings**:
- `tcheck` (missing from schema)

### NAMSTATSDUMP Namelist *(defined in modstatsdump.f90)*

🟠 **Warnings**:
- `anymore` (missing from schema)
- `is` (missing from schema)
- `khigh` (missing from schema)
- `klow` (missing from schema)
- `lmintdump` (missing from schema)
- `ltdump` (missing from schema)
- `ltkedump` (missing from schema)
- `ltreedump` (missing from schema)
- `lxydump` (missing from schema)
- `lxytdump` (missing from schema)
- `lydump` (missing from schema)
- `lytdump` (missing from schema)
- `maybe` (missing from schema)
- `namstatsdump` (missing from schema)
- `removed` (missing from schema)
- `tsample` (missing from schema)
- `tstatsdump` (missing from schema)

### NAMSUBGRID Namelist *(defined in modsubgrid.f90)*

🟠 **Warnings**:
- `c_vreman` (missing from schema)
- `cf` (missing from schema)
- `cn` (missing from schema)
- `cs` (missing from schema)
- `lbuoycorr` (missing from schema)
- `ldelta` (missing from schema)
- `lmason` (missing from schema)
- `loneeqn` (missing from schema)
- `lsmagorinsky` (missing from schema)
- `lvreman` (missing from schema)
- `nmason` (missing from schema)
- `prandtl` (missing from schema)
- `rigc` (missing from schema)

### OUTPUT Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `fieldvars`, `islice`, `jslice`, `kslice`, `lislicedump`, `ljslicedump`, `lkslicedump`, `ltkedump`, `lxydump`, `lxytdump`, `lydump`, `lytdump`, `tfielddump`, `tsample`, `tstatsdump`, `tstatstart`

🟠 **Warnings**:
- `lfielddump` (duplicate operations: Broadcast: 2x)
- `lmintdump` (duplicate operations: Broadcast: 2x)
- `ltdump` (duplicate operations: Broadcast: 2x)
- `slicevars` (missing from schema)

### PHYSICS Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `dpdx`, `ifixuinf`, `igrw_damp`, `lbuoyancy`, `lcoriol`, `lmoist`, `lnudge`, `lnudgevel`, `lprofforc`, `ltempeq`, `luoutflowr`, `luvolflowr`, `lvinf`, `lvoutflowr`, `lvvolflowr`, `nnudge`, `ntimedeplw`, `ntimedepnudge`, `ntimedepsurf`, `ntimedepsw`, `ps`, `tnudge`, `tscale`, `uflowrate`, `vflowrate`

🟠 **Warnings**:
- `ltimedeplw` (duplicate operations: Broadcast: 2x)
- `ltimedepnudge` (duplicate operations: Broadcast: 2x)
- `ltimedepsurf` (duplicate operations: Broadcast: 2x)
- `ltimedepsw` (duplicate operations: Broadcast: 2x)

### PURIFS Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `epu`, `lpurif`, `qpu`

🟠 **Warnings**:
- `npurif` (duplicate operations: Broadcast: 2x)

### RUN Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `author`, `courant`, `diffnr`, `dtmax`, `iexpnr`, `irandom`, `krand`, `ladaptive`, `libm`, `lles`, `lper2inout`, `lrandomize`, `lreadmean`, `lstratstart`, `lwalldist`, `lwarmstart`, `nprocx`, `nprocy`, `randqt`, `randthl`, `randu`, `runtime`, `startfile`, `trestart`

🟠 **Warnings**:
- `runmode` (missing from schema)

### SCALARS Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `lreadscal`, `lscasrc`, `lscasrcl`, `lscasrcr`, `nscasrc`, `nscasrcl`, `nsv`

### TREES Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `cd`, `dec`, `dqdt`, `lad`, `lsize`, `ltrees`, `qstar`, `r_s`, `ud`

🟠 **Warnings**:
- `ltreedump` (duplicate operations: Broadcast: 2x)
- `ntrees` (duplicate operations: Broadcast: 2x)

### WALLS Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `dtfac`, `fkar`, `iwallmoist`, `iwallmom`, `iwallscal`, `iwalltemp`, `lbottom`, `lnorec`, `lwritefac`, `nblocks`, `nbndpts_c`, `nbndpts_u`, `nbndpts_v`, `nbndpts_w`, `nfcts`, `nfctsecs_c`, `nfctsecs_u`, `nfctsecs_v`, `nfctsecs_w`, `nsolpts_c`, `nsolpts_u`, `nsolpts_v`, `nsolpts_w`, `prandtlturb`

## 📋 Recommendations

### Medium Priority
3. **Update JSON schema** to include missing namelist variables
4. **Review schema variables** that don't correspond to namelists

### General
5. **Improve JSON coverage**: Currently 245/251 namelist variables support JSON
6. **Focus on namelist variables**: 0 namelist variables lack JSON support despite being broadcast