# uDALES Variable Status Report

This report shows the status of variables with color coding:

- 🟢 **GREEN**: Full support (Namelist + JSON + Broadcast + Schema)
- 🟠 **ORANGE**: Warnings (Schema mismatches, duplicates)
- 🔴 **RED**: Errors (JSON/Broadcast issues)

## Executive Summary

- **Total Namelist Variables**: 249
- **JSON Support**: 220 variables
- **Namelist Broadcasts**: 199 variables
- **Non-namelist Broadcasts**: 106 variables (internal/computed)
- **Schema Coverage**: 225 variables

**Status Distribution:**
- 🟢 **170 variables** with full support
- 🟠 **47 variables** with warnings
- 🔴 **42 variables** with errors

## 🚨 Critical Configuration Issue

**JSON Coverage Gap**: 220 JSON reads vs 199 namelist broadcasts
**Missing JSON Support**: 21 namelist variables

This means many namelist variables are broadcast to all MPI processes but cannot be configured via JSON.
Additionally, 106 internal/computed variables are also broadcast (this is normal).

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

🟠 **Warnings**:
- `idriver` (duplicate operations: Broadcast: 2x)

🔴 **Errors**:
- `chunkread_size` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `driverjobnr` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `driverstore` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `dtdriver` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `iangledeg` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `iplane` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `lchunkread` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `tdriverstart` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast

### DYNAMICS Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `iadv_mom`, `iadv_qt`, `iadv_thl`, `iadv_tke`, `ipoiss`, `lqlnr`

🟠 **Warnings**:
- `iadv_sv` (duplicate operations: Broadcast: 2x)

### ENERGYBALANCE Namelist *(defined in modstartup.f90)*

🔴 **Errors**:
- `bldt` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `dteb` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `flrt` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `fraction` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `grlai` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `lconstw` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `leb` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `lfactlyrs` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `lperiodicebcorr` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `lvfsparse` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `lwriteebfiles` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `nfaclyrs` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `nnz` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `rsmin` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `sinkbase` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `skylw` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `wfc` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `wgrmax` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `wsoil` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `wwilt` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast

### HEATPUMP Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `lfan_hp`, `lheatpump`, `nhppoints`, `q_dot_hp`, `qh_dot_hp`

### INFO Namelist *(defined in modinlet.f90)*

🟠 **Warnings**:
- `dtin` (missing from schema)
- `inl` (missing from schema)
- `jgtotinl` (missing from schema)
- `kmaxin` (missing from schema)
- `namezinlet` (missing from schema)
- `nprocsinl` (missing from schema)
- `totalreadu` (missing from schema)
- `wtop` (missing from schema)
- `zgrid` (missing from schema)

### INLET Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `di`, `dti`, `inletav`, `lfixinlet`, `lfixutauin`, `linletra`, `lreadminl`, `lstoreplane`, `lwallfunc`, `uinf`, `vinf`

### NAMCHECKSIM Namelist *(defined in modchecksim.f90)*

🟠 **Warnings**:
- `tcheck` (missing from schema)

### NAMSTATSDUMP Namelist *(defined in modstatsdump.f90)*

🟠 **Warnings**:
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

🟢 **Full Support**: `dpdx`, `ifixuinf`, `igrw_damp`, `lbuoyancy`, `lcoriol`, `lmoist`, `lprofforc`, `ltempeq`, `ltimedeplw`, `ltimedepnudge`, `ltimedepsurf`, `ltimedepsw`, `lvinf`, `ps`, `tscale`

🔴 **Errors**:
- `lnudge` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `lnudgevel` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `luoutflowr` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `luvolflowr` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `lvoutflowr` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `lvvolflowr` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `nnudge` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `ntimedeplw` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `ntimedepnudge` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `ntimedepsurf` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `ntimedepsw` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `tnudge` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `uflowrate` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast
- `vflowrate` (JSON: ✓, Broadcast: ✗, Schema: ✓) - No MPI broadcast

### PURIFS Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `epu`, `lpurif`, `qpu`

🟠 **Warnings**:
- `npurif` (duplicate operations: Broadcast: 2x)

### RUN Namelist *(defined in modstartup.f90)*

🟢 **Full Support**: `author`, `courant`, `diffnr`, `dtmax`, `iexpnr`, `irandom`, `krand`, `ladaptive`, `libm`, `lles`, `lper2inout`, `lrandomize`, `lreadmean`, `lstratstart`, `lwalldist`, `lwarmstart`, `nprocx`, `nprocy`, `randqt`, `randthl`, `randu`, `runmode`, `runtime`, `startfile`, `trestart`

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

### High Priority
1. **Add JSON reading support** for variables in namelists
2. **Add MPI broadcast calls** for JSON-read variables

### Medium Priority
3. **Update JSON schema** to include missing namelist variables
4. **Review schema variables** that don't correspond to namelists

### General
5. **Improve JSON coverage**: Currently 220/249 namelist variables support JSON
6. **Focus on namelist variables**: 21 namelist variables lack JSON support despite being broadcast