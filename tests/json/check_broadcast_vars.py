#!/usr/bin/env python3
"""
Check which variables from namelists are missing from broadcast_config_parameters.
"""

import re

# Variables defined in namelists (extracted from modstartup.f90)
namelist_variables = {
    "RUN": [
        "iexpnr", "lwarmstart", "lstratstart", "startfile",
        "runtime", "dtmax", "trestart", "ladaptive",
        "irandom", "randu", "randthl", "randqt", "krand",
        "courant", "diffnr", "author",
        "libm", "lles",
        "lper2inout", "lwalldist",
        "lreadmean",
        "nprocx", "nprocy",
        "lrandomize", "runmode"
    ],
    "DOMAIN": [
        "itot", "jtot", "ktot", "xlen", "ylen",
        "xlat", "xlon", "xday", "xtime", "ksp"
    ],
    "PHYSICS": [
        "ps", "igrw_damp", "lmoist", "lcoriol", "lbuoyancy", "ltempeq",
        "lprofforc", "ifixuinf", "lvinf", "tscale", "dpdx",
        "luoutflowr", "lvoutflowr", "luvolflowr", "lvvolflowr",
        "uflowrate", "vflowrate",
        "lnudge", "lnudgevel", "tnudge", "nnudge",
        "ltimedepsurf", "ntimedepsurf", "ltimedepnudge", "ntimedepnudge",
        "ltimedeplw", "ntimedeplw", "ltimedepsw", "ntimedepsw"
    ],
    "DYNAMICS": [
        "lqlnr", "ipoiss",
        "iadv_mom", "iadv_tke", "iadv_thl", "iadv_qt", "iadv_sv"
    ],
    "BC": [
        "BCxm", "BCxT", "BCxq", "BCxs",
        "BCym", "BCyT", "BCyq", "BCys",
        "BCtopm", "BCtopT", "BCtopq", "BCtops",
        "BCbotm", "BCbotT", "BCbotq", "BCbots",
        "bctfxm", "bctfxp", "bctfym", "bctfyp", "bctfz",
        "bcqfxm", "bcqfxp", "bcqfym", "bcqfyp", "bcqfz",
        "wttop", "thl_top", "qt_top", "qts", "wsvsurfdum", "wsvtopdum",
        "wtsurf", "wqsurf", "thls", "z0", "z0h", "BCzp", "ds"
    ],
    "INLET": [
        "Uinf", "Vinf", "di", "dti", "inletav", "linletRA",
        "lstoreplane", "lreadminl", "lfixinlet", "lfixutauin",
        "lwallfunc"
    ],
    "DRIVER": [
        "idriver", "tdriverstart", "driverjobnr", "dtdriver",
        "driverstore", "iplane", "iangledeg",
        "lchunkread", "chunkread_size"
    ],
    "WALLS": [
        "nblocks", "nfcts", "iwallmom", "iwalltemp", "iwallmoist", "iwallscal",
        "nsolpts_u", "nsolpts_v", "nsolpts_w", "nsolpts_c",
        "nbndpts_u", "nbndpts_v", "nbndpts_w", "nbndpts_c",
        "nfctsecs_u", "nfctsecs_v", "nfctsecs_w", "nfctsecs_c", "lbottom", "lnorec",
        "prandtlturb", "fkar", "lwritefac", "dtfac"
    ],
    "ENERGYBALANCE": [
        "lEB", "lwriteEBfiles", "lperiodicEBcorr", "sinkbase", "lconstW", "dtEB", "bldT", "flrT", "wsoil", "wgrmax", "wwilt", "wfc",
        "skyLW", "GRLAI", "rsmin", "nfaclyrs", "lfacTlyrs", "lvfsparse", "nnz", "fraction"
    ],
    "SCALARS": [
        "lreadscal", "lscasrc", "lscasrcl", "lscasrcr",
        "nsv", "nscasrc", "nscasrcl"
    ],
    "CHEMISTRY": [
        "lchem", "k1", "JNO2"
    ],
    "OUTPUT": [
        "lfielddump", "tfielddump", "fieldvars",
        "ltdump", "lydump", "lytdump", "lxydump", "lxytdump", "lmintdump", "ltkedump",
        "slicevars", "lkslicedump", "kslice", "lislicedump", "islice", "ljslicedump", "jslice",
        "tstatsdump", "tsample", "tstatstart"
    ],
    "TREES": [
        "ltrees", "ntrees", "cd", "dec", "ud", "lad", "Qstar", "dQdt", "lsize", "r_s", "ltreedump"
    ],
    "PURIFS": [
        "lpurif", "npurif", "Qpu", "epu"
    ],
    "HEATPUMP": [
        "lheatpump", "lfan_hp", "nhppoints", "Q_dot_hp", "QH_dot_hp"
    ]
}

# Variables currently being broadcast (extracted from broadcast_config_parameters)
broadcast_variables = [
    "itot", "jtot", "ktot", "nprocx", "nprocy", "iexpnr", "runmode",
    "runtime", "dtmax", "lwarmstart", "lfielddump", "ltempeq", "lmoist",
    "lstratstart", "lreadscal", "lscasrc", "lscasrcl", "lscasrcr", "lbuoyancy",
    "lper2inout", "libm", "lnudge", "lnudgevel", "nnudge", "tnudge",
    "ltimedepsurf", "ltimedepnudge", "ltimedeplw", "ltimedepsw",
    "ntimedepsurf", "ntimedepnudge", "ntimedeplw", "ntimedepsw",
    "lwalldist", "lles", "linletRA", "lfixinlet", "lfixutauin",
    "idriver", "tdriverstart", "driverjobnr", "dtdriver", "driverstore",
    "lchunkread", "chunkread_size",
    "BCxm", "BCxT", "BCxq", "BCxs", "BCym", "BCyT", "BCyq", "BCys",
    "BCtopm", "BCtopT", "BCtopq", "BCtops", "BCbotm", "BCbotT", "BCbotq", "BCbots",
    "BCzp", "ds",
    "trestart", "tfielddump", "tsample", "tstatsdump", "tstatstart", "nsv",
    "nscasrc", "nscasrcl", "fieldvars", "slicevars",
    "xlen", "ylen", "xlat", "xlon", "xday", "xtime", "ladaptive", "courant", "diffnr", "lrandomize"
]

# Convert to set for easier comparison
broadcast_set = set(broadcast_variables)

print("=== MISSING VARIABLES FROM BROADCAST ===\n")

total_missing = 0
for section, variables in namelist_variables.items():
    missing = [var for var in variables if var not in broadcast_set]
    if missing:
        print(f"{section}:")
        for var in missing:
            print(f"  - {var}")
            total_missing += 1
        print()

print(f"=== SUMMARY ===")
print(f"Total namelist variables: {sum(len(vars) for vars in namelist_variables.values())}")
print(f"Total broadcast variables: {len(broadcast_variables)}")
print(f"Missing from broadcast: {total_missing}")

# Check for variables in broadcast that might not be in namelists
all_namelist_vars = set()
for variables in namelist_variables.values():
    all_namelist_vars.update(variables)

extra_broadcast = [var for var in broadcast_variables if var not in all_namelist_vars]
if extra_broadcast:
    print(f"\n=== EXTRA VARIABLES IN BROADCAST (not in namelists) ===")
    for var in extra_broadcast:
        print(f"  - {var}")

print(f"\nMissing variables that should be added to broadcast_config_parameters:")
all_missing = []
for section, variables in namelist_variables.items():
    missing = [var for var in variables if var not in broadcast_set]
    all_missing.extend(missing)

for var in sorted(all_missing):
    print(f"call MPI_BCAST({var}, 1, MPI_TYPE, 0, comm3d, mpierr)")