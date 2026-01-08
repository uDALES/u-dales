from __future__ import annotations

from typing import Any, Callable, Dict

import numpy as np

_SKIP = object()


class UDPrep:
    """
    Preprocessing config derived from UDBase with defaults from preprocessing.m.
    """

    DEFAULTS: Dict[str, Any | Callable[["UDPrep"], Any]] = {
        # &RUN
        "ltrees": 0,
        "ltreesfile": 0,
        "tree_dz": lambda self: 0 if self.ltrees and not self.ltreesfile else _SKIP,
        "tree_dx": lambda self: 0 if self.ltrees and not self.ltreesfile else _SKIP,
        "tree_dy": lambda self: 0 if self.ltrees and not self.ltreesfile else _SKIP,
        "tree_h": lambda self: 0 if self.ltrees and not self.ltreesfile else _SKIP,
        "tree_w": lambda self: 0 if self.ltrees and not self.ltreesfile else _SKIP,
        "tree_b": lambda self: 0 if self.ltrees and not self.ltreesfile else _SKIP,
        "nrows": lambda self: 0 if self.ltrees and not self.ltreesfile else _SKIP,
        "treesfile": lambda self: "" if self.ltrees and self.ltreesfile else _SKIP,
        "lpurif": 0,
        "luoutflowr": 0,
        "lvoutflowr": 0,
        "luvolflowr": 0,
        "lvvolflowr": 0,
        # &DOMAIN
        "itot": 64,
        "xlen": 64,
        "jtot": 64,
        "ylen": 64,
        "ktot": 96,
        "dx": lambda self: self.xlen / self.itot,
        "dy": lambda self: self.ylen / self.jtot,
        # BCs
        "BCxm": 1,
        "BCym": 1,
        # &ENERGYBALANCE
        "lEB": 0,
        "lfacTlyrs": 0,
        # &WALLS
        "iwallmom": 3,
        "iwalltemp": 1,
        "lbottom": 0,
        "lwritefac": 0,
        # &PHYSICS
        "ltempeq": 0,
        "lmoist": 0,
        "lchem": 0,
        "lprofforc": 0,
        "lcoriol": 0,
        "idriver": 0,
        "ldp": lambda self: (
            1
            if (
                not self.luoutflowr
                and not self.lvoutflowr
                and not self.luvolflowr
                and not self.lvvolflowr
                and not self.lprofforc
                and not self.lcoriol
                and self.idriver != 2
            )
            else 0
        ),
        # &INPS
        "zsize": 96,
        "lzstretch": 0,
        "stl_file": "",
        "gen_geom": True,
        "geom_path": "",
        "diag_neighbs": True,
        "stl_ground": True,
        "stretchconst": lambda self: 0.01 if self.lzstretch else _SKIP,
        "lstretchexp": lambda self: 0 if self.lzstretch else _SKIP,
        "lstretchexpcheck": lambda self: 0 if self.lzstretch else _SKIP,
        "lstretchtanh": lambda self: 0 if self.lzstretch else _SKIP,
        "lstretch2tanh": lambda self: 0 if self.lzstretch else _SKIP,
        "hlin": lambda self: 0 if self.lzstretch else _SKIP,
        "dzlin": lambda self: 0 if self.lzstretch else _SKIP,
        "dz": lambda self: self.dzlin if self.lzstretch else (self.zsize / self.ktot),
        "maxlen": lambda self: 10 if self.lEB else np.inf,
        "u0": 0,
        "v0": 0,
        "tke": 0,
        "dpdx": 0,
        "dpdy": 0,
        "thl0": 288,
        "qt0": 0,
        "nsv": 0,
        "sv10": lambda self: 0 if self.nsv > 0 else _SKIP,
        "sv20": lambda self: 0 if self.nsv > 0 else _SKIP,
        "sv30": lambda self: 0 if self.nsv > 0 else _SKIP,
        "sv40": lambda self: 0 if self.nsv > 0 else _SKIP,
        "sv50": lambda self: 0 if self.nsv > 0 else _SKIP,
        "lscasrc": lambda self: 0 if self.nsv > 0 else _SKIP,
        "lscasrcl": lambda self: 0 if self.nsv > 0 else _SKIP,
        "lscasrcr": lambda self: 0 if self.nsv > 0 else _SKIP,
        "xS": lambda self: -1 if self.nsv > 0 else _SKIP,
        "yS": lambda self: -1 if self.nsv > 0 else _SKIP,
        "zS": lambda self: -1 if self.nsv > 0 else _SKIP,
        "SSp": lambda self: -1 if self.nsv > 0 else _SKIP,
        "sigSp": lambda self: -1 if self.nsv > 0 else _SKIP,
        "nscasrc": lambda self: 0 if self.nsv > 0 else _SKIP,
        "xSb": lambda self: -1 if self.nsv > 0 else _SKIP,
        "ySb": lambda self: -1 if self.nsv > 0 else _SKIP,
        "zSb": lambda self: -1 if self.nsv > 0 else _SKIP,
        "xSe": lambda self: -1 if self.nsv > 0 else _SKIP,
        "ySe": lambda self: -1 if self.nsv > 0 else _SKIP,
        "zSe": lambda self: -1 if self.nsv > 0 else _SKIP,
        "SSl": lambda self: -1 if self.nsv > 0 else _SKIP,
        "sigSl": lambda self: -1 if self.nsv > 0 else _SKIP,
        "nscasrcl": lambda self: 0 if self.nsv > 0 else _SKIP,
        "lapse": 0,
        "w_s": 0,
        "R": 0,
        "libm": 1,
        "isolid_bound": 1,
        "ifacsec": 1,
        "read_types": 0,
        "types_path": lambda self: 0 if self.read_types else _SKIP,
        "xazimuth": lambda self: 90 if self.lEB else _SKIP,
        "ltimedepsw": lambda self: 0 if self.lEB else _SKIP,
        "ishortwave": lambda self: 1 if self.lEB else _SKIP,
        "isolar": lambda self: 1 if self.lEB else _SKIP,
        "runtime": lambda self: 0 if self.lEB else _SKIP,
        "dtEB": lambda self: 10.0 if self.lEB else _SKIP,
        "dtSP": lambda self: (self.dtEB if self.lEB else _SKIP),
        "solarazimuth": lambda self: 135 if (self.lEB and self.isolar == 1) else _SKIP,
        "solarzenith": lambda self: 28.4066 if (self.lEB and self.isolar == 1) else _SKIP,
        "I": lambda self: 800 if (self.lEB and self.isolar == 1) else _SKIP,
        "Dsky": lambda self: 418.8041 if (self.lEB and self.isolar == 1) else _SKIP,
        "longitude": lambda self: -0.13 if (self.lEB and self.isolar == 2) else _SKIP,
        "latitude": lambda self: 51.5 if (self.lEB and self.isolar == 2) else _SKIP,
        "timezone": lambda self: 0 if (self.lEB and self.isolar == 2) else _SKIP,
        "elevation": lambda self: 0 if (self.lEB and self.isolar == 2) else _SKIP,
        "hour": lambda self: 6 if (self.lEB and self.isolar == 2) else (0 if (self.lEB and self.isolar == 3) else _SKIP),
        "minute": lambda self: 0 if (self.lEB and self.isolar in (2, 3)) else _SKIP,
        "second": lambda self: 0 if (self.lEB and self.isolar in (2, 3)) else _SKIP,
        "year": lambda self: 2011 if (self.lEB and self.isolar == 2) else (0 if (self.lEB and self.isolar == 3) else _SKIP),
        "month": lambda self: 9 if (self.lEB and self.isolar == 2) else (6 if (self.lEB and self.isolar == 3) else _SKIP),
        "day": lambda self: 30 if (self.lEB and self.isolar == 2) else (1 if (self.lEB and self.isolar == 3) else _SKIP),
        "weatherfname": lambda self: "" if (self.lEB and self.isolar == 3) else _SKIP,
        "psc_res": lambda self: 0.1 if self.lEB else _SKIP,
        "lvfsparse": lambda self: False if self.lEB else _SKIP,
        "calc_vf": lambda self: True if self.lEB else _SKIP,
        "maxD": lambda self: np.inf if self.lEB else _SKIP,
        "vf_path": lambda self: "" if (self.lEB and not self.calc_vf) else _SKIP,
        "view3d_out": lambda self: 0 if self.lEB else _SKIP,
        "facT": 288.0,
        "nfaclyrs": 3,
        "nfcts": 0,
        "facT_file": "",
        "factypes": lambda self: _SKIP,
    }

    DEFAULT_ORDER = [
        "ltrees",
        "ltreesfile",
        "lpurif",
        "luoutflowr",
        "lvoutflowr",
        "luvolflowr",
        "lvvolflowr",
        "itot",
        "xlen",
        "jtot",
        "ylen",
        "ktot",
        "dx",
        "dy",
        "BCxm",
        "BCym",
        "lEB",
        "lfacTlyrs",
        "iwallmom",
        "iwalltemp",
        "lbottom",
        "lwritefac",
        "ltempeq",
        "lmoist",
        "lchem",
        "lprofforc",
        "lcoriol",
        "idriver",
        "ldp",
        "zsize",
        "lzstretch",
        "stl_file",
        "gen_geom",
        "geom_path",
        "diag_neighbs",
        "stl_ground",
        "stretchconst",
        "lstretchexp",
        "lstretchexpcheck",
        "lstretchtanh",
        "lstretch2tanh",
        "hlin",
        "dzlin",
        "dz",
        "maxlen",
        "u0",
        "v0",
        "tke",
        "dpdx",
        "dpdy",
        "thl0",
        "qt0",
        "nsv",
        "sv10",
        "sv20",
        "sv30",
        "sv40",
        "sv50",
        "lscasrc",
        "lscasrcl",
        "lscasrcr",
        "xS",
        "yS",
        "zS",
        "SSp",
        "sigSp",
        "nscasrc",
        "xSb",
        "ySb",
        "zSb",
        "xSe",
        "ySe",
        "zSe",
        "SSl",
        "sigSl",
        "nscasrcl",
        "lapse",
        "w_s",
        "R",
        "libm",
        "isolid_bound",
        "ifacsec",
        "read_types",
        "types_path",
        "xazimuth",
        "ltimedepsw",
        "ishortwave",
        "isolar",
        "runtime",
        "dtEB",
        "dtSP",
        "solarazimuth",
        "solarzenith",
        "I",
        "Dsky",
        "longitude",
        "latitude",
        "timezone",
        "elevation",
        "hour",
        "minute",
        "second",
        "year",
        "month",
        "day",
        "weatherfname",
        "psc_res",
        "lvfsparse",
        "calc_vf",
        "maxD",
        "vf_path",
        "view3d_out",
        "facT",
        "nfaclyrs",
        "nfcts",
        "facT_file",
        "factypes",
        "tree_dz",
        "tree_dx",
        "tree_dy",
        "tree_h",
        "tree_w",
        "tree_b",
        "nrows",
        "treesfile",
    ]

    def __init__(self, sim):
        self.sim = sim
        self.expnr = getattr(sim, "expnr", None)
        self.path = getattr(sim, "path", None)

        for key in self.DEFAULT_ORDER:
            if hasattr(sim, key):
                setattr(self, key, getattr(sim, key))
                continue
            default = self._get_default(key)
            if default is _SKIP:
                continue
            setattr(self, key, default)

        if self.ltempeq == 0 or (self.iwalltemp == 1 and self.iwallmom == 2):
            self.iwallmom = 3

    def _get_default(self, key: str) -> Any:
        default = self.DEFAULTS.get(key, _SKIP)
        return default(self) if callable(default) else default

    def addvar(self, name: str, value: Any) -> None:
        if not hasattr(self, name):
            setattr(self, name, value)

    def __repr__(self) -> str:
        lines = ["UDPrep:"]
        header_keys = ["expnr", "path"]
        for key in header_keys:
            if hasattr(self, key):
                lines.append(f"  {key}: {getattr(self, key)}")
        for key in self.DEFAULT_ORDER:
            if hasattr(self, key):
                lines.append(f"  {key}: {getattr(self, key)}")
        return "\n".join(lines)

    def __str__(self) -> str:
        return self.__repr__()
