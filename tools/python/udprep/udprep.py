from __future__ import annotations

from typing import Any, List

from .udprep_common import SKIP, Section
from .udprep_ic import SPEC as IC_SPEC
from .udprep_ibm import SPEC as IBM_SPEC
from .udprep_radiation import SPEC as RADIATION_SPEC
from .udprep_scalars import SPEC as SCALARS_SPEC
from .udprep_seb import SPEC as SEB_SPEC
from .udprep_vegetation import SPEC as VEGETATION_SPEC


class UDPrep:
    """
    Preprocessing config derived from UDBase with defaults from preprocessing.m.
    """

    SECTION_SPECS = [
        SEB_SPEC,
        IBM_SPEC,
        IC_SPEC,
        VEGETATION_SPEC,
        SCALARS_SPEC,
        RADIATION_SPEC,
    ]

    def __init__(self, expnr, path=None, load_geometry: bool = True):
        from udbase import UDBase

        if isinstance(expnr, UDBase):
            sim = expnr
        else:
            sim = UDBase(expnr, path, load_geometry=load_geometry)

        self.sim = sim
        self.expnr = getattr(sim, "expnr", None)
        self.path = getattr(sim, "path", None)
        self._section_spec_map = {spec.name: spec for spec in self.SECTION_SPECS}

        for spec in self.SECTION_SPECS:
            for key in spec.fields:
                if hasattr(sim, key):
                    setattr(self, key, getattr(sim, key))
                    continue
                default = spec.get_default(key, self)
                if default is SKIP:
                    continue
                setattr(self, key, default)

        for spec in self.SECTION_SPECS:
            values = {key: getattr(self, key) for key in spec.fields if hasattr(self, key)}
            setattr(self, spec.name, spec.section_cls(spec.name, values, sim=self.sim))

    def addvar(self, name: str, value: Any, section: str | None = None) -> None:
        if section is not None and hasattr(self, section):
            section_values = self._section_values(section)
            section_values[name] = value
            spec = self._section_spec_map.get(section)
            section_cls = spec.section_cls if spec is not None else Section
            setattr(self, section, section_cls(section, section_values, sim=self.sim))
        if not hasattr(self, name):
            setattr(self, name, value)

    def _section_values(self, section: str) -> dict[str, Any]:
        section_obj = getattr(self, section, None)
        if section_obj is None:
            return {}
        return {
            key: val
            for key, val in section_obj.__dict__.items()
            if key != "_name"
        }

    def __repr__(self) -> str:
        lines = ["UDPrep:"]
        header_keys = ["expnr", "path"]
        for key in header_keys:
            if hasattr(self, key):
                lines.append(f"  {key}: {getattr(self, key)}")
        for spec in self.SECTION_SPECS:
            section_obj = getattr(self, spec.name, None)
            if section_obj is None:
                continue
            lines.append(f"  [{spec.name}]")
            for key, val in sorted(section_obj.__dict__.items()):
                if key == "_name":
                    continue
                lines.append(f"    {key}: {val}")
        return "\n".join(lines)

    def run_all(self) -> None:
        """Run preprocessing in a write_inputs.m-style sequence."""
        self.ibm.run_all()
        self.ic.run_all()
        if self.scalars.nsv > 0:
            self.scalars.run_all()
        if self.vegetation.ltrees or self.vegetation.ltreesfile:
            self.vegetation.run_all()
        if self.seb.lEB:
            self.radiation.run_all()
            self.seb.run_all()

    def __str__(self) -> str:
        return self.__repr__()
