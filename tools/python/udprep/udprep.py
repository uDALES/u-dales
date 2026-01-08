from __future__ import annotations

from dataclasses import dataclass
import time
from typing import Any, Callable, Dict, List, Tuple, Type

# Shared base utilities for preprocessing sections.
SKIP = object()


class Section:
    def __init__(self, name: str, values: Dict[str, Any], sim: Any | None = None):
        self._name = name
        self.sim = sim
        for key, val in values.items():
            setattr(self, key, val)

    def run_steps(self, label: str, steps: List[Tuple[str, Callable[[], None]]]) -> None:
        for name, func in steps:
            start = time.perf_counter()
            print(f"[{label}] {name}...")
            func()
            elapsed = time.perf_counter() - start
            print(f"[{label}] {name} done in {elapsed:.3f} s")

    def write_changed_params(self) -> None:
        """
        Compare current values against defaults and write only changed parameters.

        This is a placeholder for logic that will diff section attributes against
        their default values and emit a minimal input file update.
        """

    def __repr__(self) -> str:
        lines = [f"{self._name}:"]
        for key, val in sorted(self.__dict__.items()):
            if key in ("_name", "sim"):
                continue
            lines.append(f"  {key}: {val}")
        return "\n".join(lines)


@dataclass(frozen=True)
class SectionSpec:
    name: str
    fields: List[str]
    defaults: Dict[str, Any | Callable[[Any], Any]]
    section_cls: Type[Section] = Section

    def get_default(self, key: str, ctx: Any) -> Any:
        default = self.defaults.get(key, SKIP)
        return default(ctx) if callable(default) else default


class _DefaultContext:
    def __init__(self, sim: Any, values: Dict[str, Any]):
        self._sim = sim
        self._values = values

    def __getattr__(self, name: str) -> Any:
        if name in self._values:
            return self._values[name]
        if self._sim is not None and hasattr(self._sim, name):
            return getattr(self._sim, name)
        raise AttributeError(name)
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
        self._section_spec_map = {spec.name: spec for spec in self.SECTION_SPECS}

        for spec in self.SECTION_SPECS:
            values = {}
            ctx = _DefaultContext(sim, values)
            for key in spec.fields:
                if hasattr(sim, key):
                    values[key] = getattr(sim, key)
                else:
                    default = spec.get_default(key, ctx)
                    if default is SKIP:
                        continue
                    values[key] = default
            setattr(self, spec.name, spec.section_cls(spec.name, values, sim=self.sim))

    def addvar(self, name: str, value: Any, section: str | None = None) -> None:
        if section is None:
            raise ValueError("section is required; add variables to a specific section")
        spec = self._section_spec_map.get(section)
        if spec is None:
            raise ValueError(f"Unknown section: {section}")
        if name not in spec.fields:
            raise ValueError(f"{name} is not a valid field for section '{section}'")
        section_values = self._section_values(section)
        section_values[name] = value
        setattr(self, section, spec.section_cls(section, section_values, sim=self.sim))

    def _section_values(self, section: str) -> dict[str, Any]:
        section_obj = getattr(self, section, None)
        if section_obj is None:
            return {}
        return {
            key: val
            for key, val in section_obj.__dict__.items()
            if key not in ("_name", "sim")
        }

    def __repr__(self) -> str:
        lines = ["UDPrep:"]
        if self.sim is not None:
            sim_expnr = getattr(self.sim, "expnr", None)
            sim_path = getattr(self.sim, "path", None)
            if sim_expnr is not None:
                lines.append(f"  expnr: {sim_expnr}")
            if sim_path is not None:
                lines.append(f"  path: {sim_path}")
        for spec in self.SECTION_SPECS:
            section_obj = getattr(self, spec.name, None)
            if section_obj is None:
                continue
            lines.append(section_obj.__repr__())
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
