from __future__ import annotations

from dataclasses import dataclass
import inspect
import json
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple, Type

import numpy as np

# Shared base utilities for preprocessing sections.
SKIP = object()


class Section:
    def __init__(
        self,
        name: str,
        values: Dict[str, Any],
        sim: Any | None = None,
        defaults: Dict[str, Any] | None = None,
    ):
        self._name = name
        self.sim = sim
        self._defaults = defaults or {}
        for key, val in values.items():
            setattr(self, key, val)

    def run_steps(
        self,
        label: str,
        steps: List[Tuple[str, Callable[[], None]]],
        **kwargs: Any,
    ) -> None:
        for name, func in steps:
            start = time.perf_counter()
            print(f"[{label}] {name}...")
            if kwargs:
                sig = inspect.signature(func)
                params = sig.parameters.values()
                if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in params):
                    func(**kwargs)
                else:
                    call_kwargs = {k: v for k, v in kwargs.items() if k in sig.parameters}
                    func(**call_kwargs)
            else:
                func()
            elapsed = time.perf_counter() - start
            print(f"[{label}] {name} done in {elapsed:.3f} s")

    @staticmethod
    def load_defaults_json(path: str | None = None) -> Dict[str, Dict[str, Any]]:
        defaults_path = Path(path) if path is not None else Path(__file__).with_name("defaults.json")
        try:
            with defaults_path.open("r", encoding="ascii") as f:
                return json.load(f)
        except OSError:
            return {}

    @staticmethod
    def resolve_default(
        section: str,
        key: str,
        ctx: Any,
        defaults_json: Dict[str, Dict[str, Any]],
        fallback: Dict[str, Any | Callable[[Any], Any]],
    ) -> Any:
        if section in defaults_json and key in defaults_json[section]:
            value = defaults_json[section][key]
            if isinstance(value, str) and "/" in value:
                num, den = value.split("/", 1)
                num = num.strip()
                den = den.strip()
                if num.isidentifier() and den.isidentifier():
                    return getattr(ctx, num) / getattr(ctx, den)
            return value
        default = fallback.get(key, SKIP)
        return default(ctx) if callable(default) else default

    def write_changed_params(self) -> None:
        """
        Compare current values against defaults and write only changed parameters.

        This is a placeholder for logic that will diff section attributes against
        their default values and emit a minimal input file update.
        """
        changed = self._changed_params()
        for key, value, default in changed:
            print(f"[{self._name}] writing {key} = {value} to input file")

    def show_changed_params(self) -> None:
        """
        Show changed parameters for this section.

        This is a placeholder for logic that will diff section attributes against
        their default values and print a summary.
        """
        changed = self._changed_params()
        if not changed:
            print(f"[{self._name}] no changes")
            return
        for key, value, default in changed:
            print(f"[{self._name}] {key}: {value} (default: {default})")

    def _changed_params(self) -> List[Tuple[str, Any, Any]]:
        changed = []
        for key, default in self._defaults.items():
            if not hasattr(self, key):
                continue
            value = getattr(self, key)
            if isinstance(value, np.ndarray) or isinstance(default, np.ndarray):
                same = np.array_equal(value, default)
            else:
                same = value == default
            if not same:
                changed.append((key, value, default))
        return changed

    def __repr__(self) -> str:
        lines = [f"{self._name}:"]
        for key, val in sorted(self.__dict__.items()):
            if key in ("_name", "sim", "_defaults"):
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
    DEFAULTS_JSON = None

    def __init__(self, expnr, path=None, load_geometry: bool = True):
        from udbase import UDBase

        if isinstance(expnr, UDBase):
            sim = expnr
        else:
            sim = UDBase(expnr, path, load_geometry=load_geometry)

        self.sim = sim
        self._section_spec_map = {spec.name: spec for spec in self.SECTION_SPECS}
        if self.DEFAULTS_JSON is None:
            self.DEFAULTS_JSON = Section.load_defaults_json()

        for spec in self.SECTION_SPECS:
            defaults = {}
            ctx = _DefaultContext(sim, defaults)
            for key in spec.fields:
                default = Section.resolve_default(
                    spec.name,
                    key,
                    ctx,
                    self.DEFAULTS_JSON,
                    spec.defaults,
                )
                if default is SKIP:
                    continue
                defaults[key] = default

            values = {}
            for key in spec.fields:
                if hasattr(sim, key):
                    values[key] = getattr(sim, key)
                elif key in defaults:
                    values[key] = defaults[key]
            setattr(
                self,
                spec.name,
                spec.section_cls(spec.name, values, sim=self.sim, defaults=defaults),
            )

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
            if key not in ("_name", "sim", "_defaults")
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

    def run_all(self, **kwargs: Any) -> None:
        """Run preprocessing in a write_inputs.m-style sequence."""
        self.ibm.run_all()
        if self.vegetation.ltrees or self.vegetation.ltreesfile:
            self.vegetation.run_all()
        self.ic.run_all()
        if self.scalars.nsv > 0:
            self.scalars.run_all()
        if self.seb.lEB:
            run_all = self.radiation.run_all
            sig = inspect.signature(run_all)
            params = sig.parameters.values()
            if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in params):
                run_all(**kwargs)
            else:
                call_kwargs = {k: v for k, v in kwargs.items() if k in sig.parameters}
                run_all(**call_kwargs)
            self.seb.run_all()

    def write_changed_params(self) -> None:
        """Write changed parameters for every section."""
        for spec in self.SECTION_SPECS:
            section_obj = getattr(self, spec.name, None)
            if section_obj is None:
                continue
            section_obj.write_changed_params()

    def show_changed_params(self) -> None:
        """Show changed parameters for every section."""
        for spec in self.SECTION_SPECS:
            section_obj = getattr(self, spec.name, None)
            if section_obj is None:
                continue
            section_obj.show_changed_params()

    def __str__(self) -> str:
        return self.__repr__()
