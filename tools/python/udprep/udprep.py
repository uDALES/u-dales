"""Core preprocessing framework for uDALES.

Provides the :class:`Section` base class that all preprocessing modules
(grid, IBM, radiation, vegetation, ...) inherit from, and the
:class:`UDPrep` orchestrator that discovers, configures, and runs
sections in order.  Section defaults are loaded from ``defaults.json``
and can be overridden per-case via namoptions or keyword arguments.
"""
from __future__ import annotations

import inspect
import logging
from pathlib import Path
from typing import Any, Dict

# Section base classes live in _section (imported here for backward
# compatibility and to break the section<->orchestrator import cycle).
from ._section import SKIP, Section, SectionSpec  # noqa: F401  (re-exported)

logger = logging.getLogger(__name__)


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
from .udprep_init import validate_expnr
from .udprep_bcs import SPEC as BCS_SPEC
from .udprep_grid import SPEC as GRID_SPEC
from .udprep_forcing import SPEC as FORCING_SPEC
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
        GRID_SPEC,
        FORCING_SPEC,
        BCS_SPEC,
        SCALARS_SPEC,
        VEGETATION_SPEC,
        IBM_SPEC,
        RADIATION_SPEC,
        SEB_SPEC,
    ]
    DEFAULTS_JSON = None

    def __init__(self, expnr, path=None, load_geometry: bool = True, suppress_load_warnings: bool = False):
        from udbase import UDBase

        if isinstance(expnr, UDBase):
            sim = expnr
        else:
            if path is None:
                if isinstance(expnr, (str, Path)):
                    path = Path(expnr)
                else:
                    raise TypeError(
                        f"expnr must be a string or Path pointing to a case directory, "
                        f"got {type(expnr).__name__!r}.\n"
                        f"Example: UDPrep('examples/101') or "
                        f"         UDPrep(101, 'examples/101') or "
                        f"         UDPrep('101', path='/path/to/cases/101')"
                    )
            expnr = validate_expnr(Path(path))
            sim = UDBase(expnr, path, load_geometry=load_geometry, suppress_load_warnings=suppress_load_warnings)

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
        # For proxy-style Section, fields are on sim. Use _fields to gather them.
        try:
            fields = object.__getattribute__(section_obj, "_fields")
        except AttributeError:
            fields = None
        if fields is not None:
            result: dict[str, Any] = {}
            for key in fields:
                try:
                    result[key] = getattr(section_obj, key)
                except AttributeError:
                    pass
            # Also include any locally-stored attributes (computed outputs, etc.)
            for key, val in section_obj.__dict__.items():
                if key in Section._MANAGED_INTERNALS or key.startswith("_") or key in fields:
                    continue
                result[key] = val
            return result
        # FakeSection or other non-proxy sections: fall back to __dict__ scan.
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
        """Run preprocessing for all sections."""
        ibm_backend = kwargs.pop("ibm_backend", "f2py")
        self.grid.run_all()
        self.forcing.run_all()
        if self.scalars.nsv > 0:
            self.scalars.run_all()
        if self.vegetation.ltrees:
            self.vegetation.run_all()
        if self.ibm.libm:
            self.ibm.run_all(backend=ibm_backend)
            if self.radiation.lEB:
                run_all = self.radiation.run_all
                sig = inspect.signature(run_all)
                params = sig.parameters.values()
                if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in params):
                    run_all(**kwargs)
                else:
                    call_kwargs = {k: v for k, v in kwargs.items() if k in sig.parameters}
                    run_all(**call_kwargs)
            # P7: SEB is gated independently of radiation EB. SEBSection.run_all
            # writes Tfacinit whenever lEB OR any of iwallmom/iwalltemp/
            # iwallmoist == 2, so the orchestrator must call it under the same
            # condition — otherwise a non-EB case with e.g. iwalltemp=2 never
            # gets Tfacinit.inp from the pipeline.
            if self._seb_needs_run():
                self.seb.run_all()

    def _seb_needs_run(self) -> bool:
        """Whether SEBSection.run_all would do work for the current case.

        Mirrors the gate inside :meth:`SEBSection.run_all` (reads the same
        sim-owned fields), so the orchestrator and the section stay aligned.
        """
        sim = self.sim
        return bool(
            getattr(sim, "lEB", False)
            or getattr(sim, "iwallmom", 1) == 2
            or getattr(sim, "iwalltemp", 1) == 2
            or getattr(sim, "iwallmoist", 1) == 2
        )

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
