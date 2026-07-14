"""Core preprocessing framework for uDALES.

Provides the :class:`Section` base class that all preprocessing modules
(grid, IBM, radiation, vegetation, ...) inherit from, and the
:class:`UDPrep` orchestrator that discovers, configures, and runs
sections in order.  Section defaults are loaded from ``defaults.json``
and can be overridden per-case via namoptions or keyword arguments.
"""
from __future__ import annotations

from dataclasses import dataclass
import inspect
import json
import logging
import time
import warnings
from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple, Type

import numpy as np

# Shared base utilities for preprocessing sections.
SKIP = object()


logger = logging.getLogger(__name__)


class Section:
    # Internal attributes always stored on the instance itself, never proxied to sim.
    _MANAGED_INTERNALS: frozenset = frozenset({"_name", "_defaults", "_fields", "sim"})
    _namelist_map_cache: Dict[str, str] | None = None  # per-instance lazy cache

    def __init__(
        self,
        name: str,
        values: Dict[str, Any],
        sim: Any | None = None,
        defaults: Dict[str, Any] | None = None,
    ):
        object.__setattr__(self, "_name", name)
        object.__setattr__(self, "sim", sim)
        object.__setattr__(self, "_defaults", defaults or {})
        object.__setattr__(self, "_fields", frozenset(values.keys()))
        if sim is not None:
            # Write all resolved field values to sim.
            # Values dict is already prioritised (namoptions > defaults) by UDPrep.__init__,
            # so we can write unconditionally — sim becomes the single source of truth.
            for key, val in values.items():
                setattr(sim, key, val)
        else:
            # No sim available — store fields locally (no proxy possible).
            for key, val in values.items():
                object.__setattr__(self, key, val)

    def __getattr__(self, name: str) -> Any:
        # Called only when normal lookup (instance __dict__, class) fails.
        # Proxy to sim so that section.attr and sim.attr are the same object.
        try:
            sim = object.__getattribute__(self, "sim")
        except AttributeError:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")
        if sim is not None and hasattr(sim, name):
            return getattr(sim, name)
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def __setattr__(self, name: str, value: Any) -> None:
        if name in Section._MANAGED_INTERNALS or name.startswith("_"):
            object.__setattr__(self, name, value)
            return
        try:
            fields = object.__getattribute__(self, "_fields")
            sim = object.__getattribute__(self, "sim")
        except AttributeError:
            object.__setattr__(self, name, value)
            return
        if name in fields:
            if sim is not None:
                # Write through to sim — sim is the single source of truth for fields.
                setattr(sim, name, value)
            else:
                object.__setattr__(self, name, value)
        else:
            # Computed output or new local attribute — store on the section instance.
            object.__setattr__(self, name, value)

    def run_steps(
        self,
        label: str,
        steps: List[Tuple[str, Callable[[], None]]],
        **kwargs: Any,
    ) -> None:
        _orig_formatwarning = warnings.formatwarning
        warnings.formatwarning = lambda msg, cat, *_a, **_kw: f"{cat.__name__}: {msg}\n"
        try:
            for name, func in steps:
                start = time.perf_counter()
                logger.info("[%s] %s...", label, name)
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
                logger.info("[%s] %s done in %.3f s", label, name, elapsed)
        finally:
            warnings.formatwarning = _orig_formatwarning

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

    @staticmethod
    def _load_namelist_map() -> Dict[str, str]:
        """Load variable -> namelist block mapping from namelists.json."""
        map_path = Path(__file__).resolve().parent.parent / "namelists.json"
        if not map_path.is_file():
            return {}
        try:
            data = json.loads(map_path.read_text(encoding="ascii"))
        except json.JSONDecodeError:
            return {}
        return {k.lower(): v for k, v in data.get("variables", {}).items()}

    def save_param(self, varname: str, value: Any) -> Path:
        """Update a namelist variable in namoptions.<expnr> and sync sim in memory."""
        if self._namelist_map_cache is None:
            self._namelist_map_cache = type(self)._load_namelist_map()
        namelist = self._namelist_map_cache.get(varname.lower(), "INP")

        def _format_value(val: Any) -> str:
            if isinstance(val, (list, tuple, np.ndarray)):
                arr = np.asarray(val)
                if arr.ndim == 0:
                    return _format_value(arr.item())
                return ", ".join(_format_value(v) for v in arr.ravel())
            if isinstance(val, (np.bool_, bool)):
                return ".true." if bool(val) else ".false."
            if isinstance(val, (np.integer, int)):
                return f"{int(val):d}"
            if isinstance(val, (np.floating, float)):
                s = f"{float(val):.6g}"
                if "e" not in s and "E" not in s and "." not in s:
                    s = f"{s}."
                return s
            if isinstance(val, str):
                trimmed = val
                if (trimmed.startswith("'") and trimmed.endswith("'")) or (
                    trimmed.startswith('"') and trimmed.endswith('"')
                ):
                    trimmed = trimmed[1:-1]
                return f"'{trimmed}'"
            return str(val)

        namelist_path = Path(self.path) / f"namoptions.{self.expnr}"
        if not namelist_path.is_file():
            raise FileNotFoundError(f"Missing {namelist_path}")

        lines = namelist_path.read_text(encoding="ascii").splitlines(keepends=True)
        nml_lower = namelist.strip().lower()
        var_lower = varname.strip().lower()
        value_str = _format_value(value)

        start_idx = None
        end_idx = None
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped.lower().startswith("&") and stripped[1:].strip().lower() == nml_lower:
                start_idx = i
                break

        if start_idx is not None:
            for i in range(start_idx + 1, len(lines)):
                if lines[i].strip().startswith("/"):
                    end_idx = i
                    break
            if end_idx is None:
                raise ValueError(
                    f"Namelist block '{namelist}' in {namelist_path} has no terminator '/'"
                )

            updated = False
            for i in range(start_idx + 1, end_idx):
                line = lines[i]
                if "=" not in line:
                    continue
                left, right = line.split("=", 1)
                if left.strip().lower() != var_lower:
                    continue
                comment = ""
                if "!" in right:
                    _, comment = right.split("!", 1)
                    comment = "!" + comment.rstrip("\n")
                newline = f"{left.rstrip()} = {value_str}"
                if comment:
                    newline = f"{newline} {comment}"
                lines[i] = newline + ("\n" if line.endswith("\n") else "")
                updated = True
                break

            if not updated:
                indent = "  "
                insert_line = f"{indent}{varname} = {value_str}\n"
                lines.insert(end_idx, insert_line)
        else:
            block_name = namelist.strip().upper()
            lines.append(f"&{block_name}\n")
            lines.append(f"  {varname} = {value_str}\n")
            lines.append("/\n")

        with namelist_path.open("w", encoding="ascii", newline="\n") as f:
            f.write("".join(lines))
        return namelist_path

    def write_changed_params(self) -> None:
        """
        Compare current values against defaults and write only changed parameters.
        """
        if self.sim is None:
            return
        changed = self._changed_params()
        for key, value, default in changed:
            if isinstance(value, np.ndarray):
                continue
            if isinstance(value, (list, tuple)) and any(
                isinstance(v, (list, tuple, np.ndarray, dict)) for v in value
            ):
                continue
            try:
                self.save_param(key, value)
            except Exception:
                logger.info("[%s] skipping writeback for %s = %r", self._name, key, value)

    def show_changed_params(self) -> None:
        """
        Show changed parameters for this section.

        This is a placeholder for logic that will diff section attributes against
        their default values and print a summary.
        """
        changed = self._changed_params()
        if not changed:
            logger.info("[%s] no changes", self._name)
            return
        for key, value, default in changed:
            logger.info("[%s] %s: %s (default: %s)", self._name, key, value, default)

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
        try:
            fields = object.__getattribute__(self, "_fields")
            sim = object.__getattribute__(self, "sim")
        except AttributeError:
            fields, sim = frozenset(), None
        shown: set = set()
        for key in sorted(fields):
            try:
                val = getattr(self, key)
            except AttributeError:
                continue
            lines.append(f"  {key}: {val!r}")
            shown.add(key)
        # Also show locally-stored attributes (computed outputs, etc.)
        for key, val in sorted(self.__dict__.items()):
            if key in Section._MANAGED_INTERNALS or key.startswith("_") or key in shown:
                continue
            lines.append(f"  {key}: {val!r}")
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
