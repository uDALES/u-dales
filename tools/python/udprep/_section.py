"""Base classes for uDALES preprocessing sections.

Kept separate from udprep.py so section modules (udprep_radiation, ...) can
import Section/SectionSpec without importing the UDPrep orchestrator, which
imports every section back — that mutual import was a latent cycle exposed by
making the package __init__ lazy.
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
            except FileNotFoundError as exc:
                # Benign: the case's namoptions.<expnr> is not on disk (e.g.
                # writeback invoked for a section whose namelist was never
                # written). Surface it visibly instead of hiding it at INFO, but
                # do not abort the rest of the writeback. Any other failure
                # (PermissionError, a malformed namelist raising ValueError, or a
                # genuine bug such as a TypeError) is a real problem that would
                # otherwise leave namoptions silently un-updated, so it
                # propagates.
                warnings.warn(
                    f"[{self._name}] skipping writeback for {key} = {value!r}: {exc}"
                )

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
