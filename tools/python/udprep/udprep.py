from __future__ import annotations

from dataclasses import dataclass
import json
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Dict, List, Tuple, Type, Union
import warnings

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

    def run_steps(self, label: str, steps: List[Tuple[str, Callable[[], None]]]) -> None:
        for name, func in steps:
            start = time.perf_counter()
            print(f"[{label}] {name}...")
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
            # Handle expression strings like "xlen/itot" or "0.1*zsize"
            if isinstance(value, str):
                # Division: "xlen/itot"
                if "/" in value:
                    num, den = value.split("/", 1)
                    num = num.strip()
                    den = den.strip()
                    if num.isidentifier() and den.isidentifier():
                        num_val = getattr(ctx, num)
                        den_val = getattr(ctx, den)
                        return num_val / den_val
                # Multiplication: "0.1*zsize" or "zsize*0.1"
                elif "*" in value:
                    parts = value.split("*", 1)
                    left = parts[0].strip()
                    right = parts[1].strip()
                    # Try to evaluate both sides (could be number or identifier)
                    try:
                        left_val = float(left) if not left.isidentifier() else getattr(ctx, left)
                        right_val = float(right) if not right.isidentifier() else getattr(ctx, right)
                        return left_val * right_val
                    except (ValueError, AttributeError):
                        pass  # Return string as-is if evaluation fails
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
        # Check sim FIRST (namoptions values take precedence)
        if self._sim is not None and hasattr(self._sim, name):
            return getattr(self._sim, name)
        # Then check defaults being built up
        if name in self._values:
            return self._values[name]
        raise AttributeError(name)
from .udprep_ic import SPEC as IC_SPEC
from .udprep_ibm import SPEC as IBM_SPEC
from .udprep_grid import SPEC as GRID_SPEC
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
        GRID_SPEC,
        IC_SPEC,
        VEGETATION_SPEC,
        SCALARS_SPEC,
        RADIATION_SPEC,
    ]
    DEFAULTS_JSON = None

    def __init__(self, expdir: Union[Path, Any], load_geometry: bool = True, 
                 use_udbase: bool = False):
        """
        Initialize UDPrep for preprocessing.
        
        Parameters
        ----------
        expdir : Path or UDBase
            - If Path: experiment directory containing namoptions file
              (directory name must match experiment number)
            - If UDBase: existing UDBase object (backward compatibility)
        load_geometry : bool, default True
            Whether to load STL geometry
        use_udbase : bool, default False
            If True, use UDBase for initialization (legacy mode)
            If False, use lightweight namoptions parser (recommended)
            
        Examples
        --------
        # Recommended: lightweight mode
        >>> prep = UDPrep(Path('experiments/001'))
        
        # Legacy: use UDBase
        >>> prep = UDPrep(Path('experiments/001'), use_udbase=True)
        
        # Backward compatible: pass UDBase object
        >>> from udbase import UDBase
        >>> sim = UDBase('001', 'experiments/001')
        >>> prep = UDPrep(sim)
        """
        # Import here to avoid circular dependency
        try:
            from udbase import UDBase
            UDBASE_AVAILABLE = True
        except ImportError:
            UDBASE_AVAILABLE = False
            UDBase = type(None)  # Placeholder for type checking
        
        # Option 1: UDBase object passed directly (backward compatibility)
        if isinstance(expdir, UDBase):
            sim = expdir
        
        # Option 2: Use UDBase (legacy mode)
        elif use_udbase:
            if not UDBASE_AVAILABLE:
                raise ImportError("UDBase not available. Set use_udbase=False to use lightweight parser.")
            if not isinstance(expdir, Path):
                expdir = Path(expdir)
            # Extract expnr from directory for UDBase
            from .udprep_init import setup_paths_from_config
            expnr = setup_paths_from_config(expdir)
            sim = UDBase(expnr, str(expdir), load_geometry=load_geometry)
        
        # Option 3: Lightweight parser (recommended default)
        else:
            sim = self._create_sim_from_namoptions(expdir, load_geometry)
        
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

    def _create_sim_from_namoptions(self, expdir: Union[Path, str], 
                                      load_geometry: bool) -> SimpleNamespace:
        """
        Create lightweight sim object from namoptions file.
        
        This replaces UDBase for preprocessing, reading only namoptions and geometry.
        
        Parameters
        ----------
        expdir : Path or str
            Experiment directory containing namoptions file
        load_geometry : bool
            Whether to load STL geometry
        """
        from .udprep_readnamelist import read_namoptions
        from .udprep_init import setup_paths_from_config
        
        # Convert to Path and validate directory structure
        expdir = Path(expdir)
        expnr_str = setup_paths_from_config(expdir)
        
        # Read namoptions
        namoptions_path = expdir / f"namoptions.{expnr_str}"
        params = read_namoptions(namoptions_path)
        
        # Create sim object with all namoptions as attributes
        sim = SimpleNamespace(**params)
        sim.expnr = expnr_str
        sim.path = str(expdir)
        
        # Load geometry if STL file specified and loading requested
        if load_geometry and hasattr(sim, 'stl_file') and sim.stl_file:
            try:
                from udgeom import UDGeom
                stl_path = expdir / sim.stl_file
                if stl_path.exists():
                    sim.geom = UDGeom(path=expdir)
                    sim.geom.load(sim.stl_file)
                else:
                    warnings.warn(f"STL file not found: {stl_path}")
                    sim.geom = None
            except ImportError:
                warnings.warn("UDGeom not available. Geometry will not be loaded.")
                sim.geom = None
        else:
            sim.geom = None
        
        return sim
    
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
