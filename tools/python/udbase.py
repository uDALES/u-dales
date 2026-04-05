from __future__ import annotations

"""
uDALES Post-Processing Module

Python implementation of the MATLAB udbase class for analyzing
uDALES simulation outputs.

Aug 2024, Maarten van Reeuwijk, Jingzi Huang. first version
Dec 2024, Maarten van Reeuwijk, Chris Wilson. added facet functionality.
Oct 2025, Maarten van Reeuwijk, Jingzi Huang, Dipanjan Majumdar. Major upgrade of facet functionality and tree handling.

Copyright (C) 2016- the uDALES Team.
"""

import numpy as np
import xarray as xr
import json
from pathlib import Path
from typing import Optional, Union, Dict, Any, List
import importlib.util
import warnings
import sys

# Import UDGeom from the udgeom package
try:
    from udgeom import UDGeom
except ImportError:
    from .udgeom import UDGeom

try:
    from udvis import UDVis
except ImportError:
    from .udvis import UDVis


def _file_has_data(path: Path, skiprows: int = 0) -> bool:
    try:
        with path.open("r", encoding="ascii", errors="ignore") as f:
            for _ in range(skiprows):
                next(f, None)
            for line in f:
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    return True
    except OSError:
        return False
    return False


def _empty_2d(dtype: np.dtype) -> np.ndarray:
    return np.empty((0, 0), dtype=dtype)


class UDBase:
    """
    Post-processing class for uDALES simulations.
    
    This class provides methods to load and analyze uDALES simulation outputs,
    including field data, statistics, facet data, and geometry.
    
    Parameters
    ----------
    expnr : int or str
        Experiment number
    path : str or Path, optional
        Path to experiment directory. Defaults to current directory.
    load_geometry : bool, optional
        If True, load STL geometry when available.
    
    Examples
    --------
    >>> sim = UDBase(expnr=65, path='experiments/065')
    >>> u_field = sim.load_field('u')
    >>> stats = sim.load_stat_xyt('u')
    """
    
    def __init__(
        self,
        expnr: Union[int, str],
        path: Optional[Union[str, Path]] = None,
        load_geometry: bool = True,
        suppress_load_warnings: bool = False,
    ):
        """
        Initialize a UDBase instance.

        Parameters
        ----------
        expnr : int or str
            Experiment number. Converted to a zero-padded 3-digit string.
        path : str or Path, optional
            Path to the experiment directory. Defaults to current working directory.
        load_geometry : bool, optional
            If True, load STL geometry when available.
        suppress_load_warnings : bool, optional
            If True, suppress missing-file warnings during constructor loading.

        Attributes
        ----------
        expnr : str
            Experiment number as a zero-padded string.
        path : Path
            Base path to the experiment directory.
        geom : UDGeom or None
            Loaded geometry instance when STL is present and load_geometry is True.
        facs : dict
            Facet data (e.g., area, typeid, normals) if available.
        factypes : dict
            Facet type properties if available.
        facsec : dict
            Facet section data if available.
        trees : ndarray or None
            Tree index data if available.

        Raises
        ------
        ImportError
            If trimesh is not installed.
        FileNotFoundError
            If the namoptions file is missing.
        RuntimeError
            If geometry loading fails when requested and STL is present.
        """

        if importlib.util.find_spec("trimesh") is None:
            raise ImportError("trimesh is required for UDBase. Install with: pip install trimesh")

        # Store experiment number
        self.expnr = f"{int(expnr):03d}"
        
        # Set paths
        self.cpath = Path.cwd()
        if callable(path):
            path = path()
        self.path = Path(path) if path else self.cpath

        # Standard filenames and prefixes used across loaders, kept aligned
        # with the legacy MATLAB udbase defaults.
        self.fnamoptions = "namoptions"
        self.fprof = "prof.inp"
        self.fxytdump = "xytdump"
        self.ftdump = "tdump"
        self.ffielddump = "fielddump"
        self.fislicedump = "islicedump"
        self.fjslicedump = "jslicedump"
        self.fkslicedump = "kslicedump"
        self.fsolid = "solid"
        self.ffacEB = "facEB"
        self.ffacT = "facT"
        self.ffac = "fac"
        self.ffacets = "facets.inp"
        self.ffactypes = "factypes.inp"
        self.ffacetarea = "facetarea.inp"
        self.ffluid_boundary = "fluid_boundary"
        self.ffacet_sections = "facet_sections"
        self.ftreedump = "treedump"
        self.ftrees = "trees.inp"

        # Load-time warning control
        self._suppress_load_warnings = suppress_load_warnings

        # File presence flags
        self._lfprof = True
        self._lfsolid = True
        self._lffacets = True
        self._lffactypes = True
        self._lffacetarea = True
        self._lffacet_sections = True
        self._lfgeom = True
        self._lftrees = True
        
        # Read namoptions file
        self._read_namoptions()
        self._namelist_map = self._load_namelist_map()
        
        # Load grid
        self._load_grid()
        
        # Load geometry if present (can be disabled for faster startup)
        if load_geometry:
            self._load_geometry()
        
        # Load solid masks if present
        self._load_solid_masks()
        
        # Load facet data if present
        self._load_facet_data()
        
        # Load tree data if present
        self._load_tree_data()

        # Load vegetation data if present
        self._load_veg_data()

        # Visualization facade. `UDBase` owns the simulation state; `self.vis`
        # provides plotting methods on top of that state.
        self.vis = UDVis(self)
    
    def _read_namoptions(self):
        """
        Read namoptions file and store all parameters as attributes.
        
        Parses the namoptions.XXX file and converts values:
        - .true./.false. -> bool
        - Numbers -> int or float
        - Strings -> str (with quotes removed)
        """
        filepath = self.path / f"namoptions.{self.expnr}"
        
        if not filepath.exists():
            raise FileNotFoundError(f"namoptions.{self.expnr} not found in {self.path}")
        
        with open(filepath, 'r') as f:
            for line in f:
                # Skip comments and namelist headers
                if line.strip().startswith('&') or line.strip().startswith('!') or '=' not in line:
                    continue
                
                # Parse key=value pairs
                if '=' in line:
                    parts = line.split('=', 1)
                    key = parts[0].strip()
                    val = parts[1].split('!')[0].strip()  # Remove inline comments
                    
                    # Convert value types
                    if val.lower() == '.true.':
                        val = True
                    elif val.lower() == '.false.':
                        val = False
                    else:
                        # Try to parse as number
                        try:
                            if '.' in val or 'e' in val.lower():
                                val = float(val)
                            else:
                                val = int(val)
                        except ValueError:
                            # String value - remove quotes
                            val = val.strip("'\"")
                    
                    # Store as attribute
                    setattr(self, key, val)
        
        # Compute derived values
        if hasattr(self, 'xlen') and hasattr(self, 'itot'):
            self.dx = self.xlen / self.itot
        if hasattr(self, 'ylen') and hasattr(self, 'jtot'):
            self.dy = self.ylen / self.jtot

    def _load_namelist_map(self) -> Dict[str, str]:
        """Load variable->namelist mapping from namelists.json."""
        map_path = Path(__file__).resolve().parent / "namelists.json"
        if not map_path.is_file():
            return {}
        try:
            data = json.loads(map_path.read_text(encoding="ascii"))
        except json.JSONDecodeError:
            return {}
        return {k.lower(): v for k, v in data.get("variables", {}).items()}

    def _load_sparse_file(
        self,
        path: Path,
        *,
        skiprows: int = 0,
        dtype: np.dtype = float,
        min_cols: int | None = None,
        zero_based_cols: list[int] | None = None,
    ) -> np.ndarray:
        """Load a sparse text file with comments, returning a 2D array."""
        if not path.exists():
            return _empty_2d(dtype)
        rows = []
        with path.open("r", encoding="ascii", errors="ignore") as f:
            for _ in range(skiprows):
                next(f, None)
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = stripped.split()
                if min_cols is not None and len(parts) < min_cols:
                    raise ValueError(f"Expected at least {min_cols} columns in {path}")
                rows.append(parts)
        if not rows:
            return _empty_2d(dtype)
        arr = np.asarray(rows, dtype=dtype)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        if min_cols is not None and arr.shape[1] < min_cols:
            raise ValueError(f"Expected at least {min_cols} columns in {path}")
        if zero_based_cols:
            arr[:, zero_based_cols] -= 1
        return arr
    
    def _load_grid(self):
        """
        Load grid coordinates.
        
        Attempts to load from prof.inp.XXX file. If not present,
        generates uniform grid based on domain size and number of cells.
        """
        prof_file = self.path / f"prof.inp.{self.expnr}"
        
        if prof_file.exists():
            try:
                # Load prof.inp - skip header lines
                data = np.loadtxt(prof_file, skiprows=1)
                
                # Cell centers (zt)
                self.zt = data[:, 0]
                
                # Cell edges (zm)
                zm_full = np.concatenate([[0], 0.5 * (self.zt[:-1] + self.zt[1:]), [self.zsize]])
                self.zm = zm_full[:-1]  # Match dimensions
                
                # Grid spacing
                self.dzt = np.diff(np.concatenate([self.zm, [self.zsize]]))
                self.dzm = np.concatenate([[2 * self.zt[0]], np.diff(self.zt)])
                
            except Exception as e:
                self._warn_load(
                    f"Error loading prof.inp.{self.expnr}: {e}. Using uniform grid."
                )
                self._lfprof = False
                self._generate_uniform_zgrid()
        else:
            self._warn_load(f"prof.inp.{self.expnr} not found. Assuming equidistant grid.")
            self._lfprof = False
            self._generate_uniform_zgrid()
        
        # Generate x and y grids
        self.xm = np.arange(self.itot) * self.dx
        self.xt = self.xm + 0.5 * self.dx
        
        self.ym = np.arange(self.jtot) * self.dy
        self.yt = self.ym + 0.5 * self.dy
    
    def _generate_uniform_zgrid(self):
        """Generate uniform z-grid when prof.inp is not available."""
        if hasattr(self, 'zsize') and hasattr(self, 'ktot'):
            dz = self.zsize / self.ktot
            
            self.zm = np.arange(self.ktot) * dz
            self.zt = self.zm + 0.5 * dz
            
            self.dzm = np.full(self.ktot, dz)
            self.dzt = np.full(self.ktot, dz)
        else:
            self._warn_load("Cannot generate z-grid: zsize or ktot not found in namoptions")

    def _warn_load(self, message: str) -> None:
        if self._suppress_load_warnings:
            return
        print("=" * 67, file=sys.stderr)
        print(f"WARNING: {message}", file=sys.stderr)
        print("=" * 67, file=sys.stderr)
    
    def _load_geometry(self):
        """Load STL geometry file if present."""
        if hasattr(self, 'stl_file'):
            stl_path = self.path / self.stl_file
            if stl_path.exists():
                try:
                    self.geom = UDGeom(self.path)
                    self.geom.load(self.stl_file)
                except Exception as e:
                    raise RuntimeError(f"Error loading geometry from {stl_path}: {e}") from e
            else:
                raise FileNotFoundError(f"STL file not found: {stl_path}")
        else:
            self.geom = None
    
    def _load_solid_masks(self):
        """Load solid point arrays for u, v, w, c grids."""
        for grid_type in ['u', 'v', 'w', 'c']:
            solid_file = self.path / f"solid_{grid_type}.txt"
            
            if solid_file.exists():
                try:
                    indices = self._load_sparse_file(
                        solid_file,
                        skiprows=1,
                        dtype=int,
                        min_cols=3,
                        zero_based_cols=[0, 1, 2],
                    )

                    mask = np.zeros((self.itot, self.jtot, self.ktot), dtype=bool)
                    if indices.size:
                        mask[indices[:, 0], indices[:, 1], indices[:, 2]] = True
                    
                    # Store as attribute (Su, Sv, Sw, Sc)
                    setattr(self, f'S{grid_type}', mask)
                    
                except Exception as e:
                    self._warn_load(f"Error loading solid_{grid_type}.txt: {e}")
                    setattr(self, f'S{grid_type}', None)
                    self._lfsolid = False
            else:
                setattr(self, f'S{grid_type}', None)
    
    def _load_facet_data(self):
        """Load facet information if present."""
        self.facs = {}
        self.factypes = {}
        self._lffacetarea = True
        self._lffacets = True
        self._lffactypes = True
        self._lffacet_sections = True
        
        # Load facet areas
        facetarea_file = self.path / f"facetarea.inp.{self.expnr}"
        if facetarea_file.exists():
            try:
                self.facs['area'] = np.loadtxt(facetarea_file, skiprows=1)
            except Exception as e:
                self._warn_load(f"Error loading facetarea.inp.{self.expnr}: {e}")
                self._lffacetarea = False
        else:
            self._lffacetarea = False
        
        # Load facet types
        facets_file = self.path / f"facets.inp.{self.expnr}"
        if facets_file.exists():
            try:
                data = np.loadtxt(facets_file, skiprows=1)
                self.facs['typeid'] = data[:, 0].astype(int)
                self.facs['normals'] = data[:, 1:4]  # Surface normals
            except Exception as e:
                self._warn_load(f"Error loading facets.inp.{self.expnr}: {e}")
                self._lffacets = False
        else:
            self._lffacets = False
        
        # Load facet type properties
        factypes_file = self.path / f"factypes.inp.{self.expnr}"
        if factypes_file.exists():
            try:
                data = np.loadtxt(factypes_file, skiprows=3)
                
                if not hasattr(self, 'nfaclyrs'):
                    self.nfaclyrs = 3  # Default
                
                self.factypes['id'] = data[:, 0].astype(int)
                self.factypes['lGR'] = data[:, 1].astype(bool)
                self.factypes['z0'] = data[:, 2]
                self.factypes['z0h'] = data[:, 3]
                self.factypes['al'] = data[:, 4]  # Albedo
                self.factypes['em'] = data[:, 5]  # Emissivity
                self.factypes['d'] = data[:, 6:6 + self.nfaclyrs]  # Layer thicknesses
                self.factypes['C'] = data[:, 6 + self.nfaclyrs:6 + 2 * self.nfaclyrs]  # Heat capacity
                self.factypes['lam'] = data[:, 6 + 2 * self.nfaclyrs:6 + 3 * self.nfaclyrs + 1]  # Conductivity
                
                # Add names for common wall types
                wall_names = {
                    0: "Default dummy",
                    -1: "Asphalt floor",
                    -101: "Concrete bounding wall",
                    1: "Concrete",
                    2: "Brick",
                    3: "Stone",
                    4: "Painted wood",
                    11: "Green 1",
                    12: "Green 2"
                }
                
                self.factypes['name'] = [
                    wall_names.get(int(id_), f"Custom walltype {int(id_)}")
                    for id_ in self.factypes['id']
                ]
                
            except Exception as e:
                self._warn_load(f"Error loading factypes.inp.{self.expnr}: {e}")
                self._lffactypes = False
        else:
            self._lffactypes = False
        
        # Load facet sections
        self._load_facet_sections()
    
    def _load_facet_sections(self):
        """Load facet section information for u, v, w, c grids."""
        self.facsec = {}

        for grid_type in ['u', 'v', 'w', 'c']:
            facsec_file = self.path / f"facet_sections_{grid_type}.txt"
            fluid_boundary_file = self.path / f"fluid_boundary_{grid_type}.txt"
            
            if facsec_file.exists() and fluid_boundary_file.exists():
                try:
                    facsec_data = self._load_sparse_file(
                        facsec_file,
                        skiprows=1,
                        dtype=float,
                        min_cols=4,
                    )
                    fluid_boundary = self._load_sparse_file(
                        fluid_boundary_file,
                        skiprows=1,
                        dtype=int,
                        min_cols=3,
                        zero_based_cols=[0, 1, 2],
                    )
                    if facsec_data.size == 0 or fluid_boundary.size == 0:
                        self._lffacet_sections = False
                        continue
                    
                    # Store in structure
                    self.facsec[grid_type] = {
                        'facid': facsec_data[:, 0].astype(int) - 1,
                        'area': facsec_data[:, 1],
                        'locs': fluid_boundary[facsec_data[:, 2].astype(int) - 1, :].astype(int),
                        'distance': facsec_data[:, 3]
                    }
                    
                except Exception as e:
                    self._warn_load(f"Error loading facet_sections_{grid_type}.txt: {e}")
                    self._lffacet_sections = False
            else:
                self._lffacet_sections = False
    
    def _load_tree_data(self):
        """Load tree bounding box information if present."""
        trees_file = self.path / f"trees.inp.{self.expnr}"
        
        if trees_file.exists():
            try:
                # Load tree data (skip first 2 header lines)
                self.trees = np.loadtxt(trees_file, skiprows=2, dtype=int)
                if self.trees.ndim == 1:
                    self.trees = self.trees.reshape(1, -1)
                self.trees = self.trees.astype(int) - 1
            except Exception as e:
                self._warn_load(f"Error loading trees.inp.{self.expnr}: {e}")
                self._lftrees = False
                self.trees = None
        else:
            self._lftrees = False
            self.trees = None

    def _load_veg_data(self):
        """Load vegetation sparse data if present."""
        points_path = self.path / f"veg.inp.{self.expnr}"
        params_path = self.path / f"veg_params.inp.{self.expnr}"

        if not points_path.exists() or not params_path.exists():
            self.veg = None
            return

        try:
            self.load_veg(zero_based=True, cache=True)
        except Exception as exc:
            if not self._suppress_load_warnings:
                warnings.warn(f"Error loading vegetation data: {exc}")
            self.veg = None

    def save_param(self, varname: str, value: Any) -> Path:
        """Update a namelist variable in namoptions.<id> using a lookup map."""
        if not hasattr(self, "_namelist_map") or self._namelist_map is None:
            self._namelist_map = self._load_namelist_map()
        namelist = self._namelist_map.get(varname.lower(), "INP")

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
                raise ValueError(f"Namelist block '{namelist}' in {namelist_path} has no terminator '/'")

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
        setattr(self, varname.strip(), value)
        return namelist_path

    def save_veg(
        self,
        points: np.ndarray,
        ids: np.ndarray,
        lad_values: np.ndarray,
        cd: float,
        ud: float,
        dec: float,
        lsize: float,
        r_s: float,
        *,
        write_ids: bool = False,
    ) -> Dict[str, Path]:
        """Write vegetation sparse inputs (veg.inp, veg_params, optional veg_id)."""
        if points.ndim != 2 or points.shape[1] < 3:
            raise ValueError("points must be an (n, 3) array of 1-based indices")
        if len(points) != len(lad_values) or len(points) != len(ids):
            raise ValueError("points, ids, and lad_values must be the same length")

        sim_dir = Path(self.path)
        sim_id = self.expnr
        out_paths: Dict[str, Path] = {}

        veg_path = sim_dir / f"veg.inp.{sim_id}"
        with veg_path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# position (i,j,k)\n")
            for i, j, k in points:
                f.write(f"{int(i):7d} {int(j):7d} {int(k):7d}\n")
        out_paths["veg"] = veg_path

        params_path = sim_dir / f"veg_params.inp.{sim_id}"
        with params_path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# id lad cd ud dec lsize r_s\n")
            for bid, lad_val in zip(ids, lad_values):
                f.write(
                    f"{int(bid):7d} {float(lad_val):12.6f} {cd:12.6f} {ud:12.6f} "
                    f"{dec:12.6f} {lsize:12.6f} {r_s:12.6f}\n"
                )
        out_paths["params"] = params_path

        if write_ids:
            ids_path = sim_dir / f"veg_id.inp.{sim_id}"
            with ids_path.open("w", encoding="ascii", newline="\n") as f:
                f.write("# block_id for each point in veg.inp, same order\n")
                for block_id in ids:
                    f.write(f"{int(block_id):7d}\n")
            out_paths["ids"] = ids_path

        points = np.asarray(points)
        ids = np.asarray(ids)
        lad_values = np.asarray(lad_values)
        points_zero = points[:, :3].astype(int, copy=False) - 1
        self.veg = {
            "points": points_zero,
            "params": {
                "id": ids.astype(int, copy=False),
                "lad": lad_values.astype(float, copy=False),
                "cd": np.full(len(ids), cd, dtype=float),
                "ud": np.full(len(ids), ud, dtype=float),
                "dec": np.full(len(ids), dec, dtype=float),
                "lsize": np.full(len(ids), lsize, dtype=float),
                "r_s": np.full(len(ids), r_s, dtype=float),
            },
        }
        self._lftrees = True

        return out_paths

    def save_trees(
        self,
        points: np.ndarray,
        ids: np.ndarray,
        lad_values: np.ndarray,
        cd: float,
        ud: float,
        dec: float,
        lsize: float,
        r_s: float,
        *,
        write_ids: bool = False,
    ) -> Dict[str, Path]:
        """Backward-compatible alias for save_veg."""
        return self.save_veg(points, ids, lad_values, cd, ud, dec, lsize, r_s, write_ids=write_ids)

    def load_veg(self, *, zero_based: bool = True, cache: bool = True) -> Dict[str, Any]:
        """Load vegetation sparse points and parameters."""
        if cache and hasattr(self, "veg") and self.veg is not None:
            return self.veg

        points_path = self.path / f"veg.inp.{self.expnr}"
        params_path = self.path / f"veg_params.inp.{self.expnr}"
        point_cols = [0, 1, 2] if zero_based else None
        points = self._load_sparse_file(
            points_path,
            dtype=int,
            min_cols=3,
            zero_based_cols=point_cols,
        )
        params = self._load_sparse_file(
            params_path,
            dtype=float,
            min_cols=7,
        )

        if points.size == 0 and params.size == 0:
            veg = {
                "points": np.empty((0, 3), dtype=int),
                "params": {
                    "id": np.empty((0,), dtype=int),
                    "lad": np.empty((0,), dtype=float),
                    "cd": np.empty((0,), dtype=float),
                    "ud": np.empty((0,), dtype=float),
                    "dec": np.empty((0,), dtype=float),
                    "lsize": np.empty((0,), dtype=float),
                    "r_s": np.empty((0,), dtype=float),
                },
            }
            if cache:
                self.veg = veg
            return veg

        if points.size == 0 or params.size == 0:
            raise ValueError("veg.inp and veg_params.inp must both be present and non-empty")

        points = points[:, :3].astype(int, copy=False)
        params = params[:, :7]
        if len(points) != len(params):
            raise ValueError(
                f"veg.inp.{self.expnr} has {len(points)} points but veg_params has {len(params)} rows"
            )

        veg = {
            "points": points,
            "params": {
                "id": params[:, 0].astype(int),
                "lad": params[:, 1],
                "cd": params[:, 2],
                "ud": params[:, 3],
                "dec": params[:, 4],
                "lsize": params[:, 5],
                "r_s": params[:, 6],
            },
        }
        if cache:
            self.veg = veg
        return veg
    
    def __repr__(self):
        """String representation of UDBase object."""
        info = [
            f"UDBase(expnr='{self.expnr}')",
            f"  path: {self.path}",
        ]
        
        if hasattr(self, 'itot'):
            info.append(f"  grid: {self.itot} x {self.jtot} x {self.ktot}")
        
        if hasattr(self, 'xlen'):
            info.append(f"  domain: {self.xlen} x {self.ylen} x {self.zsize}")
        
        if self.geom is not None:
            info.append(f"  geometry: loaded")
        
        if hasattr(self, 'nfcts'):
            info.append(f"  facets: {self.nfcts}")

        scalar_types = (int, float, bool, str, np.integer, np.floating)
        if self.__dict__:
            for key, val in sorted(self.__dict__.items()):
                if isinstance(val, scalar_types):
                    info.append(f"  {key}: {val}")
                elif isinstance(val, np.ndarray):
                    info.append(f"  {key}: ndarray[{val.dtype}] shape={val.shape}")
                else:
                    info.append(f"  {key}: {type(val).__name__}")

        return "\n".join(info)
    
    # ===== Data Loading Methods =====
    
    def load_field(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load 3D instantaneous field data.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_field()  # Display available variables
        >>> u = sim.load_field('u')  # Returns numpy array
        >>> print(u.shape)
        """
        filename = self.path / f"fielddump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_stat_xyt(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load slab-averaged (xy-plane) statistics.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_stat_xyt()  # Display available variables
        >>> u_avg = sim.load_stat_xyt('u')  # Returns numpy array
        """
        filename = self.path / f"xytdump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_stat_t(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load time-averaged 3D statistics.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_stat_t()  # Display available variables
        >>> u_tavg = sim.load_stat_t('u')  # Returns numpy array
        """
        filename = self.path / f"tdump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_stat_tree(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load time-averaged statistics of tree source terms.
        
        Retrieves tree drag, heat, and moisture source terms from the treedump file.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_stat_tree()  # Display available variables
        >>> tree_drag = sim.load_stat_tree('tree_drag_u')  # Returns numpy array
        """
        filename = self.path / f"treedump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_slice(self, plane: str, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load 2D slice data.
        
        Parameters
        ----------
        plane : str
            Slice plane: 'i', 'j', or 'k' for x, y, or z slices
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_slice('k')  # Display available variables in horizontal slice
        >>> u_slice = sim.load_slice('k', 'u')  # Returns numpy array
        """
        if plane not in ['i', 'j', 'k']:
            raise ValueError("plane must be 'i', 'j', or 'k'")
        
        filename = self.path / f"{plane}slicedump.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def _load_ncdata(self, filename: Path, var: Optional[str]) -> Union[xr.Dataset, np.ndarray]:
        """
        Helper method to load NetCDF data using xarray.
        
        Automatically reverses dimension order to match MATLAB conventions.
        When a specific variable is requested, returns it as a numpy array for
        memory efficiency with large datasets.
        
        NetCDF files use C-style (row-major) dimension ordering, while MATLAB
        uses Fortran-style (column-major) ordering. Due to this fundamental
        difference, MATLAB automatically reverses dimension order when reading
        NetCDF files. This method reverses dimensions to match MATLAB's behavior.
        
        Examples of transformations:
        - (time, zt) -> (zt, time)
        - (time, zt, yt, xm) -> (xm, yt, zt, time)
        - (time, fct) -> (fct, time)
        - (time, lyr, fct) -> (fct, lyr, time)
        
        Parameters
        ----------
        filename : Path
            Path to NetCDF file
        var : str, optional
            Variable to extract. If None, displays available variables.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), numpy array if var is specified.
        """
        if not filename.exists():
            raise FileNotFoundError(f"File not found: {filename}")
        
        ds = xr.open_dataset(filename)
        
        # Transpose all data variables to match MATLAB's column-major convention
        # Reverse dimension order for all variables with 2+ dimensions
        transposed_vars = {}
        for var_name in ds.data_vars:
            data_var = ds[var_name]
            if len(data_var.dims) >= 2:
                transposed_vars[var_name] = data_var.transpose(*reversed(data_var.dims))
            else:
                transposed_vars[var_name] = data_var
        
        # Create new dataset with transposed variables
        ds_transposed = xr.Dataset(transposed_vars, coords=ds.coords, attrs=ds.attrs)
        
        if var is None:
            # Display file contents
            self._display_ncinfo(ds_transposed, filename.name)
            return ds_transposed
        else:
            if var not in ds_transposed:
                raise KeyError(f"Variable '{var}' not found in {filename.name}")
            
            # Return as numpy array for consistency with MATLAB and memory efficiency
            return ds_transposed[var].values
    
    def _display_ncinfo(self, ds: xr.Dataset, filename: str):
        """Display information about NetCDF dataset."""
        print(f"\nContents of {filename}:")
        print(f"{'Name':<20} {'Dimensions':<30} {'Shape':<20}")
        print("-" * 70)
        
        for var_name in sorted(ds.data_vars):
            var = ds[var_name]
            dims = ', '.join(var.dims)
            shape = ' x '.join(map(str, var.shape))
            print(f"{var_name:<20} {dims:<30} {shape:<20}")
    
    # ===== Facet Data Loading Methods =====
    
    def load_fac_momentum(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load facet momentum data (pressure and shear stress).
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_fac_momentum()  # Display available variables
        >>> pressure = sim.load_fac_momentum('pres')  # Returns numpy array
        """
        filename = self.path / f"fac.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_fac_eb(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load facet surface energy balance data.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_fac_eb()  # Display available variables
        >>> H = sim.load_fac_eb('hf')  # Sensible heat flux - Returns numpy array
        >>> K = sim.load_fac_eb('netsw')  # Net shortwave - Returns numpy array
        """
        filename = self.path / f"facEB.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_fac_temperature(self, var: Optional[str] = None) -> Union[xr.Dataset, np.ndarray]:
        """
        Load facet temperature data.
        
        Parameters
        ----------
        var : str, optional
            Variable name to load. If None, displays available variables and returns full dataset.
        
        Returns
        -------
        xarray.Dataset or numpy.ndarray
            Full dataset if var is None (for browsing), otherwise numpy array of the requested variable.
        
        Examples
        --------
        >>> sim.load_fac_temperature()  # Display available variables
        >>> T = sim.load_fac_temperature('T')  # Temperature in layers - Returns numpy array
        >>> dTdz = sim.load_fac_temperature('dTdz')  # Temperature gradient - Returns numpy array
        """
        filename = self.path / f"facT.{self.expnr}.nc"
        return self._load_ncdata(filename, var)
    
    def load_seb(self) -> Dict[str, np.ndarray]:
        """
        Load all surface energy balance terms.
        
        Returns
        -------
        dict
            Dictionary containing all SEB terms:
            - Kstar: Net shortwave radiation
            - Lstar: Net longwave radiation
            - Lin: Incoming longwave
            - Lout: Outgoing longwave
            - H: Sensible heat flux
            - E: Latent heat flux
            - G: Ground heat flux
            - Tsurf: Surface temperature
            - t: Time array
        
        Examples
        --------
        >>> seb = sim.load_seb()
        >>> print(seb['Kstar'].shape)
        >>> seb_avg = sim.area_average_seb(seb)
        """
        # Load energy balance terms
        try:
            t = self.load_fac_eb('t')
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Surface energy balance file facEB.{self.expnr}.nc not found in {self.path}"
            )
        K = self.load_fac_eb('netsw')
        Lin = self.load_fac_eb('LWin')
        Lout = self.load_fac_eb('LWout')
        H = self.load_fac_eb('hf')
        E = self.load_fac_eb('ef')
        
        # Load temperature data for ground heat flux calculation
        T = self.load_fac_temperature('T')
        dTdz = self.load_fac_temperature('dTdz')
        lam = self.assign_prop_to_fac('lam')
        
        # Calculate derived quantities
        L = Lin - Lout
        # Data now comes as (fct, lyr, time) from _load_ncdata, transposed to match MATLAB
        # lam is (n_facets, n_layers), dTdz is (n_facets, n_layers, n_time)
        G = -lam[:, 0, np.newaxis] * dTdz[:, 0, :]  # Ground heat flux: (n_facets, n_time)
        Tsurf = T[:, 0, :]  # Surface temperature: (n_facets, n_time)
        
        # All data is already in MATLAB convention: (n_facets, n_time)
        # _load_ncdata transposes from NetCDF (time, fct) to (fct, time)
        return {
            'Kstar': K,
            'Lstar': L,
            'Lin': Lin,
            'Lout': Lout,
            'H': -H,  # Sign convention
            'E': -E,  # Sign convention
            'G': G,
            'Tsurf': Tsurf,
            't': t
        }
    
    # ===== Facet Analysis Methods =====
    
    def assign_prop_to_fac(self, prop: str) -> np.ndarray:
        """
        Assign facet type property to individual facets.
        
        Parameters
        ----------
        prop : str
            Property name from factypes (e.g., 'al', 'em', 'z0', 'lam')
        
        Returns
        -------
        ndarray
            Property values for each facet
        
        Examples
        --------
        >>> albedo = sim.assign_prop_to_fac('al')
        >>> emissivity = sim.assign_prop_to_fac('em')
        """
        if not self._lffactypes or not self._lffacets:
            # Preprocessing may have generated facet metadata after UDBase was
            # initialized. Only fall back to re-reading from disk when the
            # required data is currently missing in memory.
            self._load_facet_data()

        if not self._lffactypes:
            raise ValueError(f"factypes.inp.{self.expnr} required for this method")
        
        if not self._lffacets:
            raise ValueError(f"facets.inp.{self.expnr} required for this method")
        
        if prop not in self.factypes:
            raise KeyError(f"Property '{prop}' not found in factypes")
        
        # Get property values and type IDs
        prop_vals = self.factypes[prop]
        type_ids = self.factypes['id']
        facet_types = self.facs['typeid']
        
        # Map type IDs to property values
        if prop_vals.ndim == 1:
            # 1D property (e.g., albedo)
            result = np.full_like(facet_types, np.nan, dtype=float)
            for i, tid in enumerate(type_ids):
                result[facet_types == tid] = prop_vals[i]
        else:
            # 2D property (e.g., layer thicknesses)
            result = np.full((len(facet_types), prop_vals.shape[1]), np.nan)
            for i, tid in enumerate(type_ids):
                result[facet_types == tid, :] = prop_vals[i, :]
        
        return result
    
    def area_average_fac(self, var: np.ndarray, sel: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Area-weighted average of facet data.
        
        Parameters
        ----------
        var : ndarray
            Facet variable to average. Can be:
            - 1D array with shape (n_facets,)
            - 2D array with shape (n_facets, n_times) or (n_times, n_facets)
            Method auto-detects which dimension is facets based on size.
        sel : ndarray, optional
            Boolean mask or indices to select subset of facets
        
        Returns
        -------
        ndarray
            Area-averaged values. Same shape as input except facet dimension is removed.
        
        Examples
        --------
        >>> H_avg = sim.area_average_fac(H)
        >>> # Average over specific facets
        >>> roof_facets = sim.facs['typeid'] == 1
        >>> H_avg_roof = sim.area_average_fac(H, roof_facets)
        """
        if not self._lffacetarea:
            raise ValueError(f"facetarea.inp.{self.expnr} required for this method")
        
        areas = self.facs['area']
        n_facets = len(areas)
        
        if sel is None:
            sel = slice(None)
            sel_areas = areas
        else:
            sel_areas = areas[sel]
        
        # Handle different array dimensions
        if var.ndim == 1:
            # Simple 1D case
            total_var = np.sum(var[sel] * sel_areas)
            total_area = np.sum(sel_areas)
            return total_var / total_area
        elif var.ndim == 2:
            # Determine which axis is facets
            if var.shape[0] == n_facets:
                # Shape (n_facets, n_other) - facets in first dimension
                if isinstance(sel, slice):
                    total_var = np.sum(var * areas[:, np.newaxis], axis=0)
                else:
                    total_var = np.sum(var[sel, :] * sel_areas[:, np.newaxis], axis=0)
                total_area = np.sum(sel_areas)
                return total_var / total_area
            elif var.shape[1] == n_facets:
                # Shape (n_other, n_facets) - facets in second dimension  
                if isinstance(sel, slice):
                    total_var = np.sum(var * areas[np.newaxis, :], axis=1)
                else:
                    total_var = np.sum(var[:, sel] * sel_areas[np.newaxis, :], axis=1)
                total_area = np.sum(sel_areas)
                return total_var / total_area
            else:
                raise ValueError(f"Neither dimension of var {var.shape} matches number of facets {n_facets}")
        else:
            raise ValueError(f"area_average_fac only supports 1D and 2D arrays, got shape {var.shape}")
    
    def area_average_seb(self, seb: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Area-average all surface energy balance terms.
        
        Parameters
        ----------
        seb : dict
            Dictionary from load_seb()
        
        Returns
        -------
        dict
            Area-averaged SEB terms
        
        Examples
        --------
        >>> seb = sim.load_seb()
        >>> seb_avg = sim.area_average_seb(seb)
        >>> import matplotlib.pyplot as plt
        >>> plt.plot(seb_avg['t'], seb_avg['Kstar'])
        """
        return {
            'Kstar': self.area_average_fac(seb['Kstar']),
            'Lstar': self.area_average_fac(seb['Lstar']),
            'Lin': self.area_average_fac(seb['Lin']),
            'Lout': self.area_average_fac(seb['Lout']),
            'H': self.area_average_fac(seb['H']),
            'E': self.area_average_fac(seb['E']),
            'G': self.area_average_fac(seb['G']),
            't': seb['t']
        }
    
    @staticmethod
    def time_average(var: np.ndarray, other: Optional[np.ndarray] = None):
        """
        Time-average variables.
        
        Parameters
        ----------
        var : ndarray
            Variable to average. Time must be the last dimension.
        other : ndarray, optional
            Second variable (same shape as ``var``) to compute covariance with.
        
        Returns
        -------
        tuple
            - (mean, variance) when only ``var`` is provided
            - (x_mean, y_mean, covariance) when ``other`` is provided
        
        Notes
        -----
        Variance/covariance are computed by delegating to ``merge_stat`` with
        zero instantaneous contributions.
        """
        X = np.asarray(var)
        n = X.shape[-1]
        
        if other is None:
            zeros = np.zeros_like(X)
            return UDBase.merge_stat(X, zeros, n)
        
        Y = np.asarray(other)
        zeros = np.zeros_like(X)
        return UDBase.merge_stat(X, Y, zeros, n)

    @staticmethod
    def merge_stat(X: np.ndarray, *args, Y: Optional[np.ndarray] = None,
                   XpXp: Optional[np.ndarray] = None,
                   XpYp: Optional[np.ndarray] = None):
        """
        Merge short-term statistics into longer windows.

        Supported call patterns:
          - merge_stat(X, n)
          - merge_stat(X, XpXp, n)                  
          - merge_stat(X, Y, XpYp, n)               
          - merge_stat(X, n, XpXp=array)            # single variable mean/variance
          - merge_stat(X, n, Y=array, XpYp=array)   # two variables mean/cov

        Parameters
        ----------
        X : ndarray
            First variable. Final dimension is time or short-term windows.
        args : tuple
            Positional parameters parsed according to the patterns above.
        Y : ndarray, optional
            Second variable (same trailing dimension as ``X``).
        XpXp : ndarray, optional
            Variance contribution for ``X`` aligned with ``X``.
        XpYp : ndarray, optional
            Covariance contribution for ``X`` and ``Y`` aligned with ``X``.

        Returns
        -------
        tuple | ndarray
            ``Xmean`` if only ``X`` provided;
            ``Xmean, var`` if ``XpXp`` given;
            ``Xmean, Ymean, cov`` if ``Y`` provided (and optionally ``XpYp``).

        Notes
        -----
        Discards the oldest samples that do not fill a complete window so the
        most recent data is retained. Variance/covariance combine the mean of
        short-window contributions with the variance/covariance of the short
        means inside each merged window.
        """
        # Parse positional arguments to support both MATLAB and Python styles
        n = None
        X = np.asarray(X)
        if len(args) == 1:
            n = int(args[0])
        elif len(args) == 2 and Y is None:
            # MATLAB style: merge_stat(X, XpXp, n)
            XpXp = np.asarray(args[0])
            n = int(args[1])
        elif len(args) == 3:
            # MATLAB style: merge_stat(X, Y, XpYp, n)
            Y = np.asarray(args[0])
            XpYp = np.asarray(args[1])
            n = int(args[2])
        else:
            raise ValueError("merge_stat expects 1, 2, or 3 positional arguments after X")

        if n <= 0:
            raise ValueError("n must be positive")
        if X.shape[-1] < n:
            raise ValueError("Not enough samples to form a single merged window")

        nwin = X.shape[-1] // n
        start = X.shape[-1] - nwin * n  # discard oldest incomplete window
        X_use = X[..., start:]
        X_group = X_use.reshape(*X.shape[:-1], nwin, n)
        Xmean = X_group.mean(axis=-1)

        if Y is None:
            if XpXp is None:
                return Xmean

            XpXp = np.asarray(XpXp)
            if XpXp.shape[-1] != X.shape[-1]:
                raise ValueError("XpXp must match X shape in the last dimension")
            XpXp_use = XpXp[..., start:]
            XpXp_group = XpXp_use.reshape(*XpXp.shape[:-1], nwin, n)
            within = XpXp_group.mean(axis=-1)
            between = ((X_group - Xmean[..., None]) ** 2).mean(axis=-1)
            return Xmean, within + between

        Y = np.asarray(Y)
        if Y.shape[-1] < n:
            raise ValueError("Y does not have enough samples to form a merged window")
        if Y.shape[-1] != X.shape[-1]:
            raise ValueError("X and Y must share the same length in the last dimension")

        Y_use = Y[..., start:]
        Y_group = Y_use.reshape(*Y.shape[:-1], nwin, n)
        Ymean = Y_group.mean(axis=-1)

        if XpYp is None:
            cov_within = ((X_group - Xmean[..., None]) * (Y_group - Ymean[..., None])).mean(axis=-1)
            between = 0.0
        else:
            XpYp = np.asarray(XpYp)
            if XpYp.shape[-1] != X.shape[-1]:
                raise ValueError("XpYp must match X and Y in the last dimension")
            XpYp_use = XpYp[..., start:]
            XpYp_group = XpYp_use.reshape(*XpYp.shape[:-1], nwin, n)
            cov_within = XpYp_group.mean(axis=-1)
            between = ((X_group - Xmean[..., None]) * (Y_group - Ymean[..., None])).mean(axis=-1)

        cov = cov_within + between
        return Xmean, Ymean, cov

    @staticmethod
    def coarsegrain_field(var: np.ndarray, Lflt: np.ndarray,
                          xm: np.ndarray, ym: np.ndarray) -> np.ndarray:
        """
        Apply 2D periodic box filters to a 3D field.

        Parameters
        ----------
        var : ndarray, shape (nx, ny, nz)
            Field data (time or height on the last axis).
        Lflt : array-like
            Filter lengths (meters). Multiple lengths are allowed.
        xm, ym : array-like
            Grid coordinates in meters; must be 1D and roughly uniform.

        Returns
        -------
        ndarray
            Filtered field with shape (nx, ny, nz, n_filters).
        """
        var = np.asarray(var)
        if var.ndim != 3:
            raise ValueError("var must be 3D with shape (nx, ny, nz)")

        xm = np.asarray(xm).ravel()
        ym = np.asarray(ym).ravel()
        if xm.size < 2 or ym.size < 2:
            raise ValueError("xm and ym must contain at least two points")

        dx = float(np.mean(np.diff(xm)))
        dy = float(np.mean(np.diff(ym)))
        if dx <= 0 or dy <= 0:
            raise ValueError("Grid spacings must be positive")

        Lflt_arr = np.atleast_1d(Lflt)
        nx, ny, nz = var.shape
        out = np.empty((nx, ny, nz, len(Lflt_arr)))

        # Build filters using half-width convention (matches MATLAB implementation)
        ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
        di = np.minimum(ii, nx - ii)  # periodic distance in i
        dj = np.minimum(jj, ny - jj)  # periodic distance in j

        for i, L in enumerate(Lflt_arr):
            ngx = max(int(round((L / dx) / 2.0)), 1)
            ngy = max(int(round((L / dy) / 2.0)), 1)
            mask = (di <= ngx) & (dj <= ngy)
            kernel = mask.astype(float)
            kernel /= kernel.sum()
            k_hat = np.fft.fftn(kernel)

            for k in range(nz):
                v_hat = np.fft.fftn(var[:, :, k])
                out[:, :, k, i] = np.real(np.fft.ifftn(v_hat * k_hat))

        return out
    
    def plot_veg(self, veg: Optional[Dict[str, Any]] = None, show: bool = False):
        """Plot vegetation points on top of the geometry using the visualization facade."""
        return self.vis.plot_veg(veg=veg, show=show)

    def plot_trees(self, show: bool = False):
        """Backward-compatible alias for plot_veg."""
        return self.vis.plot_trees(show=show)
    
    def plot_fac(self, var: np.ndarray, building_ids: Optional[np.ndarray] = None, show: bool = True):
        """Plot facet data as a 3D surface."""
        return self.vis.plot_fac(var=var, building_ids=building_ids, show=show)
    
    def _create_colored_mesh(self, var: np.ndarray, building_ids: Optional[np.ndarray] = None):
        """Create a colored trimesh object from facet data, optionally filtered by building IDs."""
        return self.vis._create_colored_mesh(var=var, building_ids=building_ids)
    
    def _render_scene(self, mesh, show_outlines: bool = True, angle_threshold: float = 45.0, building_ids: Optional[np.ndarray] = None, custom_edges: Optional[List[tuple]] = None, show: bool = True):
        """Render the mesh scene using trimesh/plotly. Returns figure handle if available."""
        return self.vis._render_scene(
            mesh=mesh,
            show_outlines=show_outlines,
            angle_threshold=angle_threshold,
            building_ids=building_ids,
            custom_edges=custom_edges,
            show=show,
        )
    
    def _render_plotly(self, meshes, outline_edges, show: bool = True):
        """Render using plotly for notebook display. Returns the figure object."""
        return self.vis._render_plotly(meshes=meshes, outline_edges=outline_edges, show=show)
    
    def _render_trimesh(self, scene, num_outline_edges, show: bool = True):
        """Render using trimesh viewer."""
        return self.vis._render_trimesh(scene=scene, num_outline_edges=num_outline_edges, show=show)
    
    def _add_building_outlines_to_scene(self, building_ids=None):
        """Add building outlines to the current scene (placeholder for future implementation)."""
        return self.vis._add_building_outlines_to_scene(building_ids=building_ids)
    
    def _add_building_outlines(self, ax, building_ids=None, angle_threshold: float = 45.0):
        """Helper method to add building outline edges to current 3D matplotlib plot."""
        return self.vis._add_building_outlines(
            ax=ax,
            building_ids=building_ids,
            angle_threshold=angle_threshold,
        )
    
    def plot_fac_type(self, building_ids: Optional[np.ndarray] = None, 
                      show_outlines: bool = True, angle_threshold: float = 45.0, show: bool = True):
        """Plot the different surface types in the geometry using the visualization facade."""
        return self.vis.plot_fac_type(
            building_ids=building_ids,
            show_outlines=show_outlines,
            angle_threshold=angle_threshold,
            show=show,
        )
    
    def convert_fac_to_field(self, var: np.ndarray, facsec: Optional[Dict] = None,
                            dz: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Convert a facet variable to a 3D field density.
        
        This method converts facet-based data to a 3D volumetric field by
        distributing facet areas into grid cells. Useful for computing volume
        integrals of surface quantities.
        
        Parameters
        ----------
        var : ndarray, shape (n_facets,)
            Facet variable to convert to field
        facsec : dict, optional
            Facet sections dictionary. If None, uses self.facsec['c']
        dz : ndarray, optional
            Vertical grid spacing. If None, uses self.dzt
            
        Returns
        -------
        fld : ndarray, shape (itot, jtot, ktot)
            3D field density [var_units / m]
            
        Raises
        ------
        ValueError
            If facet section data is not available
            
        Notes
        -----
        The conversion computes a density field where each cell contains:
        fld[i,j,k] = sum(var[facid] * area / (dx * dy * dz))
        
        This is used internally by calculate_frontal_properties() to compute
        frontal area densities from projected facet areas.
        
        Examples
        --------
        Convert surface temperature to a field:
        >>> Ts = sim.load_fac_temperature('Ts')
        >>> T_field = sim.convert_fac_to_field(Ts[:, 0])
        
        Use custom facet sections (e.g., for u-grid):
        >>> var = np.random.randn(sim.geom.n_faces)
        >>> fld = sim.convert_fac_to_field(var, facsec=sim.facsec['u'])
        """
        # Check that facet sections are available
        if not hasattr(self, 'facsec') or self.facsec is None:
            raise ValueError(
                "This method requires facet section data. "
                f"Ensure {self.ffacet_sections}_(u,v,w,c).{self.expnr} and "
                f"{self.ffluid_boundary}_(u,v,w,c).{self.expnr} files exist."
            )
        
        # Use defaults if not provided
        if facsec is None:
            facsec = self.facsec['c']
        
        if dz is None:
            dz = self.dzt
        
        # Initialize field
        fld = np.zeros((self.itot, self.jtot, self.ktot), dtype=np.float32)
        
        # Get facet section data
        facids = facsec['facid']
        areas = facsec['area']
        locs = facsec['locs']  # (i, j, k) locations (0-based)
        
        i_idx = locs[:, 0].astype(int)
        j_idx = locs[:, 1].astype(int)
        k_idx = locs[:, 2].astype(int)
        
        # Loop over all facet sections and create density field
        for m in range(len(areas)):
            facid = facids[m]
            i, j, k = i_idx[m], j_idx[m], k_idx[m]
            
            # Add contribution to cell
            cell_volume = self.dx * self.dy * dz[k]
            fld[i, j, k] += var[facid] * areas[m] / cell_volume
        
        return fld
    
    def convert_facvar_to_field(self, var: np.ndarray, facsec: Dict,
                                building_ids: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Transfer a facet variable onto the grid.
        
        Matches MATLAB implementation: convert_facvar_to_field(obj, var, facsec, building_ids)
        
        Parameters
        ----------
        var : ndarray
            Facet variable (e.g., from load_fac_eb, load_fac_temperature)
        facsec : dict
            Facet section structure (e.g., obj.facsec['u'])
        building_ids : array-like, optional
            Array of building IDs to include. If None, all buildings included.
            
        Returns
        -------
        ndarray
            Variable on the grid (itot x jtot x ktot)
            
        Examples
        --------
        >>> # Convert all facets
        >>> fld = sim.convert_facvar_to_field(var, sim.facsec['c'])
        
        >>> # Convert only specific buildings
        >>> fld = sim.convert_facvar_to_field(var, sim.facsec['c'], [1, 5, 10])
        """
        # Normalize by converting ones, then divide actual values by normalization
        norm = self.convert_facflx_to_field(np.ones_like(var), facsec, self.dzt, building_ids)
        norm[norm == 0] = 1  # Avoid NaNs
        fld = self.convert_facflx_to_field(var, facsec, self.dzt, building_ids) / norm
        return fld
    
    def convert_facflx_to_field(self, var: np.ndarray, facsec: Dict, dz: np.ndarray,
                                building_ids: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Convert a facet flux variable to a density in a 3D field.
        
        Matches MATLAB implementation: convert_facflx_to_field(obj, var, facsec, dz, building_ids)
        
        Parameters
        ----------
        var : ndarray
            Facet flux variable (e.g., from load_fac_eb)
        facsec : dict
            Facet section structure (e.g., obj.facsec['u'])
        dz : ndarray
            Vertical grid spacing at cell centers (obj.dzt)
        building_ids : array-like, optional
            Array of building IDs to include. If None, all buildings included.
            
        Returns
        -------
        ndarray
            3D field density (itot x jtot x ktot)
            
        Raises
        ------
        ValueError
            If facet section data is not loaded
            
        Examples
        --------
        >>> # Convert all facets
        >>> fld = sim.convert_facflx_to_field(var, sim.facsec['c'], sim.dzt)
        
        >>> # Convert only specific buildings
        >>> fld = sim.convert_facflx_to_field(var, sim.facsec['c'], sim.dzt, [1, 5, 10])
        """
        if not hasattr(self, 'facsec') or self.facsec is None:
            raise ValueError("This method requires facet section data. "
                           "Ensure facet_sections_*.txt and fluid_boundary_*.txt files exist.")
        
        # Initialize output field
        fld = np.zeros((self.itot, self.jtot, self.ktot), dtype=np.float32)
        
        # Get facet section data
        facids = facsec['facid']
        areas = facsec['area']
        locs = facsec['locs']  # (i, j, k) locations (0-based)
        
        # If building IDs specified, get face mask for filtering
        face_mask = None
        if building_ids is not None:
            # Get face to building mapping
            face_to_building = self.geom.get_face_to_building_map()
            building_ids = np.asarray(building_ids)
            
            # Create face mask for requested building IDs
            face_mask = np.isin(face_to_building, building_ids)
        
        # Loop over all facet sections and create density field
        for m in range(len(facids)):
            facid = int(facids[m])
            
            # If building filtering is active, check if this face should be included
            if face_mask is not None and not face_mask[facid]:
                continue
            
            # Get grid location (convert from 1-indexed to 0-indexed)
            i, j, k = int(locs[m, 0]), int(locs[m, 1]), int(locs[m, 2])
            
            # Add contribution to cell
            cell_volume = self.dx * self.dy * dz[k]
            fld[i, j, k] += var[facid] * areas[m] / cell_volume
        
        return fld
    
    def load_prof(self) -> np.ndarray:
        """
        Load information from prof.inp file.
        
        Matches MATLAB implementation: load_prof(obj)
        
        Returns
        -------
        ndarray
            Profile data from prof.inp file
            
        Raises
        ------
        FileNotFoundError
            If prof.inp file does not exist
            
        Examples
        --------
        >>> prof_data = sim.load_prof()
        """
        fname = f"{self.fprof}.{self.expnr}"
        fpath = self.path / fname
        
        if not fpath.exists():
            raise FileNotFoundError(f"Profile file not found: {fpath}")
        
        # Read data starting from line 3 (skip 2 header lines)
        data = np.loadtxt(fpath, skiprows=2)
        return data
    
    def load_facsec(self, var: str) -> Dict[str, np.ndarray]:
        """
        Load facet section information for a specific variable.
        
        Matches MATLAB implementation: load_facsec(obj, strvar)
        
        Parameters
        ----------
        var : str
            Variable name ('u', 'v', 'w', or 'c')
            
        Returns
        -------
        dict
            Dictionary with keys:
            - 'facid': Facet IDs (0-based)
            - 'area': Facet section areas
            - 'locs': Fluid boundary point locations (i, j, k), 0-based
            - 'distance': Distances from facet to fluid point
            
        Raises
        ------
        FileNotFoundError
            If facet section or fluid boundary files don't exist
            
        Examples
        --------
        >>> facsec_u = sim.load_facsec('u')
        >>> print(facsec_u.keys())
        """
        # Load facet section data
        fname_sec = self.path / f"{self.ffacet_sections}_{var}.txt"
        if not fname_sec.exists():
            raise FileNotFoundError(f"Facet section file not found: {fname_sec}")
        
        facsecs = self._load_sparse_file(fname_sec, skiprows=1, dtype=float, min_cols=4)
        
        # Load fluid boundary points
        fname_fluid = self.path / f"{self.ffluid_boundary}_{var}.txt"
        if not fname_fluid.exists():
            raise FileNotFoundError(f"Fluid boundary file not found: {fname_fluid}")
        
        fluid_boundary = self._load_sparse_file(
            fname_fluid,
            skiprows=1,
            dtype=int,
            min_cols=3,
            zero_based_cols=[0, 1, 2],
        )

        # Create structure matching MATLAB
        data = {
            'facid': facsecs[:, 0].astype(int) - 1,
            'area': facsecs[:, 1],
            'locs': fluid_boundary[facsecs[:, 2].astype(int) - 1, :].astype(int),
            'distance': facsecs[:, 3]
        }
        
        return data
    
    def calculate_frontal_properties(self) -> Dict[str, Any]:
        """
        Calculate skyline, frontal areas, and blockage ratios.
        
        Computes geometric properties of the urban canopy including:
        - Skyline profiles in x and y directions
        - Total frontal areas perpendicular to x and y
        - Blockage ratios (fraction of frontal area blocked)
        
        Returns
        -------
        res : dict
            Dictionary containing:
            - 'skylinex' : ndarray, shape (jtot, ktot)
                Binary indicator of blocked cells in x-direction (y-z plane)
            - 'skyliney' : ndarray, shape (itot, ktot)
                Binary indicator of blocked cells in y-direction (x-z plane)
            - 'Afx' : float
                Total frontal area perpendicular to x-direction [m²]
            - 'Afy' : float
                Total frontal area perpendicular to y-direction [m²]
            - 'brx' : float
                Blockage ratio in x-direction (dimensionless, 0-1)
            - 'bry' : float
                Blockage ratio in y-direction (dimensionless, 0-1)
                
        Raises
        ------
        ValueError
            If geometry or facet section data is not available
            
        Notes
        -----
        The frontal area is the projected area of all surfaces onto a plane
        perpendicular to the flow direction. The blockage ratio is the fraction
        of the domain cross-section that is blocked by buildings.
        
        These quantities are important for understanding drag and flow resistance
        in urban canopies.
        
        Examples
        --------
        >>> props = sim.calculate_frontal_properties()
        >>> print(f"Frontal area (x): {props['Afx']:.1f} m²")
        >>> print(f"Blockage ratio (x): {props['brx']:.3f}")
        >>> 
        >>> # Visualize skyline
        >>> import matplotlib.pyplot as plt
        >>> plt.imshow(props['skylinex'].T, origin='lower', aspect='auto')
        >>> plt.xlabel('j (grid points)')
        >>> plt.ylabel('k (grid points)')
        >>> plt.title('Skyline in x-direction')
        >>> plt.colorbar(label='Blocked (1) / Open (0)')
        >>> plt.show()
        """
        # Check required data
        if self.geom is None:
            raise ValueError("This method requires a geometry (STL) file. "
                           "Ensure stl_file is specified in namoptions.")
        
        if not hasattr(self, 'facsec') or self.facsec is None:
            raise ValueError(
                "This method requires facet section data. "
                f"Ensure {self.ffacet_sections}_(u,v,w,c).{self.expnr} and "
                f"{self.ffluid_boundary}_(u,v,w,c).{self.expnr} files exist."
            )
        
        # Get face normals
        norms = self.geom.face_normals
        
        # Create surface quantity for projected area in x and y directions
        # Project onto planes perpendicular to x and y axes
        # Only count outward-facing (negative dot product with axis)
        phix = -np.minimum(np.dot(norms, np.array([1, 0, 0])), 0)
        phiy = -np.minimum(np.dot(norms, np.array([0, 1, 0])), 0)
        
        # Convert to density fields using convert_facflx_to_field (matches MATLAB)
        rhoLx = self.convert_facflx_to_field(phix, self.facsec['c'], self.dzt)
        rhoLy = self.convert_facflx_to_field(phiy, self.facsec['c'], self.dzt)
        
        # Calculate indicator functions for blockage
        # Ibx[j,k] = 1 if any cell along x-direction at (j,k) is blocked
        Ibx = (np.sum(rhoLx, axis=0) > 0).astype(float)
        # Iby[i,k] = 1 if any cell along y-direction at (i,k) is blocked
        Iby = (np.sum(rhoLy, axis=1) > 0).astype(float)
        
        # Integrate to get frontal areas and blockage ratios
        Afx = 0.0
        Afy = 0.0
        brx = 0.0
        bry = 0.0
        
        for k in range(self.ktot):
            # Frontal areas (integrate density over volume)
            Afx += np.sum(rhoLx[:, :, k]) * self.dx * self.dy * self.dzt[k]
            Afy += np.sum(rhoLy[:, :, k]) * self.dx * self.dy * self.dzt[k]
            
            # Blockage ratios (integrate indicator over cross-section)
            brx += np.sum(Ibx[:, k]) * self.dy * self.dzt[k]
            bry += np.sum(Iby[:, k]) * self.dx * self.dzt[k]
        
        # Normalize blockage ratios by total cross-sectional area
        brx /= (self.ylen * self.zsize)
        bry /= (self.xlen * self.zsize)
        
        # Print results
        print(f"x-direction: frontal area = {Afx:8.1f} m², blockage ratio = {brx:8.3f}")
        print(f"y-direction: frontal area = {Afy:8.1f} m², blockage ratio = {bry:8.3f}")
        
        return {
            'skylinex': Ibx,
            'skyliney': Iby,
            'Afx': Afx,
            'Afy': Afy,
            'brx': brx,
            'bry': bry
        }
    
    def plot_building_ids(self, show: bool = True):
        """Plot building IDs from above (x,y view) with distinct colors."""
        return self.vis.plot_building_ids(show=show)
    
    def plot_2dmap(self, val: Union[float, np.ndarray], 
                   labels: Optional[Union[str, list]] = None,
                   show: bool = True):
        """Plot a 2D map of buildings colored by a value per building."""
        return self.vis.plot_2dmap(val=val, labels=labels, show=show)
    
    def __str__(self):
        """User-friendly string representation."""
        return self.__repr__()


if __name__ == "__main__":
    # Basic test
    print("UDBase module loaded successfully")
    print("\nUsage:")
    print("  from udbase import UDBase")
    print("  sim = UDBase(expnr=65, path='experiments/065')")
