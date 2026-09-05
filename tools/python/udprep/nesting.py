"""Preprocessing for one-way nesting (``nesting.inp.<expnr>.nc``).

Implements work package W1-PY of the nesting feature:

* **conservative flux interpolation** of a parent velocity field onto the child
  faces (design ``docs/udales-nesting-design.md`` §1.3);
* the offline **divergence correction** that makes the net volume flux through
  the child's lateral boundary vanish on every stored time level (§3.1, §3.3);
* the **writer/reader** for the file format of ``docs/udales-nesting-spec.md`` §5,
  with schema validation and a refinement-ratio guard;
* an **analytic-field generator** used by the in-solver unit tests (U15--U22).

Interpolation
-------------
The prolongation is the tensor product of

* piecewise-constant distribution in the two directions **tangential** to the
  face, so that a parent face flux is shared over the child faces it contains
  in exact proportion to their areas (the defining property, design §1.3), and
* linear interpolation in the direction **normal** to the face, between the two
  bracketing parent faces.

A child face that is coplanar with a parent face therefore takes the parent
value exactly, and the sum of the child fluxes over a parent face equals the
parent face flux to round-off.  Because the normal direction is interpolated
linearly rather than injected, the child target additionally reproduces the
parent's discrete divergence *cell by cell*, so a discretely solenoidal parent
gives a discretely solenoidal child target.  The design only requires the
weaker statement that the divergence integrated over a parent cell is exact;
this scheme satisfies both, and is second-order rather than first-order in the
normal direction.

Analytic field
--------------
:func:`analytic_field` evaluates, at each variable's own stagger,

    f(x, y, z, t) = sin(a x + b y) cos(c z) (1 + d t) + e x y z

with a distinct coefficient set per velocity component
(:data:`ANALYTIC_COEFFS`).  It is non-separable in (x, y) and differs between
components, so a transposed, swapped or off-by-one index cannot reproduce it.
The Fortran side (`TEST_NESTING_IO`) must use the same formula and constants.
"""

from __future__ import annotations

import getpass
import json
import os
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

from exceptions import ConfigurationError, DataFormatError, DependencyError

__all__ = [
    "ANALYTIC_COEFFS",
    "COMPONENTS",
    "FACES",
    "MAX_SPATIAL_REFINEMENT",
    "MAX_TEMPORAL_REFINEMENT",
    "REQUIRED_GLOBAL_ATTRIBUTES",
    "SCHEMA_VERSION",
    "SLAB_VARIABLES",
    "STAGGER",
    "TOOL_VERSION",
    "FaceMasks",
    "NestGrid",
    "NestingData",
    "NestingRefinementError",
    "NestingSchemaError",
    "analytic_field",
    "analytic_slabs",
    "apply_divergence_correction",
    "boundary_faces",
    "check_refinement",
    "conservative_interpolate",
    "fluid_face_area",
    "interpolate_child_fields",
    "net_volume_flux",
    "nesting_data_from_parent",
    "nesting_filename",
    "read_nesting_file",
    "refinement_ratios",
    "slab_coordinates",
    "slab_dimensions",
    "slab_indices",
    "slab_shape",
    "slabs_from_fields",
    "slabs_from_parent",
    "validate_nesting_file",
    "write_analytic_nesting_file",
    "write_nesting_file",
]

# --------------------------------------------------------------------------- #
# Contract constants (docs/udales-nesting-spec.md §5)
# --------------------------------------------------------------------------- #

SCHEMA_VERSION = 1
TOOL_VERSION = "udprep.nesting/1.0"

FACES = ("west", "east", "south", "north")
COMPONENTS = ("u", "v", "w")
SLAB_VARIABLES = tuple(f"{c}_{f}" for f in FACES for c in COMPONENTS)

#: ``stagger`` attribute required on every slab variable.
STAGGER = {"u": "xh yf zf", "v": "xf yh zf", "w": "xf yf zh"}

#: Per component, per axis (x, y, z): True where the component lives on a face.
_IS_FACE = {"u": (True, False, False),
            "v": (False, True, False),
            "w": (False, False, True)}

#: Axis (0 = x, 1 = y) normal to each lateral boundary face.
_FACE_AXIS = {"west": 0, "east": 0, "south": 1, "north": 1}

#: True for the faces at the high-index end of the domain.
_FACE_UPPER = {"west": False, "east": True, "south": False, "north": True}

#: Outward unit normal component of each face, in the axis's own direction.
_FACE_SIGN = {"west": -1.0, "east": +1.0, "south": -1.0, "north": +1.0}

#: Velocity component whose faces coincide with each lateral boundary.
_FACE_NORMAL_COMPONENT = {"west": "u", "east": "u", "south": "v", "north": "v"}

MAX_SPATIAL_REFINEMENT = 4.0
MAX_TEMPORAL_REFINEMENT = 30.0

REQUIRED_GLOBAL_ATTRIBUTES = (
    "Conventions",
    "udales_nesting_schema",
    "divergence_corrected",
    "itot",
    "jtot",
    "ktot",
    "nzone",
    "xlen",
    "ylen",
    "parent_model",
    "parent_dx",
    "parent_dt",
    "child_origin_x",
    "child_origin_y",
    "rotation_deg",
    "created",
    "creator",
    "tool_version",
)

_COORD_VARIABLES = ("xf", "xh", "yf", "yh", "zf", "zh")

#: Snapping tolerance for "this child face is coplanar with a parent face",
#: relative to the smallest parent spacing.
_ALIGN_RTOL = 1.0e-9


class NestingSchemaError(DataFormatError):
    """A nesting file does not match the schema of the implementation contract."""


class NestingRefinementError(ConfigurationError):
    """The parent/child refinement ratio exceeds the supported range."""


def _import_dataset():
    """Return ``netCDF4.Dataset``, with the repo's standard dependency error."""
    try:
        from netCDF4 import Dataset
    except ImportError as exc:  # pragma: no cover - netCDF4 is a hard dependency
        raise DependencyError(
            "netCDF4 is required to read or write nesting.inp.<expnr>.nc"
        ) from exc
    return Dataset


def nesting_filename(expnr: Any) -> str:
    """Return the default nesting input file name for experiment ``expnr``."""
    if isinstance(expnr, str):
        return f"nesting.inp.{expnr}.nc"
    return f"nesting.inp.{int(expnr):03d}.nc"


# --------------------------------------------------------------------------- #
# Grid
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class NestGrid:
    """A staggered Cartesian grid, in the child's (or parent's) own metres.

    ``xh``/``yh``/``zh`` are the cell faces (``itot+1``, ``jtot+1``, ``ktot+1``
    long); ``xf``/``yf``/``zf`` the cell centres.  The vertical may be stretched;
    the horizontal is assumed uniform, as everywhere else in uDALES.
    """

    xh: np.ndarray
    yh: np.ndarray
    zh: np.ndarray
    xf: np.ndarray
    yf: np.ndarray
    zf: np.ndarray

    def __post_init__(self) -> None:
        for name in ("xh", "yh", "zh", "xf", "yf", "zf"):
            arr = np.ascontiguousarray(np.asarray(getattr(self, name), dtype=np.float64))
            object.__setattr__(self, name, arr)
        for name in ("xh", "yh", "zh"):
            edges = getattr(self, name)
            if edges.ndim != 1 or edges.size < 2:
                raise ConfigurationError(f"NestGrid.{name} must be a 1-D array of length >= 2")
            if not np.all(np.diff(edges) > 0.0):
                raise ConfigurationError(f"NestGrid.{name} must be strictly increasing")
        for face, centre in (("xh", "xf"), ("yh", "yf"), ("zh", "zf")):
            if getattr(self, centre).size != getattr(self, face).size - 1:
                raise ConfigurationError(
                    f"NestGrid.{centre} must have one element fewer than NestGrid.{face}"
                )

    @classmethod
    def uniform(
        cls,
        itot: int,
        jtot: int,
        ktot: int,
        xlen: float,
        ylen: float,
        zsize: float,
        x0: float = 0.0,
        y0: float = 0.0,
        z0: float = 0.0,
    ) -> "NestGrid":
        """Build a uniform grid of the given size and origin."""
        xh = x0 + np.linspace(0.0, xlen, itot + 1)
        yh = y0 + np.linspace(0.0, ylen, jtot + 1)
        zh = z0 + np.linspace(0.0, zsize, ktot + 1)
        return cls.from_faces(xh, yh, zh)

    @classmethod
    def from_faces(cls, xh, yh, zh) -> "NestGrid":
        """Build a grid from its face coordinates, centres at the midpoints."""
        xh = np.asarray(xh, dtype=np.float64)
        yh = np.asarray(yh, dtype=np.float64)
        zh = np.asarray(zh, dtype=np.float64)
        return cls(
            xh=xh, yh=yh, zh=zh,
            xf=0.5 * (xh[:-1] + xh[1:]),
            yf=0.5 * (yh[:-1] + yh[1:]),
            zf=0.5 * (zh[:-1] + zh[1:]),
        )

    # -- sizes ------------------------------------------------------------- #

    @property
    def itot(self) -> int:
        return int(self.xf.size)

    @property
    def jtot(self) -> int:
        return int(self.yf.size)

    @property
    def ktot(self) -> int:
        return int(self.zf.size)

    @property
    def xlen(self) -> float:
        return float(self.xh[-1] - self.xh[0])

    @property
    def ylen(self) -> float:
        return float(self.yh[-1] - self.yh[0])

    @property
    def zsize(self) -> float:
        return float(self.zh[-1] - self.zh[0])

    # -- spacings ---------------------------------------------------------- #

    @property
    def dx(self) -> np.ndarray:
        return np.diff(self.xh)

    @property
    def dy(self) -> np.ndarray:
        return np.diff(self.yh)

    @property
    def dzf(self) -> np.ndarray:
        return np.diff(self.zh)

    def edges(self, axis: int) -> np.ndarray:
        """Face coordinates along ``axis`` (0 = x, 1 = y, 2 = z)."""
        return (self.xh, self.yh, self.zh)[axis]

    def centres(self, axis: int) -> np.ndarray:
        """Cell-centre coordinates along ``axis`` (0 = x, 1 = y, 2 = z)."""
        return (self.xf, self.yf, self.zf)[axis]

    def component_coords(self, component: str, axis: int) -> np.ndarray:
        """Coordinates along ``axis`` at the stagger of ``component``."""
        if _IS_FACE[component][axis]:
            return self.edges(axis)
        return self.centres(axis)

    def component_shape(self, component: str) -> Tuple[int, int, int]:
        """Array shape of ``component`` on this grid, as ``(nx, ny, nz)``."""
        return tuple(self.component_coords(component, ax).size for ax in range(3))


# --------------------------------------------------------------------------- #
# Conservative flux interpolation (design §1.3)
# --------------------------------------------------------------------------- #


def _normal_map(parent_faces: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Linear map along the face-normal direction.

    Returns ``(idx, wlo)`` such that the interpolant is
    ``wlo * f[idx] + (1 - wlo) * f[idx + 1]``.  Targets that coincide with a
    parent face (to within :data:`_ALIGN_RTOL` of the smallest parent spacing)
    are snapped so that they reproduce the parent value bit for bit -- this is
    what makes the flux identity exact.
    """
    target = np.asarray(target, dtype=np.float64)
    tol = _ALIGN_RTOL * float(np.min(np.diff(parent_faces)))
    if np.any(target < parent_faces[0] - tol) or np.any(target > parent_faces[-1] + tol):
        raise ConfigurationError(
            "child grid extends outside the parent grid: target range "
            f"[{target.min():g}, {target.max():g}] vs parent "
            f"[{parent_faces[0]:g}, {parent_faces[-1]:g}]"
        )
    idx = np.searchsorted(parent_faces, target, side="right") - 1
    idx = np.clip(idx, 0, parent_faces.size - 2)
    lo = parent_faces[idx]
    hi = parent_faces[idx + 1]
    wlo = (hi - target) / (hi - lo)
    wlo = np.where(np.abs(target - lo) <= tol, 1.0, wlo)
    wlo = np.where(np.abs(target - hi) <= tol, 0.0, wlo)
    return idx, wlo


def _tangential_map(parent_faces: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Piecewise-constant map: index of the parent cell containing each target."""
    target = np.asarray(target, dtype=np.float64)
    tol = _ALIGN_RTOL * float(np.min(np.diff(parent_faces)))
    if np.any(target < parent_faces[0] - tol) or np.any(target > parent_faces[-1] + tol):
        raise ConfigurationError(
            "child grid extends outside the parent grid: target range "
            f"[{target.min():g}, {target.max():g}] vs parent "
            f"[{parent_faces[0]:g}, {parent_faces[-1]:g}]"
        )
    idx = np.searchsorted(parent_faces, target, side="right") - 1
    return np.clip(idx, 0, parent_faces.size - 2)


def _apply_normal(arr: np.ndarray, axis: int, idx: np.ndarray, wlo: np.ndarray) -> np.ndarray:
    shape = [1] * arr.ndim
    shape[axis] = -1
    w = wlo.reshape(shape)
    lo = np.take(arr, idx, axis=axis)
    hi = np.take(arr, idx + 1, axis=axis)
    return lo * w + hi * (1.0 - w)


def conservative_interpolate(
    parent: NestGrid,
    field: np.ndarray,
    component: str,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
) -> np.ndarray:
    """Interpolate one parent velocity component onto child faces conservatively.

    Parameters
    ----------
    parent
        The parent grid.
    field
        The parent field at the stagger of ``component``, shape
        ``parent.component_shape(component)``.
    component
        ``'u'``, ``'v'`` or ``'w'``.
    x, y, z
        Target coordinates, each already at ``component``'s own stagger: face
        coordinates in the component's normal direction, cell-centre
        coordinates in the other two.

    Returns
    -------
    numpy.ndarray
        The interpolated field, shape ``(x.size, y.size, z.size)``.

    Notes
    -----
    Piecewise-constant in the two tangential directions (so each parent face
    flux is shared over the child faces it contains in proportion to their
    areas) and linear in the normal direction (so the child target reproduces
    the parent's discrete divergence).  See the module docstring.
    """
    if component not in COMPONENTS:
        raise ConfigurationError(f"unknown velocity component {component!r}")
    field = np.asarray(field, dtype=np.float64)
    expected = parent.component_shape(component)
    if field.shape != expected:
        raise ConfigurationError(
            f"parent field for {component!r} has shape {field.shape}, expected {expected}"
        )
    out = field
    for axis, target in enumerate((x, y, z)):
        target = np.asarray(target, dtype=np.float64)
        faces = parent.edges(axis)
        if _IS_FACE[component][axis]:
            idx, wlo = _normal_map(faces, target)
            out = _apply_normal(out, axis, idx, wlo)
        else:
            out = np.take(out, _tangential_map(faces, target), axis=axis)
    return out


def interpolate_child_fields(
    parent: NestGrid,
    parent_u: np.ndarray,
    parent_v: np.ndarray,
    parent_w: np.ndarray,
    child: NestGrid,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate a full parent velocity field onto the whole child grid.

    Convenience wrapper around :func:`conservative_interpolate`; the writer
    itself only ever interpolates the zone slabs (:func:`slabs_from_parent`).
    """
    out = []
    for component, pf in zip(COMPONENTS, (parent_u, parent_v, parent_w)):
        coords = [child.component_coords(component, ax) for ax in range(3)]
        out.append(conservative_interpolate(parent, pf, component, *coords))
    return tuple(out)


# --------------------------------------------------------------------------- #
# Slab geometry (docs/udales-nesting-spec.md §5, "Index conventions")
# --------------------------------------------------------------------------- #


def slab_dimensions(face: str, component: str) -> Tuple[str, str, str]:
    """CDL dimension names of slab variable ``<component>_<face>``, after ``time``."""
    _check_face(face)
    if component not in COMPONENTS:
        raise ConfigurationError(f"unknown velocity component {component!r}")
    axis = _FACE_AXIS[face]
    zone = "nzh" if _IS_FACE[component][axis] else "nz"
    zdim = "zh" if component == "w" else "zf"
    if axis == 0:  # west/east: decomposed index is y
        span = "yh" if component == "v" else "yf"
    else:  # south/north: decomposed index is x
        span = "xh" if component == "u" else "xf"
    return (span, zdim, zone)


def slab_indices(grid: NestGrid, nzone: int, face: str, component: str) -> np.ndarray:
    """Zero-based global indices, along the face normal, of a slab's zone points.

    Slab index ``m`` (0-based) maps to global index ``m`` on the west/south
    slabs and to ``itot - nzone + m`` / ``jtot - nzone + m`` on the east/north
    slabs, for both cell centres and faces -- i.e. the first slab index is
    always the lowest global index, as the contract requires.
    """
    _check_face(face)
    axis = _FACE_AXIS[face]
    n = nzone + (1 if _IS_FACE[component][axis] else 0)
    if _FACE_UPPER[face]:
        ntot = grid.itot if axis == 0 else grid.jtot
        start = ntot - nzone
    else:
        start = 0
    return start + np.arange(n)


def slab_coordinates(
    grid: NestGrid, nzone: int, face: str, component: str
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Physical ``(x, y, z)`` coordinates of a slab's points, at its own stagger."""
    axis = _FACE_AXIS[face]
    coords = [grid.component_coords(component, ax) for ax in range(3)]
    coords[axis] = coords[axis][slab_indices(grid, nzone, face, component)]
    return tuple(coords)


def _slab_transpose(face: str) -> Tuple[int, int, int]:
    """Axis permutation from ``(x, y, z)`` order to the stored slab order."""
    # west/east store (span_y, z, zone_x); south/north store (span_x, z, zone_y).
    return (1, 2, 0) if _FACE_AXIS[face] == 0 else (0, 2, 1)


def _check_face(face: str) -> None:
    if face not in FACES:
        raise ConfigurationError(f"unknown lateral face {face!r}; expected one of {FACES}")


def slab_shape(grid: NestGrid, nzone: int, face: str, component: str) -> Tuple[int, int, int]:
    """Stored shape (without the time dimension) of ``<component>_<face>``."""
    sizes = {
        "xf": grid.itot, "xh": grid.itot + 1,
        "yf": grid.jtot, "yh": grid.jtot + 1,
        "zf": grid.ktot, "zh": grid.ktot + 1,
        "nz": nzone, "nzh": nzone + 1,
    }
    return tuple(sizes[d] for d in slab_dimensions(face, component))


def slabs_from_fields(
    grid: NestGrid,
    nzone: int,
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
) -> Dict[str, np.ndarray]:
    """Cut the twelve zone slabs out of full child fields, for one time level."""
    _check_nzone(grid, nzone)
    fields = {"u": np.asarray(u, dtype=np.float64),
              "v": np.asarray(v, dtype=np.float64),
              "w": np.asarray(w, dtype=np.float64)}
    for component, arr in fields.items():
        expected = grid.component_shape(component)
        if arr.shape != expected:
            raise ConfigurationError(
                f"child field {component!r} has shape {arr.shape}, expected {expected}"
            )
    slabs: Dict[str, np.ndarray] = {}
    for face in FACES:
        axis = _FACE_AXIS[face]
        for component in COMPONENTS:
            idx = slab_indices(grid, nzone, face, component)
            block = np.take(fields[component], idx, axis=axis)
            slabs[f"{component}_{face}"] = np.ascontiguousarray(
                block.transpose(_slab_transpose(face))
            )
    return slabs


def slabs_from_parent(
    parent: NestGrid,
    parent_u: np.ndarray,
    parent_v: np.ndarray,
    parent_w: np.ndarray,
    child: NestGrid,
    nzone: int,
) -> Dict[str, np.ndarray]:
    """Interpolate the parent onto the twelve zone slabs, for one time level.

    Only the slab points are interpolated, so the cost and the memory are
    proportional to the zone rather than to the child domain.
    """
    _check_nzone(child, nzone)
    parent_fields = dict(zip(COMPONENTS, (parent_u, parent_v, parent_w)))
    slabs: Dict[str, np.ndarray] = {}
    for face in FACES:
        for component in COMPONENTS:
            coords = slab_coordinates(child, nzone, face, component)
            block = conservative_interpolate(
                parent, parent_fields[component], component, *coords
            )
            slabs[f"{component}_{face}"] = np.ascontiguousarray(
                block.transpose(_slab_transpose(face))
            )
    return slabs


def _check_nzone(grid: NestGrid, nzone: int) -> None:
    if nzone < 1:
        raise ConfigurationError(f"nzone must be >= 1, got {nzone}")
    if nzone > min(grid.itot, grid.jtot):
        raise ConfigurationError(
            f"nzone = {nzone} exceeds the child domain ({grid.itot} x {grid.jtot} cells)"
        )


# --------------------------------------------------------------------------- #
# The data container
# --------------------------------------------------------------------------- #


@dataclass
class FaceMasks:
    """Fluid/solid masks on the four lateral boundary faces.

    ``True`` marks a **fluid** face; ``None`` means "all fluid".  The west and
    east masks have shape ``(jtot, ktot)``, south and north ``(itot, ktot)`` --
    the shape of the boundary face itself, matching ``IIu``/``IIv`` in the
    solver.
    """

    west: Optional[np.ndarray] = None
    east: Optional[np.ndarray] = None
    south: Optional[np.ndarray] = None
    north: Optional[np.ndarray] = None

    def get(self, grid: NestGrid, face: str) -> np.ndarray:
        """Return the mask for ``face`` as a float array of 1.0 (fluid) / 0.0."""
        _check_face(face)
        expected = (grid.jtot, grid.ktot) if _FACE_AXIS[face] == 0 else (grid.itot, grid.ktot)
        mask = getattr(self, face)
        if mask is None:
            return np.ones(expected, dtype=np.float64)
        mask = np.asarray(mask)
        if mask.shape != expected:
            raise ConfigurationError(
                f"{face} face mask has shape {mask.shape}, expected {expected}"
            )
        return mask.astype(np.float64)


_ALL_FLUID = FaceMasks()


@dataclass
class NestingData:
    """Everything that goes into ``nesting.inp.<expnr>.nc``.

    ``slabs`` maps each of the twelve variable names (``u_west`` ...
    ``w_north``) to an array whose first dimension is time and whose remaining
    dimensions are those of :func:`slab_dimensions`.
    """

    grid: NestGrid
    nzone: int
    times: np.ndarray
    slabs: Dict[str, np.ndarray]
    rhobf: Optional[np.ndarray] = None
    rhobh: Optional[np.ndarray] = None
    net_volume_flux: Optional[np.ndarray] = None
    divergence_corrected: bool = False
    parent_model: str = "unknown"
    parent_dx: float = 0.0
    parent_dt: float = 0.0
    child_origin_x: float = 0.0
    child_origin_y: float = 0.0
    child_dt: Optional[float] = None
    rotation_deg: float = 0.0
    created: str = ""
    creator: str = ""
    tool_version: str = TOOL_VERSION

    def __post_init__(self) -> None:
        self.nzone = int(self.nzone)
        _check_nzone(self.grid, self.nzone)
        self.times = np.asarray(self.times, dtype=np.float64).reshape(-1)
        if self.rhobf is None:
            self.rhobf = np.ones(self.grid.ktot, dtype=np.float64)
        if self.rhobh is None:
            self.rhobh = np.ones(self.grid.ktot + 1, dtype=np.float64)
        self.rhobf = np.asarray(self.rhobf, dtype=np.float64).reshape(-1)
        self.rhobh = np.asarray(self.rhobh, dtype=np.float64).reshape(-1)
        if self.rhobf.size != self.grid.ktot or self.rhobh.size != self.grid.ktot + 1:
            raise ConfigurationError(
                f"rhobf/rhobh must have {self.grid.ktot}/{self.grid.ktot + 1} elements"
            )
        if self.net_volume_flux is None:
            self.net_volume_flux = np.zeros(self.ntime, dtype=np.float64)
        self.net_volume_flux = np.asarray(self.net_volume_flux, dtype=np.float64).reshape(-1)
        if self.net_volume_flux.size != self.ntime:
            raise ConfigurationError("net_volume_flux must have one value per time level")
        missing = [name for name in SLAB_VARIABLES if name not in self.slabs]
        if missing:
            raise ConfigurationError(f"missing slab variables: {', '.join(missing)}")
        for face in FACES:
            for component in COMPONENTS:
                name = f"{component}_{face}"
                arr = np.ascontiguousarray(np.asarray(self.slabs[name], dtype=np.float64))
                expected = (self.ntime,) + slab_shape(self.grid, self.nzone, face, component)
                if arr.shape != expected:
                    raise ConfigurationError(
                        f"{name} has shape {arr.shape}, expected {expected}"
                    )
                self.slabs[name] = arr

    @property
    def ntime(self) -> int:
        return int(self.times.size)

    def copy(self) -> "NestingData":
        """Deep copy, so a correction can be applied without touching the original."""
        return NestingData(
            grid=self.grid,
            nzone=self.nzone,
            times=self.times.copy(),
            slabs={k: v.copy() for k, v in self.slabs.items()},
            rhobf=self.rhobf.copy(),
            rhobh=self.rhobh.copy(),
            net_volume_flux=self.net_volume_flux.copy(),
            divergence_corrected=self.divergence_corrected,
            parent_model=self.parent_model,
            parent_dx=self.parent_dx,
            parent_dt=self.parent_dt,
            child_origin_x=self.child_origin_x,
            child_origin_y=self.child_origin_y,
            child_dt=self.child_dt,
            rotation_deg=self.rotation_deg,
            created=self.created,
            creator=self.creator,
            tool_version=self.tool_version,
        )


def nesting_data_from_parent(
    parent: NestGrid,
    child: NestGrid,
    nzone: int,
    times: Sequence[float],
    fields: Sequence[Tuple[np.ndarray, np.ndarray, np.ndarray]],
    **kwargs: Any,
) -> NestingData:
    """Build a :class:`NestingData` by conservative interpolation of a parent.

    ``fields`` is a sequence, one entry per time level, of ``(u, v, w)`` parent
    arrays at their own staggers.  Extra keyword arguments are passed to
    :class:`NestingData` (provenance attributes, ``rhobf``, ...).
    """
    times = np.asarray(times, dtype=np.float64).reshape(-1)
    if len(fields) != times.size:
        raise ConfigurationError(
            f"got {len(fields)} field sets for {times.size} times"
        )
    per_time = [slabs_from_parent(parent, *fields[n], child=child, nzone=nzone)
                for n in range(times.size)]
    slabs = {name: np.stack([s[name] for s in per_time], axis=0) for name in SLAB_VARIABLES}
    kwargs.setdefault("parent_dx", float(np.min(np.diff(parent.xh))))
    kwargs.setdefault("child_origin_x", float(child.xh[0]))
    kwargs.setdefault("child_origin_y", float(child.yh[0]))
    if times.size > 1:
        kwargs.setdefault("parent_dt", float(np.min(np.diff(times))))
    return NestingData(grid=child, nzone=nzone, times=times, slabs=slabs, **kwargs)


# --------------------------------------------------------------------------- #
# Divergence correction (design §3.1, §3.3)
# --------------------------------------------------------------------------- #


def _boundary_slab_index(data: NestingData, face: str) -> int:
    """Slab index of the domain-boundary face within its own slab."""
    return data.nzone if _FACE_UPPER[face] else 0


def boundary_faces(data: NestingData) -> Dict[str, np.ndarray]:
    """Views of the normal velocity on the four lateral **domain boundary** faces.

    Each entry has shape ``(ntime, n_span, ktot)``: ``(ntime, jtot, ktot)`` for
    west/east and ``(ntime, itot, ktot)`` for south/north.  These are views into
    ``data.slabs``, so writing to them corrects the file in place.
    """
    out = {}
    for face in FACES:
        name = f"{_FACE_NORMAL_COMPONENT[face]}_{face}"
        out[face] = data.slabs[name][:, :, :, _boundary_slab_index(data, face)]
    return out


def _face_weights(data: NestingData, face: str) -> np.ndarray:
    """rho * dA on one lateral boundary face, shape ``(n_span, ktot)``."""
    grid = data.grid
    span = grid.dy if _FACE_AXIS[face] == 0 else grid.dx
    return span[:, None] * (grid.dzf * data.rhobf)[None, :]


def net_volume_flux(data: NestingData, masks: Optional[FaceMasks] = None) -> np.ndarray:
    """Net volume flux :math:`\\Phi` through the lateral boundary, per time level.

    .. math::
        \\Phi(t) = \\sum_{\\rm fluid\\ lateral\\ faces} \\rho\\, u_n\\, dA,

    outward positive, summed over **fluid** faces only (design §3.1).  With
    ``rhobf = 1`` (the Boussinesq default) this is exactly the expression of
    design §3.1.  Returns an array of shape ``(ntime,)``.
    """
    masks = masks or _ALL_FLUID
    faces = boundary_faces(data)
    phi = np.zeros(data.ntime, dtype=np.float64)
    for face in FACES:
        weight = _face_weights(data, face) * masks.get(data.grid, face)
        phi += _FACE_SIGN[face] * np.tensordot(faces[face], weight, axes=([1, 2], [0, 1]))
    return phi


def fluid_face_area(data: NestingData, masks: Optional[FaceMasks] = None) -> float:
    """Total ``rho``-weighted **fluid** area of the four lateral boundary faces.

    This is the normalisation of the correction: adding a uniform outward
    normal-velocity increment ``delta`` changes :math:`\\Phi` by
    ``delta * fluid_face_area``.
    """
    masks = masks or _ALL_FLUID
    return float(sum(
        np.sum(_face_weights(data, face) * masks.get(data.grid, face)) for face in FACES
    ))


def apply_divergence_correction(
    data: NestingData, masks: Optional[FaceMasks] = None
) -> np.ndarray:
    """Make :math:`\\Phi = 0` on every time level, in place.

    A uniform outward normal-velocity increment is added over the fluid lateral
    faces -- so the correction *flux* is distributed in proportion to fluid face
    area, and differences between faces are untouched (design §3.3, DALES's
    ``openboundary_divcorr``).  Solid faces are left alone.

    The **pre-correction** residual is stored in ``data.net_volume_flux`` and
    returned; ``data.divergence_corrected`` is set.
    """
    masks = masks or _ALL_FLUID
    area = fluid_face_area(data, masks)
    if area <= 0.0:
        raise ConfigurationError(
            "the lateral boundary has no fluid area; cannot correct the volume flux"
        )
    residual = net_volume_flux(data, masks)
    delta = -residual / area
    faces = boundary_faces(data)
    for face in FACES:
        mask = masks.get(data.grid, face)
        faces[face] += _FACE_SIGN[face] * delta[:, None, None] * mask[None, :, :]
    data.net_volume_flux = residual
    data.divergence_corrected = True
    return residual


# --------------------------------------------------------------------------- #
# Refinement guard
# --------------------------------------------------------------------------- #


def refinement_ratios(data: NestingData) -> Tuple[Optional[float], Optional[float]]:
    """Return ``(spatial, temporal)`` parent-to-child refinement ratios.

    The spatial ratio is ``parent_dx / min(dx_child, dy_child)``; the temporal
    ratio is ``parent_dt / child_dt``.  Either is ``None`` when the file does
    not carry enough information to compute it (``parent_dx``/``parent_dt``
    unset, or ``child_dt`` not supplied by the caller).
    """
    dmin = min(float(np.min(data.grid.dx)), float(np.min(data.grid.dy)))
    spatial = data.parent_dx / dmin if data.parent_dx and dmin > 0.0 else None
    temporal = None
    if data.parent_dt and data.child_dt:
        temporal = data.parent_dt / data.child_dt
    return spatial, temporal


def check_refinement(
    data: NestingData, override: bool = False
) -> Tuple[Optional[float], Optional[float]]:
    """Refuse ratios beyond the validated range (design §10.4 V5) unless overridden.

    Raises :class:`NestingRefinementError` when the spatial refinement ratio
    exceeds :data:`MAX_SPATIAL_REFINEMENT` or the temporal one exceeds
    :data:`MAX_TEMPORAL_REFINEMENT`, unless ``override`` is ``True``.
    """
    spatial, temporal = refinement_ratios(data)
    if override:
        return spatial, temporal
    tol = 1.0 + 1.0e-12
    if spatial is not None and spatial > MAX_SPATIAL_REFINEMENT * tol:
        raise NestingRefinementError(
            f"spatial refinement ratio {spatial:.3g} exceeds the supported maximum of "
            f"{MAX_SPATIAL_REFINEMENT:g} (parent_dx = {data.parent_dx:g} m); "
            "pass override=True to write anyway"
        )
    if temporal is not None and temporal > MAX_TEMPORAL_REFINEMENT * tol:
        raise NestingRefinementError(
            f"temporal refinement ratio {temporal:.3g} exceeds the supported maximum of "
            f"{MAX_TEMPORAL_REFINEMENT:g} (parent_dt = {data.parent_dt:g} s, "
            f"child_dt = {data.child_dt:g} s); pass override=True to write anyway"
        )
    return spatial, temporal


# --------------------------------------------------------------------------- #
# Writer / reader (docs/udales-nesting-spec.md §5)
# --------------------------------------------------------------------------- #


def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _default_creator() -> str:
    try:
        return getpass.getuser()
    except Exception:  # pragma: no cover - no passwd entry (batch nodes)
        return os.environ.get("USER", "unknown")


def _global_attributes(data: NestingData) -> Dict[str, Any]:
    grid = data.grid
    attrs: Dict[str, Any] = {
        "Conventions": "CF-1.8",
        "udales_nesting_schema": np.int32(SCHEMA_VERSION),
        "divergence_corrected": np.int32(1 if data.divergence_corrected else 0),
        "itot": np.int32(grid.itot),
        "jtot": np.int32(grid.jtot),
        "ktot": np.int32(grid.ktot),
        "nzone": np.int32(data.nzone),
        "xlen": np.float64(grid.xlen),
        "ylen": np.float64(grid.ylen),
        "parent_model": str(data.parent_model),
        "parent_dx": np.float64(data.parent_dx),
        "parent_dt": np.float64(data.parent_dt),
        "child_origin_x": np.float64(data.child_origin_x),
        "child_origin_y": np.float64(data.child_origin_y),
        "rotation_deg": np.float64(data.rotation_deg),
        "created": data.created or _utc_now(),
        "creator": data.creator or _default_creator(),
        "tool_version": data.tool_version or TOOL_VERSION,
    }
    if data.child_dt is not None:
        # Optional extension: lets the refinement guard run on a round-tripped
        # file.  Ignored by readers that do not know about it.
        attrs["child_dt"] = np.float64(data.child_dt)
    return attrs


def write_nesting_file(
    path: os.PathLike | str,
    data: NestingData,
    override: bool = False,
    backend: Optional[str] = None,
) -> Path:
    """Write ``data`` to ``nesting.inp.<expnr>.nc``, exactly per the contract.

    ``backend`` is ``'netcdf'`` (default, or implied by a ``.nc`` suffix) or
    ``'raw'`` (a flat stream of the slab arrays plus a JSON sidecar, design
    §6.3).  Both back-ends store bit-identical values.

    Refuses to write when the refinement ratios are outside the validated range
    (see :func:`check_refinement`) unless ``override`` is ``True``.
    """
    path = Path(path)
    backend = backend or ("raw" if path.suffix in (".dat", ".bin") else "netcdf")
    check_refinement(data, override=override)
    for name, arr in data.slabs.items():
        if not np.all(np.isfinite(arr)):
            raise NestingSchemaError(f"{name} contains non-finite values; NaN is an error")
    if backend == "raw":
        return _write_raw(path, data)
    if backend != "netcdf":
        raise ConfigurationError(f"unknown nesting file backend {backend!r}")
    return _write_netcdf(path, data)


def _write_netcdf(path: Path, data: NestingData) -> Path:
    Dataset = _import_dataset()
    grid = data.grid
    with Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("time", None)
        ds.createDimension("zf", grid.ktot)
        ds.createDimension("zh", grid.ktot + 1)
        ds.createDimension("xf", grid.itot)
        ds.createDimension("xh", grid.itot + 1)
        ds.createDimension("yf", grid.jtot)
        ds.createDimension("yh", grid.jtot + 1)
        ds.createDimension("nz", data.nzone)
        ds.createDimension("nzh", data.nzone + 1)

        var = ds.createVariable("time", "f8", ("time",))
        var.units = "s"
        var.long_name = "parent time level"
        var[:] = data.times
        for name in _COORD_VARIABLES:
            var = ds.createVariable(name, "f8", (name,))
            var.units = "m"
            var[:] = getattr(grid, name)
        for name, dim, values in (("rhobf", "zf", data.rhobf), ("rhobh", "zh", data.rhobh)):
            var = ds.createVariable(name, "f8", (dim,))
            var.units = "kg m-3"
            var[:] = values
        var = ds.createVariable("net_volume_flux", "f8", ("time",))
        var.units = "kg s-1"
        var.long_name = "net volume flux through the lateral boundary before correction"
        var[:] = data.net_volume_flux

        for face in FACES:
            for component in COMPONENTS:
                name = f"{component}_{face}"
                dims = ("time",) + slab_dimensions(face, component)
                var = ds.createVariable(name, "f8", dims)
                var.stagger = STAGGER[component]
                var.units = "m s-1"
                var[:] = data.slabs[name]

        for key, value in _global_attributes(data).items():
            ds.setncattr(key, value)
    return path


def _raw_sidecar_path(path: Path) -> Path:
    return path.with_suffix(path.suffix + ".json") if path.suffix else path.with_suffix(".json")


def _write_raw(path: Path, data: NestingData) -> Path:
    """Raw stream back-end: one contiguous ``<f8`` blob plus a JSON sidecar."""
    grid = data.grid
    attrs = {k: (int(v) if isinstance(v, np.integer)
                 else float(v) if isinstance(v, np.floating) else v)
             for k, v in _global_attributes(data).items()}
    sidecar = {
        "udales_nesting_schema": SCHEMA_VERSION,
        "attributes": attrs,
        "coordinates": {name: getattr(grid, name).tolist() for name in _COORD_VARIABLES},
        "time": data.times.tolist(),
        "rhobf": data.rhobf.tolist(),
        "rhobh": data.rhobh.tolist(),
        "net_volume_flux": data.net_volume_flux.tolist(),
        "dtype": "<f8",
        "variables": [
            {
                "name": f"{component}_{face}",
                "dimensions": ("time",) + slab_dimensions(face, component),
                "shape": list(data.slabs[f"{component}_{face}"].shape),
                "stagger": STAGGER[component],
            }
            for face in FACES for component in COMPONENTS
        ],
    }
    _raw_sidecar_path(path).write_text(json.dumps(sidecar, indent=1))
    with open(path, "wb") as stream:
        for face in FACES:
            for component in COMPONENTS:
                arr = data.slabs[f"{component}_{face}"]
                stream.write(np.ascontiguousarray(arr, dtype="<f8").tobytes(order="C"))
    return path


def _read_raw(path: Path) -> NestingData:
    sidecar = json.loads(_raw_sidecar_path(path).read_text())
    if int(sidecar.get("udales_nesting_schema", -1)) != SCHEMA_VERSION:
        raise NestingSchemaError(
            f"{path}: sidecar schema {sidecar.get('udales_nesting_schema')} != {SCHEMA_VERSION}"
        )
    coords = sidecar["coordinates"]
    grid = NestGrid(**{name: np.asarray(coords[name], dtype=np.float64)
                       for name in _COORD_VARIABLES})
    blob = np.fromfile(path, dtype="<f8")
    slabs: Dict[str, np.ndarray] = {}
    offset = 0
    for spec in sidecar["variables"]:
        shape = tuple(spec["shape"])
        size = int(np.prod(shape))
        slabs[spec["name"]] = np.ascontiguousarray(
            blob[offset:offset + size].reshape(shape).astype(np.float64)
        )
        offset += size
    if offset != blob.size:
        raise NestingSchemaError(f"{path}: {blob.size - offset} trailing values in the stream")
    attrs = sidecar["attributes"]
    return NestingData(
        grid=grid,
        nzone=int(attrs["nzone"]),
        times=np.asarray(sidecar["time"], dtype=np.float64),
        slabs=slabs,
        rhobf=np.asarray(sidecar["rhobf"], dtype=np.float64),
        rhobh=np.asarray(sidecar["rhobh"], dtype=np.float64),
        net_volume_flux=np.asarray(sidecar["net_volume_flux"], dtype=np.float64),
        divergence_corrected=bool(int(attrs["divergence_corrected"])),
        parent_model=str(attrs["parent_model"]),
        parent_dx=float(attrs["parent_dx"]),
        parent_dt=float(attrs["parent_dt"]),
        child_origin_x=float(attrs["child_origin_x"]),
        child_origin_y=float(attrs["child_origin_y"]),
        child_dt=float(attrs["child_dt"]) if "child_dt" in attrs else None,
        rotation_deg=float(attrs["rotation_deg"]),
        created=str(attrs["created"]),
        creator=str(attrs["creator"]),
        tool_version=str(attrs["tool_version"]),
    )


def validate_nesting_file(path: os.PathLike | str) -> Dict[str, Any]:
    """Check a NetCDF nesting file against the schema of the contract.

    Raises :class:`NestingSchemaError` naming the first offending item; returns
    the global attributes on success.  Checked: schema version, every required
    dimension and its size, every required variable with its exact dimension
    order, the ``stagger`` tag of each slab variable, every required global
    attribute, ``rotation_deg == 0``, monotone coordinates, and the absence of
    NaN/missing values.
    """
    Dataset = _import_dataset()
    path = Path(path)
    with Dataset(path, "r") as ds:
        attrs = {key: ds.getncattr(key) for key in ds.ncattrs()}
        missing = [key for key in REQUIRED_GLOBAL_ATTRIBUTES if key not in attrs]
        if missing:
            raise NestingSchemaError(
                f"{path.name}: missing required global attribute(s): {', '.join(missing)}"
            )
        schema = int(attrs["udales_nesting_schema"])
        if schema != SCHEMA_VERSION:
            raise NestingSchemaError(
                f"{path.name}: udales_nesting_schema = {schema}, expected {SCHEMA_VERSION}"
            )
        if float(attrs["rotation_deg"]) != 0.0:
            raise NestingSchemaError(
                f"{path.name}: rotation_deg = {attrs['rotation_deg']}, only 0 is supported in v1"
            )
        itot, jtot = int(attrs["itot"]), int(attrs["jtot"])
        ktot, nzone = int(attrs["ktot"]), int(attrs["nzone"])
        expected_dims = {
            "zf": ktot, "zh": ktot + 1, "xf": itot, "xh": itot + 1,
            "yf": jtot, "yh": jtot + 1, "nz": nzone, "nzh": nzone + 1,
        }
        if "time" not in ds.dimensions:
            raise NestingSchemaError(f"{path.name}: missing dimension 'time'")
        if not ds.dimensions["time"].isunlimited():
            raise NestingSchemaError(f"{path.name}: dimension 'time' must be UNLIMITED")
        for name, size in expected_dims.items():
            if name not in ds.dimensions:
                raise NestingSchemaError(f"{path.name}: missing dimension {name!r}")
            actual = len(ds.dimensions[name])
            if actual != size:
                raise NestingSchemaError(
                    f"{path.name}: dimension {name} = {actual}, expected {size}"
                )

        expected_vars: Dict[str, Tuple[str, ...]] = {"time": ("time",),
                                                     "net_volume_flux": ("time",),
                                                     "rhobf": ("zf",), "rhobh": ("zh",)}
        for name in _COORD_VARIABLES:
            expected_vars[name] = (name,)
        for face in FACES:
            for component in COMPONENTS:
                expected_vars[f"{component}_{face}"] = ("time",) + slab_dimensions(face, component)
        for name, dims in expected_vars.items():
            if name not in ds.variables:
                raise NestingSchemaError(f"{path.name}: missing variable {name!r}")
            var = ds.variables[name]
            if tuple(var.dimensions) != dims:
                raise NestingSchemaError(
                    f"{path.name}: variable {name} has dimensions {tuple(var.dimensions)}, "
                    f"expected {dims}"
                )
            values = np.asarray(var[:], dtype=np.float64)
            if not np.all(np.isfinite(values)):
                raise NestingSchemaError(
                    f"{path.name}: variable {name} has missing or NaN values, which are an error"
                )
        for face in FACES:
            for component in COMPONENTS:
                name = f"{component}_{face}"
                var = ds.variables[name]
                tag = getattr(var, "stagger", None)
                if tag != STAGGER[component]:
                    raise NestingSchemaError(
                        f"{path.name}: variable {name} has stagger {tag!r}, "
                        f"expected {STAGGER[component]!r}"
                    )
        for name in ("xh", "yh", "zh", "xf", "yf", "zf"):
            values = np.asarray(ds.variables[name][:], dtype=np.float64)
            if values.size > 1 and not np.all(np.diff(values) > 0.0):
                raise NestingSchemaError(f"{path.name}: coordinate {name} is not increasing")
        xlen = float(np.asarray(ds.variables["xh"][:])[-1] - np.asarray(ds.variables["xh"][:])[0])
        ylen = float(np.asarray(ds.variables["yh"][:])[-1] - np.asarray(ds.variables["yh"][:])[0])
        for label, from_attr, from_coord in (("xlen", float(attrs["xlen"]), xlen),
                                             ("ylen", float(attrs["ylen"]), ylen)):
            if not np.isclose(from_attr, from_coord, rtol=1e-12, atol=1e-12):
                raise NestingSchemaError(
                    f"{path.name}: attribute {label} = {from_attr!r} disagrees with the "
                    f"coordinate variable ({from_coord!r})"
                )
    return attrs


def read_nesting_file(
    path: os.PathLike | str, validate: bool = True, backend: Optional[str] = None
) -> NestingData:
    """Read a nesting file back into a :class:`NestingData`.

    Validates against the schema first unless ``validate`` is ``False``.
    """
    path = Path(path)
    backend = backend or ("raw" if path.suffix in (".dat", ".bin") else "netcdf")
    if backend == "raw":
        return _read_raw(path)
    if backend != "netcdf":
        raise ConfigurationError(f"unknown nesting file backend {backend!r}")
    if validate:
        validate_nesting_file(path)
    Dataset = _import_dataset()
    with Dataset(path, "r") as ds:
        attrs = {key: ds.getncattr(key) for key in ds.ncattrs()}
        grid = NestGrid(**{name: np.asarray(ds.variables[name][:], dtype=np.float64)
                           for name in _COORD_VARIABLES})
        slabs = {name: np.ascontiguousarray(np.asarray(ds.variables[name][:], dtype=np.float64))
                 for name in SLAB_VARIABLES}
        return NestingData(
            grid=grid,
            nzone=int(attrs["nzone"]),
            times=np.asarray(ds.variables["time"][:], dtype=np.float64),
            slabs=slabs,
            rhobf=np.asarray(ds.variables["rhobf"][:], dtype=np.float64),
            rhobh=np.asarray(ds.variables["rhobh"][:], dtype=np.float64),
            net_volume_flux=np.asarray(ds.variables["net_volume_flux"][:], dtype=np.float64),
            divergence_corrected=bool(int(attrs["divergence_corrected"])),
            parent_model=str(attrs["parent_model"]),
            parent_dx=float(attrs["parent_dx"]),
            parent_dt=float(attrs["parent_dt"]),
            child_origin_x=float(attrs["child_origin_x"]),
            child_origin_y=float(attrs["child_origin_y"]),
            child_dt=float(attrs["child_dt"]) if "child_dt" in attrs else None,
            rotation_deg=float(attrs["rotation_deg"]),
            created=str(attrs["created"]),
            creator=str(attrs["creator"]),
            tool_version=str(attrs["tool_version"]),
        )


# --------------------------------------------------------------------------- #
# Analytic field generator (design §10.5; drives the Fortran tests U15--U22)
# --------------------------------------------------------------------------- #

#: Coefficients ``(a, b, c, d, e)`` of the analytic field, per velocity
#: component, in SI units (1/m, 1/m, 1/m, 1/s, 1/(m^2 s)).  They differ between
#: components and ``a != b``, so a swapped or transposed index cannot reproduce
#: the field.  The Fortran unit test must use exactly these values.
ANALYTIC_COEFFS: Dict[str, Tuple[float, float, float, float, float]] = {
    "u": (0.017, 0.011, 0.023, 0.05, 1.0e-5),
    "v": (0.013, 0.019, 0.029, 0.03, -7.0e-6),
    "w": (0.023, 0.007, 0.013, 0.07, 4.0e-6),
}


def analytic_field(
    component: str,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: float,
    coeffs: Optional[Mapping[str, Sequence[float]]] = None,
) -> np.ndarray:
    """Evaluate the reference analytic velocity field.

    .. math::
        f(x, y, z, t) = \\sin(a x + b y)\\,\\cos(c z)\\,(1 + d t) + e\\,x y z

    with ``(a, b, c, d, e)`` taken from ``coeffs[component]`` (default
    :data:`ANALYTIC_COEFFS`).  ``x``, ``y`` and ``z`` are broadcast against one
    another -- pass ``*np.ix_(xs, ys, zs)`` for a 3-D block.  The field is
    non-separable in ``(x, y)`` and differs between components, which is what
    lets the solver-side test detect a transposed or off-by-one index.
    """
    table = dict(ANALYTIC_COEFFS)
    if coeffs is not None:
        table.update({k: tuple(v) for k, v in coeffs.items()})
    if component not in table:
        raise ConfigurationError(f"no analytic coefficients for component {component!r}")
    a, b, c, d, e = (float(v) for v in table[component])
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    return np.sin(a * x + b * y) * np.cos(c * z) * (1.0 + d * float(t)) + e * (x * y * z)


def analytic_slabs(
    grid: NestGrid,
    nzone: int,
    times: Sequence[float],
    coeffs: Optional[Mapping[str, Sequence[float]]] = None,
) -> Dict[str, np.ndarray]:
    """Evaluate :func:`analytic_field` on the twelve zone slabs, at their staggers."""
    times = np.asarray(times, dtype=np.float64).reshape(-1)
    _check_nzone(grid, nzone)
    slabs: Dict[str, np.ndarray] = {}
    for face in FACES:
        perm = _slab_transpose(face)
        for component in COMPONENTS:
            xs, ys, zs = slab_coordinates(grid, nzone, face, component)
            block = np.stack(
                [analytic_field(component, *np.ix_(xs, ys, zs), t, coeffs).transpose(perm)
                 for t in times],
                axis=0,
            )
            slabs[f"{component}_{face}"] = np.ascontiguousarray(block)
    return slabs


def write_analytic_nesting_file(
    path: os.PathLike | str,
    grid: NestGrid,
    times: Sequence[float],
    nzone: int,
    coeffs: Optional[Mapping[str, Sequence[float]]] = None,
    correct_divergence: bool = False,
    override: bool = True,
    backend: Optional[str] = None,
    **attributes: Any,
) -> NestingData:
    """Write a nesting file whose velocities are the analytic field of §10.5.

    Used as the fixture for the in-solver tests ``TEST_NESTING_IO`` (U15--U22)
    and, uncorrected, for the flux assertion test U27.  ``correct_divergence``
    is off by default so that every stored value is exactly
    :func:`analytic_field` -- switching it on perturbs the four boundary faces.
    The refinement guard is off by default (``override=True``) because the
    fixture grids carry no meaningful parent metadata.

    Returns the :class:`NestingData` that was written.
    """
    attributes.setdefault("parent_model", "analytic")
    attributes.setdefault("child_origin_x", float(grid.xh[0]))
    attributes.setdefault("child_origin_y", float(grid.yh[0]))
    data = NestingData(
        grid=grid,
        nzone=nzone,
        times=np.asarray(times, dtype=np.float64),
        slabs=analytic_slabs(grid, nzone, times, coeffs),
        **attributes,
    )
    if correct_divergence:
        apply_divergence_correction(data)
    else:
        data.net_volume_flux = net_volume_flux(data)
    write_nesting_file(path, data, override=override, backend=backend)
    return data
