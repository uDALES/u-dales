"""Preprocessing for one-way nesting (``nesting.inp.<expnr>.nc``).

Implements work package W1-PY of the nesting feature:

* **conservative flux interpolation** of a parent velocity field onto the child
  faces (design ``docs/udales-nesting-design.md`` §1.3);
* the offline **divergence correction** that makes the net volume flux through
  the child's lateral boundary vanish on every stored time level (§3.1, §3.3);
* the **writer/reader** for the file format of ``docs/udales-nesting-spec.md`` §5,
  with schema validation and a refinement-ratio guard;
* the OPTIONAL **full-domain initial condition** of schema 2 and its
  **projection** onto the discretely solenoidal subspace, so that a cold start
  with ``nest_linitfromparent`` begins divergence free (design §10.6 item 4);
* an **analytic-field generator** used by the in-solver unit tests (U15--U22).

Schema versions
---------------
Version 1 is the original file: zone slabs, a grid, provenance, and the
*pre*-correction ``net_volume_flux``.  Version 2 adds

* ``flux_residual(time)`` -- the residual of the data **as stored**, so the
  solver can validate every stored level without re-reading the boundary slabs
  (design §10.6 item 3), together with the ``fluid_lateral_area`` the residual
  was summed over, so a mask mismatch between writer and solver is caught rather
  than trusted;
* ``u_init``/``v_init``/``w_init`` -- the optional full-domain initial condition
  at ``times[0]``, flagged by ``has_initial_condition``.

Both versions are read, by this module and by ``src/nesting_read.f90``.

Interpolation
-------------
The prolongation is the tensor product of

* a conservative reconstruction in the two directions **tangential** to the
  face -- by default piecewise **linear** within each parent cell, with the
  slope from central differences of the neighbouring parent cells (one-sided
  at the ends and next to a solid cell, zero when both neighbours are solid),
  no limiter; ``prolongation="constant"`` gives the piecewise-constant
  distribution of the original scheme.  Either way the child face values
  integrate to the parent face flux exactly, because the child cells nest in
  the parent cells (:func:`check_alignment`) and a linear term integrates to
  zero over the parent cell -- the defining property, design §1.3; and
* linear interpolation in the direction **normal** to the face, between the two
  bracketing parent faces.

The sum of the child fluxes over a parent face equals the parent face flux to
round-off, so :math:`\Phi` and the stored ``flux_residual`` are the same for
both reconstructions, and the divergence integrated over a parent cell is
preserved, which is what the design requires.  The **constant** scheme in
addition reproduces the parent's discrete divergence *cell by cell* and the
parent value exactly on a coplanar child face; the **linear** scheme gives
those up in exchange for removing the staircase a piecewise-constant
prolongation leaves in a sheared profile (a two-level sawtooth of
:math:`\pm\,\partial_z u\,\Delta z_c/2` at ratio 2, which the V0
refinement run saw as :math:`\pm 0.1 u_*` in the child's interior mean).
The interpolant stays linear in the parent data, as the time interpolant
must (spec §4, ``nest_timeinterp``).

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
import logging
import os
from dataclasses import dataclass, field
from datetime import datetime, timezone
from types import SimpleNamespace
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

from exceptions import ConfigurationError, DataFormatError, DependencyError

from ._section import Section, SectionSpec

logger = logging.getLogger(__name__)

__all__ = [
    "ANALYTIC_COEFFS",
    "CADENCE_COURANT_MAX",
    "COMPONENTS",
    "CORRECTION_WARN_FRACTION",
    "FACES",
    "FLUX_UNITS",
    "INIT_VARIABLES",
    "MAX_SPATIAL_REFINEMENT",
    "MAX_TEMPORAL_REFINEMENT",
    "OPTIONAL_GLOBAL_ATTRIBUTES",
    "PARENT_DT_RTOL",
    "PROLONGATIONS",
    "DEFAULT_PROLONGATION",
    "REQUIRED_GLOBAL_ATTRIBUTES",
    "REQUIRED_GLOBAL_ATTRIBUTES_V2",
    "SCHEMA_VERSION",
    "SLAB_VARIABLES",
    "SPEC",
    "STAGGER",
    "SUPPORTED_SCHEMA_VERSIONS",
    "TOOL_VERSION",
    "FaceMasks",
    "NestGrid",
    "NestingAlignmentError",
    "NestingData",
    "NestingRefinementError",
    "NestingSchemaError",
    "NestingWriter",
    "analytic_field",
    "analytic_initial_fields",
    "analytic_slabs",
    "apply_divergence_correction",
    "boundary_faces",
    "cadence_courant",
    "check_alignment",
    "check_refinement",
    "check_time_axis",
    "conservative_interpolate",
    "correction_report",
    "discrete_divergence",
    "face_masks_from_ibm",
    "fluid_face_area",
    "fluid_lateral_area",
    "init_dimensions",
    "initial_fields_from_fields",
    "initial_fields_from_parent",
    "interpolate_child_fields",
    "net_volume_flux",
    "nesting_data_from_parent",
    "nesting_diagnostics",
    "nesting_filename",
    "project_initial_condition",
    "read_nesting_file",
    "refinement_ratios",
    "refinement_ratios_by_axis",
    "refinement_verdict",
    "slab_coordinates",
    "slab_dimensions",
    "slab_indices",
    "slab_shape",
    "slabs_from_fields",
    "stagger_masks_from_ibm",
    "slabs_from_parent",
    "stored_coordinates",
    "sync_initial_condition",
    "validate_nesting_file",
    "verify_stored_residual",
    "write_analytic_nesting_file",
    "write_nesting_file",
]

# --------------------------------------------------------------------------- #
# Contract constants (docs/udales-nesting-spec.md §5)
# --------------------------------------------------------------------------- #

#: Schema this writer emits.  Version 2 adds the per-time-level
#: post-correction ``flux_residual`` and the OPTIONAL full-domain initial
#: condition (``u_init``/``v_init``/``w_init``); version 1 files have neither
#: and are still read, by both this module and ``src/nesting_read.f90``.
SCHEMA_VERSION = 2

#: Schema versions this module can read.
SUPPORTED_SCHEMA_VERSIONS = (1, 2)

TOOL_VERSION = "udprep.nesting/2.0"

FACES = ("west", "east", "south", "north")
COMPONENTS = ("u", "v", "w")
SLAB_VARIABLES = tuple(f"{c}_{f}" for f in FACES for c in COMPONENTS)

#: Full-domain initial-condition variables (schema 2, optional).
INIT_VARIABLES = tuple(f"{c}_init" for c in COMPONENTS)

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

#: The largest spatial refinement ratio the validation campaign has exercised,
#: **not** a limit of the scheme or of the interpolation.  V0 ran r = 2 and 4
#: end to end; nothing in `conservative_interpolate` assumes a bound, and its
#: unit tests cover larger ratios.  The writer refuses more than this so that a
#: production case cannot quietly run outside what has been measured; raise it
#: together with evidence at the new ratio, not on its own.
MAX_SPATIAL_REFINEMENT = 4.0

#: Withdrawn as a physical criterion by the C0 cadence study: what matters is
#: the dump Courant number `C_dump` (see `CADENCE_COURANT_MAX`), not the ratio
#: of the boundary cadence to the child's timestep.  Kept only as a coarse
#: sanity bound on absurd inputs.
MAX_TEMPORAL_REFINEMENT = 30.0

#: Boundary-cadence criterion: the writer warns when
#: ``C_dump = max|u_n| parent_dt / parent_dx`` exceeds this -- a feature then
#: crosses more than two parent cells between stored levels, and the linear
#: (or Hermite) time interpolation of the boundary cannot represent it.
CADENCE_COURANT_MAX = 2.0

#: The divergence correction is reported at WARNING level when the uniform
#: normal-velocity increment it adds exceeds this fraction of the boundary
#: velocity scale (the rho-weighted rms of the normal velocity over the fluid
#: lateral faces) on any time level.
CORRECTION_WARN_FRACTION = 0.05

#: Units of ``net_volume_flux`` and ``flux_residual``.  They are
#: :math:`\sum \rho u_n dA` with ``rhobf == 1`` always in uDALES (design
#: finding F1), i.e. a volume flux; the density weighting is kept for parity
#: with DALES's ``openboundary_divcorr`` but carries no dimension here.
FLUX_UNITS = "m3 s-1"

#: Optional global attributes this writer emits when it knows them; readers
#: that do not know them ignore them.
OPTIONAL_GLOBAL_ATTRIBUTES = ("child_dt", "parent_dy", "parent_dz")

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

#: Additional global attributes required by schema 2.
REQUIRED_GLOBAL_ATTRIBUTES_V2 = REQUIRED_GLOBAL_ATTRIBUTES + (
    "has_initial_condition",
    "fluid_lateral_area",
)

_COORD_VARIABLES = ("xf", "xh", "yf", "yh", "zf", "zh")

#: ``parent_dt`` must agree with the median spacing of the time axis to this
#: relative tolerance.  Loose enough for a parent whose dump times carry
#: round-off, tight enough to catch a cadence that is simply wrong.
PARENT_DT_RTOL = 1.0e-2

#: Tangential reconstructions :func:`conservative_interpolate` offers.
#:
#: ``"constant"`` is **divergence-preserving**: within a parent cell each of
#: du/dx, dv/dy and dw/dz equals the parent's, so a discretely solenoidal
#: parent gives an exactly solenoidal child target (design section 1.3).
#: ``"linear"`` reconstructs the tangential directions with unlimited central
#: slopes.  It removes the mean-profile staircase the constant scheme leaves --
#: on a log profile at r = 2 the interior mean error falls from 0.110 to 0.006
#: u* -- but it is **not** divergence-preserving: measured on the V0 coarse
#: arm's own 4 m LES field, the child target's divmax goes from the parent's
#: own 1.3e-7 to 1.7e-1, about 5 % of u/dx, which the child's projection then
#: has to remove inside the zone at every timestep.
#:
#: The default is ``"linear"``, **settled by measurement** in V0b (job 4000813,
#: design section 10.5): both schemes run at r = 2 and 4 with the final cadence
#: (0.5 s) and interpolant (Catmull-Rom), against the same fine truth.  The
#: linear scheme's lost solenoidality is real but turns out to be cheap, and
#: what it buys is not:
#:
#:   quantity (r = 2 / r = 4)     constant          linear
#:   pre-projection child divmax  9.5e-8 / 2.3e-8   0.363 / 0.287
#:   criterion A [u*]             0.128 / 0.280     0.098 / 0.121
#:   staircase RMS [u*]           0.044 / 0.055     0.005 / 0.009
#:   TKE deficit above z/h = 2    -6.5 / -16.9 %    -2.6 / -10.6 %
#:
#: The mean-flow error falls by up to 2.3x and the staircase by 6-9x, and
#: post-projection divmax and Phi are unchanged at round-off.  That is the
#: whole of the case for linear, and it is an empirical one for the cases
#: tested.
#:
#: **What is NOT established.**  An earlier version of this comment claimed the
#: projection absorbs the extra divergence "locally, at 1-3 % on the pressure
#: response".  That rested on two pressure samples per run -- V0b inherited the
#: default diagnostic interval, so only the compulsory first and last reports
#: exist, and the first is a startup transient outside the averaging window.
#: Two endpoints are not a time mean, and at r = 4 the interior norm rose more
#: than the zone norm, so locality is not shown either.  The claim is withdrawn
#: pending a run with a real diagnostic interval.
#:
#: ``"constant"`` remains available and is the only scheme satisfying design
#: section 1.3's exact identity; the identity tests name it explicitly.  Caveat:
#: at r = 4 linear overshoots the parent-resolved band and generates less
#: sub-filter energy, so it is not uniformly better.
PROLONGATIONS = ("constant", "linear")
DEFAULT_PROLONGATION = "linear"

#: Snapping tolerance for "this child face is coplanar with a parent face",
#: relative to the smallest parent spacing.  The same number is the alignment
#: tolerance of :func:`check_alignment`.
_ALIGN_RTOL = 1.0e-9


# --------------------------------------------------------------------------- #
# UDPrep section (the &NESTING namelist, docs/udales-nesting-spec.md section 4)
# --------------------------------------------------------------------------- #

#: Namelist defaults, from ``defaults.json`` like every other section.  The
#: section carries the &NESTING switches through ``UDPrep`` (``prep.nesting``)
#: so that ``save_param``/``write_changed_params`` see them; the file itself
#: is written by :func:`write_nesting_file` or :class:`NestingWriter`, not by
#: a section step, because its input -- a parent record -- is not part of a
#: case directory.
DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("nesting", {})
FIELDS = list(DEFAULTS.keys())

SPEC = SectionSpec(name="nesting", fields=FIELDS, defaults=DEFAULTS, section_cls=Section)


class NestingSchemaError(DataFormatError):
    """A nesting file does not match the schema of the implementation contract."""


class NestingRefinementError(ConfigurationError):
    """The parent/child refinement ratio exceeds the supported range."""


class NestingAlignmentError(ConfigurationError):
    """A child cell straddles a parent face, so the prolongation is not conservative."""


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


def _axis_misalignment(parent_faces: np.ndarray, child_faces: np.ndarray) -> Tuple[float, float]:
    """How far the child cells along one axis are from being nested in the parent's.

    Returns ``(worst, tol)`` in metres: ``worst`` is the largest distance from a
    parent face that lies inside the child's extent to the nearest child face,
    i.e. by how much a child cell straddles a parent face (0 when every parent
    face inside the child is also a child face); ``tol`` is the alignment
    tolerance, :data:`_ALIGN_RTOL` times the smallest parent spacing.
    """
    parent_faces = np.asarray(parent_faces, dtype=np.float64)
    child_faces = np.asarray(child_faces, dtype=np.float64)
    tol = _ALIGN_RTOL * float(np.min(np.diff(parent_faces)))
    inside = parent_faces[(parent_faces > child_faces[0] + tol)
                          & (parent_faces < child_faces[-1] - tol)]
    if inside.size == 0:
        return 0.0, tol
    pos = np.searchsorted(child_faces, inside)
    lo = child_faces[np.clip(pos - 1, 0, child_faces.size - 1)]
    hi = child_faces[np.clip(pos, 0, child_faces.size - 1)]
    return float(np.max(np.minimum(np.abs(inside - lo), np.abs(hi - inside)))), tol


def check_alignment(
    parent: NestGrid, child: NestGrid, allow_misaligned: bool = False
) -> Dict[str, float]:
    """Require every child cell to lie inside exactly one parent cell.

    The prolongation is conservative -- the child fluxes over a parent face sum
    to the parent flux, and the divergence integrated over a parent cell is
    preserved -- **only** when the grids nest: every parent face within the
    child's extent must coincide with a child face, along all three axes, to
    :data:`_ALIGN_RTOL` of the smallest parent spacing.  A refinement ratio of
    1.5, or a stretched parent vertical under a uniform child, breaks this
    silently and leaves the child target with a divergence of order half the
    velocity scale over a parent spacing (review 2026-09-06 §3 item 1).

    Returns the worst misalignment per axis, in metres, as
    ``{"x": ..., "y": ..., "z": ..., "tolerance": ...}``.  Raises
    :class:`NestingAlignmentError` when any exceeds the tolerance, unless
    ``allow_misaligned`` is set, in which case the violation is logged at
    WARNING level and the numbers are returned for the caller to record.
    """
    report: Dict[str, float] = {}
    worst_axis = None
    for name, axis in (("x", 0), ("y", 1), ("z", 2)):
        worst, tol = _axis_misalignment(parent.edges(axis), child.edges(axis))
        report[name] = worst
        report["tolerance"] = max(report.get("tolerance", 0.0), tol)
        if worst > tol and (worst_axis is None or worst > report[worst_axis]):
            worst_axis = name
    if worst_axis is not None:
        message = (
            "child cells are not nested in parent cells: a parent face lies "
            f"{report[worst_axis]:.3g} m inside a child cell along {worst_axis} "
            f"(tolerance {report['tolerance']:.3g} m; per axis x = {report['x']:.3g}, "
            f"y = {report['y']:.3g}, z = {report['z']:.3g} m). The prolongation is "
            "conservative only for nested grids: use an integer refinement ratio "
            "and a child vertical whose faces contain the parent's"
        )
        if not allow_misaligned:
            raise NestingAlignmentError(message + "; pass allow_misaligned=True to proceed anyway")
        logger.warning("udprep.nesting: %s (allow_misaligned=True, proceeding)", message)
    return report


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


def _along(axis: int, ndim: int, values: np.ndarray) -> np.ndarray:
    """Reshape a 1-D array so it broadcasts along ``axis`` of an ``ndim`` array."""
    shape = [1] * ndim
    shape[axis] = -1
    return np.asarray(values, dtype=np.float64).reshape(shape)


def _tangential_slopes(
    arr: np.ndarray, axis: int, centres: np.ndarray, mask: Optional[np.ndarray]
) -> np.ndarray:
    """Per-parent-cell slope along a cell-centred axis, for the linear reconstruction.

    Central differences of the neighbouring cell values over the distance
    between their centres; one-sided at the two ends.  With ``mask`` (True
    where the parent cell is fluid) a solid neighbour is not used: the slope
    is one-sided from the fluid neighbour, zero when both are solid, and zero
    in a solid cell itself -- so the near-zero in-building velocity of an IBM
    parent never leaks into the fluid cell beside it, and the reconstruction
    there degrades to first order rather than to a wrong slope.  No limiter:
    the result is linear in ``arr``.
    """
    n = arr.shape[axis]
    slope = np.zeros_like(arr)
    if n < 2:
        return slope
    c = _along(axis, arr.ndim, centres)
    fwd = np.diff(arr, axis=axis) / np.diff(c, axis=axis)      # between k and k+1
    interior = [slice(None)] * arr.ndim
    first = [slice(None)] * arr.ndim
    last = [slice(None)] * arr.ndim
    interior[axis] = slice(1, n - 1)
    first[axis] = slice(0, 1)
    last[axis] = slice(n - 1, n)
    head = [slice(None)] * arr.ndim
    tail = [slice(None)] * arr.ndim
    head[axis] = slice(0, n - 2)          # fwd[k-1] for interior k
    tail[axis] = slice(1, n - 1)          # fwd[k]   for interior k
    if n > 2:
        span = (np.take(c, np.arange(2, n), axis=axis) - np.take(c, np.arange(0, n - 2), axis=axis))
        central = (np.take(arr, np.arange(2, n), axis=axis)
                   - np.take(arr, np.arange(0, n - 2), axis=axis)) / span
        slope[tuple(interior)] = central
    slope[tuple(first)] = fwd[tuple(first)]
    # With exactly two parent cells there is a single forward difference
    # (fwd has length 1 on this axis); it is the correct slope at both
    # endpoints, so re-use fwd[first] rather than the out-of-range
    # fwd[last] (which would slice an empty range and fail to broadcast).
    slope[tuple(last)] = fwd[tuple(first)] if n == 2 else np.take(fwd, [n - 2], axis=axis)
    if mask is not None:
        fluid = np.asarray(mask, dtype=bool)
        left = np.zeros_like(fluid)
        right = np.zeros_like(fluid)
        left[tuple(interior)] = np.take(fluid, np.arange(0, n - 2), axis=axis)
        left[tuple(last)] = np.take(fluid, [n - 2], axis=axis)
        right[tuple(interior)] = np.take(fluid, np.arange(2, n), axis=axis)
        right[tuple(first)] = np.take(fluid, [1], axis=axis)
        bwd = np.zeros_like(arr)
        fwd_full = np.zeros_like(arr)
        bwd[tuple(interior)] = fwd[tuple(head)]
        bwd[tuple(last)] = np.take(fwd, [n - 2], axis=axis)
        fwd_full[tuple(interior)] = fwd[tuple(tail)]
        fwd_full[tuple(first)] = fwd[tuple(first)]
        slope = np.where(left & right, slope,
                         np.where(left, bwd, np.where(right, fwd_full, 0.0)))
        slope = np.where(fluid, slope, 0.0)
    return slope


def _apply_tangential(
    arr: np.ndarray,
    axis: int,
    idx: np.ndarray,
    target: np.ndarray,
    centres: np.ndarray,
    prolongation: str,
    mask: Optional[np.ndarray],
) -> np.ndarray:
    """Distribute parent cell values onto child points along a tangential axis."""
    out = np.take(arr, idx, axis=axis)
    if prolongation == "linear":
        slope = _tangential_slopes(arr, axis, centres, mask)
        offset = _along(axis, arr.ndim, target - centres[idx])
        out = out + np.take(slope, idx, axis=axis) * offset
    return out


def conservative_interpolate(
    parent: NestGrid,
    field: np.ndarray,
    component: str,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    prolongation: str = DEFAULT_PROLONGATION,
    parent_mask: Optional[np.ndarray] = None,
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
        coordinates in the other two.  The tangential targets must be the
        child cell centres (midpoints), which is what makes the linear
        reconstruction conservative.
    prolongation
        ``"linear"`` (default) or ``"constant"``: the tangential
        reconstruction, see the module docstring.
    parent_mask
        Optional fluid mask of the parent at ``field``'s stagger (True =
        fluid), shape of ``field``; steers the slopes of the linear
        reconstruction away from solid cells (:func:`_tangential_slopes`).
        Ignored by the constant scheme.

    Returns
    -------
    numpy.ndarray
        The interpolated field, shape ``(x.size, y.size, z.size)``.

    Notes
    -----
    Conservative in the two tangential directions (each parent face flux is
    shared over the child faces it contains, to round-off) and linear in the
    normal direction.  The tangential axes are applied first so the parent
    mask stays on the parent stagger while it is needed; the three 1-D
    operators act on separate axes, so the order does not change the result.
    """
    if component not in COMPONENTS:
        raise ConfigurationError(f"unknown velocity component {component!r}")
    if prolongation not in PROLONGATIONS:
        raise ConfigurationError(
            f"unknown prolongation {prolongation!r}; expected one of {PROLONGATIONS}"
        )
    field = np.asarray(field, dtype=np.float64)
    expected = parent.component_shape(component)
    if field.shape != expected:
        raise ConfigurationError(
            f"parent field for {component!r} has shape {field.shape}, expected {expected}"
        )
    mask = None
    if parent_mask is not None and prolongation == "linear":
        mask = np.asarray(parent_mask, dtype=bool)
        if mask.shape != expected:
            raise ConfigurationError(
                f"parent mask for {component!r} has shape {mask.shape}, expected {expected}"
            )
    targets = [np.asarray(t, dtype=np.float64) for t in (x, y, z)]
    out = field
    tangential = [ax for ax in range(3) if not _IS_FACE[component][ax]]
    normal = [ax for ax in range(3) if _IS_FACE[component][ax]]
    for axis in tangential:
        idx = _tangential_map(parent.edges(axis), targets[axis])
        out = _apply_tangential(out, axis, idx, targets[axis], parent.centres(axis),
                                prolongation, mask)
        if mask is not None:
            mask = np.take(mask, idx, axis=axis)
    for axis in normal:
        idx, wlo = _normal_map(parent.edges(axis), targets[axis])
        out = _apply_normal(out, axis, idx, wlo)
    return out


def stagger_masks_from_ibm(fluid: np.ndarray) -> Dict[str, np.ndarray]:
    """Fluid masks at the three velocity staggers from a cell-centred IBM mask.

    ``fluid`` has shape ``(itot, jtot, ktot)``, True where the cell is fluid.
    A face is fluid when both cells it separates are fluid, and a domain
    boundary face when the one cell behind it is -- the solver's ``IIu``,
    ``IIv``, ``IIw``.  Use the result as ``parent_masks`` of
    :func:`slabs_from_parent` and friends.
    """
    fluid = np.asarray(fluid, dtype=bool)
    if fluid.ndim != 3:
        raise ConfigurationError(
            f"the IBM mask must be a 3-D (itot, jtot, ktot) array, got shape {fluid.shape}"
        )
    out: Dict[str, np.ndarray] = {}
    for component, axis in zip(COMPONENTS, range(3)):
        pad = [(0, 0)] * 3
        pad[axis] = (1, 1)
        padded = np.pad(fluid, pad, mode="edge")
        lo = [slice(None)] * 3
        hi = [slice(None)] * 3
        lo[axis] = slice(0, fluid.shape[axis] + 1)
        hi[axis] = slice(1, fluid.shape[axis] + 2)
        out[component] = padded[tuple(lo)] & padded[tuple(hi)]
    return out


def interpolate_child_fields(
    parent: NestGrid,
    parent_u: np.ndarray,
    parent_v: np.ndarray,
    parent_w: np.ndarray,
    child: NestGrid,
    allow_misaligned: Optional[bool] = False,
    prolongation: str = DEFAULT_PROLONGATION,
    parent_masks: Optional[Mapping[str, np.ndarray]] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate a full parent velocity field onto the whole child grid.

    Convenience wrapper around :func:`conservative_interpolate`; the writer
    itself only ever interpolates the zone slabs (:func:`slabs_from_parent`).
    The grids must nest (:func:`check_alignment`); ``allow_misaligned=True``
    logs the violation and proceeds, ``None`` means the caller has already
    checked.  ``prolongation`` and ``parent_masks`` (per component, see
    :func:`stagger_masks_from_ibm`) are passed on.
    """
    if allow_misaligned is not None:
        check_alignment(parent, child, allow_misaligned=allow_misaligned)
    out = []
    for component, pf in zip(COMPONENTS, (parent_u, parent_v, parent_w)):
        coords = [child.component_coords(component, ax) for ax in range(3)]
        mask = None if parent_masks is None else parent_masks.get(component)
        out.append(conservative_interpolate(parent, pf, component, *coords,
                                            prolongation=prolongation, parent_mask=mask))
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


def init_dimensions(component: str) -> Tuple[str, str, str]:
    """CDL dimension names of the full-domain initial-condition variable.

    ``u_init(xh, yf, zf)``, ``v_init(xf, yh, zf)``, ``w_init(xf, yf, zh)`` --
    each component on its own stagger, in ``(x, y, z)`` CDL order, so Fortran
    sees ``(z, y, x)`` and a rank's ``(i, j)`` block is contiguous in ``z``.
    """
    if component not in COMPONENTS:
        raise ConfigurationError(f"unknown velocity component {component!r}")
    return (
        "xh" if component == "u" else "xf",
        "yh" if component == "v" else "yf",
        "zh" if component == "w" else "zf",
    )


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
    allow_misaligned: Optional[bool] = False,
    prolongation: str = DEFAULT_PROLONGATION,
    parent_masks: Optional[Mapping[str, np.ndarray]] = None,
) -> Dict[str, np.ndarray]:
    """Interpolate the parent onto the twelve zone slabs, for one time level.

    Only the slab points are interpolated, so the cost and the memory are
    proportional to the zone rather than to the child domain.  The grids must
    nest (:func:`check_alignment`); ``allow_misaligned=True`` logs the
    violation and proceeds, ``None`` means the caller has already checked.
    ``prolongation`` and ``parent_masks`` (per component, see
    :func:`stagger_masks_from_ibm`) are passed on to
    :func:`conservative_interpolate`.
    """
    _check_nzone(child, nzone)
    if allow_misaligned is not None:
        check_alignment(parent, child, allow_misaligned=allow_misaligned)
    parent_fields = dict(zip(COMPONENTS, (parent_u, parent_v, parent_w)))
    slabs: Dict[str, np.ndarray] = {}
    for face in FACES:
        for component in COMPONENTS:
            coords = slab_coordinates(child, nzone, face, component)
            mask = None if parent_masks is None else parent_masks.get(component)
            block = conservative_interpolate(
                parent, parent_fields[component], component, *coords,
                prolongation=prolongation, parent_mask=mask,
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
# Full-domain initial condition (schema 2)
# --------------------------------------------------------------------------- #


def initial_fields_from_fields(
    grid: NestGrid, u: np.ndarray, v: np.ndarray, w: np.ndarray
) -> Dict[str, np.ndarray]:
    """Validate ``(u, v, w)`` on the child grid as an initial-condition block."""
    out: Dict[str, np.ndarray] = {}
    for component, arr in zip(COMPONENTS, (u, v, w)):
        arr = np.ascontiguousarray(np.asarray(arr, dtype=np.float64))
        expected = grid.component_shape(component)
        if arr.shape != expected:
            raise ConfigurationError(
                f"initial field {component!r} has shape {arr.shape}, expected {expected}"
            )
        out[component] = arr
    return out


def initial_fields_from_parent(
    parent: NestGrid,
    parent_u: np.ndarray,
    parent_v: np.ndarray,
    parent_w: np.ndarray,
    child: NestGrid,
    allow_misaligned: Optional[bool] = False,
    prolongation: str = DEFAULT_PROLONGATION,
    parent_masks: Optional[Mapping[str, np.ndarray]] = None,
) -> Dict[str, np.ndarray]:
    """Conservatively interpolate a parent field onto the whole child grid."""
    return initial_fields_from_fields(
        child, *interpolate_child_fields(parent, parent_u, parent_v, parent_w, child,
                                         allow_misaligned=allow_misaligned,
                                         prolongation=prolongation, parent_masks=parent_masks)
    )


def discrete_divergence(
    grid: NestGrid, u: np.ndarray, v: np.ndarray, w: np.ndarray
) -> np.ndarray:
    """Cell-centred discrete divergence, the operator ``fillps`` applies.

    ``(u[i+1] - u[i])/dx + (v[j+1] - v[j])/dy + (w[k+1] - w[k])/dzf`` on the
    child grid (``modpois.f90:966-973``), shape ``(itot, jtot, ktot)``.
    """
    u = np.asarray(u, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    return (np.diff(u, axis=0) / grid.dx[:, None, None]
            + np.diff(v, axis=1) / grid.dy[None, :, None]
            + np.diff(w, axis=2) / grid.dzf[None, None, :])


def _dzh(grid: NestGrid) -> np.ndarray:
    """Spacing between successive cell centres, ``dzh[k] = zf[k] - zf[k-1]``."""
    out = np.empty(grid.ktot + 1, dtype=np.float64)
    out[0] = 2.0 * (grid.zf[0] - grid.zh[0])
    out[1:grid.ktot] = np.diff(grid.zf)
    out[grid.ktot] = 2.0 * (grid.zh[-1] - grid.zf[-1])
    return out


def _neumann_poisson(grid: NestGrid, rhs: np.ndarray) -> np.ndarray:
    """Solve ``D G p = rhs`` with homogeneous Neumann on all six faces.

    Diagonalised by a type-II DCT in ``x`` and ``y`` (the cosine transform
    ``modpois`` uses for non-periodic laterals, ``FFTW_REDFT10``) and solved by
    a tridiagonal sweep in the stretched vertical.  The constant mode is
    singular; it is fixed by ``p = 0`` in the bottom layer, which is legitimate
    exactly when ``rhs`` is compatible -- i.e. when the net boundary flux
    vanishes, which is what :func:`apply_divergence_correction` enforces.
    """
    try:
        from scipy.fft import dctn, idctn
    except ImportError as exc:  # pragma: no cover - scipy is a core dependency
        raise DependencyError(
            "scipy is required to project the nesting initial condition"
        ) from exc

    itot, jtot, ktot = rhs.shape
    dx = float(grid.dx[0])
    dy = float(grid.dy[0])
    if not (np.allclose(grid.dx, dx, rtol=1e-12, atol=0.0)
            and np.allclose(grid.dy, dy, rtol=1e-12, atol=0.0)):
        raise ConfigurationError("the initial-condition projection needs a uniform x-y grid")

    rhat = dctn(rhs, type=2, norm="ortho", axes=(0, 1))

    kx = np.arange(itot)
    ky = np.arange(jtot)
    lam = ((-4.0 / dx**2) * np.sin(np.pi * kx / (2 * itot))**2)[:, None] \
        + ((-4.0 / dy**2) * np.sin(np.pi * ky / (2 * jtot))**2)[None, :]

    dzf = grid.dzf
    dzh = _dzh(grid)
    a = np.zeros(ktot)
    c = np.zeros(ktot)
    a[1:] = 1.0 / (dzf[1:] * dzh[1:ktot])
    c[:-1] = 1.0 / (dzf[:-1] * dzh[1:ktot])
    b = -(a + c)

    phat = np.empty_like(rhat)

    # Singular (constant) mode: integrate the flux form directly.
    flux = np.concatenate(([0.0], np.cumsum(dzf * rhat[0, 0, :])))
    pcol = np.empty(ktot)
    pcol[0] = 0.0
    for k in range(1, ktot):
        pcol[k] = pcol[k - 1] + flux[k] * dzh[k]
    phat[0, 0, :] = pcol

    # Every other mode: vectorised Thomas sweep over the horizontal modes.
    sel = np.ones((itot, jtot), dtype=bool)
    sel[0, 0] = False
    lam_s = lam[sel]
    rhs_s = rhat[sel]
    n = lam_s.size
    cp = np.empty((n, ktot))
    dp = np.empty((n, ktot))
    beta = b[0] + lam_s
    cp[:, 0] = c[0] / beta
    dp[:, 0] = rhs_s[:, 0] / beta
    for k in range(1, ktot):
        beta = (b[k] + lam_s) - a[k] * cp[:, k - 1]
        cp[:, k] = c[k] / beta
        dp[:, k] = (rhs_s[:, k] - a[k] * dp[:, k - 1]) / beta
    sol = np.empty((n, ktot))
    sol[:, -1] = dp[:, -1]
    for k in range(ktot - 2, -1, -1):
        sol[:, k] = dp[:, k] - cp[:, k] * sol[:, k + 1]
    phat[sel] = sol

    return idctn(phat, type=2, norm="ortho", axes=(0, 1))


def project_initial_condition(
    grid: NestGrid, fields: Mapping[str, np.ndarray], rtol: float = 1.0e-9
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    """Project ``(u, v, w)`` onto the discretely solenoidal subspace, in place.

    The boundary-normal velocities on all six faces are **not** touched: the
    pressure carries homogeneous Neumann conditions, exactly as ``bcp`` applies
    them in the solver (design finding F2), so the imposed lateral data and the
    lid survive the projection.  That also means the projection is only solvable
    when the net boundary flux already vanishes; a residual larger than
    ``rtol`` times the boundary flux scale is an error, not something to absorb.

    Returns ``(u, v, w, div_before, div_after)``; the arrays are the same
    objects that were passed in.
    """
    u = fields["u"]
    v = fields["v"]
    w = fields["w"]

    div = discrete_divergence(grid, u, v, w)
    before = float(np.max(np.abs(div)))

    # Compatibility: the net boundary flux telescopes out of the divergence sum
    # and must vanish, or the pure-Neumann system has no solution.  Judge it
    # against the flux a typical velocity would carry through the whole boundary,
    # which is the only scale that stays meaningful when the field is already
    # solenoidal and both sides are at round-off.
    vol = (grid.dx[:, None, None] * grid.dy[None, :, None] * grid.dzf[None, None, :])
    net = float(np.sum(vol * div))
    area = 2.0 * (grid.xlen + grid.ylen) * grid.zsize + 2.0 * grid.xlen * grid.ylen
    speed = max(float(np.max(np.abs(u))), float(np.max(np.abs(v))),
                float(np.max(np.abs(w))), np.finfo(np.float64).tiny)
    if abs(net) > rtol * speed * area:
        raise ConfigurationError(
            "the initial condition is not compatible with Neumann pressure: its net "
            f"boundary flux is {net:g} {FLUX_UNITS}, against {rtol * speed * area:g} allowed. "
            "Run apply_divergence_correction (or sync_initial_condition) first."
        )

    p = _neumann_poisson(grid, div)

    # Interior faces only, so every boundary-normal velocity is left alone.  The
    # horizontal spacings are the cell widths because _neumann_poisson has
    # already required a uniform x-y grid, where they equal the centre spacings;
    # the vertical uses the centre spacing explicitly, since z may be stretched.
    u[1:-1, :, :] -= np.diff(p, axis=0) / grid.dx[1:, None, None]
    v[:, 1:-1, :] -= np.diff(p, axis=1) / grid.dy[None, 1:, None]
    w[:, :, 1:-1] -= np.diff(p, axis=2) / _dzh(grid)[None, None, 1:grid.ktot]

    after = float(np.max(np.abs(discrete_divergence(grid, u, v, w))))
    return u, v, w, before, after


# --------------------------------------------------------------------------- #
# Time axis
# --------------------------------------------------------------------------- #


def check_time_axis(times: np.ndarray, parent_dt: float, context: str = "") -> float:
    """Validate the stored time axis and ``parent_dt`` against each other.

    The contract (docs/udales-nesting-spec.md section 5) is that ``time`` is
    the **child's** clock: it starts at exactly 0, since the solver positions
    its parent-time buffer on ``timee`` and a record starting later would
    leave the child frozen on the first level with no message; it is strictly
    increasing; and ``parent_dt`` -- the cadence the temporal-refinement guard
    and the cadence criterion use -- is positive and agrees with the median
    spacing to :data:`PARENT_DT_RTOL`.  A producer working from absolute
    parent times subtracts ``times[0]`` first, as ``make_child_case`` does.

    A single time level has no spacing, so only ``parent_dt >= 0`` is asked
    of it.  Returns the median spacing (``parent_dt`` when there is only one
    level).  Raises :class:`ConfigurationError` naming the violation.
    """
    times = np.asarray(times, dtype=np.float64).reshape(-1)
    where = f"{context}: " if context else ""
    if times.size == 0:
        raise ConfigurationError(f"{where}the time axis is empty")
    if not np.all(np.isfinite(times)):
        raise ConfigurationError(f"{where}the time axis contains non-finite values")
    if times[0] != 0.0:
        raise ConfigurationError(
            f"{where}time[0] = {times[0]:g} s, but the stored axis is the child's own "
            "clock and must start at exactly 0; subtract the first parent time "
            "(the time-origin convention make_child_case applies)"
        )
    if times.size > 1:
        spacing = np.diff(times)
        if not np.all(spacing > 0.0):
            bad = int(np.argmin(spacing > 0.0))
            raise ConfigurationError(
                f"{where}the time axis is not strictly increasing: time[{bad}] = "
                f"{times[bad]:g} s is followed by time[{bad + 1}] = {times[bad + 1]:g} s"
            )
        median = float(np.median(spacing))
        if not parent_dt > 0.0:
            raise ConfigurationError(
                f"{where}parent_dt = {parent_dt:g} s must be positive; the time axis "
                f"has a median spacing of {median:g} s"
            )
        if abs(parent_dt - median) > PARENT_DT_RTOL * median:
            raise ConfigurationError(
                f"{where}parent_dt = {parent_dt:g} s disagrees with the median spacing "
                f"of the time axis, {median:g} s (tolerance {PARENT_DT_RTOL:g} relative)"
            )
        return median
    if parent_dt < 0.0:
        raise ConfigurationError(f"{where}parent_dt = {parent_dt:g} s must not be negative")
    return float(parent_dt)


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


def face_masks_from_ibm(fluid: np.ndarray, solid: bool = False) -> FaceMasks:
    """Derive the four lateral :class:`FaceMasks` from the child's IBM cell mask.

    ``fluid`` is the child's cell-centred mask, shape ``(itot, jtot, ktot)``,
    ``True`` where the cell is fluid (pass ``solid=True`` for the opposite
    convention, e.g. a mask built straight from ``solid_c.txt``).  A boundary
    face is fluid exactly when the boundary cell behind it is: that is what
    the solver's ``IIu``/``IIv`` say at ``ib``/``ie+1`` and ``jb``/``je+1``,
    and what ``fluid_lateral_area`` on the solver side sums over, so the
    writer's residual and the solver's agree.
    """
    fluid = np.asarray(fluid, dtype=bool)
    if fluid.ndim != 3:
        raise ConfigurationError(
            f"the IBM mask must be a 3-D (itot, jtot, ktot) array, got shape {fluid.shape}"
        )
    if solid:
        fluid = ~fluid
    return FaceMasks(
        west=np.ascontiguousarray(fluid[0, :, :]),
        east=np.ascontiguousarray(fluid[-1, :, :]),
        south=np.ascontiguousarray(fluid[:, 0, :]),
        north=np.ascontiguousarray(fluid[:, -1, :]),
    )


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
    #: Post-correction residual of :func:`net_volume_flux`, one value per time
    #: level (schema 2).  ``None`` means "recompute it at write time".
    flux_residual: Optional[np.ndarray] = None
    #: Geometric fluid area of the four lateral boundary faces the residual was
    #: summed over.  The solver compares it against its own before trusting
    #: ``flux_residual`` (docs/udales-nesting-spec.md section 5).
    fluid_lateral_area: Optional[float] = None
    #: OPTIONAL full-domain initial condition at ``times[0]``, keyed ``u``/
    #: ``v``/``w`` at each component's own stagger (schema 2).
    initial_fields: Optional[Dict[str, np.ndarray]] = None
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
    #: Parent spacings in ``y`` and ``z`` (smallest), when known; they enter
    #: the spatial refinement ratio next to ``parent_dx`` and are written as
    #: optional attributes.  ``None`` means unknown.
    parent_dy: Optional[float] = None
    parent_dz: Optional[float] = None
    #: The fluid/solid masks the correction used and the residual was summed
    #: over.  Not stored in the file (the solver derives its own from the IBM);
    #: kept here so that the writer can recompute ``flux_residual`` and
    #: ``fluid_lateral_area`` and verify the cached values rather than trust
    #: them.  ``None`` means all fluid.
    masks: Optional[FaceMasks] = None
    #: What :func:`apply_divergence_correction` did, per :func:`correction_report`;
    #: ``None`` until it has run.  Not stored in the file.
    correction: Optional[Dict[str, Any]] = None

    def __post_init__(self) -> None:
        self.nzone = int(self.nzone)
        _check_nzone(self.grid, self.nzone)
        self.times = np.asarray(self.times, dtype=np.float64).reshape(-1)
        self.parent_dt = float(self.parent_dt)
        check_time_axis(self.times, self.parent_dt, context="NestingData")
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
        if self.flux_residual is not None:
            self.flux_residual = np.asarray(self.flux_residual, dtype=np.float64).reshape(-1)
            if self.flux_residual.size != self.ntime:
                raise ConfigurationError("flux_residual must have one value per time level")
        if self.initial_fields is not None:
            self.initial_fields = initial_fields_from_fields(
                self.grid,
                self.initial_fields["u"],
                self.initial_fields["v"],
                self.initial_fields["w"],
            )
            # The initial condition is projected with the solver's OWN divergence
            # operator, which carries no density (design finding F1: `fillps`
            # differs from the Laplacian in exactly this respect).  The slab flux
            # correction, in contrast, is rho-weighted, to match DALES's
            # `openboundary_divcorr`.  The two agree only while rho == 1 -- which
            # in uDALES is always, since `rhobf`/`rhobh` are set to 1 and never
            # assigned anywhere else.  Refuse rather than silently produce an
            # initial condition whose boundary flux does not close.
            if not (np.allclose(self.rhobf, 1.0, rtol=0.0, atol=1e-12)
                    and np.allclose(self.rhobh, 1.0, rtol=0.0, atol=1e-12)):
                raise ConfigurationError(
                    "a full-domain initial condition needs rhobf == rhobh == 1: it is "
                    "projected with the solver's density-free divergence operator "
                    "(design finding F1), which cannot be reconciled with the "
                    "rho-weighted lateral flux correction"
                )
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
            flux_residual=None if self.flux_residual is None else self.flux_residual.copy(),
            fluid_lateral_area=self.fluid_lateral_area,
            initial_fields=(None if self.initial_fields is None
                            else {k: v.copy() for k, v in self.initial_fields.items()}),
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
            parent_dy=self.parent_dy,
            parent_dz=self.parent_dz,
            masks=self.masks,
            correction=None if self.correction is None else dict(self.correction),
        )


def nesting_data_from_parent(
    parent: NestGrid,
    child: NestGrid,
    nzone: int,
    times: Sequence[float],
    fields: Sequence[Tuple[np.ndarray, np.ndarray, np.ndarray]],
    initial: bool = False,
    allow_misaligned: bool = False,
    prolongation: str = DEFAULT_PROLONGATION,
    parent_masks: Optional[Mapping[str, np.ndarray]] = None,
    **kwargs: Any,
) -> NestingData:
    """Build a :class:`NestingData` by conservative interpolation of a parent.

    ``fields`` is a sequence, one entry per time level, of ``(u, v, w)`` parent
    arrays at their own staggers.  With ``initial=True`` the first time level is
    additionally interpolated onto the **whole** child grid and carried as the
    schema-2 initial-condition block, for a cold start with
    ``nest_linitfromparent``.  Extra keyword arguments are passed to
    :class:`NestingData` (provenance attributes, ``rhobf``, ...).

    The grids must nest (:func:`check_alignment`, checked once here);
    ``allow_misaligned=True`` logs the violation and proceeds.

    This holds every level in memory; for a long parent record use
    :class:`NestingWriter` and append level by level instead.
    """
    times = np.asarray(times, dtype=np.float64).reshape(-1)
    if len(fields) != times.size:
        raise ConfigurationError(
            f"got {len(fields)} field sets for {times.size} times"
        )
    check_alignment(parent, child, allow_misaligned=allow_misaligned)
    per_time = [slabs_from_parent(parent, *fields[n], child=child, nzone=nzone,
                                  allow_misaligned=None, prolongation=prolongation,
                                  parent_masks=parent_masks)
                for n in range(times.size)]
    slabs = {name: np.stack([s[name] for s in per_time], axis=0) for name in SLAB_VARIABLES}
    if initial:
        kwargs.setdefault(
            "initial_fields",
            initial_fields_from_parent(parent, *fields[0], child=child, allow_misaligned=None,
                                       prolongation=prolongation, parent_masks=parent_masks),
        )
    kwargs.setdefault("parent_dx", float(np.min(np.diff(parent.xh))))
    kwargs.setdefault("parent_dy", float(np.min(np.diff(parent.yh))))
    kwargs.setdefault("parent_dz", float(np.min(np.diff(parent.zh))))
    kwargs.setdefault("child_origin_x", float(child.xh[0]))
    kwargs.setdefault("child_origin_y", float(child.yh[0]))
    if times.size > 1:
        kwargs.setdefault("parent_dt", float(np.median(np.diff(times))))
    return NestingData(grid=child, nzone=nzone, times=times, slabs=slabs, **kwargs)


# --------------------------------------------------------------------------- #
# Divergence correction (design §3.1, §3.3)
# --------------------------------------------------------------------------- #


@dataclass
class _LevelView:
    """One time level of slab data, with a time axis of length 1.

    The flux functions below only ever touch ``grid``, ``nzone``, ``rhobf``
    and ``slabs`` and are vectorised over a leading time axis, so a single
    level wrapped this way goes through exactly the code a whole
    :class:`NestingData` does -- which is how :class:`NestingWriter` corrects
    level by level without a second implementation.
    """

    grid: NestGrid
    nzone: int
    rhobf: np.ndarray
    slabs: Dict[str, np.ndarray]

    @property
    def ntime(self) -> int:
        return 1

    @classmethod
    def wrap(cls, grid: NestGrid, nzone: int, rhobf: np.ndarray,
             slabs: Mapping[str, np.ndarray]) -> "_LevelView":
        """Wrap 3-D per-level slabs as views with a leading time axis."""
        return cls(grid=grid, nzone=nzone, rhobf=rhobf,
                   slabs={name: arr[None] for name, arr in slabs.items()})


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


def fluid_lateral_area(data: NestingData, masks: Optional[FaceMasks] = None) -> float:
    """Total **geometric** fluid area of the four lateral boundary faces [m2].

    Unlike :func:`fluid_face_area` this carries no density, because it is
    written to the file for the solver to compare against its own IIu/IIv
    boundary area before it trusts the stored ``flux_residual``.  The two agree
    whenever ``rhobf == 1``, which in uDALES is always (design finding F1).
    """
    masks = masks or _ALL_FLUID
    grid = data.grid
    total = 0.0
    for face in FACES:
        span = grid.dy if _FACE_AXIS[face] == 0 else grid.dx
        area = span[:, None] * grid.dzf[None, :]
        total += float(np.sum(area * masks.get(grid, face)))
    return total


def _boundary_velocity_scale(data: Any, masks: FaceMasks) -> np.ndarray:
    """rho-weighted rms of the normal velocity over the fluid lateral faces, per level."""
    faces = boundary_faces(data)
    num = np.zeros(data.ntime, dtype=np.float64)
    den = 0.0
    for face in FACES:
        weight = _face_weights(data, face) * masks.get(data.grid, face)
        num += np.tensordot(faces[face] ** 2, weight, axes=([1, 2], [0, 1]))
        den += float(np.sum(weight))
    return np.sqrt(num / den) if den > 0.0 else num


def _face_mean_outward_velocity(data: Any, masks: FaceMasks) -> Dict[str, np.ndarray]:
    """Area-weighted mean **outward** normal velocity on each face, per level.

    Negative means the face is, on the mean, an inflow face at that level.
    A face with no fluid area gets 0.
    """
    faces = boundary_faces(data)
    out: Dict[str, np.ndarray] = {}
    for face in FACES:
        weight = _face_weights(data, face) * masks.get(data.grid, face)
        total = float(np.sum(weight))
        if total > 0.0:
            out[face] = _FACE_SIGN[face] * np.tensordot(
                faces[face], weight, axes=([1, 2], [0, 1])) / total
        else:
            out[face] = np.zeros(data.ntime, dtype=np.float64)
    return out


class _CorrectionTally:
    """Accumulates, over time levels, what the divergence correction did.

    Fed level by level (or all at once) with the pre-correction residual, the
    increment ``delta``, the boundary velocity scale and each face's mean
    outward normal velocity; :meth:`report` returns the dictionary that
    :func:`correction_report` documents.
    """

    def __init__(self, data: Any, masks: FaceMasks) -> None:
        total = fluid_face_area(data, masks)
        self.area_fraction = {
            face: float(np.sum(_face_weights(data, face) * masks.get(data.grid, face))) / total
            for face in FACES
        }
        self.ntime = 0
        self.residual_max = 0.0
        self.delta_max = 0.0
        self.fraction_max = 0.0
        self.fraction_level = -1
        self.scale_at_max = 0.0
        self.delta_at_max = 0.0
        self.outward_sum = {face: 0.0 for face in FACES}

    def add(self, residual: np.ndarray, delta: np.ndarray, scale: np.ndarray,
            outward: Mapping[str, np.ndarray], first_level: int) -> None:
        residual = np.atleast_1d(residual)
        delta = np.atleast_1d(delta)
        scale = np.atleast_1d(scale)
        fraction = np.abs(delta) / np.where(scale > 0.0, scale, np.inf)
        n = int(np.argmax(fraction))
        if fraction[n] > self.fraction_max or self.fraction_level < 0:
            self.fraction_max = float(fraction[n])
            self.fraction_level = first_level + n
            self.scale_at_max = float(scale[n])
            self.delta_at_max = float(delta[n])
        self.residual_max = max(self.residual_max, float(np.max(np.abs(residual))))
        self.delta_max = max(self.delta_max, float(np.max(np.abs(delta))))
        for face in FACES:
            self.outward_sum[face] += float(np.sum(outward[face]))
        self.ntime += residual.size

    def report(self) -> Dict[str, Any]:
        n = max(self.ntime, 1)
        faces = {}
        for face in FACES:
            mean = self.outward_sum[face] / n
            faces[face] = {
                "area_fraction": self.area_fraction[face],
                "mean_outward_normal_velocity": mean,
                "inflow": bool(mean < 0.0),
            }
        return {
            "ntime": self.ntime,
            "residual_max_abs": self.residual_max,
            "delta_max_abs": self.delta_max,
            "delta_fraction_max": self.fraction_max,
            "delta_fraction_level": self.fraction_level,
            "velocity_scale_at_max": self.scale_at_max,
            "delta_at_max": self.delta_at_max,
            "warn_fraction": CORRECTION_WARN_FRACTION,
            "exceeded": bool(self.fraction_max > CORRECTION_WARN_FRACTION),
            "faces": faces,
            "inflow_faces": [face for face in FACES if faces[face]["inflow"]],
        }

    def log(self) -> Dict[str, Any]:
        rep = self.report()
        split = ", ".join(
            f"{face} {100.0 * rep['faces'][face]['area_fraction']:.0f} %"
            f"{' (inflow)' if rep['faces'][face]['inflow'] else ''}"
            for face in FACES
        )
        summary = (
            f"divergence correction over {rep['ntime']} level(s): max |delta| = "
            f"{rep['delta_max_abs']:.3g} m/s; largest relative to the boundary velocity "
            f"scale {100.0 * rep['delta_fraction_max']:.1f} % (level "
            f"{rep['delta_fraction_level']}, delta = {rep['delta_at_max']:.3g} m/s against "
            f"{rep['velocity_scale_at_max']:.3g} m/s); the net flux (max |Phi| = "
            f"{rep['residual_max_abs']:.3g} {FLUX_UNITS}) is spread over the fluid lateral "
            f"faces by area: {split}"
        )
        if rep["exceeded"]:
            inflow = ", ".join(rep["inflow_faces"]) or "none"
            logger.warning(
                "udprep.nesting: %s -- this exceeds %.0f %%: the lid flux the parent "
                "carried is being pushed through the lateral faces, inflow faces (%s) "
                "included; check the parent's top boundary and the child's lid",
                summary, 100.0 * CORRECTION_WARN_FRACTION, inflow,
            )
        else:
            logger.info("udprep.nesting: %s", summary)
        return rep


def _max_normal_speed(data: Any, masks: FaceMasks) -> np.ndarray:
    """Largest |u_n| over the fluid lateral boundary faces, per level."""
    faces = boundary_faces(data)
    out = np.zeros(data.ntime, dtype=np.float64)
    for face in FACES:
        mask = masks.get(data.grid, face)
        out = np.maximum(out, np.max(np.abs(faces[face]) * mask[None, :, :], axis=(1, 2)))
    return out


class _CadenceTally:
    """Running ``C_dump = max|u_n| parent_dt / parent_dx`` over the stored levels."""

    def __init__(self, parent_dt: float, parent_dx: float) -> None:
        self.parent_dt = float(parent_dt)
        self.parent_dx = float(parent_dx)
        self.speed = 0.0
        self.level = -1
        self.ntime = 0

    def add(self, speed: np.ndarray, first_level: int) -> None:
        speed = np.atleast_1d(speed)
        n = int(np.argmax(speed))
        if speed[n] > self.speed or self.level < 0:
            self.speed = float(speed[n])
            self.level = first_level + n
        self.ntime += speed.size

    def report(self) -> Dict[str, Any]:
        known = self.parent_dt > 0.0 and self.parent_dx > 0.0
        c_dump = self.speed * self.parent_dt / self.parent_dx if known else None
        return {
            "C_dump": c_dump,
            "max_normal_speed": self.speed,
            "level": self.level,
            "parent_dt": self.parent_dt,
            "parent_dx": self.parent_dx,
            "limit": CADENCE_COURANT_MAX,
            "exceeded": bool(c_dump is not None and c_dump > CADENCE_COURANT_MAX),
        }

    def log(self) -> Dict[str, Any]:
        rep = self.report()
        if rep["C_dump"] is None:
            logger.info("udprep.nesting: cadence criterion not evaluated (parent_dt = %g s, "
                        "parent_dx = %g m)", self.parent_dt, self.parent_dx)
        elif rep["exceeded"]:
            logger.warning(
                "udprep.nesting: boundary cadence C_dump = max|u_n| dt_P / dx_P = %.3g x %g / %g "
                "= %.2f exceeds %g (level %d): between stored levels a feature crosses more "
                "than %g parent cells and the time interpolation of the boundary cannot "
                "follow it; store the parent more often or coarsen it",
                rep["max_normal_speed"], self.parent_dt, self.parent_dx, rep["C_dump"],
                CADENCE_COURANT_MAX, rep["level"], CADENCE_COURANT_MAX,
            )
        else:
            logger.info("udprep.nesting: boundary cadence C_dump = %.2f (max|u_n| = %.3g m/s, "
                        "dt_P = %g s, dx_P = %g m), within %g",
                        rep["C_dump"], rep["max_normal_speed"], self.parent_dt,
                        self.parent_dx, CADENCE_COURANT_MAX)
        return rep


def cadence_courant(data: NestingData, masks: Optional[FaceMasks] = None,
                    log: bool = True) -> Dict[str, Any]:
    """The boundary-cadence criterion ``C_dump = max|u_n| parent_dt / parent_dx``.

    ``max|u_n|`` is over the fluid lateral boundary faces and every stored
    level.  Returns ``{"C_dump", "max_normal_speed", "level", "parent_dt",
    "parent_dx", "limit", "exceeded"}``; ``C_dump`` is ``None`` when
    ``parent_dt`` or ``parent_dx`` is unknown.  With ``log`` the verdict is
    logged, at WARNING level above :data:`CADENCE_COURANT_MAX`.
    """
    masks = masks if masks is not None else (data.masks or _ALL_FLUID)
    tally = _CadenceTally(data.parent_dt, data.parent_dx)
    tally.add(_max_normal_speed(data, masks), 0)
    return tally.log() if log else tally.report()


def nesting_diagnostics(data: NestingData, masks: Optional[FaceMasks] = None) -> Dict[str, Any]:
    """Everything a manifest should record about a nesting file's data.

    ``{"cadence": cadence_courant(...), "correction": correction_report(...),
    "refinement": refinement_verdict(...)}`` -- computed without logging and
    without changing ``data``.
    """
    return {
        "cadence": cadence_courant(data, masks, log=False),
        "correction": correction_report(data, masks),
        "refinement": refinement_verdict(data),
    }


def correction_report(data: NestingData, masks: Optional[FaceMasks] = None) -> Dict[str, Any]:
    """What the divergence correction of ``data`` did, or would do.

    Keys: ``residual_max_abs`` (largest pre-correction |Phi|, m3 s-1),
    ``delta_max_abs`` (largest uniform increment, m/s), ``delta_fraction_max``
    (largest |delta| relative to the boundary velocity scale -- the
    rho-weighted rms normal velocity over the fluid faces -- with the level it
    occurs on, ``delta_fraction_level``), ``exceeded`` (above
    :data:`CORRECTION_WARN_FRACTION`), and per face under ``faces`` the
    ``area_fraction`` of the correction flux it carries and the time-mean
    outward normal velocity with an ``inflow`` flag, so a correction that
    pushes the lid flux through an inflow face is on record; ``inflow_faces``
    lists them.

    On corrected data (``data.correction`` set) this is what was applied;
    otherwise it is computed from the current slabs without changing them.
    """
    if data.correction is not None:
        return dict(data.correction)
    masks = masks if masks is not None else (data.masks or _ALL_FLUID)
    area = fluid_face_area(data, masks)
    residual = net_volume_flux(data, masks)
    tally = _CorrectionTally(data, masks)
    tally.add(residual, -residual / area, _boundary_velocity_scale(data, masks),
              _face_mean_outward_velocity(data, masks), 0)
    return tally.report()


def _correct_in_place(data: Any, masks: FaceMasks, area: float) -> Tuple[np.ndarray, np.ndarray]:
    """Add the uniform outward increment that zeroes Phi on every level of ``data``.

    Returns ``(residual, delta)``, both per level.  Works on a
    :class:`NestingData` and on a :class:`_LevelView` alike.
    """
    residual = net_volume_flux(data, masks)
    delta = -residual / area
    faces = boundary_faces(data)
    for face in FACES:
        mask = masks.get(data.grid, face)
        faces[face] += _FACE_SIGN[face] * delta[:, None, None] * mask[None, :, :]
    return residual, delta


def sync_initial_condition(
    data: "NestingData", masks: Optional[FaceMasks] = None, close_lid: bool = True
) -> Tuple[float, float]:
    """Make the initial condition consistent with the boundary data, then project.

    Three steps, in this order:

    1. the boundary-normal velocity on the four lateral faces is taken from the
       (corrected) slabs at ``times[0]``, so the first substep does not see a
       jump between the stored initial condition and the imposed boundary;
    2. ``w`` on the floor and the lid is set to zero -- the rigid-lid case A of
       design section 3.2, which is what v1 supports and what ``boundary``
       imposes on the child regardless of what the parent did;
    3. the field is projected (:func:`project_initial_condition`), which leaves
       every boundary-normal velocity untouched.

    Returns the peak discrete divergence before and after the projection.
    """
    if data.initial_fields is None:
        raise ConfigurationError("this NestingData carries no initial condition")
    masks = masks or _ALL_FLUID
    fields = data.initial_fields
    faces = boundary_faces(data)
    for face in FACES:
        component = _FACE_NORMAL_COMPONENT[face]
        index = -1 if _FACE_UPPER[face] else 0
        if _FACE_AXIS[face] == 0:
            fields[component][index, :, :] = faces[face][0]
        else:
            fields[component][:, index, :] = faces[face][0]
    if close_lid:
        fields["w"][:, :, 0] = 0.0
        fields["w"][:, :, -1] = 0.0
    _, _, _, before, after = project_initial_condition(data.grid, fields)
    return before, after


def apply_divergence_correction(
    data: NestingData, masks: Optional[FaceMasks] = None, project_initial: bool = True
) -> np.ndarray:
    """Make :math:`\\Phi = 0` on every time level, in place.

    A uniform outward normal-velocity increment is added over the fluid lateral
    faces -- so the correction *flux* is distributed in proportion to fluid face
    area, and differences between faces are untouched (design §3.3, DALES's
    ``openboundary_divcorr``).  Solid faces are left alone.

    The **pre-correction** residual is stored in ``data.net_volume_flux`` and
    returned; the **post-correction** one -- which is what the solver validates
    against at initialisation under schema 2 -- in ``data.flux_residual``, next
    to the fluid lateral area it was summed over.  ``data.divergence_corrected``
    is set.

    What was done is recorded in ``data.correction`` (:func:`correction_report`)
    and logged: at WARNING level when the increment exceeds
    :data:`CORRECTION_WARN_FRACTION` of the boundary velocity scale on any
    level -- the net flux the parent pushed through the child's lid is then
    being forced through the lateral faces, inflow faces included, and the
    message names them.

    When ``data.initial_fields`` is present and ``project_initial`` is true, the
    initial condition is made consistent with the corrected boundary data and
    then projected onto the discretely solenoidal subspace, so that a cold start
    from it begins divergence free (:func:`sync_initial_condition`, design
    section 10.6 item 4).
    """
    masks = masks or _ALL_FLUID
    area = fluid_face_area(data, masks)
    if area <= 0.0:
        raise ConfigurationError(
            "the lateral boundary has no fluid area; cannot correct the volume flux"
        )
    tally = _CorrectionTally(data, masks)
    scale = _boundary_velocity_scale(data, masks)
    outward = _face_mean_outward_velocity(data, masks)
    residual, delta = _correct_in_place(data, masks, area)
    tally.add(residual, delta, scale, outward, 0)
    data.correction = tally.log()
    data.net_volume_flux = residual
    data.flux_residual = net_volume_flux(data, masks)
    data.fluid_lateral_area = fluid_lateral_area(data, masks)
    data.masks = None if masks is _ALL_FLUID else masks
    data.divergence_corrected = True
    if data.initial_fields is not None and project_initial:
        sync_initial_condition(data, masks)
    return residual


def verify_stored_residual(data: NestingData, masks: Optional[FaceMasks] = None) -> None:
    """Recompute ``flux_residual``/``fluid_lateral_area`` and check any cached value.

    ``masks`` defaults to ``data.masks``.  Unset fields are filled in; a set
    field that disagrees with the recomputation -- beyond round-off on the
    flux a typical boundary velocity carries -- raises
    :class:`ConfigurationError`.  A ``fluid_lateral_area`` that disagrees
    while no masks are known means the data were corrected with masks that
    were not passed on, which is the same error.
    """
    masks = masks if masks is not None else data.masks
    area = fluid_lateral_area(data, masks)
    if data.fluid_lateral_area is not None and not np.isclose(
            data.fluid_lateral_area, area, rtol=1e-12, atol=0.0):
        raise ConfigurationError(
            f"cached fluid_lateral_area = {data.fluid_lateral_area:.15g} m2 disagrees with "
            f"the {area:.15g} m2 of the {'given' if masks is not None else 'all-fluid'} "
            "masks; pass masks= (or set data.masks) to the ones the correction used"
        )
    residual = net_volume_flux(data, masks)
    if data.flux_residual is not None:
        speed = max(float(np.max(np.abs(v))) for v in boundary_faces(data).values())
        tol = 1e-12 * max(fluid_face_area(data, masks) * speed, np.finfo(np.float64).tiny)
        worst = int(np.argmax(np.abs(data.flux_residual - residual)))
        if abs(data.flux_residual[worst] - residual[worst]) > tol:
            raise ConfigurationError(
                f"cached flux_residual[{worst}] = {data.flux_residual[worst]:.6g} {FLUX_UNITS} "
                f"disagrees with the residual of the slabs as stored, {residual[worst]:.6g}; "
                "the slabs were changed after the correction, or the masks differ"
            )
    data.flux_residual = residual
    data.fluid_lateral_area = area
    data.masks = masks


# --------------------------------------------------------------------------- #
# Refinement guard
# --------------------------------------------------------------------------- #


def refinement_ratios_by_axis(data: NestingData) -> Dict[str, Optional[float]]:
    """Parent-to-child refinement ratio per axis, ``None`` where unknown.

    ``x`` is ``parent_dx / min(dx_child)``; ``y`` uses ``parent_dy`` when the
    file carries it and falls back to ``parent_dx`` (uDALES parents are
    horizontally isotropic) otherwise; ``z`` uses ``parent_dz`` against the
    smallest child ``dzf`` and is ``None`` without it; ``t`` is
    ``parent_dt / child_dt``.
    """
    grid = data.grid
    out: Dict[str, Optional[float]] = {"x": None, "y": None, "z": None, "t": None}
    if data.parent_dx:
        out["x"] = float(data.parent_dx) / float(np.min(grid.dx))
    parent_dy = data.parent_dy if data.parent_dy else data.parent_dx
    if parent_dy:
        out["y"] = float(parent_dy) / float(np.min(grid.dy))
    if data.parent_dz:
        out["z"] = float(data.parent_dz) / float(np.min(grid.dzf))
    if data.parent_dt and data.child_dt:
        out["t"] = float(data.parent_dt) / float(data.child_dt)
    return out


def refinement_ratios(data: NestingData) -> Tuple[Optional[float], Optional[float]]:
    """Return ``(spatial, temporal)`` parent-to-child refinement ratios.

    The spatial ratio is the **largest** of the per-axis ratios of
    :func:`refinement_ratios_by_axis` that are known (``x``, ``y`` and, when
    ``parent_dz`` is set, ``z``); the temporal ratio is ``parent_dt /
    child_dt``.  Either is ``None`` when the file does not carry enough
    information to compute it.
    """
    by_axis = refinement_ratios_by_axis(data)
    known = [by_axis[a] for a in ("x", "y", "z") if by_axis[a] is not None]
    spatial = max(known) if known else None
    return spatial, by_axis["t"]


def refinement_verdict(data: NestingData) -> Dict[str, Any]:
    """The refinement guard's verdict, without raising or logging.

    ``{"spatial", "temporal", "by_axis", "spatial_max", "temporal_max",
    "violations", "within_limits"}``: the ratios of :func:`refinement_ratios`
    and :func:`refinement_ratios_by_axis`, the limits of design section 10.4
    (V5), the list of violation messages (empty when within limits or when a
    ratio is unknown), and the boolean summary.
    """
    spatial, temporal = refinement_ratios(data)
    tol = 1.0 + 1.0e-12
    violations = []
    if spatial is not None and spatial > MAX_SPATIAL_REFINEMENT * tol:
        violations.append(
            f"spatial refinement ratio {spatial:.3g} exceeds the supported maximum of "
            f"{MAX_SPATIAL_REFINEMENT:g} (parent_dx = {data.parent_dx:g} m)"
        )
    if temporal is not None and temporal > MAX_TEMPORAL_REFINEMENT * tol:
        violations.append(
            f"temporal refinement ratio {temporal:.3g} exceeds the supported maximum of "
            f"{MAX_TEMPORAL_REFINEMENT:g} (parent_dt = {data.parent_dt:g} s, "
            f"child_dt = {data.child_dt:g} s)"
        )
    return {
        "spatial": spatial,
        "temporal": temporal,
        "by_axis": refinement_ratios_by_axis(data),
        "spatial_max": MAX_SPATIAL_REFINEMENT,
        "temporal_max": MAX_TEMPORAL_REFINEMENT,
        "violations": violations,
        "within_limits": not violations,
    }


def check_refinement(
    data: NestingData,
    allow_refinement_violation: bool = False,
    reason: Optional[str] = None,
    override: Optional[bool] = None,
) -> Dict[str, Any]:
    """Refuse ratios beyond the validated range (design §10.4 V5) unless allowed.

    Raises :class:`NestingRefinementError` when the spatial refinement ratio
    exceeds :data:`MAX_SPATIAL_REFINEMENT` or the temporal one exceeds
    :data:`MAX_TEMPORAL_REFINEMENT`.  With ``allow_refinement_violation`` the
    violation is logged at WARNING level together with ``reason`` -- the
    caller's stated grounds, which is what makes the choice explicit -- and
    the verdict of :func:`refinement_verdict` is returned with ``"allowed"``
    and ``"reason"`` added, for the caller to record.  ``override`` is the
    former name of ``allow_refinement_violation`` and is accepted as an alias.
    """
    if override is not None:
        allow_refinement_violation = bool(override)
    verdict = refinement_verdict(data)
    verdict["allowed"] = bool(allow_refinement_violation)
    verdict["reason"] = reason
    if verdict["violations"]:
        message = "; ".join(verdict["violations"])
        if not allow_refinement_violation:
            raise NestingRefinementError(
                message + "; pass allow_refinement_violation=True, with a reason, to write anyway"
            )
        logger.warning("udprep.nesting: %s -- allowed by the caller%s", message,
                       f": {reason}" if reason else " (no reason given)")
    return verdict


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


def stored_coordinates(data: NestingData) -> Dict[str, np.ndarray]:
    """The six coordinate variables **as written to the file**: child-relative.

    The solver validates ``xh``/``yh`` against its own grid, which starts at
    0, to ``nestio_tol = 1e-10`` of ``xlen`` (``nesting_read.f90``,
    ``nestio_validate``); so a child grid built in the parent's coordinates
    (``nesting_data_from_parent`` does that) is shifted to its own origin here
    and the offset is carried by the ``child_origin_x``/``child_origin_y``
    attributes.  ``z`` is not shifted: the child's vertical is the run's.

    A grid that is neither at 0 nor at ``child_origin_*`` is ambiguous and
    refused rather than guessed.
    """
    grid = data.grid
    out: Dict[str, np.ndarray] = {}
    for axis, origin_attr in (("x", "child_origin_x"), ("y", "child_origin_y")):
        faces = getattr(grid, f"{axis}h")
        centres = getattr(grid, f"{axis}f")
        first = float(faces[0])
        origin = float(getattr(data, origin_attr))
        tol = 1.0e-10 * max(1.0, abs(faces[-1] - faces[0]))
        if abs(first) > tol and abs(first - origin) > tol:
            raise ConfigurationError(
                f"the child grid starts at {axis} = {first:g} m but {origin_attr} = "
                f"{origin:g} m; build the grid at its own origin (0) and carry the "
                f"offset in {origin_attr}, or set {origin_attr} = {first:g}"
            )
        out[f"{axis}h"] = faces - first
        out[f"{axis}f"] = centres - first
    out["zh"] = grid.zh
    out["zf"] = grid.zf
    return out


def _global_attributes(data: NestingData, schema: int = SCHEMA_VERSION) -> Dict[str, Any]:
    grid = data.grid
    attrs: Dict[str, Any] = {
        "Conventions": "CF-1.8",
        "udales_nesting_schema": np.int32(schema),
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
    if schema >= 2:
        attrs["has_initial_condition"] = np.int32(1 if data.initial_fields is not None else 0)
        area = data.fluid_lateral_area
        if area is None:
            area = fluid_lateral_area(data)
        attrs["fluid_lateral_area"] = np.float64(area)
    # Optional extensions (OPTIONAL_GLOBAL_ATTRIBUTES): let the refinement
    # guard run on a round-tripped file.  Ignored by readers that do not know
    # about them.
    for name in OPTIONAL_GLOBAL_ATTRIBUTES:
        value = getattr(data, name)
        if value is not None:
            attrs[name] = np.float64(value)
    return attrs


#: Metadata keyword arguments :class:`NestingWriter` accepts, with their
#: defaults -- the provenance fields of :class:`NestingData`.
_WRITER_METADATA: Dict[str, Any] = {
    "rhobf": None, "rhobh": None,
    "parent_model": "unknown", "parent_dx": 0.0, "parent_dt": 0.0,
    "child_origin_x": 0.0, "child_origin_y": 0.0, "child_dt": None,
    "rotation_deg": 0.0, "created": "", "creator": "", "tool_version": TOOL_VERSION,
    "parent_dy": None, "parent_dz": None,
}


class NestingWriter:
    """Per-level append writer for ``nesting.inp.<expnr>.nc``.

    The whole-record route -- :func:`nesting_data_from_parent` into
    :func:`write_nesting_file` -- holds every parent level in memory, which
    the production record cannot afford (21 600 levels of 9.8 MB is 212 GB).
    This writer takes the header once and then one level at a time::

        with NestingWriter(path, grid, nzone, parent_dt=..., parent_dx=..., ...) as w:
            for t, (pu, pv, pw) in parent_levels():
                w.append_level(t - t0, slabs_from_parent(pgrid, pu, pv, pw, grid, nzone))
        w.diagnostics   # cadence, correction, refinement, as nesting_diagnostics()

    Each level is divergence-corrected on arrival (``correct_divergence``,
    the per-level operation :func:`apply_divergence_correction` performs on a
    whole record; the correction of one level never depends on another), its
    ``net_volume_flux`` and ``flux_residual`` are stored, and the optional
    schema-2 initial condition is synced to the corrected first level and
    projected (``project_initial``) when it is passed with level 0.  The
    cadence and correction diagnostics are accumulated level by level and
    logged at :meth:`close`; the refinement guard runs at :meth:`open`,
    before anything is written.  A failure at any point removes the partial
    file.

    ``masks`` are the child's lateral fluid masks (:func:`face_masks_from_ibm`).
    Metadata keywords are those of :class:`NestingData`: ``rhobf``, ``rhobh``,
    ``parent_model``, ``parent_dx``, ``parent_dy``, ``parent_dz``,
    ``parent_dt``, ``child_origin_x``, ``child_origin_y``, ``child_dt``,
    ``rotation_deg``, ``created``, ``creator``, ``tool_version``.  The
    stored time axis is the child's clock (:func:`check_time_axis`): pass
    ``t - times[0]``.  NetCDF only; the raw back-end is whole-record.
    """

    def __init__(
        self,
        path: os.PathLike | str,
        grid: NestGrid,
        nzone: int,
        masks: Optional[FaceMasks] = None,
        correct_divergence: bool = True,
        project_initial: bool = True,
        schema: Optional[int] = None,
        allow_refinement_violation: bool = False,
        refinement_reason: Optional[str] = None,
        divergence_corrected: Optional[bool] = None,
        **metadata: Any,
    ) -> None:
        unknown = sorted(set(metadata) - set(_WRITER_METADATA))
        if unknown:
            raise ConfigurationError(
                f"NestingWriter: unknown metadata {', '.join(unknown)}; accepted: "
                f"{', '.join(_WRITER_METADATA)}"
            )
        self.path = Path(path)
        self.grid = grid
        self.nzone = int(nzone)
        _check_nzone(grid, self.nzone)
        self.masks = masks if masks is not None else _ALL_FLUID
        self.correct_divergence = bool(correct_divergence)
        self.project_initial = bool(project_initial)
        self.schema = SCHEMA_VERSION if schema is None else int(schema)
        if self.schema not in SUPPORTED_SCHEMA_VERSIONS:
            raise ConfigurationError(
                f"unsupported nesting schema {self.schema}; this writer emits "
                f"{SUPPORTED_SCHEMA_VERSIONS}"
            )
        self.allow_refinement_violation = bool(allow_refinement_violation)
        self.refinement_reason = refinement_reason
        self.divergence_corrected = (self.correct_divergence if divergence_corrected is None
                                     else bool(divergence_corrected))
        meta = dict(_WRITER_METADATA)
        meta.update(metadata)
        rhobf = meta.pop("rhobf")
        rhobh = meta.pop("rhobh")
        self.rhobf = (np.ones(grid.ktot) if rhobf is None
                      else np.asarray(rhobf, dtype=np.float64).reshape(-1))
        self.rhobh = (np.ones(grid.ktot + 1) if rhobh is None
                      else np.asarray(rhobh, dtype=np.float64).reshape(-1))
        if self.rhobf.size != grid.ktot or self.rhobh.size != grid.ktot + 1:
            raise ConfigurationError(
                f"rhobf/rhobh must have {grid.ktot}/{grid.ktot + 1} elements"
            )
        # The header, with the attribute names _global_attributes,
        # stored_coordinates and refinement_verdict read off a NestingData.
        self.header = SimpleNamespace(
            grid=grid, nzone=self.nzone, rhobf=self.rhobf, rhobh=self.rhobh,
            divergence_corrected=self.divergence_corrected, initial_fields=None,
            fluid_lateral_area=fluid_lateral_area(self, self.masks), **meta,
        )
        self.header.parent_dt = float(self.header.parent_dt)
        self.header.parent_dx = float(self.header.parent_dx)
        self._ds = None
        self._times: list = []
        self._area = fluid_face_area(self, self.masks)
        self._correction = _CorrectionTally(self, self.masks)
        self._cadence = _CadenceTally(self.header.parent_dt, self.header.parent_dx)
        self._has_initial = False
        self.refinement: Optional[Dict[str, Any]] = None
        self.diagnostics: Optional[Dict[str, Any]] = None

    # -- what the flux functions need to see this object as --------------- #

    @property
    def ntime(self) -> int:
        return len(self._times)

    @classmethod
    def from_data(cls, path: os.PathLike | str, data: NestingData, **options: Any
                  ) -> "NestingWriter":
        """A writer carrying ``data``'s header, for writing its levels through."""
        options.setdefault("masks", data.masks)
        options.setdefault("divergence_corrected", data.divergence_corrected)
        metadata = {name: getattr(data, name) for name in _WRITER_METADATA}
        return cls(path, data.grid, data.nzone, **options, **metadata)

    # -- lifecycle --------------------------------------------------------- #

    def open(self) -> "NestingWriter":
        """Run the refinement guard, create the file and write the header."""
        if self._ds is not None:
            raise ConfigurationError(f"{self.path.name} is already open")
        if self._area <= 0.0:
            raise ConfigurationError(
                "the lateral boundary has no fluid area; cannot correct the volume flux"
            )
        self.refinement = check_refinement(
            self.header, allow_refinement_violation=self.allow_refinement_violation,
            reason=self.refinement_reason,
        )
        coords = stored_coordinates(self.header)
        Dataset = _import_dataset()
        grid = self.grid
        try:
            ds = Dataset(self.path, "w", format="NETCDF4")
            self._ds = ds
            ds.createDimension("time", None)
            ds.createDimension("zf", grid.ktot)
            ds.createDimension("zh", grid.ktot + 1)
            ds.createDimension("xf", grid.itot)
            ds.createDimension("xh", grid.itot + 1)
            ds.createDimension("yf", grid.jtot)
            ds.createDimension("yh", grid.jtot + 1)
            ds.createDimension("nz", self.nzone)
            ds.createDimension("nzh", self.nzone + 1)

            var = ds.createVariable("time", "f8", ("time",))
            var.units = "s"
            var.long_name = "parent time level, on the child's clock (starts at 0)"
            for name in _COORD_VARIABLES:
                var = ds.createVariable(name, "f8", (name,))
                var.units = "m"
                var.long_name = (f"child-relative {name}" if name[0] in "xy" else name)
                var[:] = coords[name]
            for name, dim, values in (("rhobf", "zf", self.rhobf), ("rhobh", "zh", self.rhobh)):
                var = ds.createVariable(name, "f8", (dim,))
                var.units = "kg m-3"
                var[:] = values
            var = ds.createVariable("net_volume_flux", "f8", ("time",))
            var.units = FLUX_UNITS
            var.long_name = "net volume flux through the lateral boundary before correction"
            if self.schema >= 2:
                var = ds.createVariable("flux_residual", "f8", ("time",))
                var.units = FLUX_UNITS
                var.long_name = ("net volume flux through the lateral boundary as stored, "
                                 "i.e. after any divergence correction")
            for face in FACES:
                for component in COMPONENTS:
                    name = f"{component}_{face}"
                    dims = ("time",) + slab_dimensions(face, component)
                    var = ds.createVariable(name, "f8", dims)
                    var.stagger = STAGGER[component]
                    var.units = "m s-1"
        except Exception:
            self._abandon()
            raise
        return self

    def append_level(
        self,
        time: float,
        slabs: Mapping[str, np.ndarray],
        initial_fields: Optional[Mapping[str, np.ndarray]] = None,
        net_volume_flux_before: Optional[float] = None,
        flux_residual_verified: Optional[float] = None,
    ) -> Dict[str, float]:
        """Correct (unless ``correct_divergence`` is off) and write one level.

        ``slabs`` maps the twelve variable names to 3-D arrays of
        :func:`slab_shape`; the caller's arrays are not modified.
        ``initial_fields`` (``u``/``v``/``w`` on the whole child grid) may
        only come with level 0.  The last two are for the pass-through of
        :func:`write_nesting_file`, whose levels were corrected elsewhere:
        ``net_volume_flux_before`` is the pre-correction residual to store
        (otherwise: this level's residual on arrival) and
        ``flux_residual_verified`` the residual of the stored data as
        :func:`verify_stored_residual` recomputed it, so that both back-ends
        store the same bits (otherwise: recomputed here from the level).

        Returns ``{"level", "time", "net_volume_flux", "flux_residual", "delta"}``.
        """
        if self._ds is None:
            raise ConfigurationError(f"{self.path.name}: append_level before open()")
        try:
            return self._append(float(time), slabs, initial_fields,
                                net_volume_flux_before, flux_residual_verified)
        except Exception:
            self._abandon()
            raise

    def _append(self, time, slabs, initial_fields, net_before, residual_verified
                ) -> Dict[str, float]:
        n = self.ntime
        if n == 0 and time != 0.0:
            raise ConfigurationError(
                f"the first stored time must be exactly 0 (the child's clock), got "
                f"{time:g} s; subtract the first parent time"
            )
        if n > 0 and time <= self._times[-1]:
            raise ConfigurationError(
                f"time level {n} at {time:g} s does not follow level {n - 1} at "
                f"{self._times[-1]:g} s: the time axis must be strictly increasing"
            )
        missing = [name for name in SLAB_VARIABLES if name not in slabs]
        if missing:
            raise ConfigurationError(f"missing slab variables: {', '.join(missing)}")
        level: Dict[str, np.ndarray] = {}
        for face in FACES:
            for component in COMPONENTS:
                name = f"{component}_{face}"
                arr = np.asarray(slabs[name], dtype=np.float64)
                expected = slab_shape(self.grid, self.nzone, face, component)
                if arr.shape != expected:
                    raise ConfigurationError(
                        f"{name} at level {n} has shape {arr.shape}, expected {expected}"
                    )
                if not np.all(np.isfinite(arr)):
                    raise NestingSchemaError(
                        f"{name} at level {n} contains non-finite values; NaN is an error"
                    )
                # the correction writes into the boundary-normal arrays: copy those
                if self.correct_divergence and component == _FACE_NORMAL_COMPONENT[face]:
                    arr = arr.copy()
                level[name] = np.ascontiguousarray(arr)
        view = _LevelView.wrap(self.grid, self.nzone, self.rhobf, level)
        if self.correct_divergence:
            scale = _boundary_velocity_scale(view, self.masks)
            outward = _face_mean_outward_velocity(view, self.masks)
            residual, delta = _correct_in_place(view, self.masks, self._area)
            self._correction.add(residual, delta, scale, outward, n)
            before = float(residual[0]) if net_before is None else float(net_before)
            delta = float(delta[0])
        else:
            residual = net_volume_flux(view, self.masks)
            before = float(residual[0]) if net_before is None else float(net_before)
            delta = 0.0
        after = (float(net_volume_flux(view, self.masks)[0]) if residual_verified is None
                 else float(residual_verified))
        self._cadence.add(_max_normal_speed(view, self.masks), n)

        if initial_fields is not None:
            if n != 0:
                raise ConfigurationError("the initial condition belongs with time level 0")
            if self.schema < 2:
                raise ConfigurationError(
                    "a full-domain initial condition needs schema 2; schema 1 has no place for it"
                )
            self._write_initial(view, initial_fields)

        ds = self._ds
        ds.variables["time"][n] = time
        ds.variables["net_volume_flux"][n] = before
        if self.schema >= 2:
            ds.variables["flux_residual"][n] = after
        for name, arr in level.items():
            ds.variables[name][n] = arr
        self._times.append(time)
        return {"level": n, "time": time, "net_volume_flux": before,
                "flux_residual": after, "delta": delta}

    def _write_initial(self, view: _LevelView, initial_fields: Mapping[str, np.ndarray]) -> None:
        fields = initial_fields_from_fields(
            self.grid, initial_fields["u"], initial_fields["v"], initial_fields["w"]
        )
        if not (np.allclose(self.rhobf, 1.0, rtol=0.0, atol=1e-12)
                and np.allclose(self.rhobh, 1.0, rtol=0.0, atol=1e-12)):
            raise ConfigurationError(
                "a full-domain initial condition needs rhobf == rhobh == 1 (design finding F1)"
            )
        if self.correct_divergence and self.project_initial:
            fields = {c: a.copy() for c, a in fields.items()}
            holder = SimpleNamespace(grid=self.grid, nzone=self.nzone, rhobf=self.rhobf,
                                     slabs=view.slabs, initial_fields=fields, ntime=1)
            sync_initial_condition(holder, self.masks)
        for name, arr in fields.items():
            if not np.all(np.isfinite(arr)):
                raise NestingSchemaError(f"{name}_init contains non-finite values")
        ds = self._ds
        for component in COMPONENTS:
            var = ds.createVariable(f"{component}_init", "f8", init_dimensions(component))
            var.stagger = STAGGER[component]
            var.units = "m s-1"
            var.long_name = ("full-domain initial condition at the first stored time, "
                             f"velocity component {component}")
            var[:] = fields[component]
        self._has_initial = True

    def close(self) -> Path:
        """Validate the time axis, write the global attributes, close the file.

        A time axis that fails :func:`check_time_axis` (``parent_dt`` against
        the median spacing) removes the file and raises.  Logs the cadence
        and, when correcting, the correction summary; ``self.diagnostics``
        holds both plus the refinement verdict afterwards.
        """
        if self._ds is None:
            raise ConfigurationError(f"{self.path.name}: close() before open()")
        try:
            if self.ntime == 0:
                raise ConfigurationError(f"{self.path.name}: no time level was appended")
            times = np.asarray(self._times, dtype=np.float64)
            check_time_axis(times, self.header.parent_dt, context=self.path.name)
            self.header.initial_fields = True if self._has_initial else None
            for key, value in _global_attributes(self.header, self.schema).items():
                self._ds.setncattr(key, value)
        except Exception:
            self._abandon()
            raise
        self._ds.close()
        self._ds = None
        self.diagnostics = {
            "ntime": self.ntime,
            "cadence": self._cadence.log(),
            "correction": self._correction.log() if self.correct_divergence else None,
            "refinement": self.refinement,
        }
        return self.path

    def _abandon(self) -> None:
        """Close and delete the partial file after a failure."""
        if self._ds is not None:
            try:
                self._ds.close()
            except Exception:  # pragma: no cover - already broken
                pass
            self._ds = None
        try:
            self.path.unlink()
        except OSError:
            pass

    def __enter__(self) -> "NestingWriter":
        return self.open()

    def __exit__(self, exc_type, exc, tb) -> None:
        if exc_type is not None:
            self._abandon()
            return
        self.close()


def write_nesting_file(
    path: os.PathLike | str,
    data: NestingData,
    allow_refinement_violation: bool = False,
    backend: Optional[str] = None,
    schema: Optional[int] = None,
    masks: Optional[FaceMasks] = None,
    refinement_reason: Optional[str] = None,
    override: Optional[bool] = None,
) -> Path:
    """Write ``data`` to ``nesting.inp.<expnr>.nc``, exactly per the contract.

    A thin wrapper over :class:`NestingWriter`: the levels of ``data`` are
    passed through one at a time, as they are (a correction is applied with
    :func:`apply_divergence_correction` beforehand, not here).  For a record
    that does not fit in memory use the writer directly.

    ``backend`` is ``'netcdf'`` (default, or implied by a ``.nc`` suffix) or
    ``'raw'`` (a flat stream of the slab arrays plus a JSON sidecar, design
    §6.3).  Both back-ends store bit-identical values.

    ``schema`` selects the file version: 2 (the default) writes
    ``flux_residual`` and, when ``data.initial_fields`` is set, the full-domain
    initial condition; 1 writes neither, which is what the pre-v2 writer
    produced and what the back-compatibility tests need.

    ``masks`` are the lateral fluid masks the residual is summed over; they
    default to ``data.masks``, which :func:`apply_divergence_correction`
    records.  The stored ``flux_residual`` and ``fluid_lateral_area`` are
    recomputed from the slabs here and a cached value that disagrees is an
    error (:func:`verify_stored_residual`), so a file can never claim a
    residual its data do not have.

    Refuses to write when the refinement ratios are outside the validated range
    (see :func:`check_refinement`) unless ``allow_refinement_violation`` is
    ``True``, in which case the violation and ``refinement_reason`` are logged
    (``override`` is the former name and still accepted).  The boundary
    cadence criterion :func:`cadence_courant` is evaluated and logged, at
    WARNING level above :data:`CADENCE_COURANT_MAX`.
    """
    path = Path(path)
    backend = backend or ("raw" if path.suffix in (".dat", ".bin") else "netcdf")
    schema = SCHEMA_VERSION if schema is None else int(schema)
    if schema not in SUPPORTED_SCHEMA_VERSIONS:
        raise ConfigurationError(
            f"unsupported nesting schema {schema}; this writer emits "
            f"{SUPPORTED_SCHEMA_VERSIONS}"
        )
    if schema < 2 and data.initial_fields is not None:
        raise ConfigurationError(
            "a full-domain initial condition needs schema 2; schema 1 has no place to put it"
        )
    if override is not None:
        allow_refinement_violation = bool(override)
    verify_stored_residual(data, masks)
    if backend == "raw":
        check_refinement(data, allow_refinement_violation=allow_refinement_violation,
                         reason=refinement_reason)
        arrays = dict(data.slabs)
        if data.initial_fields is not None:
            arrays.update({f"{c}_init": a for c, a in data.initial_fields.items()})
        for name, arr in arrays.items():
            if not np.all(np.isfinite(arr)):
                raise NestingSchemaError(f"{name} contains non-finite values; NaN is an error")
        cadence_courant(data, data.masks)
        return _write_raw(path, data, schema)
    if backend != "netcdf":
        raise ConfigurationError(f"unknown nesting file backend {backend!r}")
    writer = NestingWriter.from_data(
        path, data, schema=schema, correct_divergence=False, project_initial=False,
        allow_refinement_violation=allow_refinement_violation,
        refinement_reason=refinement_reason,
    )
    with writer:
        for n in range(data.ntime):
            writer.append_level(
                data.times[n],
                {name: arr[n] for name, arr in data.slabs.items()},
                initial_fields=data.initial_fields if n == 0 else None,
                net_volume_flux_before=float(data.net_volume_flux[n]),
                flux_residual_verified=float(data.flux_residual[n]),
            )
    return path


def _raw_sidecar_path(path: Path) -> Path:
    return path.with_suffix(path.suffix + ".json") if path.suffix else path.with_suffix(".json")


def _write_raw(path: Path, data: NestingData, schema: int = SCHEMA_VERSION) -> Path:
    """Raw stream back-end: one contiguous ``<f8`` blob plus a JSON sidecar."""
    grid = data.grid
    attrs = {k: (int(v) if isinstance(v, np.integer)
                 else float(v) if isinstance(v, np.floating) else v)
             for k, v in _global_attributes(data, schema).items()}
    sidecar = {
        "udales_nesting_schema": schema,
        "attributes": attrs,
        "coordinates": {name: arr.tolist() for name, arr in stored_coordinates(data).items()},
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
    written = [data.slabs[f"{component}_{face}"]
               for face in FACES for component in COMPONENTS]
    if schema >= 2:
        sidecar["flux_residual"] = data.flux_residual.tolist()
        if data.initial_fields is not None:
            for component in COMPONENTS:
                arr = data.initial_fields[component]
                sidecar["variables"].append({
                    "name": f"{component}_init",
                    "dimensions": init_dimensions(component),
                    "shape": list(arr.shape),
                    "stagger": STAGGER[component],
                })
                written.append(arr)
    _raw_sidecar_path(path).write_text(json.dumps(sidecar, indent=1))
    with open(path, "wb") as stream:
        for arr in written:
            stream.write(np.ascontiguousarray(arr, dtype="<f8").tobytes(order="C"))
    return path


def _read_raw(path: Path) -> NestingData:
    sidecar = json.loads(_raw_sidecar_path(path).read_text())
    schema = int(sidecar.get("udales_nesting_schema", -1))
    if schema not in SUPPORTED_SCHEMA_VERSIONS:
        raise NestingSchemaError(
            f"{path}: sidecar schema {sidecar.get('udales_nesting_schema')} is not one of "
            f"{SUPPORTED_SCHEMA_VERSIONS}"
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
    initial = {c: slabs.pop(f"{c}_init") for c in COMPONENTS} \
        if f"{COMPONENTS[0]}_init" in slabs else None
    residual = (np.asarray(sidecar["flux_residual"], dtype=np.float64)
                if "flux_residual" in sidecar else None)
    return NestingData(
        grid=grid,
        nzone=int(attrs["nzone"]),
        times=np.asarray(sidecar["time"], dtype=np.float64),
        slabs=slabs,
        rhobf=np.asarray(sidecar["rhobf"], dtype=np.float64),
        rhobh=np.asarray(sidecar["rhobh"], dtype=np.float64),
        net_volume_flux=np.asarray(sidecar["net_volume_flux"], dtype=np.float64),
        flux_residual=residual,
        fluid_lateral_area=(float(attrs["fluid_lateral_area"])
                            if "fluid_lateral_area" in attrs else None),
        initial_fields=initial,
        divergence_corrected=bool(int(attrs["divergence_corrected"])),
        parent_model=str(attrs["parent_model"]),
        parent_dx=float(attrs["parent_dx"]),
        parent_dt=float(attrs["parent_dt"]),
        child_origin_x=float(attrs["child_origin_x"]),
        child_origin_y=float(attrs["child_origin_y"]),
        child_dt=float(attrs["child_dt"]) if "child_dt" in attrs else None,
        rotation_deg=float(attrs["rotation_deg"]),
        parent_dy=float(attrs["parent_dy"]) if "parent_dy" in attrs else None,
        parent_dz=float(attrs["parent_dz"]) if "parent_dz" in attrs else None,
        created=str(attrs["created"]),
        creator=str(attrs["creator"]),
        tool_version=str(attrs["tool_version"]),
    )


def validate_nesting_file(path: os.PathLike | str) -> Dict[str, Any]:
    """Check a NetCDF nesting file against the schema of the contract.

    Raises :class:`NestingSchemaError` naming the first offending item; returns
    the global attributes on success.  Checked: schema version, every required
    dimension and its size, every required variable with its exact dimension
    order, the ``stagger`` tag of each slab variable and of the optional
    initial-condition block, every required global attribute, ``rotation_deg ==
    0``, monotone coordinates, and the absence of NaN/missing values.

    Both schema 1 and schema 2 are accepted; the schema-2 items are required
    only of a schema-2 file.
    """
    Dataset = _import_dataset()
    path = Path(path)
    with Dataset(path, "r") as ds:
        attrs = {key: ds.getncattr(key) for key in ds.ncattrs()}
        if "udales_nesting_schema" not in attrs:
            raise NestingSchemaError(
                f"{path.name}: missing required global attribute(s): udales_nesting_schema"
            )
        schema = int(attrs["udales_nesting_schema"])
        if schema not in SUPPORTED_SCHEMA_VERSIONS:
            raise NestingSchemaError(
                f"{path.name}: udales_nesting_schema = {schema}, expected one of "
                f"{SUPPORTED_SCHEMA_VERSIONS}"
            )
        required = REQUIRED_GLOBAL_ATTRIBUTES_V2 if schema >= 2 else REQUIRED_GLOBAL_ATTRIBUTES
        missing = [key for key in required if key not in attrs]
        if missing:
            raise NestingSchemaError(
                f"{path.name}: missing required global attribute(s): {', '.join(missing)}"
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
        if schema >= 2:
            expected_vars["flux_residual"] = ("time",)
        for name in _COORD_VARIABLES:
            expected_vars[name] = (name,)
        for face in FACES:
            for component in COMPONENTS:
                expected_vars[f"{component}_{face}"] = ("time",) + slab_dimensions(face, component)
        has_init = schema >= 2 and int(attrs["has_initial_condition"]) != 0
        if has_init:
            for component in COMPONENTS:
                expected_vars[f"{component}_init"] = init_dimensions(component)
        elif any(f"{c}_init" in ds.variables for c in COMPONENTS):
            raise NestingSchemaError(
                f"{path.name}: an initial-condition variable is present but "
                "has_initial_condition is 0"
            )
        ntime = len(ds.dimensions["time"])
        for name, dims in expected_vars.items():
            if name not in ds.variables:
                raise NestingSchemaError(f"{path.name}: missing variable {name!r}")
            var = ds.variables[name]
            if tuple(var.dimensions) != dims:
                raise NestingSchemaError(
                    f"{path.name}: variable {name} has dimensions {tuple(var.dimensions)}, "
                    f"expected {dims}"
                )
            # One time level at a time: the production file does not fit in
            # memory, and a NaN check needs no more than a level.
            chunks = ((var[n] for n in range(ntime)) if dims[0] == "time" and len(dims) > 1
                      else (var[:],))
            for n, chunk in enumerate(chunks):
                if not np.all(np.isfinite(np.asarray(chunk, dtype=np.float64))):
                    where = f" at time level {n}" if dims[0] == "time" and len(dims) > 1 else ""
                    raise NestingSchemaError(
                        f"{path.name}: variable {name} has missing or NaN values{where}, "
                        "which are an error"
                    )
        staggered = [f"{c}_{f}" for f in FACES for c in COMPONENTS]
        if has_init:
            staggered += list(INIT_VARIABLES)
        for name in staggered:
            component = name.split("_")[0]
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
        try:
            check_time_axis(np.asarray(ds.variables["time"][:], dtype=np.float64),
                            float(attrs["parent_dt"]), context=path.name)
        except ConfigurationError as exc:
            raise NestingSchemaError(str(exc)) from exc
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
        # Filled level by level into preallocated arrays: no transient second
        # copy of the record, which is the difference between fitting and not.
        ntime = len(ds.dimensions["time"])
        slabs = {}
        for name in SLAB_VARIABLES:
            var = ds.variables[name]
            out = np.empty(var.shape, dtype=np.float64)
            for n in range(ntime):
                out[n] = var[n]
            slabs[name] = out
        initial = None
        if all(f"{c}_init" in ds.variables for c in COMPONENTS):
            initial = {c: np.ascontiguousarray(
                np.asarray(ds.variables[f"{c}_init"][:], dtype=np.float64))
                for c in COMPONENTS}
        residual = (np.asarray(ds.variables["flux_residual"][:], dtype=np.float64)
                    if "flux_residual" in ds.variables else None)
        return NestingData(
            grid=grid,
            nzone=int(attrs["nzone"]),
            times=np.asarray(ds.variables["time"][:], dtype=np.float64),
            slabs=slabs,
            rhobf=np.asarray(ds.variables["rhobf"][:], dtype=np.float64),
            rhobh=np.asarray(ds.variables["rhobh"][:], dtype=np.float64),
            net_volume_flux=np.asarray(ds.variables["net_volume_flux"][:], dtype=np.float64),
            flux_residual=residual,
            fluid_lateral_area=(float(attrs["fluid_lateral_area"])
                                if "fluid_lateral_area" in attrs else None),
            initial_fields=initial,
            divergence_corrected=bool(int(attrs["divergence_corrected"])),
            parent_model=str(attrs["parent_model"]),
            parent_dx=float(attrs["parent_dx"]),
            parent_dt=float(attrs["parent_dt"]),
            child_origin_x=float(attrs["child_origin_x"]),
            child_origin_y=float(attrs["child_origin_y"]),
            child_dt=float(attrs["child_dt"]) if "child_dt" in attrs else None,
            rotation_deg=float(attrs["rotation_deg"]),
            parent_dy=float(attrs["parent_dy"]) if "parent_dy" in attrs else None,
            parent_dz=float(attrs["parent_dz"]) if "parent_dz" in attrs else None,
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


def analytic_initial_fields(
    grid: NestGrid,
    t: float = 0.0,
    coeffs: Optional[Mapping[str, Sequence[float]]] = None,
) -> Dict[str, np.ndarray]:
    """:func:`analytic_field` on the whole child grid, at each component's stagger."""
    return {
        component: analytic_field(
            component, *np.ix_(*(grid.component_coords(component, ax) for ax in range(3))),
            t, coeffs,
        )
        for component in COMPONENTS
    }


def write_analytic_nesting_file(
    path: os.PathLike | str,
    grid: NestGrid,
    times: Sequence[float],
    nzone: int,
    coeffs: Optional[Mapping[str, Sequence[float]]] = None,
    correct_divergence: bool = False,
    allow_refinement_violation: bool = True,
    backend: Optional[str] = None,
    initial: bool = False,
    schema: Optional[int] = None,
    override: Optional[bool] = None,
    **attributes: Any,
) -> NestingData:
    """Write a nesting file whose velocities are the analytic field of §10.5.

    Used as the fixture for the in-solver tests ``TEST_NESTING_IO`` (U15--U22)
    and, uncorrected, for the flux assertion test U27.  ``correct_divergence``
    is off by default so that every stored value is exactly
    :func:`analytic_field` -- switching it on perturbs the four boundary faces.
    ``initial=True`` adds the full-domain initial-condition block, also exactly
    the analytic field when ``correct_divergence`` is off, which is what the
    cold-start tests read back point by point.
    The refinement guard is off by default (``allow_refinement_violation=True``;
    ``override`` is the former name) because the fixture grids carry no
    meaningful parent metadata.

    Returns the :class:`NestingData` that was written.
    """
    attributes.setdefault("parent_model", "analytic")
    attributes.setdefault("child_origin_x", float(grid.xh[0]))
    attributes.setdefault("child_origin_y", float(grid.yh[0]))
    times_arr = np.asarray(times, dtype=np.float64).reshape(-1)
    if times_arr.size > 1:
        attributes.setdefault("parent_dt", float(np.median(np.diff(times_arr))))
    if initial:
        attributes.setdefault(
            "initial_fields",
            analytic_initial_fields(grid, float(np.asarray(times).reshape(-1)[0]), coeffs),
        )
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
        data.flux_residual = data.net_volume_flux.copy()
    if override is not None:
        allow_refinement_violation = bool(override)
    write_nesting_file(path, data, allow_refinement_violation=allow_refinement_violation,
                       backend=backend, schema=schema,
                       refinement_reason="analytic fixture: no parent metadata")
    return data
