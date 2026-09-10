"""Python replacement for the MATLAB synthetic-inflow preprocessing script.

``tools/syntheticInflow/write_Reynolds_stress.m`` prepares the four input
files read by the standalone Fortran synthetic-inflow generator
(``tools/syntheticInflow/modSyntheticInflow.f90``), from a precursor run's
``tdump.<expnr>.nc`` statistics. This module reproduces that step in Python
(:func:`profile_from_tdump`, :func:`write_reynolds_stress_and_scales`), and
adds a second path that derives the same inputs directly from a precursor
run's own ``*driver_*`` plane records (:func:`profile_from_driver_files`,
:func:`scales_from_driver_files`) instead of from ``tdump.nc``.

Only the velocity files are produced (temperature/moisture variants,
``ltempeq``/``lmoist``, are out of scope here). The generator's input
contract (read in ``modSyntheticInflow.f90``'s ``read_synInflow_inputs``,
around lines 480-600) is:

``Reynolds_stress_profiles_velocity.txt``
    One header line (skipped with Fortran ``'(a80)'``), then ``ktot + 1``
    rows ``k = 0..nz`` of ``z umean R11 R21 R22 R31 R32 R33`` at cell
    *edges* from ``z = 0`` to ``z = zsize`` (``R11 = u'u'``, ``R21 = u'v'``,
    ``R22 = v'v'``, ``R31 = u'w'``, ``R32 = v'w'``, ``R33 = w'w'``, all in
    m^2/s^2; ``umean`` in m/s).

``length_time_scales_{u,v,w}.txt``
    One header line, then ``ktot + 1`` rows ``z nl_y nl_z t_scale``. ``nl_y``
    and ``nl_z`` are **integer grid-point counts** (not metres) -- the
    Fortran declares them ``INTEGER`` and later does
    ``NUY = 2*nluy`` etc. to size filter stencils in grid points. ``t_scale``
    is a time in **seconds** (the Lund/Xie-Castro AR(1) filter coefficients
    use ``EXP(-(pi/2)*(dt_sig/t_scale))`` with ``dt_sig = dtmax`` in
    seconds). These files are only read at all when the namelist flag
    ``lcalc_time_and_length_scale`` (block ``&STG``) is ``.FALSE.``; when it
    is ``.TRUE.`` (the Fortran default, set before ``read_namelist`` even
    runs) the generator computes them itself in
    ``calc_time_and_length_scale`` from the local grid spacing and the bulk
    velocity, and never opens these files. This module can still produce
    them (for the ``.FALSE.`` case, or simply for inspection/QA).

Both input paths write cell-edge profiles ``k = 0..nz`` (``nz = ktot``,
``ktot + 1`` rows), matching what ``read_synInflow_inputs`` expects.
"""

from __future__ import annotations

import argparse
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np

try:
    from exceptions import ConfigurationError, DataFormatError, DependencyError
except ImportError:  # pragma: no cover - exceptions module always ships with tools/python
    class ConfigurationError(ValueError):
        pass

    class DataFormatError(ValueError):
        pass

    class DependencyError(ImportError):
        pass


_TINY = 1.0e-12

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class ReynoldsStressProfile:
    """Cell-edge profiles for ``Reynolds_stress_profiles_velocity.txt``.

    All arrays have length ``nz + 1`` (``k = 0..nz``), ``z`` running from
    ``0`` to ``zsize`` at cell edges, matching ``obtain_profile`` in
    ``write_Reynolds_stress.m`` and the ``dummyz(k), umean(k), R11(k), ...``
    read loop in ``read_synInflow_inputs``.
    """

    z: np.ndarray
    umean: np.ndarray
    R11: np.ndarray
    R21: np.ndarray
    R22: np.ndarray
    R31: np.ndarray
    R32: np.ndarray
    R33: np.ndarray

    def validate(self) -> None:
        """Check the shape/physical-realizability constraints the task asks for.

        Raises :class:`DataFormatError` with a specific message on failure.
        """
        arrays = {
            "z": self.z,
            "umean": self.umean,
            "R11": self.R11,
            "R21": self.R21,
            "R22": self.R22,
            "R31": self.R31,
            "R32": self.R32,
            "R33": self.R33,
        }
        n = len(self.z)
        if n < 2:
            raise DataFormatError(
                f"Reynolds stress profile must have at least 2 rows (k=0..nz), got {n}"
            )
        for name, arr in arrays.items():
            if arr.shape != (n,):
                raise DataFormatError(
                    f"Reynolds stress profile field '{name}' has shape {arr.shape}, "
                    f"expected ({n},) to match z"
                )
            if not np.all(np.isfinite(arr)):
                raise DataFormatError(f"Reynolds stress profile field '{name}' has non-finite values")

        if abs(self.z[0]) > 1.0e-9:
            raise DataFormatError(f"z[0] must be 0 (bottom edge), got {self.z[0]!r}")
        if np.any(np.diff(self.z) <= 0.0):
            raise DataFormatError("z must be strictly monotonically increasing")

        for name in ("R11", "R22", "R33"):
            arr = arrays[name]
            if np.any(arr < -1.0e-9):
                raise DataFormatError(f"{name} must be >= 0 everywhere (found negative variance)")

        # Cauchy-Schwarz: |R_ab| <= sqrt(Raa * Rbb)
        def _check_cs(cross: np.ndarray, a: np.ndarray, b: np.ndarray, label: str) -> None:
            bound = np.sqrt(np.clip(a, 0.0, None) * np.clip(b, 0.0, None))
            # small numerical slack: profiles measured/estimated independently
            # need not satisfy the identity to machine precision
            slack = 1.0e-6 + 1.0e-3 * bound
            if np.any(np.abs(cross) > bound + slack):
                raise DataFormatError(
                    f"{label} violates |{label}| <= sqrt(Raa*Rbb) at at least one level"
                )

        _check_cs(self.R21, self.R11, self.R22, "R21")
        _check_cs(self.R31, self.R11, self.R33, "R31")
        _check_cs(self.R32, self.R22, self.R33, "R32")


@dataclass
class LengthTimeScales:
    """Cell-edge length/time scales for one velocity component.

    ``nl_y``/``nl_z`` are integer grid-point counts; ``t_scale`` is in
    seconds. See the module docstring for how the Fortran interprets these.
    """

    z: np.ndarray
    nl_y: np.ndarray
    nl_z: np.ndarray
    t_scale: np.ndarray

    def validate(self) -> None:
        n = len(self.z)
        if n < 2:
            raise DataFormatError(f"length/time scale profile must have at least 2 rows, got {n}")
        for name, arr, dtype_check in (
            ("nl_y", self.nl_y, True),
            ("nl_z", self.nl_z, True),
            ("t_scale", self.t_scale, False),
        ):
            if arr.shape != (n,):
                raise DataFormatError(
                    f"length/time scale field '{name}' has shape {arr.shape}, expected ({n},)"
                )
            if not np.all(np.isfinite(arr)):
                raise DataFormatError(f"length/time scale field '{name}' has non-finite values")
        if np.any(self.nl_y < 1) or np.any(self.nl_z < 1):
            raise DataFormatError("nl_y and nl_z must be >= 1 grid point")
        if np.any(self.t_scale <= 0.0):
            raise DataFormatError("t_scale must be > 0 (seconds)")
        if abs(self.z[0]) > 1.0e-9:
            raise DataFormatError(f"z[0] must be 0 (bottom edge), got {self.z[0]!r}")
        if np.any(np.diff(self.z) <= 0.0):
            raise DataFormatError("z must be strictly monotonically increasing")


# ---------------------------------------------------------------------------
# File I/O: the four generator input files
# ---------------------------------------------------------------------------

_REYNOLDS_HEADER = "z u u'u' u'v' v'v' u'w' v'w' w'w'"
_SCALE_HEADER = "z nLy nLz T"


def write_reynolds_stress_file(path: Path | str, profile: ReynoldsStressProfile) -> None:
    """Write ``Reynolds_stress_profiles_velocity.txt``.

    Replicates the ``fprintf(filename,'%15.10f ...')`` block of
    ``write_Reynolds_stress.m`` (the ``profile_vel`` write). The header line
    is skipped by the Fortran with ``READ (80, '(a80)') chmess`` so its exact
    text does not matter, only that it is a single line.
    """
    profile.validate()
    path = Path(path)
    with path.open("w") as f:
        f.write(_REYNOLDS_HEADER + "\n")
        for k in range(len(profile.z)):
            f.write(
                "%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n"
                % (
                    profile.z[k],
                    profile.umean[k],
                    profile.R11[k],
                    profile.R21[k],
                    profile.R22[k],
                    profile.R31[k],
                    profile.R32[k],
                    profile.R33[k],
                )
            )


def read_reynolds_stress_file(path: Path | str) -> ReynoldsStressProfile:
    """Parse a previously-written ``Reynolds_stress_profiles_velocity.txt``."""
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != 8:
        raise DataFormatError(
            f"{path}: expected 8 columns (z u R11 R21 R22 R31 R32 R33), got {data.shape[1]}"
        )
    profile = ReynoldsStressProfile(
        z=data[:, 0],
        umean=data[:, 1],
        R11=data[:, 2],
        R21=data[:, 3],
        R22=data[:, 4],
        R31=data[:, 5],
        R32=data[:, 6],
        R33=data[:, 7],
    )
    profile.validate()
    return profile


def write_length_time_scales_file(path: Path | str, scales: LengthTimeScales) -> None:
    """Write one ``length_time_scales_{u,v,w}.txt`` file.

    Replicates the ``fprintf(fileID,'%15.10f %8d %8d %15.10f\\n', ...)``
    block of ``write_Reynolds_stress.m``.
    """
    scales.validate()
    path = Path(path)
    with path.open("w") as f:
        f.write(_SCALE_HEADER + "\n")
        for k in range(len(scales.z)):
            f.write(
                "%15.10f %8d %8d %15.10f\n"
                % (scales.z[k], int(scales.nl_y[k]), int(scales.nl_z[k]), scales.t_scale[k])
            )


def read_length_time_scales_file(path: Path | str) -> LengthTimeScales:
    """Parse a previously-written ``length_time_scales_{u,v,w}.txt``."""
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != 4:
        raise DataFormatError(f"{path}: expected 4 columns (z nLy nLz T), got {data.shape[1]}")
    scales = LengthTimeScales(
        z=data[:, 0],
        nl_y=data[:, 1].astype(int),
        nl_z=data[:, 2].astype(int),
        t_scale=data[:, 3],
    )
    scales.validate()
    return scales


def write_reynolds_stress_and_scales(
    directory: Path | str,
    profile: ReynoldsStressProfile,
    scales_u: LengthTimeScales,
    scales_v: LengthTimeScales,
    scales_w: LengthTimeScales,
) -> None:
    """Write all four velocity input files into ``directory``.

    ``directory`` is normally ``experiments/<expnr>/syntheticInflow_inputs``.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    write_reynolds_stress_file(directory / "Reynolds_stress_profiles_velocity.txt", profile)
    write_length_time_scales_file(directory / "length_time_scales_u.txt", scales_u)
    write_length_time_scales_file(directory / "length_time_scales_v.txt", scales_v)
    write_length_time_scales_file(directory / "length_time_scales_w.txt", scales_w)


# ---------------------------------------------------------------------------
# Shared edge-grid helpers
# ---------------------------------------------------------------------------


def uniform_edges(zsize: float, ktot: int) -> np.ndarray:
    """Cell edges ``z = 0..zsize`` for a uniform grid, length ``ktot + 1``.

    Matches ``zdriver(0:nz)`` in ``modSyntheticInflow.f90`` for a uniform
    ``prof.inp`` (``zh(kb)=0``, ``zh(k+1)=zh(k)+2*(zf(k)-zh(k))``, plus the
    extrapolated top edge), which reduces to ``k*dz`` when ``dz`` is
    constant.
    """
    if ktot < 1:
        raise ConfigurationError(f"ktot must be >= 1, got {ktot}")
    return np.linspace(0.0, zsize, ktot + 1)


def uniform_zf_zh(zsize: float, ktot: int) -> Tuple[np.ndarray, np.ndarray]:
    """Cell centres ``zf`` (length ``ktot``) and the bottom ``ktot`` edges ``zh``.

    ``zh[k]`` is the bottom face of cell ``k`` (``zh[0] = 0``); it excludes
    the domain top edge ``zsize``, matching the Fortran ``zh(kb:ktot)`` array
    in ``read_namelist`` (as opposed to ``zdriver(0:nz)``, which is one
    element longer and does include ``zsize``). See :func:`uniform_edges` for
    the full ``ktot + 1`` edge set used by the *output* profiles.
    """
    dz = zsize / ktot
    zf = (np.arange(ktot) + 0.5) * dz
    zh = np.arange(ktot) * dz
    return zf, zh


def _interp_to_edges_like_matlab(z_source: np.ndarray, values: np.ndarray, z_target: np.ndarray) -> np.ndarray:
    """Interpolate a profile from its native z grid onto ``z_target`` edges.

    Reproduces ``obtain_profile`` in ``write_Reynolds_stress.m``:

    - the bottom edge (``z_target[0]``, always 0) is left at ``0`` rather
      than interpolated/extrapolated (``profile(1)`` in the MATLAB is never
      assigned and stays at its ``zeros(len,1)`` initial value) -- a no-slip
      / zero-variance-at-the-wall convention applied uniformly to the mean
      and all stress components;
    - interior edges are linearly interpolated from ``(z_source, values)``;
    - the top edge is linearly *extrapolated* from the last two source
      points: ``vec(end) + (vec(end)-vec(end-1))/(z(end)-z(end-1)) *
      (z_target[-1]-z(end))``, matching ``profile(len) = vec(nvec) + ((...))``.
    """
    z_source = np.asarray(z_source, dtype=float)
    values = np.asarray(values, dtype=float)
    if z_source.shape != values.shape:
        raise DataFormatError("z_source and values must have the same shape")
    if len(z_source) < 2:
        raise DataFormatError("need at least 2 source points to interpolate/extrapolate")

    out = np.empty(len(z_target), dtype=float)
    out[0] = 0.0
    if len(z_target) > 2:
        out[1:-1] = np.interp(z_target[1:-1], z_source, values)
    slope = (values[-1] - values[-2]) / (z_source[-1] - z_source[-2])
    out[-1] = values[-1] + slope * (z_target[-1] - z_source[-1])
    return out


# ---------------------------------------------------------------------------
# Path (a): from a tdump.<expnr>.nc precursor statistics file
# ---------------------------------------------------------------------------

# uDALES tdump.nc dimension-name fragments used to auto-detect axes,
# independent of axis order (netCDF4-python preserves the on-disk dimension
# order, which need not match MATLAB's ncread convention).
_X_DIM_HINTS = ("xt", "xm")
_Y_DIM_HINTS = ("yt", "ym")
_T_DIM_HINTS = ("time", "t")


def _find_axis(dimnames: Sequence[str], hints: Sequence[str]) -> Optional[int]:
    for i, name in enumerate(dimnames):
        if name.lower() in hints:
            return i
    return None


def _profile_at_iplane(dataset, varname: str, iplane: int) -> Tuple[np.ndarray, np.ndarray]:
    """y- and time-averaged 1-D profile of ``varname`` at x-index ``iplane``.

    Reproduces the ``var_y_m = squeeze(mean(squeeze(mean(var(:,:,:,column1:
    column2),4)),2)); vec = var_y_m(iplane,:)`` lines of ``obtain_profile``
    in ``write_Reynolds_stress.m`` (generalised to average over *all*
    available time records rather than a hard-coded column range, since a
    general tdump.nc need not share that particular site's time indexing).

    Returns ``(z_source, profile)`` where ``z_source`` is whichever of
    ``zt``/``zm`` the variable is defined on.
    """
    var = dataset.variables[varname]
    dimnames = list(var.dimensions)
    data = np.asarray(var[:])

    x_axis = _find_axis(dimnames, _X_DIM_HINTS)
    y_axis = _find_axis(dimnames, _Y_DIM_HINTS)
    t_axis = _find_axis(dimnames, _T_DIM_HINTS)
    z_axis = next((i for i in range(len(dimnames)) if i not in (x_axis, y_axis, t_axis)), None)
    if x_axis is None or z_axis is None:
        raise DataFormatError(f"could not identify x/z axes of '{varname}' from dims {dimnames}")
    z_dim = dimnames[z_axis]
    z_source = np.asarray(dataset.variables[z_dim][:])

    # Transpose to a canonical (x, [y], [t], z) axis order, then average away
    # whichever of y/t are present: the result is always (x, z), regardless
    # of the variable's on-disk axis order.
    order = [x_axis] + [a for a in (y_axis, t_axis) if a is not None] + [z_axis]
    data = np.transpose(data, order)
    if len(order) > 2:
        data = data.mean(axis=tuple(range(1, len(order) - 1)))

    profile = np.asarray(data[iplane, :]).reshape(-1)
    if profile.shape != z_source.shape:
        raise DataFormatError(
            f"'{varname}': extracted profile length {profile.shape} does not match "
            f"'{z_dim}' length {z_source.shape}"
        )
    return z_source, profile


def profile_from_tdump(
    nc_path: Path | str,
    iplane: int,
    zsize: float,
    ktot: int,
) -> ReynoldsStressProfile:
    """Build a :class:`ReynoldsStressProfile` from a tdump.<expnr>.nc file.

    Replaces ``write_Reynolds_stress.m``'s main body: reads ``ut``,
    ``upuptc``, ``vpvptc``, ``wpwptc``, ``upvpt``, ``upwpt``, ``vpwpt``,
    y-and-time-averages each at x-index ``iplane`` (:func:`_profile_from_at_iplane`),
    then interpolates each onto the ``ktot + 1`` target edges
    (:func:`_interp_to_edges_like_matlab`), exactly as ``obtain_profile``
    does for ``profile_vel``.

    Parameters
    ----------
    nc_path:
        Path to ``tdump.<expnr>.nc``.
    iplane:
        Zero-based x-index of the precursor plane to extract (MATLAB's
        ``iplane`` is 1-based; adjust by one when porting a value from an
        existing ``.m`` script).
    zsize, ktot:
        Target simulation's domain height and number of z cells; the output
        profile is on the ``ktot + 1``-point edge grid from
        :func:`uniform_edges`.
    """
    try:
        from netCDF4 import Dataset
    except ImportError as exc:  # pragma: no cover - exercised only without netCDF4
        raise DependencyError("profile_from_tdump requires the 'netCDF4' package") from exc

    z_target = uniform_edges(zsize, ktot)

    with Dataset(nc_path, "r") as ds:
        z_u, u = _profile_at_iplane(ds, "ut", iplane)
        _, upup = _profile_at_iplane(ds, "upuptc", iplane)
        _, vpvp = _profile_at_iplane(ds, "vpvptc", iplane)
        _, wpwp = _profile_at_iplane(ds, "wpwptc", iplane)
        _, upvp = _profile_at_iplane(ds, "upvpt", iplane)
        z_w, upwp = _profile_at_iplane(ds, "upwpt", iplane)
        _, vpwp = _profile_at_iplane(ds, "vpwpt", iplane)

    umean = _interp_to_edges_like_matlab(z_u, u, z_target)
    R11 = _interp_to_edges_like_matlab(z_u, upup, z_target)
    R22 = _interp_to_edges_like_matlab(z_u, vpvp, z_target)
    R33 = _interp_to_edges_like_matlab(z_u, wpwp, z_target)
    R21 = _interp_to_edges_like_matlab(z_u, upvp, z_target)
    R31 = _interp_to_edges_like_matlab(z_w, upwp, z_target)
    R32 = _interp_to_edges_like_matlab(z_w, vpwp, z_target)

    profile = ReynoldsStressProfile(
        z=z_target, umean=umean, R11=R11, R21=R21, R22=R22, R31=R31, R32=R32, R33=R33
    )
    profile.validate()
    return profile


# ---------------------------------------------------------------------------
# Path (b): from the solver's own *driver plane records
# ---------------------------------------------------------------------------


def jh_kh_for_advection(iadv_mom: int, iadv_tke: int = -1, iadv_thl: int = -1, iadv_qt: int = -1) -> Tuple[int, int]:
    """Replicate the ``jh``/``kh`` halo-width selection in ``read_namelist``.

    ``modSyntheticInflow.f90`` (and ``modglobal.f90``) size the halo from the
    widest advection stencil in use: ``iadv_kappa`` (7) needs ``jh=2, kh=1``;
    ``iadv_cd2`` (2) needs ``jh=1, kh=1``; ``iadv_fluxlimiter`` (5) needs
    ``jh=2, kh=2``. Unset (``-1``) secondary schemes fall back to
    ``iadv_mom``, as the Fortran does before building ``advarr``.
    """
    iadv_cd2, iadv_fluxlimiter, iadv_kappa = 2, 5, 7
    if iadv_tke < 0:
        iadv_tke = iadv_mom
    if iadv_thl < 0:
        iadv_thl = iadv_mom
    if iadv_qt < 0:
        iadv_qt = iadv_mom
    advarr = {iadv_mom, iadv_tke, iadv_thl, iadv_qt}
    if iadv_kappa in advarr:
        return 2, 1
    if iadv_cd2 in advarr:
        return 1, 1
    if iadv_fluxlimiter in advarr:
        return 2, 2
    raise ConfigurationError(
        f"no recognised advection scheme in {sorted(advarr)}; expected one of "
        f"iadv_cd2={iadv_cd2}, iadv_fluxlimiter={iadv_fluxlimiter}, iadv_kappa={iadv_kappa}"
    )


def read_driver_time(directory: Path | str, expnr: int, ipy: int = 0) -> np.ndarray:
    """Read ``tdriver_<ipy>.<expnr>``: one ``real(8)`` timestamp per record."""
    path = Path(directory) / f"tdriver_{ipy:03d}.{expnr:03d}"
    raw = path.read_bytes()
    if len(raw) % 8 != 0:
        raise DataFormatError(f"{path}: size {len(raw)} is not a multiple of 8 bytes (real(8))")
    return np.frombuffer(raw, dtype="<f8").copy()


def read_driver_component(
    directory: Path | str,
    expnr: int,
    component: str,
    ipy: int,
    jtot_local: int,
    ktot: int,
    jh: int,
    kh: int,
) -> np.ndarray:
    """Read one rank's ``{u,v,w}driver_<ipy>.<expnr>`` plane records.

    Each record holds ``u0(irecydriver, jb-jh:je+jh, kb-kh:ke+kh)`` (see
    ``writedriverfile`` in ``src/moddriver.f90``, and the equivalent
    ``write_driver`` in the generator itself) written by an implied
    Fortran do-loop, i.e. with the ``j`` index varying fastest. Returns the
    **physical** (halo-stripped) block, shape ``(nt, jtot_local, ktot)``.
    """
    if component not in ("u", "v", "w"):
        raise ConfigurationError(f"component must be one of 'u', 'v', 'w', got {component!r}")
    path = Path(directory) / f"{component}driver_{ipy:03d}.{expnr:03d}"
    n_j = jtot_local + 2 * jh
    n_k = ktot + 2 * kh
    record_reals = n_j * n_k
    record_bytes = record_reals * 8
    raw = path.read_bytes()
    if len(raw) % record_bytes != 0:
        raise DataFormatError(
            f"{path}: size {len(raw)} is not a multiple of the expected record size "
            f"{record_bytes} = ({jtot_local}+2*{jh})*({ktot}+2*{kh})*8"
        )
    nt = len(raw) // record_bytes
    # Each record is a flat run of record_reals reals with j (the leftmost
    # Fortran subscript) varying fastest, i.e. Fortran/column-major order;
    # reshape each record as (n_j, n_k) with order='F' to recover that.
    flat = np.frombuffer(raw, dtype="<f8").reshape(nt, record_reals)
    field = np.stack([rec.reshape(n_j, n_k, order="F") for rec in flat], axis=0)
    return field[:, jh : jh + jtot_local, kh : kh + ktot]


def assemble_driver_plane(
    directory: Path | str,
    expnr: int,
    component: str,
    nprocy: int,
    jtot: int,
    ktot: int,
    jh: int,
    kh: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Concatenate one component's per-rank driver files into ``(t, field)``.

    ``field`` has shape ``(nt, jtot, ktot)``: halos dropped, rank columns
    concatenated along y in rank order (``driverid = 0..nprocy-1`` maps to
    contiguous y-blocks of size ``jtot/nprocy``, per ``write_driver``'s
    ``y_driver_len*driverid+1-jh : y_driver_len*(driverid+1)+jh``).
    """
    if jtot % nprocy != 0:
        raise ConfigurationError(f"jtot={jtot} must be divisible by nprocy={nprocy}")
    jtot_local = jtot // nprocy
    t = read_driver_time(directory, expnr, ipy=0)
    blocks = [
        read_driver_component(directory, expnr, component, ipy, jtot_local, ktot, jh, kh)
        for ipy in range(nprocy)
    ]
    nt = min(b.shape[0] for b in blocks + [t.reshape(-1, 1)])
    field = np.concatenate([b[:nt] for b in blocks], axis=1)
    return t[:nt], field


# --- turbulence-scale estimators -------------------------------------------


def _autocorr_1d(x: np.ndarray) -> np.ndarray:
    """Normalised autocorrelation ``rho[lag] = <x'(n) x'(n+lag)> / <x'^2>``.

    ``x`` is a single (1-D) realisation; ``rho`` has the same length as
    ``x`` (only the first half is meaningful for an integral-scale estimate).
    """
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = len(x)
    var = np.dot(x, x) / n
    if var < _TINY:
        return np.zeros(n)
    # full autocorrelation via FFT (fast, exact for the direct/linear case
    # padded with zeros -- we only use the first half of lags anyway)
    size = 1
    while size < 2 * n:
        size *= 2
    fx = np.fft.rfft(x, n=size)
    acf = np.fft.irfft(fx * np.conjugate(fx), n=size)[:n]
    acf /= n * var
    # correct for the shrinking overlap count at each lag (unbiased-ish)
    acf *= n / (n - np.arange(n))
    return acf


def _integral_scale_to_first_zero(rho: np.ndarray, spacing: float, max_lag: Optional[int] = None) -> float:
    """Trapezoidal integral of ``rho`` from lag 0 to its first zero crossing.

    This is the standard "integrate the autocorrelation out to where it
    first crosses zero" definition of an integral time/length scale.
    """
    if max_lag is not None:
        rho = rho[: max_lag + 1]
    negative = np.nonzero(rho <= 0.0)[0]
    stop = negative[0] if negative.size else len(rho) - 1
    stop = max(stop, 1)
    trapezoid = getattr(np, "trapezoid", np.trapz)
    return float(trapezoid(rho[: stop + 1], dx=spacing))


def integral_time_scale(series_ty: np.ndarray, dt: float) -> float:
    """Integral time scale of a fluctuating field, averaged over realisations.

    ``series_ty`` has shape ``(nt, n_realisations)`` (e.g. all y-columns at
    one z level); the autocorrelation is computed per realisation and
    averaged before integrating to the first zero crossing.
    """
    series_ty = np.asarray(series_ty, dtype=float)
    if series_ty.ndim == 1:
        series_ty = series_ty[:, None]
    nt = series_ty.shape[0]
    acc = np.zeros(nt)
    for j in range(series_ty.shape[1]):
        acc += _autocorr_1d(series_ty[:, j])
    acc /= series_ty.shape[1]
    max_lag = max(nt // 2, 1)
    return _integral_scale_to_first_zero(acc, dt, max_lag=max_lag)


def integral_length_scale_periodic(field_ty: np.ndarray, dy: float) -> float:
    """Integral length scale of a field periodic in y (e.g. across the domain).

    Two-point correlation as a function of periodic y-lag, averaged over
    time (and, implicitly, all starting points via the circular
    autocorrelation), integrated to the first zero crossing.
    """
    field_ty = np.asarray(field_ty, dtype=float)
    if field_ty.ndim == 1:
        field_ty = field_ty[None, :]
    ny = field_ty.shape[1]
    acc = np.zeros(ny)
    for row in field_ty:
        row = row - row.mean()
        var = np.dot(row, row) / ny
        if var < _TINY:
            continue
        f = np.fft.rfft(row)
        acf = np.fft.irfft(f * np.conjugate(f), n=ny)
        acc += acf / (ny * var)
    acc /= field_ty.shape[0]
    max_lag = ny // 2
    return _integral_scale_to_first_zero(acc, dy, max_lag=max_lag)


def integral_length_scale_nonperiodic(field_ty_allz: np.ndarray, k0: int, dz: float) -> float:
    """Integral length scale in z (non-periodic) around level ``k0``.

    ``field_ty_allz`` has shape ``(nt, ny, nz)``; the two-point correlation
    is built from the available ``k0 + lag`` offsets only (no wraparound),
    averaged over time and y, and integrated to the first zero crossing (or
    to the domain edge, whichever comes first).
    """
    field_ty_allz = np.asarray(field_ty_allz, dtype=float)
    nz = field_ty_allz.shape[2]
    anom = field_ty_allz - field_ty_allz.mean(axis=(0, 1), keepdims=True)
    ref = anom[:, :, k0]
    var0 = np.mean(ref * ref)
    if var0 < _TINY:
        return dz  # flat field: fall back to one grid spacing
    max_lag = nz - 1 - k0
    rho = np.empty(max_lag + 1)
    for lag in range(max_lag + 1):
        rho[lag] = np.mean(ref * anom[:, :, k0 + lag]) / var0
    return _integral_scale_to_first_zero(rho, dz)


def profile_from_driver_files(
    directory: Path | str,
    expnr: int,
    nprocy: int,
    jtot: int,
    ktot: int,
    zsize: float,
    jh: int,
    kh: int,
) -> ReynoldsStressProfile:
    """Build a :class:`ReynoldsStressProfile` from a precursor's own driver files.

    Reads ``tdriver_000.<expnr>`` and the per-rank ``{u,v,w}driver_<ipy>.<expnr>``
    files (:func:`assemble_driver_plane`), forms the y-t mean and the six
    second moments at each native z level, then interpolates onto the
    ``ktot + 1`` target edges with the same rule as the tdump path
    (:func:`_interp_to_edges_like_matlab`).

    ``u`` and ``v`` are natively at cell centres (``zf``); ``w`` is natively
    at the bottom ``ktot`` cell edges (``zh``) -- see ``bilinear_interp`` in
    ``modSyntheticInflow.f90``, which builds ``u0``/``v0`` from a
    ``dz``-weighted average of the two edge values bracketing ``zf(k)``,
    while ``w0(k) = wdriver(k-1)``, i.e. exactly the edge value.
    """
    zf, zh = uniform_zf_zh(zsize, ktot)
    z_target = uniform_edges(zsize, ktot)

    _, u = assemble_driver_plane(directory, expnr, "u", nprocy, jtot, ktot, jh, kh)
    _, v = assemble_driver_plane(directory, expnr, "v", nprocy, jtot, ktot, jh, kh)
    _, w = assemble_driver_plane(directory, expnr, "w", nprocy, jtot, ktot, jh, kh)

    u_mean_z = u.mean(axis=(0, 1))
    v_mean_z = v.mean(axis=(0, 1))
    w_mean_z = w.mean(axis=(0, 1))

    up = u - u_mean_z
    vp = v - v_mean_z
    wp = w - w_mean_z

    R11_native = (up * up).mean(axis=(0, 1))
    R22_native = (vp * vp).mean(axis=(0, 1))
    R33_native = (wp * wp).mean(axis=(0, 1))
    R21_native = (up * vp).mean(axis=(0, 1))
    R31_native = (up * wp).mean(axis=(0, 1))
    R32_native = (vp * wp).mean(axis=(0, 1))

    umean = _interp_to_edges_like_matlab(zf, u_mean_z, z_target)
    R11 = _interp_to_edges_like_matlab(zf, R11_native, z_target)
    R22 = _interp_to_edges_like_matlab(zf, R22_native, z_target)
    R21 = _interp_to_edges_like_matlab(zf, R21_native, z_target)
    R33 = _interp_to_edges_like_matlab(zh, R33_native, z_target)
    R31 = _interp_to_edges_like_matlab(zh, R31_native, z_target)
    R32 = _interp_to_edges_like_matlab(zh, R32_native, z_target)

    profile = ReynoldsStressProfile(
        z=z_target, umean=umean, R11=R11, R21=R21, R22=R22, R31=R31, R32=R32, R33=R33
    )
    profile.validate()
    return profile


def scales_from_driver_files(
    directory: Path | str,
    expnr: int,
    nprocy: int,
    jtot: int,
    ktot: int,
    zsize: float,
    ylen: float,
    dt: float,
    jh: int,
    kh: int,
) -> Tuple[LengthTimeScales, LengthTimeScales, LengthTimeScales]:
    """Estimate the per-component length/time scale files from driver files.

    For each native z level: the temporal integral scale
    (:func:`integral_time_scale`, pooling all y-columns as realisations),
    the periodic y integral length scale (:func:`integral_length_scale_periodic`),
    and the non-periodic z integral length scale
    (:func:`integral_length_scale_nonperiodic`). Lengths are converted to
    grid-point counts with ``max(int(L/dy), 1)`` / ``max(int(L/dz), 1)``,
    mirroring the truncating ``INT()`` used in ``calc_time_and_length_scale``.
    Returns ``(scales_u, scales_v, scales_w)``, each interpolated onto the
    ``ktot + 1`` target edges like the Reynolds-stress profile.
    """
    zf, zh = uniform_zf_zh(zsize, ktot)
    z_target = uniform_edges(zsize, ktot)
    dy = ylen / jtot
    dz = zsize / ktot

    def _scales_for(component: str, z_native: np.ndarray) -> LengthTimeScales:
        _, field = assemble_driver_plane(directory, expnr, component, nprocy, jtot, ktot, jh, kh)
        nz_native = field.shape[2]
        t_scale = np.empty(nz_native)
        nl_y = np.empty(nz_native, dtype=int)
        nl_z = np.empty(nz_native, dtype=int)
        for k in range(nz_native):
            t_scale[k] = integral_time_scale(field[:, :, k], dt)
            l_y = integral_length_scale_periodic(field[:, :, k], dy)
            nl_y[k] = max(int(l_y / dy), 1)
            l_z = integral_length_scale_nonperiodic(field, k, dz)
            nl_z[k] = max(int(l_z / dz), 1)

        t_edge = _interp_to_edges_like_matlab(z_native, t_scale, z_target)
        t_edge = np.clip(t_edge, dt, None)  # keep strictly positive at the extrapolated edges
        nl_y_edge = np.round(_interp_to_edges_like_matlab(z_native, nl_y.astype(float), z_target))
        nl_z_edge = np.round(_interp_to_edges_like_matlab(z_native, nl_z.astype(float), z_target))
        nl_y_edge = np.clip(nl_y_edge, 1, None).astype(int)
        nl_z_edge = np.clip(nl_z_edge, 1, None).astype(int)
        # the k=0 (ground) edge follows the same "hold nearest interior value"
        # convention the Fortran itself applies to its own self-computed
        # scales (calc_time_and_length_scale sets index 0 := index 1) rather
        # than the "force to zero" rule used for the mean/stress profiles,
        # since a zero length/time scale is not a meaningful input.
        t_edge[0] = t_scale[0]
        nl_y_edge[0] = nl_y[0]
        nl_z_edge[0] = nl_z[0]
        scales = LengthTimeScales(z=z_target, nl_y=nl_y_edge, nl_z=nl_z_edge, t_scale=t_edge)
        scales.validate()
        return scales

    return (
        _scales_for("u", zf),
        _scales_for("v", zf),
        _scales_for("w", zh),
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _cli_from_tdump(args: argparse.Namespace) -> None:
    profile = profile_from_tdump(args.tdump, args.iplane, args.zsize, args.ktot)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_reynolds_stress_file(out_dir / "Reynolds_stress_profiles_velocity.txt", profile)
    print(f"Wrote {out_dir / 'Reynolds_stress_profiles_velocity.txt'}")
    print(
        "Note: length_time_scales_{u,v,w}.txt were not written -- a tdump.nc "
        "does not carry the two-point statistics needed to estimate them. "
        "Either set '&STG lcalc_time_and_length_scale = .TRUE.' (the Fortran "
        "default) so the generator computes them itself, or use the "
        "'from-driver' subcommand, which estimates them from a precursor's "
        "own *driver_* files."
    )


def _cli_from_driver(args: argparse.Namespace) -> None:
    jh, kh = args.jh, args.kh
    profile = profile_from_driver_files(
        args.directory, args.expnr, args.nprocy, args.jtot, args.ktot, args.zsize, jh, kh
    )
    out_dir = Path(args.output)
    scales_u, scales_v, scales_w = scales_from_driver_files(
        args.directory,
        args.expnr,
        args.nprocy,
        args.jtot,
        args.ktot,
        args.zsize,
        args.ylen,
        args.dt,
        jh,
        kh,
    )
    write_reynolds_stress_and_scales(out_dir, profile, scales_u, scales_v, scales_w)
    print(f"Wrote synthetic-inflow input files to {out_dir}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m udprep.synthetic_inflow",
        description=(
            "Generate synthetic-inflow input files "
            "(Reynolds_stress_profiles_velocity.txt, length_time_scales_{u,v,w}.txt) "
            "for tools/syntheticInflow/modSyntheticInflow.f90."
        ),
    )
    sub = parser.add_subparsers(dest="mode", required=True)

    p_tdump = sub.add_parser("from-tdump", help="from a tdump.<expnr>.nc precursor statistics file")
    p_tdump.add_argument("tdump", type=Path, help="path to tdump.<expnr>.nc")
    p_tdump.add_argument("--iplane", type=int, required=True, help="zero-based x-index of the precursor plane")
    p_tdump.add_argument("--zsize", type=float, required=True)
    p_tdump.add_argument("--ktot", type=int, required=True)
    p_tdump.add_argument("--output", type=Path, required=True, help="syntheticInflow_inputs directory to write")
    p_tdump.set_defaults(func=_cli_from_tdump)

    p_driver = sub.add_parser("from-driver", help="from a precursor run's own *driver_* files")
    p_driver.add_argument("directory", type=Path, help="directory containing the *driver_* files")
    p_driver.add_argument("--expnr", type=int, required=True)
    p_driver.add_argument("--nprocy", type=int, required=True)
    p_driver.add_argument("--jtot", type=int, required=True)
    p_driver.add_argument("--ktot", type=int, required=True)
    p_driver.add_argument("--zsize", type=float, required=True)
    p_driver.add_argument("--ylen", type=float, required=True)
    p_driver.add_argument("--dt", type=float, required=True, help="dtdriver, seconds between records")
    p_driver.add_argument("--jh", type=int, default=1)
    p_driver.add_argument("--kh", type=int, default=1)
    p_driver.add_argument("--output", type=Path, required=True, help="syntheticInflow_inputs directory to write")
    p_driver.set_defaults(func=_cli_from_driver)

    ns = parser.parse_args(argv)
    ns.func(ns)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
