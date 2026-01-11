from __future__ import annotations

from typing import Any, Dict, List, Tuple

from pathlib import Path

import numpy as np
from .udprep import Section, SectionSpec
from .directshortwave import DirectShortwaveSolver
from udgeom.view3d import (
    compute_svf,
    read_view3d_output,
    resolve_view3d_exe,
    run_view3d,
    stl_to_view3d,
    write_svf,
    write_vfsparse,
)

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("radiation", {})
FIELDS: List[str] = list(DEFAULTS.keys())

class RadiationSection(Section):
    def __init__(
        self,
        name: str,
        values: Dict[str, Any],
        sim: Any | None = None,
        defaults: Dict[str, Any] | None = None,
    ) -> None:
        super().__init__(name, values, sim=sim, defaults=defaults)
        self._direct_sw_solver: DirectShortwaveSolver | None = None
        self._direct_sw_solver_cfg: Dict[str, Any] | None = None
        self._vf_cache: Any | None = None
        self._svf_cache: np.ndarray | None = None
        self._vf_cache_key: tuple | None = None

    def run_all(self) -> None:
        """Run radiation preprocessing steps."""
        steps = [
            ("compute_direct_shortwave", self.compute_direct_shortwave),
            ("compute_view_factors", self.compute_view_factors),
            ("write_vf", self.write_vf),
            ("write_svf", self.write_svf),
            ("write_netsw", self.write_netsw),
            ("write_timedepsw", self.write_timedepsw),
        ]
        self.run_steps("radiation", steps)

    def calc_direct_sw(
        self,
        nsun: np.ndarray,
        irradiance: float,
        method: str | None = None,
        **kwargs: Any,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        """
        Compute direct shortwave using the requested method.

        Parameters
        ----------
        nsun : array-like, shape (3,)
            Vector pointing toward the sun.
        irradiance : float
            Direct normal irradiance [W/m^2].
        method : {"moller", "facsec", "scanline"}, optional
            Dispatches to the corresponding implementation. Defaults to radiation.directsw_method.
            - "moller": ray casting with Moller-Trumbore triangle hits (most accurate, most expensive)
            - "facsec": ray casting with solid mask + facet-section reconstruction (accurate, cheaper)
            - "scanline": f2py scanline rasterization (no vegetation, fastest)
        kwargs : dict
            Optional arguments that define the solver configuration:
            - ray_density (float): ray grid density (default radiation.ray_density)
            - periodic_xy (bool): periodic boundaries in x/y (default radiation.periodic_xy)
            - ray_jitter (float): jitter factor for ray positions (default radiation.ray_jitter)
            - veg_data (dict): vegetation data (default: loaded from sim)
            - resolution (float): scanline resolution override (scanline only)
            The solver is cached; changing any of these options recreates it.
        """
        if method is None:
            method = getattr(self, "directsw_method", "facsec")
        method_key = method.strip().lower()
        if method_key in ("moller", "raycast", "mt"):
            method_key = "moller"
        elif method_key in ("facsec", "section"):
            method_key = "facsec"
        elif method_key in ("scanline", "f2py"):
            method_key = "scanline"
        else:
            raise ValueError(f"Unknown direct shortwave method: {method}")

        veg_data = kwargs.pop("veg_data", self._get_veg_data())
        ray_density = kwargs.pop("ray_density", getattr(self, "ray_density", 4.0))
        periodic_xy = kwargs.pop("periodic_xy", getattr(self, "periodic_xy", False))
        ray_jitter = kwargs.pop("ray_jitter", getattr(self, "ray_jitter", 1.0))
        resolution = kwargs.pop("resolution", None)
        if kwargs:
            raise ValueError(f"Unknown direct shortwave options: {', '.join(sorted(kwargs.keys()))}")

        cfg = {
            "method": method_key,
            "ray_density": float(ray_density),
            "ray_jitter": float(ray_jitter),
            "veg_key": id(veg_data) if veg_data is not None else None,
        }
        if self._direct_sw_solver is None or self._direct_sw_solver_cfg != cfg:
            sim = self._require_sim()
            self._direct_sw_solver = DirectShortwaveSolver(
                sim,
                method_key,
                ray_density=ray_density,
                ray_jitter=ray_jitter,
                veg_data=veg_data,
            )
            self._direct_sw_solver_cfg = cfg

        return self._direct_sw_solver.compute(
            nsun,
            irradiance,
            periodic_xy=periodic_xy,
            resolution=resolution,
        )

    def calc_direct_sw_moller(
        self,
        nsun: np.ndarray,
        irradiance: float,
        **kwargs: Any,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        """Most accurate and most expensive: ray casting with Moller-Trumbore hits."""
        return self.calc_direct_sw(nsun, irradiance, method="moller", **kwargs)

    def calc_direct_sw_facsec(
        self,
        nsun: np.ndarray,
        irradiance: float,
        **kwargs: Any,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        """Accurate and cheaper: ray casting on cell solids + facsec reconstruction."""
        return self.calc_direct_sw(nsun, irradiance, method="facsec", **kwargs)

    def calc_direct_sw_scanline(
        self,
        nsun: np.ndarray,
        irradiance: float,
        ray_density: float = 4.0,
        resolution: float | None = None,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
        """
        Scanline polygon rasterization via f2py-wrapped Fortran (no vegetation).
        """
        return self.calc_direct_sw(
            nsun,
            irradiance,
            method="scanline",
            ray_density=ray_density,
            resolution=resolution,
        )

    def solar_position(self) -> None:
        """Compute solar position (SPA/solarPosition.m)."""

    def ashrae_coeffs(self) -> None:
        """Load ASHRAE coefficients (ASHRAE.m)."""

    def calc_view_factors(self, maxD: float = 250.0, force: bool = False):
        """
        Export geometry, run View3D, and load view factors + sky view factors.

        Returns
        -------
        vf : scipy.sparse.spmatrix
            View factor matrix.
        svf : np.ndarray
            Sky view factor per facet.
        paths : dict
            Paths for the generated files (vs3, vf, svf, vfsparse).
        """
        sim = self._require_sim()
        if sim.geom is None or sim.geom.stl is None:
            raise RuntimeError("UDBase geometry not loaded; cannot run View3D")

        expnr = getattr(sim, "expnr", "")
        out_dir = Path(sim.path) if getattr(sim, "path", None) is not None else Path.cwd()
        out_dir.mkdir(parents=True, exist_ok=True)

        view3d_out = int(self.view3d_out)
        maxD = float(maxD)

        stl_path = Path(sim.path) / sim.stl_file
        stl_mtime = stl_path.stat().st_mtime if stl_path.exists() else None
        nfacets = sim.geom.stl.faces.shape[0]
        cache_key = (str(stl_path), stl_mtime, view3d_out, maxD, nfacets)
        if self._vf_cache is not None and self._svf_cache is not None and self._vf_cache_key == cache_key:
            return self._vf_cache, self._svf_cache, {
                "vs3": None,
                "vf": None,
                "svf": None,
                "vfsparse": None,
            }

        vf_path = None
        if view3d_out == 0:
            vf_path = out_dir / "vf.txt"
        elif view3d_out == 1:
            vf_path = out_dir / "vf.bin"
        elif view3d_out == 2:
            vf_path = out_dir / f"vfsparse.inp.{expnr}"
        else:
            raise ValueError(f"Unsupported view3d_out: {view3d_out}")

        svf_path = out_dir / f"svf.inp.{expnr}"
        if not force and vf_path.exists() and svf_path.exists():
            vf = read_view3d_output(vf_path, nfacets=nfacets, outformat=view3d_out)
            svf = np.loadtxt(svf_path)
            vfsparse_path = None
            if view3d_out in (0, 1) and bool(getattr(sim, "lvfsparse", False)):
                vfsparse_path = out_dir / f"vfsparse.inp.{expnr}"
                if not vfsparse_path.exists():
                    write_vfsparse(vfsparse_path, vf, threshold=5e-7)
            self._vf_cache = vf
            self._svf_cache = svf
            self._vf_cache_key = cache_key
            return vf, svf, {
                "vs3": None,
                "vf": vf_path,
                "svf": svf_path,
                "vfsparse": vfsparse_path,
            }

        vs3_path = out_dir / f"facets.{expnr}.vs3"
        stl_to_view3d(stl_path, vs3_path, view3d_out, maxD=maxD, row=0, col=0)

        view3d_exe = resolve_view3d_exe()
        if not view3d_exe.exists():
            raise FileNotFoundError(
                f"View3D executable not found at {view3d_exe}. "
                "Set VIEW3D_EXE or build tools/View3D to enable."
            )

        run_view3d(view3d_exe, vs3_path, vf_path, check=True)

        vf = read_view3d_output(vf_path, nfacets=nfacets, outformat=view3d_out)
        svf = compute_svf(vf)
        write_svf(svf_path, svf)

        vfsparse_path = None
        if view3d_out in (0, 1) and bool(getattr(sim, "lvfsparse", False)):
            vfsparse_path = out_dir / f"vfsparse.inp.{expnr}"
            write_vfsparse(vfsparse_path, vf, threshold=5e-7)

        self._vf_cache = vf
        self._svf_cache = svf
        self._vf_cache_key = cache_key

        return vf, svf, {
            "vs3": vs3_path,
            "vf": vf_path,
            "svf": svf_path,
            "vfsparse": vfsparse_path,
        }

    def calc_reflections_sw(
        self,
        sdir: np.ndarray,
        dsky: float,
        vf,
        svf: np.ndarray,
        albedo: np.ndarray,
        tol: float = 0.01,
        max_iter: int = 1000,
    ) -> np.ndarray:
        """
        Compute net shortwave including reflections (tools/SEB/netShortwave.m).

        Parameters
        ----------
        sdir : np.ndarray
            Direct shortwave on facets [W/m^2].
        dsky : float
            Diffuse sky irradiance [W/m^2].
        vf : array-like or sparse matrix
            View factor matrix between facets.
        svf : np.ndarray
            Sky view factor per facet.
        albedo : np.ndarray
            Facet albedo (0-1).
        tol : float
            Convergence threshold for additional absorbed energy.
        max_iter : int
            Safety cap on the iteration count.
        """
        sdir = np.asarray(sdir, dtype=float)
        svf = np.asarray(svf, dtype=float)
        albedo = np.asarray(albedo, dtype=float)
        if sdir.shape != svf.shape or sdir.shape != albedo.shape:
            raise ValueError("sdir, svf, and albedo must have matching shapes")
        if tol <= 0.0:
            raise ValueError("tol must be > 0")
        if max_iter <= 0:
            raise ValueError("max_iter must be > 0")

        kin0 = sdir + dsky * svf
        knet = (1.0 - albedo) * kin0
        kout = albedo * kin0

        for _ in range(max_iter):
            vf_kout = vf @ kout
            kadd = (1.0 - albedo) * vf_kout
            kout = albedo * vf_kout
            knet = knet + kadd

            denom = np.maximum(knet - kadd, 1.0e-12)
            if np.max(kadd / denom) < tol:
                break

        return knet

    def calc_short_wave(
        self,
        nsun: np.ndarray,
        irradiance: float,
        *,
        method: str | None = None,
        maxD: float | None = None,
        dsky: float | None = None,
        force: bool = False,
        return_vf: bool = False,
        **kwargs: Any,
    ):
        """
        Compute direct shortwave, view factors, and reflections.

        Parameters
        ----------
        nsun : array-like, shape (3,)
            Vector pointing toward the sun.
        irradiance : float
            Direct normal irradiance [W/m^2].
        method : {"moller", "facsec", "scanline"}, optional
            Direct shortwave method to use. Defaults to radiation.directsw_method.
        maxD : float, optional
            Max distance for View3D view factors. Defaults to radiation.maxD.
        dsky : float, optional
            Diffuse sky irradiance [W/m^2]. Defaults to radiation.Dsky.
        force : bool
            If True, recompute and overwrite existing outputs.
        return_vf : bool
            If True, return (vf, svf) after S_dir, K_star, S_veg.
        kwargs : dict
            Forwarded to calc_direct_sw (e.g., ray_density, periodic_xy, ray_jitter).
            View factors are cached in memory between calls when geometry/maxD match.

        Returns
        -------
        sdir : np.ndarray
            Direct shortwave on facets [W/m^2].
        K_star : np.ndarray
            Net shortwave on facets including reflections [W/m^2].
        S_veg : np.ndarray or None
            Vegetation absorption term when available; None for scanline or cached runs.
        vf : array-like or sparse matrix, optional
            View factor matrix between facets (only when return_vf=True).
        svf : np.ndarray, optional
            Sky view factor per facet (only when return_vf=True).
        """
        sim = self._require_sim()
        expnr = getattr(sim, "expnr", "")
        out_dir = Path(sim.path) if getattr(sim, "path", None) is not None else Path.cwd()
        out_dir.mkdir(parents=True, exist_ok=True)

        sdir_path = out_dir / "Sdir.txt"
        netsw_path = out_dir / f"netsw.inp.{expnr}"
        svf_path = out_dir / f"svf.inp.{expnr}"

        view3d_out = int(self.view3d_out)
        if view3d_out == 0:
            vf_path = out_dir / "vf.txt"
        elif view3d_out == 1:
            vf_path = out_dir / "vf.bin"
        else:
            vf_path = out_dir / f"vfsparse.inp.{expnr}"

        sdir = None
        s_veg = None
        if not force and sdir_path.exists():
            sdir = np.loadtxt(sdir_path)
        if sdir is None:
            sdir, s_veg, _ = self.calc_direct_sw(nsun, irradiance, method=method, **kwargs)
            np.savetxt(sdir_path, sdir, fmt="%8.2f")

        vf = None
        svf = None
        if maxD is None:
            maxD = float(getattr(self, "maxD", 250.0))

        stl_path = Path(sim.path) / sim.stl_file
        stl_mtime = stl_path.stat().st_mtime if stl_path.exists() else None
        nfacets = sim.geom.stl.faces.shape[0]
        cache_key = (str(stl_path), stl_mtime, view3d_out, float(maxD), nfacets)
        if not force and self._vf_cache is not None and self._svf_cache is not None and self._vf_cache_key == cache_key:
            vf = self._vf_cache
            svf = self._svf_cache
        if vf is None or svf is None:
            if not force and vf_path.exists() and svf_path.exists():
                vf = read_view3d_output(vf_path, nfacets=nfacets, outformat=view3d_out)
                svf = np.loadtxt(svf_path)
            if vf is None or svf is None:
                vf, svf, _ = self.calc_view_factors(maxD=maxD)
            self._vf_cache = vf
            self._svf_cache = svf
            self._vf_cache_key = cache_key

        if dsky is None:
            dsky = float(self.Dsky)
        albedo = sim.assign_prop_to_fac("al")

        k_star = None
        if not force and netsw_path.exists():
            k_star = np.loadtxt(netsw_path)
        if k_star is None:
            k_star = self.calc_reflections_sw(sdir, dsky, vf, svf, albedo)
            with netsw_path.open("w", encoding="ascii", newline="\n") as f:
                f.write("# net shortwave on facets [W/m2] (including reflections and diffusive)\n")
                np.savetxt(f, k_star, fmt="%6.4f")

        if return_vf:
            return sdir, k_star, s_veg, vf, svf
        return sdir, k_star, s_veg

    def write_timedepsw(self) -> None:
        """Write time-dependent shortwave (write_timedepsw in preprocessing.m)."""

    def _require_sim(self):
        if self.sim is None:
            raise ValueError("RadiationSection requires a UDBase instance (sim).")
        return self.sim

    def _get_veg_data(self) -> Dict[str, Any] | None:
        sim = self._require_sim()
        ltree = bool(getattr(sim, "ltrees", False))
        if not ltree:
            return None
        if hasattr(sim, "veg") and sim.veg is not None:
            return sim.veg
        if not hasattr(sim, "load_veg"):
            raise ValueError("Vegetation data not available; sim.load_veg() is required")
        return sim.load_veg(zero_based=True, cache=True)


SPEC = SectionSpec(
    name="radiation",
    fields=FIELDS,
    defaults=DEFAULTS,
    section_cls=RadiationSection,
)
