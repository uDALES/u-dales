from __future__ import annotations

from datetime import datetime, timedelta
from typing import Any, Dict, List, Tuple

from pathlib import Path

import numpy as np
from .udprep import Section, SectionSpec
from .directshortwave import DirectShortwaveSolver
from .solar import nsun_from_angles, solar_position_python, solar_state, solar_strength_ashrae
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

    def run_all(self, force: bool = False) -> None:
        """Run radiation preprocessing steps."""
        if bool(getattr(self, "ltimedepsw", False)):
            steps = [("run_short_wave_timedep", self.run_short_wave_timedep)]
        else:
            steps = [("run_short_wave", self.run_short_wave)]
        self.run_steps("radiation", steps, force=force)

    def save(self, force: bool = False) -> None:
        """
        Run radiation preprocessing and write updated parameters.
        """
        self.run_short_wave(force=force)
        self.run_short_wave_timedep(force=force)
        self.write_changed_params()

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

    def solar_position(
        self,
        time_of_day: datetime | None = None,
    ) -> Dict[str, Any]:
        """
        Compute solar position using the Python SPA port.

        Parameters
        ----------
        time_of_day : datetime, optional
            Timestamp to evaluate. Defaults to the date/time fields in the section.
        """
        if time_of_day is None:
            time_of_day = datetime(
                int(self.year),
                int(self.month),
                int(self.day),
                int(self.hour),
                int(self.minute),
                int(self.second),
            )
        return solar_position_python(
            time_of_day,
            float(self.longitude),
            float(self.latitude),
            float(self.timezone),
            float(self.elevation),
        )

    def calc_solar_state(
        self,
        *,
        backend: str = "python",
        time_of_day: datetime | None = None,
    ) -> Tuple[np.ndarray, float, float, float, float]:
        """
        Compute solar state (nsun, zenith, azimuth, I, Dsky) based on isolar.

        Returns
        -------
        nsun : np.ndarray
            Unit vector pointing toward the sun in local coordinates.
        zenith : float
            Solar zenith angle [deg].
        azimuth_local : float
            Solar azimuth angle [deg] in local coordinates (solarazimuth - xazimuth).
        I : float
            Direct normal irradiance [W/m^2].
        Dsky : float
            Diffuse sky irradiance [W/m^2].
        """
        isolar = int(getattr(self, "isolar", 1))
        xazimuth = float(getattr(self, "xazimuth", 0.0))

        if isolar == 1:
            zenith = float(self.solarzenith)
            azimuth_local = float(self.solarazimuth) - xazimuth
            nsun = nsun_from_angles(zenith, azimuth_local)
            return nsun, zenith, azimuth_local, float(self.I), float(self.Dsky)

        if isolar == 2:
            if time_of_day is None:
                time_of_day = datetime(
                    int(self.year),
                    int(self.month),
                    int(self.day),
                    int(self.hour),
                    int(self.minute),
                    int(self.second),
                )
            return solar_state(
                time_of_day,
                float(self.longitude),
                float(self.latitude),
                float(self.timezone),
                float(self.elevation),
                xazimuth=xazimuth,
                backend=backend,
            )

        if isolar == 3:
            raise NotImplementedError(
                "isolar=3 (weatherfile-based solar state) is not implemented in udprep yet."
            )

        raise ValueError(f"Unsupported isolar value: {isolar}")

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
        if view3d_out == 2 and not bool(getattr(sim, "lvfsparse", False)):
            raise ValueError("view3d_out=2 requires lvfsparse=true in the ENERGYBALANCE section")
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
            nnz = int(getattr(vf, "nnz", np.count_nonzero(vf)))
            sim.save_param("nnz", nnz)
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

        nnz = int(getattr(vf, "nnz", np.count_nonzero(vf)))
        sim.save_param("nnz", nnz)
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

    def run_short_wave(self, force: bool = False) -> None:
        """
        Compute and write single-step shortwave (Sdir + netsw) like MATLAB.
        """
        sim = self._require_sim()
        expnr = getattr(sim, "expnr", "")
        out_dir = Path(sim.path) if getattr(sim, "path", None) is not None else Path.cwd()
        out_dir.mkdir(parents=True, exist_ok=True)

        sdir_path = out_dir / "Sdir.txt"
        netsw_path = out_dir / f"netsw.inp.{expnr}"
        sveg_path = out_dir / f"sveg.inp.{expnr}"
        if not force and sdir_path.exists() and netsw_path.exists():
            return

        start = datetime(
            int(self.year),
            int(self.month),
            int(self.day),
            int(self.hour),
            int(self.minute),
            int(self.second),
        )
        nsun, _, _, irradiance, dsky = self._solar_state_time(start)

        lscatter = bool(getattr(sim, "lEB", False))
        albedo = sim.assign_prop_to_fac("al")
        face_normals = sim.geom.stl.face_normals
        fss = (1.0 + face_normals[:, 2]) * 0.5 if not lscatter else None

        method, resolution = self._shortwave_method()
        vf = None
        svf = None
        if lscatter:
            vf, svf, _ = self.calc_view_factors(maxD=float(self.maxD))

        sdir, knet, s_veg = self._compute_knet(
            nsun,
            irradiance,
            dsky,
            method,
            resolution,
            lscatter,
            albedo,
            vf,
            svf,
            fss,
        )
        self.write_netsw(knet, s_veg=s_veg)
        np.savetxt(sdir_path, sdir, fmt="%8.2f")

    def run_short_wave_timedep(self, force: bool = False) -> None:
        """
        Compute and write time-dependent shortwave (timedepsw + Sdir.nc).
        """
        if not bool(getattr(self, "ltimedepsw", False)):
            return

        sim = self._require_sim()
        out_dir = Path(sim.path) if getattr(sim, "path", None) is not None else Path.cwd()
        out_dir.mkdir(parents=True, exist_ok=True)

        expnr = getattr(sim, "expnr", "")
        sdir_nc_path = out_dir / "Sdir.nc"
        timedepsw_path = out_dir / f"timedepsw.inp.{expnr}"
        timedepsveg_path = out_dir / f"timedepsveg.inp.{expnr}"
        if (
            not force
            and sdir_nc_path.exists()
            and timedepsw_path.exists()
        ):
            if timedepsveg_path.exists() or not bool(getattr(sim, "ltrees", False)):
                return

        runtime = float(self.runtime)
        dtSP = float(self.dtSP)
        tSP = np.arange(0.0, runtime + 0.5 * dtSP, dtSP, dtype=float)
        nt = tSP.size

        lscatter = bool(getattr(sim, "lEB", False))
        albedo = sim.assign_prop_to_fac("al")
        face_normals = sim.geom.stl.face_normals
        fss = (1.0 + face_normals[:, 2]) * 0.5 if not lscatter else None

        method, resolution = self._shortwave_method()

        vf = None
        svf = None
        if lscatter:
            vf, svf, _ = self.calc_view_factors(maxD=float(self.maxD))

        sdir_all = np.zeros((albedo.size, nt), dtype=float)
        knet_all = np.zeros((albedo.size, nt), dtype=float)
        s_veg_all = np.zeros((albedo.size, nt), dtype=float)
        has_sveg = False

        start = datetime(
            int(self.year),
            int(self.month),
            int(self.day),
            int(self.hour),
            int(self.minute),
            int(self.second),
        )

        if int(getattr(self, "isolar", 1)) == 3:
            weather = self._read_weather_table(Path(self.weatherfname))
            date_val = int(start.strftime("%d%m%y"))
            rows = weather["date"] == date_val
            if not np.any(rows):
                raise ValueError(f"No weather data for date {date_val} in {self.weatherfname}")

            timedep_time = weather["TIME"][rows]
            timedep_zenith = weather["SOLAR"][rows]
            timedep_azimuth = weather["SOLAR_1"][rows] + 90.0
            timedep_I = weather["HELIOM"][rows]
            timedep_Dsky = weather["DIFSOLAR"][rows]

            shift = -start.hour
            timedep_zenith = np.roll(timedep_zenith, shift)
            timedep_azimuth = np.roll(timedep_azimuth, shift)
            timedep_I = np.roll(timedep_I, shift)
            timedep_Dsky = np.roll(timedep_Dsky, shift)

            x = np.concatenate([timedep_time, [86400.0]])
            zenith_interp = self._interp_makima(x, np.concatenate([timedep_zenith, [timedep_zenith[0]]]), tSP)
            azimuth_interp = self._interp_makima(x, np.concatenate([timedep_azimuth, [timedep_azimuth[0]]]), tSP)
            I_interp = self._interp_makima(x, np.concatenate([timedep_I, [timedep_I[0]]]), tSP)
            Dsky_interp = self._interp_makima(x, np.concatenate([timedep_Dsky, [timedep_Dsky[0]]]), tSP)

            for n, t_val in enumerate(tSP):
                solarzenith = float(zenith_interp[n])
                irradiance = float(I_interp[n])
                if solarzenith < 90.0 and irradiance > 0.0:
                    azimuth = float(azimuth_interp[n]) - float(self.xazimuth)
                    nsun = nsun_from_angles(solarzenith, azimuth)
                    dsky = float(Dsky_interp[n])
                    sdir, knet, s_veg = self._compute_knet(
                        nsun,
                        irradiance,
                        dsky,
                        method,
                        resolution,
                        lscatter,
                        albedo,
                        vf,
                        svf,
                        fss,
                    )
                    sdir_all[:, n] = sdir
                    knet_all[:, n] = knet
                    if s_veg is not None:
                        s_veg_all[:, n] = s_veg
                        has_sveg = True
        else:
            for n, t_val in enumerate(tSP):
                time_of_day = start + timedelta(seconds=float(t_val))
                nsun, solarzenith, _, irradiance, dsky = self._solar_state_time(time_of_day)
                if solarzenith < 90.0 and irradiance > 0.0:
                    sdir, knet, s_veg = self._compute_knet(
                        nsun,
                        irradiance,
                        dsky,
                        method,
                        resolution,
                        lscatter,
                        albedo,
                        vf,
                        svf,
                        fss,
                    )
                    sdir_all[:, n] = sdir
                    knet_all[:, n] = knet
                    if s_veg is not None:
                        s_veg_all[:, n] = s_veg
                        has_sveg = True

        self._write_sdir_nc(sdir_nc_path, tSP, sdir_all)
        self.write_timedepsw(tSP, knet_all)
        if has_sveg:
            self.write_timedepsveg(tSP, s_veg_all)

    def write_timedepsw(self, tSP: np.ndarray, knet: np.ndarray) -> None:
        """Write time-dependent net shortwave (timedepsw.inp.<expnr>)."""
        sim = self._require_sim()
        expnr = getattr(sim, "expnr", "")
        out_dir = Path(sim.path) if getattr(sim, "path", None) is not None else Path.cwd()
        out_dir.mkdir(parents=True, exist_ok=True)

        timedepsw_path = out_dir / f"timedepsw.inp.{expnr}"
        with timedepsw_path.open("w", encoding="ascii", newline="\n") as f:
            f.write(
                "# time-dependent net shortwave on facets [W/m2]. "
                "First line: times (1 x nt), then netsw (nfcts x nt)\n"
            )
        with timedepsw_path.open("a", encoding="ascii", newline="\n") as f:
            np.savetxt(f, tSP[None, :], fmt="%9.2f")
            np.savetxt(f, knet, fmt="%9.4f")

    def write_timedepsveg(self, tSP: np.ndarray, s_veg: np.ndarray) -> None:
        """Write time-dependent vegetation absorption (timedepsveg.inp.<expnr>)."""
        sim = self._require_sim()
        expnr = getattr(sim, "expnr", "")
        out_dir = Path(sim.path) if getattr(sim, "path", None) is not None else Path.cwd()
        out_dir.mkdir(parents=True, exist_ok=True)

        path = out_dir / f"timedepsveg.inp.{expnr}"
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write(
                "# time-dependent vegetation absorption [W/m2]. "
                "First line: times (1 x nt), then sveg (nfcts x nt)\n"
            )
        with path.open("a", encoding="ascii", newline="\n") as f:
            np.savetxt(f, tSP[None, :], fmt="%9.2f")
            np.savetxt(f, s_veg, fmt="%9.4f")

    def write_netsw(self, knet: np.ndarray, s_veg: np.ndarray | None = None) -> None:
        """Write net shortwave on facets (netsw.inp.<expnr>)."""
        sim = self._require_sim()
        expnr = getattr(sim, "expnr", "")
        out_dir = Path(sim.path) if getattr(sim, "path", None) is not None else Path.cwd()
        out_dir.mkdir(parents=True, exist_ok=True)

        netsw_path = out_dir / f"netsw.inp.{expnr}"
        with netsw_path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# net shortwave on facets [W/m2] (including reflections and diffusive)\n")
            np.savetxt(f, knet, fmt="%6.4f")
        if s_veg is not None:
            path = out_dir / f"sveg.inp.{expnr}"
            with path.open("w", encoding="ascii", newline="\n") as f:
                f.write("# vegetation absorption on facets [W/m2]\n")
                np.savetxt(f, s_veg, fmt="%6.4f")

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

    def _shortwave_method(self) -> Tuple[str, float | None]:
        ishortwave = int(getattr(self, "ishortwave", 1))
        if ishortwave == 1:
            return "scanline", float(self.psc_res)
        return getattr(self, "directsw_method", "facsec"), None

    def _solar_state_time(
        self,
        time_of_day: datetime,
    ) -> Tuple[np.ndarray, float, float, float, float]:
        isolar = int(getattr(self, "isolar", 1))
        xazimuth = float(getattr(self, "xazimuth", 0.0))

        if isolar == 1:
            zenith = float(self.solarzenith)
            azimuth = float(self.solarazimuth) - xazimuth
            nsun = nsun_from_angles(zenith, azimuth)
            return nsun, zenith, azimuth, float(self.I), float(self.Dsky)

        if isolar == 2:
            sp = solar_position_python(
                time_of_day,
                float(self.longitude),
                float(self.latitude),
                float(self.timezone),
                float(self.elevation),
            )
            zenith = float(sp["zenith"])
            azimuth = float(sp["azimuth"]) - xazimuth
            nsun = nsun_from_angles(zenith, azimuth)
            irradiance, dsky = solar_strength_ashrae(time_of_day.month, zenith)
            return nsun, zenith, azimuth, float(irradiance), float(dsky)

        if isolar == 3:
            weather = self._read_weather_table(Path(self.weatherfname))
            date_val = int(time_of_day.strftime("%d%m%y"))
            time_val = int(time_of_day.hour * 3600)
            rows = (weather["date"] == date_val) & (weather["time"] == time_val)
            if not np.any(rows):
                raise ValueError(
                    f"No weather data for date {date_val} and time {time_val} in {self.weatherfname}"
                )
            solarzenith = float(weather["SOLAR"][rows][0])
            solarazimuth = float(weather["SOLAR_1"][rows][0]) + 90.0
            irradiance = float(weather["HELIOM"][rows][0])
            dsky = float(weather["DIFSOLAR"][rows][0])
            azimuth = solarazimuth - xazimuth
            nsun = nsun_from_angles(solarzenith, azimuth)
            return nsun, solarzenith, azimuth, irradiance, dsky

        raise ValueError(f"Unsupported isolar value: {isolar}")

    def _compute_knet(
        self,
        nsun: np.ndarray,
        irradiance: float,
        dsky: float,
        method: str,
        resolution: float | None,
        lscatter: bool,
        albedo: np.ndarray,
        vf,
        svf: np.ndarray | None,
        fss: np.ndarray | None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray | None]:
        sdir, s_veg, _ = self.calc_direct_sw(
            nsun,
            irradiance,
            method=method,
            resolution=resolution,
        )
        if lscatter:
            if vf is None or svf is None:
                raise ValueError("View factors are required for shortwave reflections")
            knet = self.calc_reflections_sw(sdir, dsky, vf, svf, albedo)
        else:
            if fss is None:
                raise ValueError("Fss is required for non-scattering shortwave")
            knet = (1.0 - albedo) * (sdir + dsky * fss)
        return sdir, knet, s_veg

    @staticmethod
    def _read_weather_table(path: Path) -> Dict[str, np.ndarray]:
        if not path.exists():
            raise FileNotFoundError(f"Weather file not found: {path}")

        with path.open("r", encoding="ascii", errors="ignore") as f:
            header = ""
            for line in f:
                if line.strip() and not line.strip().startswith("#"):
                    header = line.strip()
                    break
            if not header:
                raise ValueError(f"Weather file is empty: {path}")

            delimiter = "," if "," in header else None
            names = header.split(",") if delimiter == "," else header.split()
            names = [n.strip() for n in names]

            rows = []
            for line in f:
                if not line.strip() or line.strip().startswith("#"):
                    continue
                parts = line.split(",") if delimiter == "," else line.split()
                if len(parts) != len(names):
                    continue
                rows.append([float(p) for p in parts])

        if not rows:
            raise ValueError(f"Weather file contains no data rows: {path}")

        data = np.asarray(rows, dtype=float)
        table = {name: data[:, i] for i, name in enumerate(names)}
        table["date"] = data[:, 0].astype(int)
        table["time"] = data[:, 1].astype(int)
        return table

    @staticmethod
    def _interp_makima(x: np.ndarray, y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        x_new = np.asarray(x_new, dtype=float)

        if x.size < 2:
            raise ValueError("Need at least two points for interpolation")
        if np.any(np.diff(x) <= 0):
            raise ValueError("x must be strictly increasing for interpolation")

        m = np.diff(y) / np.diff(x)
        m1 = 2.0 * m[0] - m[1]
        m2 = 2.0 * m1 - m[0]
        m_end1 = 2.0 * m[-1] - m[-2]
        m_end2 = 2.0 * m_end1 - m[-1]
        m_ext = np.concatenate(([m2, m1], m, [m_end1, m_end2]))

        d = np.zeros_like(x)
        for i in range(x.size):
            w1 = abs(m_ext[i + 3] - m_ext[i + 2]) + abs(m_ext[i + 3] + m_ext[i + 2]) * 0.5
            w2 = abs(m_ext[i + 1] - m_ext[i]) + abs(m_ext[i + 1] + m_ext[i]) * 0.5
            if w1 + w2 > 0.0:
                d[i] = (w1 * m_ext[i + 1] + w2 * m_ext[i + 2]) / (w1 + w2)
            else:
                d[i] = 0.5 * (m_ext[i + 1] + m_ext[i + 2])

        idx = np.searchsorted(x, x_new, side="right") - 1
        idx = np.clip(idx, 0, x.size - 2)
        x0 = x[idx]
        x1 = x[idx + 1]
        y0 = y[idx]
        y1 = y[idx + 1]
        d0 = d[idx]
        d1 = d[idx + 1]
        h = x1 - x0
        t = (x_new - x0) / h

        t2 = t * t
        t3 = t2 * t
        h00 = 2.0 * t3 - 3.0 * t2 + 1.0
        h10 = t3 - 2.0 * t2 + t
        h01 = -2.0 * t3 + 3.0 * t2
        h11 = t3 - t2

        return h00 * y0 + h10 * h * d0 + h01 * y1 + h11 * h * d1

    @staticmethod
    def _write_sdir_nc(path: Path, tSP: np.ndarray, sdir: np.ndarray) -> None:
        try:
            from netCDF4 import Dataset
        except ImportError as exc:
            raise ImportError("netCDF4 is required to write Sdir.nc") from exc

        with Dataset(path, "w", format="NETCDF4") as ds:
            ds.createDimension("rows", sdir.shape[0])
            ds.createDimension("columns", sdir.shape[1])
            var_t = ds.createVariable("tSP", "f4", ("columns",))
            var_s = ds.createVariable("Sdir", "f4", ("rows", "columns"))
            var_t[:] = tSP.astype(np.float32)
            var_s[:, :] = sdir.astype(np.float32)


SPEC = SectionSpec(
    name="radiation",
    fields=FIELDS,
    defaults=DEFAULTS,
    section_cls=RadiationSection,
)
