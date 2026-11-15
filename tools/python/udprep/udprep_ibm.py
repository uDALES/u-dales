"""Immersed boundary method (IBM) preprocessing for uDALES.

Generates facet geometry files, runs the IBM preprocessor (Fortran
executable or f2py module) to compute solid/fluid boundary masks and
facet-section mappings, and writes the resulting sparse input files.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
from contextlib import contextmanager
from pathlib import Path
import shutil
import subprocess
from typing import Any, Dict, List
import warnings

import numpy as np

from ._section import Section, SectionSpec

IBM_FORTRAN_DIR = Path(__file__).resolve().parents[1] / "fortran" / "ibm_preproc"
PREBUILT_IBM_EXE = Path(__file__).resolve().parents[2] / "preprocessing" / "build" / "bin" / "IBM_preproc"

logger = logging.getLogger(__name__)

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("ibm", {})
FIELDS: List[str] = list(DEFAULTS.keys())


class IBMSection(Section):
    def run_all(self, backend: str = "f2py") -> None:
        """Run IBM preprocessing steps in the standard order."""
        sim = self._require_sim()
        self.sim.calculate_facet_sections_uvw = self.iwallmom > 1
        self.sim.calculate_facet_sections_c = self.ltempeq or self.lmoist or self.lwritefac
        steps = []

        if self.libm:
            if self.gen_geom:
                if getattr(self, "lfacetuniform", False):
                    steps.append(("uniformize_facets", self.uniformize_facets))
                steps.extend(
                    [
                        ("run_ibm", lambda: self.run_ibm(backend=backend)),
                        ("write_facets", self.write_facets),
                        ("write_facetarea", self.write_facetarea),
                    ]
                )
                if self.calculate_facet_sections_c:
                    steps.extend([
                        ("write_facets_unused", self.write_facets_unused),
                    ])

                factypes_path = Path(sim.path) / f"factypes.inp.{sim.expnr}"
                if not factypes_path.exists():
                    steps.extend([
                        ("generate_factypes", self.generate_factypes),
                        ("write_factypes", self.write_factypes),
                    ])
                else:
                    self._warn_not_overwriting(factypes_path)
            else:
                steps.extend(
                    [
                        ("copy_geom_outputs", self.copy_geom_outputs),
                    ]
                )
        self.run_steps("ibm", steps)

    @staticmethod
    def _warn_not_overwriting(path: Path) -> None:
        _orig_formatwarning = warnings.formatwarning
        warnings.formatwarning = lambda msg, cat, *_a, **_kw: f"{cat.__name__}: {msg}\n"
        try:
            warnings.warn(f"{path} already exists; NOT overwriting.", stacklevel=1)
        finally:
            warnings.formatwarning = _orig_formatwarning

    def generate_factypes(self) -> None:
        """Construct default factypes array (generate_factypes in MATLAB)."""
        sim = self._require_sim()
        k = sim.nfaclyrs
        factypes = []

        def add(row_id, lgr, z0, z0h, al, em, d_vals, c_val, l_val, k_val):
            row = [row_id, lgr, z0, z0h, al, em]
            row.extend(d_vals)
            row.extend([c_val] * k)
            row.extend([l_val] * k)
            row.extend([k_val] * (k + 1))
            factypes.append(row)

        add(-101, 0, 0.0, 0.0, 0.5, 0.85, [0.0] * k, 0.0, 0.0, 0.0)
        add(-1, 0, 0.05, 0.00035, 0.5, 0.85, [0.1, 0.2, 0.2] if k == 3 else [0.5 / k] * k, 1.875e6, 0.75, 0.4e-6)
        add(0, 0, 0.0, 0.0, 0.0, 0.0, [0.3 / k] * k, 1.875e6, 0.75, 0.4e-6)
        add(1, 0, 0.05, 0.00035, 0.5, 0.85, [0.36 / k] * k, 2.5e6, 1.0, 0.4e-6)
        add(2, 0, 0.05, 0.00035, 0.5, 0.85, [0.36 / k] * k, 2.766667e6, 0.83, 0.3e-6)
        add(3, 0, 0.05, 0.00035, 0.5, 0.85, [0.36 / k] * k, 2.19e6, 2.19, 1e-6)
        add(4, 0, 0.05, 0.00035, 0.5, 0.85, [0.36 / k] * k, 1e6, 0.1, 0.1e-6)
        add(11, 1, 0.05, 0.00035, 0.25, 0.95, [0.6 / k] * k, 5e6, 2.0, 0.4e-6)
        add(12, 1, 0.05, 0.00035, 0.35, 0.90, [0.6 / k] * k, 2e6, 0.8, 0.4e-6)
        sim.factypes = np.asarray(factypes, dtype=float)

    def write_factypes(self) -> None:
        """Write factypes.inp file."""
        sim = self._require_sim()
        if not hasattr(sim, "factypes"):
            self.generate_factypes()
        k = sim.nfaclyrs
        path = Path(sim.path) / f"factypes.inp.{sim.expnr}"
        dheader = "".join(f"  d{idx} [m]" for idx in range(1, k + 1))
        cheader = "".join(f"  C{idx} [J/(K m^3)]" for idx in range(1, k + 1))
        lheader = "".join(f"  l{idx} [W/(m K)]" for idx in range(1, k + 1))
        kheader = "".join(f"  k{idx} [W/(m K)]" for idx in range(1, k + 2))
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write(f"# walltype, {k} layers per type where layer 1 is the outdoor side and layer {k} is indoor side\n")
            f.write("# 0=default dummy, -1=asphalt floors;-101=concrete bounding walls;1=concrete;2=bricks;3=stone;4=painted wood;11=GR1; 12=GR2\n")
            f.write(f"# wallid  lGR  z0 [m]  z0h [m]  al [-]  em [-]{dheader}{cheader}{lheader}{kheader}\n")
            for row in self.factypes:
                fields = [
                    f"{int(row[0]):8d}",
                    f"{int(row[1]):3d}",
                    f"{row[2]:6.2f}",
                    f"{row[3]:7.5f}",
                    f"{row[4]:6.2f}",
                    f"{row[5]:6.2f}",
                ]
                offset = 6
                fields.extend(f"{val:6.2f}" for val in row[offset:offset + k])
                offset += k
                fields.extend(f"{val:14.0f}" for val in row[offset:offset + k])
                offset += k
                fields.extend(f"{val:13.4f}" for val in row[offset:offset + k])
                offset += k
                fields.extend(f"{val:13.8f}" for val in row[offset:offset + k + 1])
                f.write("  ".join(fields) + "\n")

    # ------------------------------------------------------------------ #
    # Facet-size uniformisation
    # ------------------------------------------------------------------ #

    @staticmethod
    def _uniform_names(stl_file: str, types_path: str | None) -> Dict[str, str]:
        """Output file names derived from the input STL (and types file) names."""
        stl = Path(stl_file)
        stem = f"{stl.stem}_uniform"
        names = {
            "stl": f"{stem}{stl.suffix or '.stl'}",
            "sidecar": f"{stem}.json",
            "facemap": f"{stem}_facemap.txt",
        }
        if types_path:
            tp = Path(types_path)
            names["types"] = f"{tp.stem}_uniform{tp.suffix or '.txt'}"
        return names

    @staticmethod
    def _file_digest(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                digest.update(chunk)
        return digest.hexdigest()

    def _uniform_settings(self) -> Dict[str, Any]:
        return {
            "facet_target_edge": float(self.facet_target_edge),
            "facet_size_band": float(self.facet_size_band),
            "facet_allow_split": bool(self.facet_allow_split),
            "facet_allow_merge": bool(self.facet_allow_merge),
        }

    def _resolve_types_path(self, sim) -> Path | None:
        if not self.read_types:
            return None
        if not str(self.types_path).strip():
            raise ValueError("types_path must be specified if read_types is true in namoptions")
        types_path = Path(self.types_path)
        if not types_path.is_absolute():
            types_path = Path(sim.path) / types_path
        if not types_path.is_file():
            raise FileNotFoundError(f"Missing types file: {types_path}")
        return types_path

    def uniformize_facets(self) -> None:
        """Make the STL facet sizes uniform and continue from the new STL.

        Runs :func:`udgeom.uniformize_facets` on ``stl_file`` (splitting facets
        above the accepted size band and merging coplanar facets of the same
        type below it), writes ``<stem>_uniform.stl``, the matching facet-type
        file when ``read_types`` is set, an old-to-new face map and a JSON
        sidecar recording the source and the settings. ``stl_file`` (and
        ``types_path``) are then pointed at the new files, in memory and in
        ``namoptions``, and the geometry is reloaded, so every later step (IBM,
        facets.inp, view factors, shortwave, SEB) works on the uniform mesh.
        The original STL is never modified.

        Reruns are safe: a ``stl_file`` that is itself a uniformised geometry
        made with the current settings is accepted as is; an existing
        ``<stem>_uniform.stl`` is reused when its sidecar matches the source
        file and the settings, and refused otherwise (uDALES never overwrites
        existing inputs; delete or rename the file to regenerate it).
        """
        from udgeom import FacetCountError, uniformize_facets

        sim = self._require_sim()
        case = Path(sim.path)
        stl_name = str(sim.stl_file)
        stl_path = case / stl_name
        if not stl_path.is_file():
            raise FileNotFoundError(f"STL file not found: {stl_path}")
        settings = self._uniform_settings()
        types_path = self._resolve_types_path(sim)

        # Is stl_file already a product of this step?
        own_sidecar = case / f"{Path(stl_name).stem}.json"
        if own_sidecar.is_file():
            meta = json.loads(own_sidecar.read_text(encoding="ascii"))
            if meta.get("tool") == "udgeom.uniformize_facets" and meta.get("outputs", {}).get("stl") == stl_name:
                if meta.get("settings") == settings:
                    logger.info(
                        "[ibm] %s is already a uniformised geometry (made from %s); nothing to do",
                        stl_name, meta.get("source", {}).get("stl"),
                    )
                    return
                raise ValueError(
                    f"{stl_name} is a uniformised geometry made from {meta.get('source', {}).get('stl')} "
                    f"with settings {meta.get('settings')}, but the current settings are {settings}. "
                    f"Set stl_file back to the source STL and delete {stl_name} (and {own_sidecar.name}) "
                    "to regenerate it with the new settings."
                )

        names = self._uniform_names(stl_name, str(types_path.name) if types_path is not None else None)
        out_stl = case / names["stl"]
        sidecar = case / names["sidecar"]
        source_meta = {"stl": stl_name, "sha256": self._file_digest(stl_path)}
        if types_path is not None:
            source_meta["types"] = str(self.types_path)
            source_meta["types_sha256"] = self._file_digest(types_path)

        if out_stl.exists():
            meta = json.loads(sidecar.read_text(encoding="ascii")) if sidecar.is_file() else None
            reusable = (
                meta is not None
                and meta.get("source") == source_meta
                and meta.get("settings") == settings
                and ("types" not in names or (case / names["types"]).is_file())
            )
            if not reusable:
                raise FileExistsError(
                    f"{out_stl} already exists but was not made from {stl_name} with the current "
                    f"settings (or its sidecar {sidecar.name} is missing). uDALES does not overwrite "
                    "existing inputs: delete or rename it to regenerate."
                )
            logger.info("[ibm] reusing existing uniformised geometry %s", out_stl.name)
        else:
            face_types = None
            if types_path is not None:
                face_types = self._load_facet_types(types_path, len(sim.geom.stl.faces))
            n_in = len(sim.geom.stl.faces)
            max_faces = int(self.facet_count_max) if int(self.facet_count_max) > 0 else 4 * n_in
            try:
                result = uniformize_facets(
                    sim.geom.stl,
                    face_types=face_types,
                    target_edge=float(self.facet_target_edge) if float(self.facet_target_edge) > 0 else None,
                    size_band=float(self.facet_size_band),
                    allow_split=bool(self.facet_allow_split),
                    allow_merge=bool(self.facet_allow_merge),
                    max_faces=max_faces,
                )
            except FacetCountError as exc:
                raise FacetCountError(
                    f"{exc} (namoptions: facet_target_edge, facet_count_max)"
                ) from exc
            logger.info("[ibm] %s", result.summary())
            result.mesh.export(str(out_stl))
            if face_types is not None:
                with (case / names["types"]).open("w", encoding="ascii", newline="\n") as f:
                    f.write("# facet types of %s, in face order\n" % names["stl"])
                    for value in result.face_types:
                        f.write(f"{int(value)}\n")
            with (case / names["facemap"]).open("w", encoding="ascii", newline="\n") as f:
                f.write(f"# source face index (1-based, in {stl_name}) of every face of {names['stl']}\n")
                for value in result.face_map:
                    f.write(f"{int(value) + 1}\n")
            meta = {
                "tool": "udgeom.uniformize_facets",
                "source": source_meta,
                "settings": settings,
                "outputs": {k: v for k, v in names.items() if k != "sidecar"},
                "report": {k: v for k, v in result.report.items() if not isinstance(v, list)},
            }
            sidecar.write_text(json.dumps(meta, indent=2) + "\n", encoding="ascii")

        # Continue from the uniform geometry: in memory, in namoptions, and reload.
        self.stl_file = names["stl"]
        self.save_param("stl_file", names["stl"])
        if types_path is not None:
            self.types_path = names["types"]
            self.save_param("types_path", names["types"])
        sim._load_geometry()

    def write_facets(self) -> None:
        """Write facets.inp from the STL and facet types.

        An existing facets.inp is treated as authored input and is NOT
        overwritten (mirrors the factypes.inp protection), so hand-assigned
        materials are never destroyed by a re-run. Its facet types are loaded so
        downstream steps still see the correct state. Delete facets.inp to force
        regeneration.
        """
        sim = self._require_sim()
        stl = sim.geom.stl
        path = Path(sim.path) / f"facets.inp.{sim.expnr}"

        if path.exists():
            self._warn_not_overwriting(path)
            sim.facet_types = np.atleast_1d(
                np.loadtxt(path, skiprows=1, usecols=0, dtype=int)
            ).reshape(-1)
            if len(sim.facet_types) != len(stl.faces):
                raise ValueError(
                    f"{path} has {len(sim.facet_types)} facets but the geometry "
                    f"{sim.stl_file} has {len(stl.faces)} faces; it was written for a "
                    f"different STL. Delete {path.name} to regenerate it."
                )
            sim.nfcts = len(stl.faces)
            return

        if self.read_types:
            if not str(self.types_path).strip():
                raise ValueError("types_path must be specified if read_types is true in namoptions")
            types_path = Path(self.types_path)
            if not types_path.is_absolute():
                types_path = Path(sim.path) / types_path
            if not types_path.is_file():
                raise FileNotFoundError(f"Missing types file: {types_path}")
            facet_types = self._load_facet_types(types_path, len(stl.faces))
        else:
            facet_types = np.ones(len(stl.faces), dtype=int)
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# type, normal\n")
            for typ, normal in zip(facet_types, np.asarray(stl.face_normals, dtype=float)):
                f.write(f"{int(typ):-4d} {normal[0]:-4.4f} {normal[1]:-4.4f} {normal[2]:-4.4f}\n")
        sim.nfcts = len(stl.faces)
        sim.facet_types = facet_types

    @staticmethod
    def _load_facet_types(types_path: Path, expected_count: int) -> np.ndarray:
        attempts: list[str] = []
        for skiprows in (0, 1):
            try:
                loaded = np.loadtxt(types_path, skiprows=skiprows)
            except ValueError as exc:
                attempts.append(f"skiprows={skiprows}: {exc}")
                continue
            values = np.atleast_1d(np.asarray(loaded, dtype=int)).reshape(-1)
            if len(values) == expected_count:
                return values
            attempts.append(f"skiprows={skiprows}: found {len(values)} value(s)")
        detail = "; ".join(attempts)
        raise ValueError(
            f"Facet type file {types_path} does not contain the expected "
            f"{expected_count} facet type value(s). Tried headerless and one-line-header reads: {detail}"
        )

    def write_facets_unused(self) -> None:
        """Write facets_unused.<expnr> for facets without c-grid facet sections."""
        sim = self._require_sim()
        if not self.calculate_facet_sections_c:
            return

        facet_sections_path = Path(sim.path) / "facet_sections_c.txt"
        if not facet_sections_path.is_file():
            raise FileNotFoundError(f"Missing {facet_sections_path}")

        if self.nfcts < 1:
            raise ValueError("Cannot write facets_unused without a positive facet count")

        try:
            facet_sections = np.loadtxt(facet_sections_path, comments="#")
        except ValueError:
            facet_sections = np.empty((0, 4), dtype=float)

        if facet_sections.size == 0:
            facets_used = np.array([], dtype=int)
        else:
            facet_sections = np.atleast_2d(facet_sections)
            facets_used = np.unique(facet_sections[:, 0].astype(int))

        facets_unused = np.setdiff1d(np.arange(1, self.nfcts + 1, dtype=int), facets_used)
        out_path = Path(sim.path) / f"facets_unused.{sim.expnr}"
        with out_path.open("w", encoding="ascii", newline="\n") as f:
            if facets_unused.size:
                for facet in facets_unused:
                    f.write(f"{int(facet)}\n")

    def write_facetarea(self) -> None:
        """Write facetarea.inp file (facet areas)."""
        sim = self._require_sim()
        path = Path(sim.path) / f"facetarea.inp.{sim.expnr}"
        areas = np.asarray(sim.geom.stl.area_faces, dtype=float)
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# area of facets\n")
            for area in areas:
                f.write(f"{area:.6f}\n")

    # def compute_boundary_cells(self) -> None:
    #     """Compute fluid/solid boundary points (IBM getBoundaryCells)."""

    # def match_facets_to_cells(self) -> None:
    #     """Match facets to fluid cells (IBM matchFacetsToCells)."""

    # def classify_solid_points(self) -> None:
    #     """Classify solid points using inpolyhedron/in_mypoly methods."""

    # def compute_facet_sections(self) -> None:
    #     """Compute facet sections for u/v/w/c grids."""

    # def compute_fluid_boundary(self) -> None:
    #     """Compute fluid boundary locations for IBM grids."""

    def copy_geom_outputs(self) -> None:
        sim = self._require_sim()
        if not self.geom_path:
            raise ValueError("Need to specify geom_path when gen_geom is false")
        geom_path = Path(self.geom_path)
        if not geom_path.is_absolute():
            geom_path = Path(sim.path) / geom_path
        for pattern in ("solid_*", "fluid_boundary_*"):
            sources = list(geom_path.glob(pattern))
            if not sources:
                raise FileNotFoundError(f"No files matching '{pattern}' found in {geom_path}")
            for source in sources:
                shutil.copy2(source, Path(sim.path) / source.name)
        if self.calculate_facet_sections_uvw:
            for pattern in ("facet_sections_u*", "facet_sections_v*", "facet_sections_w*"):
                sources = list(geom_path.glob(pattern))
                if not sources:
                    raise FileNotFoundError(f"No files matching '{pattern}' found in {geom_path}")
                for source in sources:
                    shutil.copy2(source, Path(sim.path) / source.name)
        if self.calculate_facet_sections_c:
            sources = list(geom_path.glob("facet_sections_c*"))
            if not sources:
                raise FileNotFoundError(f"No files matching 'facet_sections_c*' found in {geom_path}")
            for source in sources:
                shutil.copy2(source, Path(sim.path) / source.name)
            unused_sources = sorted(geom_path.glob("facets_unused.*"))
            if unused_sources:
                shutil.copy2(unused_sources[0], Path(sim.path) / f"facets_unused.{sim.expnr}")
        for base in ("facetarea.inp", "facets.inp", "factypes.inp"):
            # Sort so that, with several matching variants, the copied file is
            # deterministic rather than filesystem-glob order (P24).
            sources = sorted(geom_path.glob(f"{base}.*"))
            if not sources:
                raise FileNotFoundError(f"No files matching '{base}.*' found in {geom_path}")
            shutil.copy2(sources[0], Path(sim.path) / f"{base}.{sim.expnr}")
        self._update_counts_from_existing_outputs()

    def run_ibm(self, backend: str = "f2py") -> None:
        sim = self._require_sim()
        backend_key = backend.strip().lower()
        if backend_key == "f2py":
            counts = self._run_ibm_via_f2py()
            self._update_counts(counts)
        elif backend_key in {"legacy", "exe", "fortran_exe"}:
            self._write_ibm_input_files()
            self._run_ibm_via_legacy()
            self._update_counts_from_info_fort()
            # These generically-named intermediates are produced only by the
            # legacy path (_write_ibm_input_files / the legacy executable), so
            # only clean them up there. Deleting them unconditionally would
            # clobber unrelated user files of the same name on the f2py path,
            # which never creates them (P23).
            for name in (
                "inmypoly_inp_info.txt",
                "faces.txt",
                "vertices.txt",
                "zfgrid.txt",
                "zhgrid.txt",
                "info_fort.txt",
            ):
                path = Path(sim.path) / name
                if path.exists():
                    path.unlink()
        else:
            raise ValueError(f"Unknown IBM backend: {backend}")
        self._reload_ibm_outputs()

    def _reload_ibm_outputs(self) -> None:
        """Reload the solid masks and facet sections just written by ``run_ibm``.

        ``UDBase`` loads ``S{u,v,w,c}`` and ``facsec`` once at construction. On a
        first-ever preprocessing run those files do not yet exist, so ``sim``
        holds empty masks/sections. ``run_ibm`` writes the real ``solid_*.txt``
        and ``facet_sections_*.txt`` files but does not refresh ``sim``; the
        downstream radiation step would then trace an empty solid volume or raise
        "Facet sections not available". Re-invoke the existing UDBase loaders so
        the fresh outputs are in memory, for both backends and for callers that
        invoke ``prep.ibm.run_ibm()`` directly.
        """
        sim = self._require_sim()
        sim._load_solid_masks()
        sim._load_facet_sections()

    def _run_ibm_via_f2py(self) -> np.ndarray:
        try:
            from . import ibm_preproc_f2py as _ibm_mod
        except ImportError as exc:
            raise RuntimeError(
                "ibm_preproc_f2py module not available; build it with "
                "tools/build_preprocessing.sh."
            ) from exc

        sim = self._require_sim()
        stl = sim.geom.stl
        vertices = self._as_fortran_real_input(stl.vertices)
        facets = np.asfortranarray(np.asarray(stl.faces, dtype=np.int32) + 1)
        incenters = self._as_fortran_real_input(sim.geom.face_incenters)
        facenormals = self._as_fortran_real_input(sim.geom.face_normals)
        zt = self._as_fortran_real_input(np.asarray(sim.zt).reshape(-1)[: int(sim.ktot)])
        zm = self._as_fortran_real_input(np.asarray(sim.zm).reshape(-1)[: int(sim.ktot)])
        ray_dir_u = self._as_fortran_real_input(self._ray_direction("ray_dir_u"))
        ray_dir_v = self._as_fortran_real_input(self._ray_direction("ray_dir_v"))
        ray_dir_w = self._as_fortran_real_input(self._ray_direction("ray_dir_w"))
        ray_dir_c = self._as_fortran_real_input(self._ray_direction("ray_dir_c"))

        with self._pushd(Path(sim.path)):
            counts = _ibm_mod.run_ibm_preproc_f2py(
                vertices,
                facets,
                incenters,
                facenormals,
                self._as_fortran_real_scalar(sim.dx),
                self._as_fortran_real_scalar(sim.dy),
                sim.itot,
                sim.jtot,
                zt,
                zm,
                ray_dir_u,
                ray_dir_v,
                ray_dir_w,
                ray_dir_c,
                self.stl_ground,
                self.diag_neighbs,
                sim.BCxm == 1,
                sim.BCym == 1,
                self.nompthreads,
                self._as_fortran_real_scalar(self.ibmtol),
            )
        return np.asarray(counts, dtype=int).reshape(-1)

    @staticmethod
    def _as_fortran_real_input(values: Any) -> np.ndarray:
        # Match the old text-file path: Python wrote %.10f, then Fortran read default real.
        return np.asfortranarray(np.round(np.asarray(values, dtype=np.float64), 10).astype(np.float32))

    @staticmethod
    def _as_fortran_real_scalar(value: Any) -> float:
        return float(np.float32(round(float(value), 10)))

    def _ray_direction(self, name: str) -> np.ndarray:
        value = getattr(self, name)
        if isinstance(value, str):
            text = value.strip().strip("[]()")
            parts = [part.strip() for part in text.split(",") if part.strip()]
            if len(parts) == 1:
                parts = [part for part in text.split() if part]
            try:
                ray = np.asarray([float(part) for part in parts], dtype=np.float64)
            except ValueError as exc:
                raise ValueError(f"{name} must contain three numeric values") from exc
        else:
            ray = np.asarray(value, dtype=np.float64).reshape(-1)

        if ray.size != 3:
            raise ValueError(f"{name} must contain exactly three values, got {ray.size}")
        if not np.all(np.isfinite(ray)):
            raise ValueError(f"{name} must contain finite numeric values")
        if np.linalg.norm(ray) == 0.0:
            raise ValueError(f"{name} must be a non-zero vector")
        return ray

    def _run_ibm_via_legacy(self) -> None:
        sim = self._require_sim()
        exe = self._resolve_ibm_preprocessor()
        subprocess.run(
            [str(exe)],
            check=True,
            cwd=Path(sim.path),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

    def _resolve_ibm_preprocessor(self) -> Path:
        if PREBUILT_IBM_EXE.exists():
            return PREBUILT_IBM_EXE

        compiler = shutil.which("gfortran")
        if compiler is None:
            raise RuntimeError(
                "gfortran is required for the legacy IBM preprocessor, or build "
                "tools/preprocessing/build/bin/IBM_preproc with tools/build_preprocessing.sh."
            )

        build_dir = Path(self._require_sim().path) / ".ibm_build"
        build_dir.mkdir(parents=True, exist_ok=True)
        exe = build_dir / "IBM_preproc.exe"
        sources = [
            "in_mypoly_functions.f90",
            "boundaryMasking.f90",
            "matchFacetsCells.f90",
            "IBM_preproc_io.f90",
            "IBM_preproc_module.f90",
            "IBM_preproc_main.f90",
        ]
        cmd = [compiler, "-O3", "-fopenmp"]
        cmd.extend(str(IBM_FORTRAN_DIR / src) for src in sources)
        cmd.extend(["-o", str(exe)])
        subprocess.run(
            cmd,
            check=True,
            cwd=build_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        return exe

    @staticmethod
    @contextmanager
    def _pushd(path: Path):
        previous = Path.cwd()
        os.chdir(path)
        try:
            yield
        finally:
            os.chdir(previous)

    def _write_ibm_input_files(self) -> None:
        sim = self._require_sim()
        stl = sim.geom.stl
        path = Path(sim.path)
        with (path / "inmypoly_inp_info.txt").open("w", encoding="ascii", newline="\n") as f:
            f.write(f"{sim.dx:15.10f} {sim.dy:15.10f}\n")
            f.write(f"{sim.itot:5d} {sim.jtot:5d} {sim.ktot:5d}\n")
            f.write(f"{self.ibmtol:15.10f}\n")
            for name in ("ray_dir_u", "ray_dir_v", "ray_dir_w", "ray_dir_c"):
                ray = self._ray_direction(name)
                f.write(f"{ray[0]:15.10f} {ray[1]:15.10f} {ray[2]:15.10f}\n")
            f.write(f"{len(stl.vertices):8d} {len(stl.faces):8d}\n")
            f.write(f"{self.nompthreads:4d}\n")
            f.write(
                f"{int(self.stl_ground):d} "
                f"{int(self.diag_neighbs):d} "
                f"{int(sim.BCxm == 1):d} "
                f"{int(sim.BCym == 1):d}\n"
            )
        np.savetxt(path / "zhgrid.txt", sim.zm, fmt="%15.10f")
        np.savetxt(path / "zfgrid.txt", sim.zt, fmt="%15.10f")
        np.savetxt(path / "vertices.txt", np.asarray(stl.vertices, dtype=float), fmt="%15.10f %15.10f %15.10f")
        faces = np.asarray(stl.faces, dtype=int) + 1 # convert to 1-based indexing for Fortran
        incenters = np.asarray(sim.geom.face_incenters, dtype=float)
        facenormals = np.asarray(sim.geom.face_normals, dtype=float)
        face_rows = np.hstack(
            [faces, incenters, facenormals]
        )
        np.savetxt(
            path / "faces.txt",
            face_rows,
            fmt="%8d %8d %8d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f",
        )

    def _update_counts_from_info_fort(self) -> None:
        sim = self._require_sim()
        info = np.loadtxt(Path(sim.path) / "info_fort.txt", skiprows=1, dtype=int).reshape(-1)
        self._update_counts(info)

    def _update_counts(self, counts: np.ndarray) -> None:
        sim = self._require_sim()
        info = np.asarray(counts, dtype=int).reshape(-1)
        names = [
            "nfcts",
            "nsolpts_u",
            "nsolpts_v",
            "nsolpts_w",
            "nsolpts_c",
            "nbndpts_u",
            "nbndpts_v",
            "nbndpts_w",
            "nbndpts_c",
            "nfctsecs_u",
            "nfctsecs_v",
            "nfctsecs_w",
            "nfctsecs_c",
        ]
        if info.size != len(names):
            raise ValueError(f"Expected {len(names)} IBM counts, got {info.size}")
        for name, value in zip(names, info):
            self.save_param(name, int(value))
            setattr(sim, name, int(value))

    def _update_counts_from_existing_outputs(self) -> None:
        sim = self._require_sim()
        self.save_param("nfcts", int(len(sim.geom.stl.faces)))
        for name in ("solid_u", "solid_v", "solid_w", "solid_c", "fluid_boundary_u", "fluid_boundary_v", "fluid_boundary_w", "fluid_boundary_c", "facet_sections_u", "facet_sections_v", "facet_sections_w", "facet_sections_c"):
            path = Path(sim.path) / f"{name}.txt"
            count = 0
            if path.exists():
                with path.open("r", encoding="ascii", errors="ignore") as f:
                    count = sum(1 for line in f if line.strip() and not line.lstrip().startswith("#"))
            target = name.replace("solid_", "nsolpts_").replace("fluid_boundary_", "nbndpts_").replace("facet_sections_", "nfctsecs_")
            self.save_param(target, int(count))
            setattr(sim, target, int(count))

    def _require_sim(self):
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        if self.geom is None and self.libm:
            raise ValueError("Geometry must be loaded for IBM preprocessing")
        return self.sim


SPEC = SectionSpec(name="ibm", fields=FIELDS, defaults=DEFAULTS, section_cls=IBMSection)
