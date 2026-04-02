from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
from typing import Any, Dict, List

import numpy as np

from .udprep import Section, SectionSpec

IBM_FORTRAN_DIR = Path(__file__).resolve().parents[2] / "IBM" / "IBM_preproc_fortran"


DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("ibm", {})
FIELDS: List[str] = list(DEFAULTS.keys())


class IBMSection(Section):
    def run_all(self) -> None:
        """Run IBM preprocessing steps in the standard order."""
        sim = self._require_sim()
        steps = [("generate_lscale", self.generate_lscale), ("write_lscale", self.write_lscale)]
        factypes_path = Path(sim.path) / f"factypes.inp.{sim.expnr}"
        if not factypes_path.exists():
            steps.extend([("generate_factypes", self.generate_factypes), ("write_factypes", self.write_factypes)])

        if bool(getattr(self, "libm", True)):
            if bool(getattr(self, "gen_geom", True)):
                steps.extend(
                    [
                        ("run_ibm_fortran", self.run_ibm_fortran),
                        ("write_facets", self.write_facets),
                        ("write_facetarea", self.write_facetarea),
                    ]
                )
            else:
                steps.extend(
                    [
                        ("copy_geom_outputs", self.copy_geom_outputs),
                        ("write_facets", self.write_facets),
                        ("write_facetarea", self.write_facetarea),
                    ]
                )
        self.run_steps("ibm", steps)

    def generate_lscale(self) -> None:
        """Compute the length scale input (generate_lscale in MATLAB)."""
        sim = self._require_sim()
        zf = np.asarray(sim.zf, dtype=float)
        ls = np.zeros((len(zf), 10), dtype=float)
        ls[:, 0] = zf
        ls[:, 5] = float(getattr(sim, "w_s", 0.0))
        ls[:, 9] = float(getattr(sim, "R", 0.0))

        forcing_flags = sum(
            bool(getattr(sim, name, False))
            for name in ("luoutflowr", "lvoutflowr", "luvolflowr", "lvvolflowr", "lprofforc", "lcoriol")
        )
        ldp = bool(getattr(sim, "ldp", forcing_flags == 0))
        if forcing_flags + int(ldp) > 1:
            raise ValueError("More than one forcing type specified")
        if bool(getattr(sim, "lprofforc", False)) or bool(getattr(sim, "lcoriol", False)):
            ls[:, 1] = float(getattr(sim, "u0", 0.0))
            ls[:, 2] = float(getattr(sim, "v0", 0.0))
        elif ldp:
            ls[:, 3] = float(getattr(sim, "dpdx", 0.0))
            ls[:, 4] = float(getattr(sim, "dpdy", 0.0))
        self.ls = ls

    def write_lscale(self) -> None:
        """Write the length scale input to disk."""
        sim = self._require_sim()
        if not hasattr(self, "ls"):
            self.generate_lscale()
        path = Path(sim.path) / f"lscale.inp.{sim.expnr}"
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# SDBL flow \n")
            f.write("# z uq vq pqx pqy wfls dqtdxls dqtdyls dqtdtls dthlrad      \n")
            for row in self.ls:
                f.write(
                    f"{row[0]:-20.15f} {row[1]:-12.6f} {row[2]:-12.6f} {row[3]:-12.9f} "
                    f"{row[4]:-12.9f} {row[5]:-15.9f} {row[6]:-12.6f} {row[7]:-12.6f} "
                    f"{row[8]:-12.6f} {row[9]:-17.12f}\n"
                )

    def generate_factypes(self) -> None:
        """Construct default factypes array (generate_factypes in MATLAB)."""
        sim = self._require_sim()
        k = int(getattr(sim, "nfaclyrs", 3))
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
        self.factypes = np.asarray(factypes, dtype=float)

    def write_factypes(self) -> None:
        """Write factypes.inp file."""
        sim = self._require_sim()
        if not hasattr(self, "factypes"):
            self.generate_factypes()
        k = int(getattr(sim, "nfaclyrs", 3))
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

    def write_facets(self) -> None:
        """Write facets.inp file from STL and facet types."""
        sim = self._require_sim()
        stl = sim.geom.stl
        path = Path(sim.path) / f"facets.inp.{sim.expnr}"
        if bool(getattr(sim, "read_types", False)):
            types_path = Path(getattr(sim, "types_path"))
            facet_types = np.loadtxt(types_path, skiprows=1 if types_path.suffix == ".txt" else 0)
            facet_types = np.asarray(facet_types, dtype=int).reshape(-1)
        else:
            facet_types = np.ones(len(stl.faces), dtype=int)
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# type, normal\n")
            for typ, normal in zip(facet_types, np.asarray(stl.face_normals, dtype=float)):
                f.write(f"{int(typ):-4d} {normal[0]:-4.4f} {normal[1]:-4.4f} {normal[2]:-4.4f}\n")
        sim.nfcts = len(stl.faces)
        self.facet_types = facet_types

    def write_facetarea(self) -> None:
        """Write facetarea.inp file (facet areas)."""
        sim = self._require_sim()
        path = Path(sim.path) / f"facetarea.inp.{sim.expnr}"
        areas = np.asarray(sim.geom.stl.area_faces, dtype=float)
        with path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# area of facets\n")
            for area in areas:
                f.write(f"{area:.6f}\n")

    def compute_boundary_cells(self) -> None:
        """Compute fluid/solid boundary points (IBM getBoundaryCells)."""

    def match_facets_to_cells(self) -> None:
        """Match facets to fluid cells (IBM matchFacetsToCells)."""

    def classify_solid_points(self) -> None:
        """Classify solid points using inpolyhedron/in_mypoly methods."""

    def compute_facet_sections(self) -> None:
        """Compute facet sections for u/v/w/c grids."""

    def compute_fluid_boundary(self) -> None:
        """Compute fluid boundary locations for IBM grids."""

    def copy_geom_outputs(self) -> None:
        sim = self._require_sim()
        geom_value = str(getattr(sim, "geom_path", "")).strip()
        if not geom_value:
            raise ValueError("Need to specify geom_path when gen_geom is false")
        geom_path = Path(geom_value)
        if not geom_path.is_absolute():
            geom_path = Path(sim.path) / geom_path
        for pattern in ("solid_*", "fluid_boundary_*", "facet_sections_*"):
            for source in geom_path.glob(pattern):
                shutil.copy2(source, Path(sim.path) / source.name)
        self._update_counts_from_existing_outputs()

    def run_ibm_fortran(self) -> None:
        sim = self._require_sim()
        compiler = shutil.which("gfortran")
        if compiler is None:
            raise RuntimeError("gfortran is required for IBM preprocessing")
        build_dir = Path(sim.path) / ".ibm_build"
        build_dir.mkdir(parents=True, exist_ok=True)
        exe = build_dir / "IBM_preproc.exe"

        sources = [
            "in_mypoly_functions.f90",
            "boundaryMasking.f90",
            "matchFacetsCells.f90",
            "IBM_preproc_io.f90",
            "IBM_preproc_main.f90",
        ]
        cmd = [compiler, "-O3", "-fopenmp"]
        cmd.extend(str(IBM_FORTRAN_DIR / src) for src in sources)
        cmd.extend(["-o", str(exe)])
        subprocess.run(cmd, check=True, cwd=build_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

        self._write_ibm_input_files()
        subprocess.run([str(exe)], check=True, cwd=Path(sim.path), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        self._update_counts_from_info_fort()
        for name in ("inmypoly_inp_info.txt", "faces.txt", "vertices.txt", "zfgrid.txt", "zhgrid.txt", "info_fort.txt"):
            path = Path(sim.path) / name
            if path.exists():
                path.unlink()

    def _write_ibm_input_files(self) -> None:
        sim = self._require_sim()
        stl = sim.geom.stl
        path = Path(sim.path)
        dx = float(sim.dx)
        dy = float(sim.dy)
        zf = np.asarray(sim.zt, dtype=float)
        zh = np.asarray(sim.zm, dtype=float)
        with (path / "inmypoly_inp_info.txt").open("w", encoding="ascii", newline="\n") as f:
            f.write(f"{dx:15.10f} {dy:15.10f}\n")
            f.write(f"{int(sim.itot):5d} {int(sim.jtot):5d} {int(sim.ktot):5d}\n")
            f.write(f"{5e-4:15.10f}\n")
            for _ in range(4):
                f.write(f"{0.0:15.10f} {0.0:15.10f} {1.0:15.10f}\n")
            f.write(f"{len(stl.vertices):8d} {len(stl.faces):8d}\n")
            f.write(f"{8:4d}\n")
            f.write(
                f"{int(bool(getattr(sim, 'stl_ground', True))):d} "
                f"{int(bool(getattr(sim, 'diag_neighbs', True))):d} "
                f"{int(int(getattr(sim, 'BCxm', 1)) == 1):d} "
                f"{int(int(getattr(sim, 'BCym', 1)) == 1):d}\n"
            )
        np.savetxt(path / "zhgrid.txt", zh, fmt="%15.10f")
        np.savetxt(path / "zfgrid.txt", zf, fmt="%15.10f")
        np.savetxt(path / "vertices.txt", np.asarray(stl.vertices, dtype=float), fmt="%15.10f %15.10f %15.10f")
        faces = np.asarray(stl.faces, dtype=int) + 1
        vertices = np.asarray(stl.vertices, dtype=float)
        incenters = np.asarray(sim.geom.face_incenters, dtype=float)
        face_rows = np.hstack(
            [faces, incenters, np.asarray(stl.face_normals, dtype=float)]
        )
        np.savetxt(
            path / "faces.txt",
            face_rows,
            fmt="%8d %8d %8d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f",
        )

    def _update_counts_from_info_fort(self) -> None:
        sim = self._require_sim()
        info = np.loadtxt(Path(sim.path) / "info_fort.txt", skiprows=1, dtype=int).reshape(-1)
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
        for name, value in zip(names, info):
            sim.save_param(name, int(value))
            setattr(sim, name, int(value))

    def _update_counts_from_existing_outputs(self) -> None:
        sim = self._require_sim()
        sim.save_param("nfcts", int(len(sim.geom.stl.faces)))
        for name in ("solid_u", "solid_v", "solid_w", "solid_c", "fluid_boundary_u", "fluid_boundary_v", "fluid_boundary_w", "fluid_boundary_c", "facet_sections_u", "facet_sections_v", "facet_sections_w", "facet_sections_c"):
            path = Path(sim.path) / f"{name}.txt"
            count = 0
            if path.exists():
                with path.open("r", encoding="ascii", errors="ignore") as f:
                    count = sum(1 for line in f if line.strip() and not line.lstrip().startswith("#"))
            target = name.replace("solid_", "nsolpts_").replace("fluid_boundary_", "nbndpts_").replace("facet_sections_", "nfctsecs_")
            sim.save_param(target, int(count))
            setattr(sim, target, int(count))

    def _require_sim(self):
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        if getattr(self.sim, "geom", None) is None and bool(getattr(self.sim, "libm", True)):
            raise ValueError("Geometry must be loaded for IBM preprocessing")
        return self.sim


SPEC = SectionSpec(name="ibm", fields=FIELDS, defaults=DEFAULTS, section_cls=IBMSection)
