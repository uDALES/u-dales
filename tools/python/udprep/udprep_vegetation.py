"""Vegetation preprocessing for uDALES.

Converts block-style tree definitions (from namoptions or STL) into
the sparse per-cell vegetation input files (veg.inp, veg_params.inp)
consumed by the solver's vegetation module.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple, TYPE_CHECKING
import warnings

from ._section import Section, SectionSpec

if TYPE_CHECKING:
    from udbase import UDBase
import numpy as np

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("vegetation", {})
FIELDS: List[str] = list(DEFAULTS.keys())

class VegetationSection(Section):
    def run_all(self) -> None:
        """Run vegetation preprocessing steps."""
        use_stl = self.treesfile.lower().endswith(".stl") if self.treesfile else False
        steps = [("load_trees_stl", self.load_stl) if use_stl else ("load_trees_block", self.load_block)]
        self.run_steps("vegetation", steps)
        self.save()

    def save(self) -> None:
        """Write veg_* to file."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        
        if not hasattr(self.sim, "veg") or self.veg is None:
            raise ValueError("No vegetation data available; run load_block/load_stl first")

        veg = self.veg
        points = np.asarray(veg.get("points", []), dtype=int)
        if points.ndim == 1 and points.size:
            points = points.reshape(1, -1)
        params = veg.get("params", {})
        ids = np.asarray(params.get("id", []), dtype=int)
        lad_values = np.asarray(params.get("lad", []), dtype=float)
        if points.size == 0 or ids.size == 0 or lad_values.size == 0:
            raise ValueError("veg data missing points/ids/lad values")
        if len(points) != len(ids) or len(points) != len(lad_values):
            raise ValueError("veg points/ids/lad lengths do not match")

        def _coerce_scalar(name: str) -> float:
            values = np.asarray(params.get(name, []), dtype=float)
            if values.size == 0:
                raise ValueError(f"veg params missing '{name}'")
            if not np.allclose(values, values.flat[0]):
                warnings.warn(
                    f"veg param '{name}' varies per point; using the first value when writing.",
                    RuntimeWarning,
                )
            return float(values.flat[0])

        cd = _coerce_scalar("cd")
        ud = _coerce_scalar("ud")
        dec = _coerce_scalar("dec")
        lsize = _coerce_scalar("lsize")
        r_s = _coerce_scalar("r_s")

        # Write veg.inp.<expnr> (1-based indices)
        points_1based = points + 1

        veg_path = Path(self.path) / f"veg.inp.{self.expnr}"
        with veg_path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# position (i,j,k)\n")
            for i, j, k in points_1based:
                f.write(f"{int(i):7d} {int(j):7d} {int(k):7d}\n")

        params_path = Path(self.path) / f"veg_params.inp.{self.expnr}"
        with params_path.open("w", encoding="ascii", newline="\n") as f:
            f.write("# id lad cd ud dec lsize r_s\n")
            for bid, lad_val in zip(ids, lad_values):
                f.write(
                    f"{int(bid):7d} {float(lad_val):12.6f} {cd:12.6f} {ud:12.6f} "
                    f"{dec:12.6f} {lsize:12.6f} {r_s:12.6f}\n"
                )

    def vegetation_block_to_veg(self, filename: str | Path | None = None):
        """Compatibility entry point: ingest trees.inp and write sparse vegetation files."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")
        
        self.load_block(filename=filename)
        self.save()
        
        return {
            "veg": Path(self.path) / f"veg.inp.{self.expnr}",
            "params": Path(self.path) / f"veg_params.inp.{self.expnr}",
        }

    def load_block(self, filename: str | Path | None = None) -> None:
        """Convert trees.inp blocks to veg.inp/veg_params.inp."""
        if self.sim is None:
            raise ValueError("UDBase instance must be provided")

        # Parse tree blocks (1-based indices) from the input file.
        resolved_filename = filename or self.treesfile
        if not resolved_filename:
            warnings.warn("treesfile is not defined in the namoptions; using default naming convention: trees.inp.<expnr>")
            tree_file = self.path / f"trees.inp.{self.expnr}"
        else:
            tree_file = Path(resolved_filename)
            if not tree_file.is_absolute():
                tree_file = self.path / tree_file
        
        if not tree_file.is_file():
            raise FileNotFoundError(f"Missing trees file: {tree_file}")

        blocks: List[Tuple[int, int, int, int, int, int, int]] = []
        next_id = 1
        with tree_file.open("r", encoding="ascii") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                values = [int(v) for v in stripped.split()]
                if len(values) == 6:
                    il, iu, jl, ju, kl, ku = values
                elif len(values) == 7:
                    _, il, iu, jl, ju, kl, ku = values
                else:
                    raise ValueError(f"Expected 6 or 7 integers, got {len(values)} in {tree_file}")

                if il > iu:
                    il, iu = iu, il
                if jl > ju:
                    jl, ju = ju, jl
                if kl > ku:
                    kl, ku = ku, kl

                if il < 1 or iu > self.itot or jl < 1 or ju > self.jtot or kl < 1 or ku > self.ktot:
                    raise ValueError(
                        f"Tree block ({il}, {iu}, {jl}, {ju}, {kl}, {ku}) is outside "
                        f"the simulation grid 1..{self.itot}, 1..{self.jtot}, 1..{self.ktot}"
                    )

                block_id = next_id
                next_id += 1
                blocks.append((block_id, il, iu, jl, ju, kl, ku))

        if not blocks:
            raise ValueError(f"No tree entries found in {tree_file}")

        # Expand blocks into sparse point lists (1-based indices).
        points: List[Tuple[int, int, int]] = []
        point_ids: List[int] = []
        for block_id, il, iu, jl, ju, kl, ku in blocks:
            for i in range(il, iu + 1):
                for j in range(jl, ju + 1):
                    for k in range(kl, ku + 1):
                        points.append((i, j, k))
                        point_ids.append(block_id)

        # Build UDBase-style veg structure (0-based points) before writing.
        points_arr = np.asarray(points, dtype=int)
        ids_arr = np.asarray(point_ids, dtype=int)
        veg = {
            "points": points_arr - 1,
            "params": {
                "id": ids_arr,
                "cd": np.full(len(points), self.cd, dtype=float),
                "lad": np.full(len(points), self.lad, dtype=float),
                "lsize": np.full(len(points), self.lsize, dtype=float),
                "ud": np.full(len(points), self.ud, dtype=float),
                "dec": np.full(len(points), self.dec, dtype=float),
                "r_s": np.full(len(points), self.r_s, dtype=float),
            },
        }
 
        if self.ntrees != len(points):
            self.ntrees = len(points)
        self.sim.veg = veg
        self.save_param("ntrees", len(points))

    def load_stl(self, filename: str | Path | None = None) -> Dict[str, np.ndarray]:
        """Convert a closed STL volume to sparse vegetation points."""
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")

        if filename is not None:
            self.treesfile = str(filename)
        elif getattr(self, "treesfile", "") != "":
            filename = getattr(self, "treesfile", "")
        else:
            warnings.warn("vegetation.treesfile is empty; STL filename is required")
            raise ValueError("Missing STL filename for load_stl")

        stl_path = Path(filename)
        if not stl_path.is_absolute():
            stl_path = Path(sim.path) / stl_path

        try:
            import trimesh
        except ImportError as exc:
            raise ImportError("trimesh is required for load_stl") from exc

        if not stl_path.is_file():
            raise FileNotFoundError(f"Missing STL file: {stl_path}")

        mesh = trimesh.load_mesh(str(stl_path))
        if isinstance(mesh, trimesh.Scene):
            mesh = trimesh.util.concatenate(
                tuple(trimesh.Trimesh(vertices=g.vertices, faces=g.faces) for g in mesh.geometry.values())
            )
        if not mesh.is_watertight:
            raise ValueError("STL mesh must be watertight (closed) for point-in-volume checks.")

        # Use voxelization to avoid rtree dependency for point-in-volume queries.
        pitch = float(min(sim.dx, sim.dy, float(np.min(sim.dzt))))
        vox = mesh.voxelized(pitch).fill()
        bounds = vox.bounds
        x_mask = (sim.xt >= bounds[0, 0]) & (sim.xt <= bounds[1, 0])
        y_mask = (sim.yt >= bounds[0, 1]) & (sim.yt <= bounds[1, 1])
        z_mask = (sim.zt >= bounds[0, 2]) & (sim.zt <= bounds[1, 2])
        x_idx = np.where(x_mask)[0]
        y_idx = np.where(y_mask)[0]
        z_idx = np.where(z_mask)[0]

        if x_idx.size == 0 or y_idx.size == 0 or z_idx.size == 0:
            raise ValueError("STL bounds do not overlap the simulation grid.")

        xx, yy = np.meshgrid(sim.xt[x_idx], sim.yt[y_idx], indexing="ij")
        xy = np.column_stack([xx.ravel(), yy.ravel()])
        n_xy = xy.shape[0]

        points_idx: List[Tuple[int, int, int]] = []
        for k in z_idx:
            z = sim.zt[k]
            pts = np.column_stack([xy, np.full(n_xy, z)])
            inside = vox.is_filled(pts)
            if np.any(inside):
                idx = np.flatnonzero(inside)
                for flat in idx:
                    i = int(x_idx[flat // len(y_idx)])
                    j = int(y_idx[flat % len(y_idx)])
                    points_idx.append((i + 1, j + 1, k + 1))

        if not points_idx:
            raise ValueError("No vegetation points found inside the STL volume.")

        points = np.asarray(points_idx, dtype=int)

        veg = {
            "points": points - 1,
            "params": {
                "id": np.ones(len(points), dtype=int),
                "cd": np.full(len(points), self.cd, dtype=float),
                "lad": np.full(len(points), self.lad, dtype=float),
                "lsize": np.full(len(points), self.lsize, dtype=float),
                "ud": np.full(len(points), self.ud, dtype=float),
                "dec": np.full(len(points), self.dec, dtype=float),
                "r_s": np.full(len(points), self.r_s, dtype=float),
            },
        }

        if self.ntrees != len(points):
            self.ntrees = len(points)
        self.sim.veg = veg
        self.save_param("ntrees", len(points))

        fig = sim.vis.plot_veg(veg, show=False)
        return fig


def vegetation_block_to_veg(prep_or_sim, filename: str | Path | None = None):
    """Compatibility helper for legacy trees.inp -> sparse vegetation conversion."""
    from .udprep import UDPrep

    prep = prep_or_sim if isinstance(prep_or_sim, UDPrep) else UDPrep(prep_or_sim)
    return prep.vegetation.vegetation_block_to_veg(filename=filename)


SPEC = SectionSpec(
    name="vegetation",
    fields=FIELDS,
    defaults=DEFAULTS,
    section_cls=VegetationSection,
)
