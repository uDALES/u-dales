from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple, TYPE_CHECKING
import warnings

from .udprep import Section, SectionSpec

if TYPE_CHECKING:
    from udbase import UDBase
import numpy as np

DEFAULTS: Dict[str, Any] = Section.load_defaults_json().get("vegetation", {})
FIELDS: List[str] = list(DEFAULTS.keys())

class VegetationSection(Section):
    def run_all(self) -> None:
        """Run vegetation preprocessing steps."""
        treefile = str(getattr(self, "treesfile", ""))
        use_stl = treefile.lower().endswith(".stl") if treefile else False
        steps = [
            ("generate_trees_from_namoptions", self.generate_trees_from_namoptions),
            ("write_trees", self.write_trees),
            ("load_stl", self.load_stl) if use_stl else ("load_block", self.load_block),
        ]
        self.run_steps("vegetation", steps)

    def save(self) -> None:
        """Write veg_* to file."""
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")
        if not hasattr(self, "veg") or self.veg is None:
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

        sim.save_veg(points + 1, ids, lad_values, cd, ud, dec, lsize, r_s)
        self.write_changed_params()

    def load_block(self, filename: str | Path | None = None) -> Path:
        """Convert trees.inp blocks to veg.inp/veg_params.inp."""
        sim = self.sim
        if sim is None:
            raise ValueError("UDBase instance must be provided")

        if filename is not None:
            self.treesfile = str(filename)
        treefile_name = getattr(self, "treesfile", "")

        # Parse tree blocks (1-based indices) from the input file.
        sim_dir = Path(sim.path)
        sim_id = sim.expnr
        if not treefile_name:
            warnings.warn("vegetation.treesfile is empty; using default trees.inp.<expnr>")
            tree_file = sim_dir / f"trees.inp.{sim_id}"
        else:
            tree_file = Path(treefile_name)
            if not tree_file.is_absolute():
                tree_file = sim_dir / tree_file
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
                    raise ValueError(f"Expected 6 or 7 integers, got {len(values)}")

                if il > iu:
                    il, iu = iu, il
                if jl > ju:
                    jl, ju = ju, jl
                if kl > ku:
                    kl, ku = ku, kl

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

        # Resolve vegetation parameters from section defaults (loaded via defaults.json).
        lad_value_resolved = float(self.lad)
        lad_source = "Vegetation section"
        cd = float(self.cd)
        ud = float(self.ud)
        dec = float(self.dec)
        lsize = float(self.lsize)
        r_s = float(self.r_s)

        # Build UDBase-style veg structure (0-based points) before writing.
        points_arr = np.asarray(points, dtype=int)
        ids_arr = np.asarray(point_ids, dtype=int)
        lad_vals = np.full(len(points), lad_value_resolved, dtype=float)
        veg = {
            "points": points_arr - 1,
            "params": {
                "id": ids_arr,
                "lad": lad_vals,
                "cd": np.full(len(points), cd, dtype=float),
                "ud": np.full(len(points), ud, dtype=float),
                "dec": np.full(len(points), dec, dtype=float),
                "lsize": np.full(len(points), lsize, dtype=float),
                "r_s": np.full(len(points), r_s, dtype=float),
            },
        }

        self.ntrees = len(points)
        self.ltrees = True
        self.veg = veg
        if getattr(sim, "veg", None) is not veg:
            warnings.warn("vegetation data is not saved to disk; call the save() method to persist veg inputs")

        fig = sim.plot_veg(veg, show=False)

        return fig

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
        ids = np.ones(len(points), dtype=int)
        lad_value_resolved = float(self.lad)
        lad_vals = np.full(len(points), lad_value_resolved, dtype=float)
        cd = float(self.cd)
        ud = float(self.ud)
        dec = float(self.dec)
        lsize = float(self.lsize)
        r_s = float(self.r_s)

        veg = {
            "points": points - 1,
            "params": {
                "id": ids,
                "lad": lad_vals,
                "cd": np.full(len(points), cd, dtype=float),
                "ud": np.full(len(points), ud, dtype=float),
                "dec": np.full(len(points), dec, dtype=float),
                "lsize": np.full(len(points), lsize, dtype=float),
                "r_s": np.full(len(points), r_s, dtype=float),
            },
        }

        self.ltrees = True
        self.veg = veg
        if getattr(sim, "veg", None) is not veg:
            warnings.warn("vegetation data is not saved to disk; call vegetation.save() to persist veg inputs")

        fig = sim.plot_veg(veg, show=False)
        return fig


SPEC = SectionSpec(
    name="vegetation",
    fields=FIELDS,
    defaults=DEFAULTS,
    section_cls=VegetationSection,
)
