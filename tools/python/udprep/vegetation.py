from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

try:
    from ..udbase import UDBase
except ImportError:
    try:
        from udbase import UDBase  # type: ignore
    except ImportError:
        UDBase = None

def _iter_block_points(il: int, iu: int, jl: int, ju: int, kl: int, ku: int) -> Iterable[Tuple[int, int, int]]:
    """Yield (i,j,k) points inside a block defined by inclusive bounds."""

    for i in range(il, iu + 1):
        for j in range(jl, ju + 1):
            for k in range(kl, ku + 1):
                yield (i, j, k)


def _read_tree_blocks(sim_dir: Path, sim_id: str) -> List[Tuple[int, int, int, int, int, int]]:
    """Read tree blocks from trees.inp.<id> and return (il,iu,jl,ju,kl,ku) tuples."""

    tree_file = sim_dir / f"trees.inp.{sim_id}"
    if not tree_file.is_file():
        raise FileNotFoundError(f"Missing trees file: {tree_file}")

    blocks: List[Tuple[int, int, int, int, int, int]] = []
    with tree_file.open("r", encoding="ascii") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            values = [int(v) for v in stripped.split()]
            if len(values) == 6:
                il, iu, jl, ju, kl, ku = values
            elif len(values) == 7:
                il, iu, jl, ju, kl, ku = values[1:]
            else:
                raise ValueError(f"Expected 6 or 7 integers, got {len(values)}")

            blocks.append((il, iu, jl, ju, kl, ku))

    if not blocks:
        raise ValueError(f"No tree entries found in {tree_file}")

    return blocks


def _blocks_to_sparse_points(blocks: Sequence[Tuple[int, int, int, int, int, int]]) -> List[Tuple[int, int, int]]:
    points: List[Tuple[int, int, int]] = []
    for il, iu, jl, ju, kl, ku in blocks:
        points.extend(_iter_block_points(il, iu, jl, ju, kl, ku))
    return points


def _write_sparse_file(out_path: Path, points: Sequence[Tuple[int, int, int]]) -> None:
    with out_path.open("w", encoding="ascii", newline="\n") as f:
        f.write("# position (i,j,k)\n")
        for i, j, k in points:
            f.write(f"{i:7d} {j:7d} {k:7d}\n")


def _write_lad_file(out_path: Path, ntree_max: int, lad_value: float) -> None:
    """Write a simple vertical LAD profile to veg_lad.inp.<id>.

    The profile is uniform and matches the Fortran expectation of one value per
    canopy face level; length is ntree_max+1.
    """

    with out_path.open("w", encoding="ascii", newline="\n") as f:
        f.write("# lad profile at canopy faces (uniform)\n")
        for _ in range(ntree_max + 1):
            f.write(f"{lad_value:.6f}\n")


def _resolve_lad_value(sim: "UDBase") -> tuple[float, str]:
    """Resolve LAD value from UDBase (namoptions lad)."""

    if sim is None:
        raise ValueError("UDBase instance must be provided")
    if not hasattr(sim, "lad"):
        raise AttributeError("UDBase loaded but missing lad in namoptions")

    return float(sim.lad), "UDBase namoptions"


def convert_block_to_sparse(sim: "UDBase") -> Path:
    """Convert trees.inp.<id> to veg.inp.<id> and emit a LAD profile.

    Parameters
    ----------
    sim : UDBase
        Populated UDBase instance; expnr and path are used to locate inputs.
    """

    if sim is None:
        raise ValueError("sim (UDBase instance) must be provided")

    sim_dir = Path(sim.path)
    sim_id = sim.expnr

    if not (sim_dir / f"namoptions.{sim_id}").is_file():
        raise FileNotFoundError(
            f"Missing namoptions.{sim_id} in expected sim dir {sim_dir}"
        )

    blocks = _read_tree_blocks(sim_dir, sim_id)
    lad_value_resolved, lad_source = _resolve_lad_value(sim)
    points = _blocks_to_sparse_points(blocks)

    out_path = sim_dir / f"veg.inp.{sim_id}"
    _write_sparse_file(out_path, points)

    # LAD profile (uniform) sized by canopy height
    k_min = min(block[4] for block in blocks)
    k_max = max(block[5] for block in blocks)
    ntree_max = k_max - k_min + 1
    lad_path = sim_dir / f"veg_lad.inp.{sim_id}"
    _write_lad_file(lad_path, ntree_max, lad_value_resolved)

    print(
        f"Loaded {len(blocks)} tree blocks from {sim_dir}, "
        f"expanded to {len(points)} grid points",
    )
    print(f"Sparse output written to {out_path}")
    print(
        f"LAD profile (uniform {lad_value_resolved} from {lad_source}) "
        f"written to {lad_path}"
    )
    return out_path


if __name__ == "__main__":
    raise SystemExit("Run via veg.py which constructs UDBase and passes it here")
