from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

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


def convert_block_to_sparse(sim_id: str, sim_path: Path | None = None) -> Path:
    """Convert trees.inp.<id> to vegetation.inp.<id> using sparse ijk format.

    Parameters
    ----------
    sim_id : str
        Experiment identifier (e.g. "525").
    sim_path : Path, required
        Directory containing namoptions.<id> and trees.inp.<id>.
    """

    if sim_path is None:
        raise ValueError("sim_path must be provided (experiment directory)")

    sim_dir = Path(sim_path)

    if not (sim_dir / f"namoptions.{sim_id}").is_file():
        raise FileNotFoundError(
            f"Missing namoptions.{sim_id} in expected sim dir {sim_dir}"
        )

    blocks = _read_tree_blocks(sim_dir, sim_id)
    points = _blocks_to_sparse_points(blocks)

    out_path = sim_dir / f"vegetation.inp.{sim_id}"
    _write_sparse_file(out_path, points)

    print(
        f"Loaded {len(blocks)} tree blocks from {sim_dir}, "
        f"expanded to {len(points)} grid points",
    )
    print(f"Sparse output written to {out_path}")
    return out_path


if __name__ == "__main__":
    convert_block_to_sparse()
