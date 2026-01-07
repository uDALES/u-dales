from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

from udbase import UDBase
from udprep.namelist import update_namelist_value

def _iter_block_points(il: int, iu: int, jl: int, ju: int, kl: int, ku: int) -> Iterable[Tuple[int, int, int]]:
    """Yield (i,j,k) points inside a block defined by inclusive bounds."""

    for i in range(il, iu + 1):
        for j in range(jl, ju + 1):
            for k in range(kl, ku + 1):
                yield (i, j, k)


def _read_tree_blocks(sim_dir: Path, sim_id: str) -> List[Tuple[int, int, int, int, int, int, int]]:
    """Read tree blocks from trees.inp.<id>, auto-numbered in file order.

    Returns tuples of (id, il, iu, jl, ju, kl, ku). Any id column in the
    source is ignored; blocks are numbered 1..N in the order encountered.
    """

    tree_file = sim_dir / f"trees.inp.{sim_id}"
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
                _, il, iu, jl, ju, kl, ku = values  # ignore provided id
            else:
                raise ValueError(f"Expected 6 or 7 integers, got {len(values)}")

            # Normalize bounds in case input orders are inverted
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

    return blocks


def _blocks_to_sparse_points(blocks: Sequence[Tuple[int, int, int, int, int, int, int]]) -> List[Tuple[int, int, int]]:
    points: List[Tuple[int, int, int]] = []
    for _, il, iu, jl, ju, kl, ku in blocks:
        points.extend(_iter_block_points(il, iu, jl, ju, kl, ku))
    return points


def _blocks_to_id_list(blocks: Sequence[Tuple[int, int, int, int, int, int, int]]) -> List[int]:
    ids: List[int] = []
    for block_id, il, iu, jl, ju, kl, ku in blocks:
        for _ in _iter_block_points(il, iu, jl, ju, kl, ku):
            ids.append(block_id)
    return ids


def _blocks_to_k_list(blocks: Sequence[Tuple[int, int, int, int, int, int, int]]) -> List[int]:
    ks: List[int] = []
    for _, il, iu, jl, ju, kl, ku in blocks:
        for _, _, k in _iter_block_points(il, iu, jl, ju, kl, ku):
            ks.append(k)
    return ks


def _write_sparse_file(out_path: Path, points: Sequence[Tuple[int, int, int]]) -> None:
    with out_path.open("w", encoding="ascii", newline="\n") as f:
        f.write("# position (i,j,k)\n")
        for i, j, k in points:
            f.write(f"{i:7d} {j:7d} {k:7d}\n")


def _write_sparse_id_file(out_path: Path, ids: Sequence[int]) -> None:
    with out_path.open("w", encoding="ascii", newline="\n") as f:
        f.write("# block_id for each point in veg.inp, same order\n")
        for block_id in ids:
            f.write(f"{block_id:7d}\n")


def _write_params_file(
    out_path: Path,
    ids: Sequence[int],
    lad_values: Sequence[float],
    cd: float,
    ud: float,
    dec: float,
    lsize: float,
    r_s: float,
) -> None:
    with out_path.open("w", encoding="ascii", newline="\n") as f:
        f.write("# id lad cd ud dec lsize r_s\n")
        for bid, lad_val in zip(ids, lad_values):
            f.write(
                f"{bid:7d} {lad_val:12.6f} {cd:12.6f} {ud:12.6f} {dec:12.6f} {lsize:12.6f} {r_s:12.6f}\n"
            )


def _lad_for_points(
    points: Sequence[Tuple[int, int, int]],
    lad_value: float,
) -> List[float]:
    """Assign constant LAD to all points.

    LAD is now constant per vegetation point, not height-dependent.
    """
    return [lad_value for _ in points]


def _get_required(sim: "UDBase", name: str) -> float:
    if not hasattr(sim, name):
        raise AttributeError(f"UDBase missing required attribute '{name}' for veg params")
    return float(getattr(sim, name))


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
    point_ids = _blocks_to_id_list(blocks)

    out_path = sim_dir / f"veg.inp.{sim_id}"
    _write_sparse_file(out_path, points)
    # Parameter file per sparse point (constant LAD for all points)
    cd = _get_required(sim, "cd")
    ud = _get_required(sim, "ud")
    dec = _get_required(sim, "dec")
    lsize = _get_required(sim, "lsize")
    try:
        r_s = _get_required(sim, "r_s")
    except AttributeError:
        r_s = _get_required(sim, "rs")

    lad_vals = _lad_for_points(points, lad_value_resolved)
    out_params_path = sim_dir / f"veg_params.inp.{sim_id}"
    _write_params_file(out_params_path, point_ids, lad_vals, cd, ud, dec, lsize, r_s)

    update_namelist_value(sim, "TREES", "ntrees", len(points))

    print(
        f"Loaded {len(blocks)} tree blocks from {sim_dir}, "
        f"expanded to {len(points)} grid points",
    )
    print(f"Sparse output written to {out_path}")
    print(f"Veg params (id, lad={lad_value_resolved} from {lad_source}, cd, ud, dec, lsize, r_s) written to {out_params_path}")
    return out_path


if __name__ == "__main__":
    raise SystemExit("Run via veg.py which constructs UDBase and passes it here")
