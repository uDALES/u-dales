"""Shared face-adjacency graph primitives for udgeom.

Several routines (independent-surface labelling, ground-face detection, building
splitting, mesh repair) each rebuilt the same face-adjacency dict and hand-rolled
the same depth-first connected-components walk. These two helpers are that common
core; callers supply their own seed order and traversal filter, and apply their
own ID/ordering policy to the returned components.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Set

import numpy as np


def build_face_adjacency(mesh) -> Dict[int, Set[int]]:
    """Return {face_index: set(neighbouring face indices)} from ``mesh.face_adjacency``."""
    adjacency: Dict[int, Set[int]] = {i: set() for i in range(len(mesh.faces))}
    for face_a, face_b in np.asarray(mesh.face_adjacency, dtype=int):
        adjacency[int(face_a)].add(int(face_b))
        adjacency[int(face_b)].add(int(face_a))
    return adjacency


def connected_components(
    adjacency: Dict[int, Set[int]],
    *,
    seeds: Optional[Iterable[int]] = None,
    keep: Optional[np.ndarray] = None,
) -> List[np.ndarray]:
    """Depth-first connected components of an adjacency graph.

    Parameters
    ----------
    adjacency : dict
        ``{node: set(neighbours)}`` (e.g. from :func:`build_face_adjacency`).
    seeds : iterable of int, optional
        Nodes to start from, in order. Defaults to ``range(len(adjacency))``,
        which yields components in ascending node order.
    keep : bool array, optional
        If given, only traverse into node ``i`` when ``keep[i]`` is True (used to
        restrict components to a subset such as horizontal faces).

    Returns
    -------
    list of ndarray
        One int array of node indices per component, in seed order.
    """
    if seeds is None:
        seeds = range(len(adjacency))
    visited: Set[int] = set()
    components: List[np.ndarray] = []
    for seed in seeds:
        seed = int(seed)
        if seed in visited:
            continue
        stack = [seed]
        visited.add(seed)
        component: List[int] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbour in adjacency[current]:
                if neighbour not in visited and (keep is None or keep[neighbour]):
                    visited.add(neighbour)
                    stack.append(neighbour)
        components.append(np.asarray(component, dtype=int))
    return components
