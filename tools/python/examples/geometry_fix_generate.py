"""
geometry_fix_generate - Create the synthetic bad STL inputs used by geometry_fix.

This script only generates the example input geometries:
- adjacent_buildings.stl
- building_below_ground.stl
- building_above_ground.stl
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import trimesh

TOOLS_PYTHON = Path(__file__).resolve().parents[1]
if str(TOOLS_PYTHON) not in sys.path:
    sys.path.insert(0, str(TOOLS_PYTHON))

from udgeom import UDGeom, add_ground


EXAMPLES_DIR = Path(__file__).resolve().parent


def save_geometry_stl(geom: UDGeom, stem: str):
    geom.stl.export(EXAMPLES_DIR / f"{stem}.stl")


def remove_legacy_outputs():
    for name in (
        "geometry_fix_before.stl",
        "geometry_fix_subground_before.stl",
    ):
        path = EXAMPLES_DIR / name
        if path.exists():
            path.unlink()


def concatenate_and_weld(meshes, digits_vertex: int = 8):
    mesh = trimesh.util.concatenate(meshes)
    mesh.merge_vertices(digits_vertex=digits_vertex)
    mesh.remove_unreferenced_vertices()
    return mesh


def make_block(x0, x1, y0, y1, z1, nx, ny, nz, include_bottom: bool = True):
    xs = np.linspace(x0, x1, nx + 1)
    ys = np.linspace(y0, y1, ny + 1)
    zs = np.linspace(0.0, z1, nz + 1)
    vertices = []
    index = {}

    def add_vertex(x, y, z):
        key = (float(x), float(y), float(z))
        if key not in index:
            index[key] = len(vertices)
            vertices.append([float(x), float(y), float(z)])
        return index[key]

    faces = []

    def add_rect(p00, p10, p11, p01):
        i00 = add_vertex(*p00)
        i10 = add_vertex(*p10)
        i11 = add_vertex(*p11)
        i01 = add_vertex(*p01)
        faces.append([i00, i10, i11])
        faces.append([i00, i11, i01])

    for i in range(nx):
        for j in range(ny):
            if include_bottom:
                add_rect((xs[i], ys[j], 0.0), (xs[i + 1], ys[j], 0.0), (xs[i + 1], ys[j + 1], 0.0), (xs[i], ys[j + 1], 0.0))
            add_rect((xs[i], ys[j], z1), (xs[i], ys[j + 1], z1), (xs[i + 1], ys[j + 1], z1), (xs[i + 1], ys[j], z1))

    for j in range(ny):
        for k in range(nz):
            add_rect((x0, ys[j], zs[k]), (x0, ys[j + 1], zs[k]), (x0, ys[j + 1], zs[k + 1]), (x0, ys[j], zs[k + 1]))
            add_rect((x1, ys[j], zs[k]), (x1, ys[j], zs[k + 1]), (x1, ys[j + 1], zs[k + 1]), (x1, ys[j + 1], zs[k]))

    for i in range(nx):
        for k in range(nz):
            add_rect((xs[i], y0, zs[k]), (xs[i + 1], y0, zs[k]), (xs[i + 1], y0, zs[k + 1]), (xs[i], y0, zs[k + 1]))
            add_rect((xs[i], y1, zs[k]), (xs[i], y1, zs[k + 1]), (xs[i + 1], y1, zs[k + 1]), (xs[i + 1], y1, zs[k]))

    return trimesh.Trimesh(vertices=np.asarray(vertices, dtype=float), faces=np.asarray(faces, dtype=int), process=False)


def make_ground(x0, x1, y0, y1, nx, ny):
    xs = np.linspace(x0, x1, nx + 1)
    ys = np.linspace(y0, y1, ny + 1)
    vertices = []
    index = {}
    faces = []

    def add_vertex(x, y):
        key = (float(x), float(y), 0.0)
        if key not in index:
            index[key] = len(vertices)
            vertices.append([float(x), float(y), 0.0])
        return index[key]

    for i in range(nx):
        for j in range(ny):
            i00 = add_vertex(xs[i], ys[j])
            i10 = add_vertex(xs[i + 1], ys[j])
            i11 = add_vertex(xs[i + 1], ys[j + 1])
            i01 = add_vertex(xs[i], ys[j + 1])
            faces.append([i00, i10, i11])
            faces.append([i00, i11, i01])

    return trimesh.Trimesh(vertices=np.asarray(vertices, dtype=float), faces=np.asarray(faces, dtype=int), process=False)


def build_adjacent_buildings_geometry() -> UDGeom:
    small = make_block(90.0, 110.0, 50.0, 130.0, 30.0, nx=1, ny=4, nz=2, include_bottom=False)
    large = make_block(110.0, 140.0, 30.0, 150.0, 60.0, nx=1, ny=4, nz=4, include_bottom=False)
    shell = concatenate_and_weld([small, large])
    geom = add_ground(shell, 240.0, 180.0, edgelength=30.0, preserve_existing_edges=True)
    geom.stl.merge_vertices(digits_vertex=8)
    geom.stl.remove_unreferenced_vertices()
    return geom


def build_building_below_ground_geometry() -> UDGeom:
    ground = make_ground(0.0, 240.0, 0.0, 240.0, nx=5, ny=5)
    building = make_block(90.0, 150.0, 90.0, 150.0, 180.0, nx=2, ny=2, nz=3, include_bottom=False)
    verts = np.asarray(building.vertices, dtype=float).copy()
    verts[:, 2] -= 60.0
    building = trimesh.Trimesh(vertices=verts, faces=np.asarray(building.faces, dtype=int), process=False)
    mesh = concatenate_and_weld([ground, building])
    return UDGeom(stl=mesh)


def build_building_above_ground_geometry() -> UDGeom:
    ground = make_ground(0.0, 240.0, 0.0, 240.0, nx=5, ny=5)
    building = make_block(90.0, 150.0, 90.0, 150.0, 120.0, nx=2, ny=2, nz=3, include_bottom=False)
    verts = np.asarray(building.vertices, dtype=float).copy()
    verts[:, 2] += 5.0
    building = trimesh.Trimesh(vertices=verts, faces=np.asarray(building.faces, dtype=int), process=False)
    mesh = concatenate_and_weld([ground, building])
    return UDGeom(stl=mesh)


def main():
    remove_legacy_outputs()

    adjacent = build_adjacent_buildings_geometry()
    below_ground = build_building_below_ground_geometry()
    above_ground = build_building_above_ground_geometry()

    save_geometry_stl(adjacent, "adjacent_buildings")
    save_geometry_stl(below_ground, "building_below_ground")
    save_geometry_stl(above_ground, "building_above_ground")

    print("Saved:")
    print(f"- {EXAMPLES_DIR / 'adjacent_buildings.stl'}")
    print(f"- {EXAMPLES_DIR / 'building_below_ground.stl'}")
    print(f"- {EXAMPLES_DIR / 'building_above_ground.stl'}")


if __name__ == "__main__":
    main()
