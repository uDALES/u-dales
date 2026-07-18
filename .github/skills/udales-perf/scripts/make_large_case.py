#!/usr/bin/env python3
"""Generate the large-tier benchmark geometry (udales-perf skill, Phase 3).

Tiles a source STL (the Xie-Castro staggered-cube unit, 160x160 m, ~1.1k
facets) nx x ny times at a given pitch and subdivides every triangle
`sub` times (each level splits a triangle into 4 at edge midpoints, changing
only the mesh resolution, not the geometry). Defaults produce ~1.13M facets
covering 1280x1280 m -- pair with a 2048x2048x256 grid for the O(1024^3)
benchmark case.

Usage:
  make_large_case.py <in.stl> <out.stl> [--nx 8] [--ny 8] [--pitch 160]
                     [--sub 2]

Reads/writes binary STL; deterministic (no randomness).
"""
import argparse
import struct
import sys

import numpy as np


def read_stl(path):
    with open(path, "rb") as f:
        header = f.read(80)
        if header[:5] == b"solid":
            sys.exit("ASCII STL not supported; convert to binary first")
        (n,) = struct.unpack("<I", f.read(4))
        data = np.fromfile(f, dtype=np.uint8, count=n * 50).reshape(n, 50)
    tri = data[:, :48].copy().view("<f4").reshape(n, 12)  # normal + 3 vertices
    return tri[:, 3:].reshape(n, 3, 3).astype(np.float64)  # drop normals


def write_stl(path, verts):
    n = verts.shape[0]
    v01 = verts[:, 1] - verts[:, 0]
    v02 = verts[:, 2] - verts[:, 0]
    nrm = np.cross(v01, v02)
    length = np.linalg.norm(nrm, axis=1, keepdims=True)
    length[length == 0.0] = 1.0
    nrm /= length
    rec = np.zeros((n, 50), dtype=np.uint8)
    rec[:, :48] = (
        np.hstack([nrm, verts.reshape(n, 9)]).astype("<f4").view(np.uint8).reshape(n, 48)
    )
    with open(path, "wb") as f:
        f.write(b"udales-perf large-tier benchmark geometry".ljust(80, b" "))
        f.write(struct.pack("<I", n))
        rec.tofile(f)


def subdivide(verts):
    """One 4-way midpoint subdivision of all triangles."""
    a, b, c = verts[:, 0], verts[:, 1], verts[:, 2]
    ab, bc, ca = (a + b) / 2, (b + c) / 2, (c + a) / 2
    return np.concatenate(
        [
            np.stack([a, ab, ca], axis=1),
            np.stack([ab, b, bc], axis=1),
            np.stack([ca, bc, c], axis=1),
            np.stack([ab, bc, ca], axis=1),
        ]
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("infile")
    p.add_argument("outfile")
    p.add_argument("--nx", type=int, default=8)
    p.add_argument("--ny", type=int, default=8)
    p.add_argument("--pitch", type=float, default=160.0)
    p.add_argument("--sub", type=int, default=2)
    args = p.parse_args()

    unit = read_stl(args.infile)
    print(f"source: {unit.shape[0]} facets")

    tiles = []
    for i in range(args.nx):
        for j in range(args.ny):
            t = unit.copy()
            t[:, :, 0] += i * args.pitch
            t[:, :, 1] += j * args.pitch
            tiles.append(t)
    verts = np.concatenate(tiles)
    print(f"tiled {args.nx}x{args.ny}: {verts.shape[0]} facets")

    for _ in range(args.sub):
        verts = subdivide(verts)
    print(f"after {args.sub} subdivision level(s): {verts.shape[0]} facets")

    write_stl(args.outfile, verts)
    ext = verts.reshape(-1, 3)
    print(f"extent: x [{ext[:,0].min():g}, {ext[:,0].max():g}]  "
          f"y [{ext[:,1].min():g}, {ext[:,1].max():g}]  "
          f"z [{ext[:,2].min():g}, {ext[:,2].max():g}]")
    print(f"wrote {args.outfile}")


if __name__ == "__main__":
    main()
