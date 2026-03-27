"""
Geometry generation helpers ported from MATLAB udgeom.

These functions mirror the MATLAB implementations in tools/matlab/+udgeom:
- createFlatSurface.m
- createCanyons.m
- createCubes.m
- createRealistic.m
"""

from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple
import warnings

import numpy as np

try:
    import trimesh
except ImportError as exc:
    raise ImportError("trimesh is required for geometry generation; install with `pip install trimesh`.") from exc

try:
    from scipy.spatial import Delaunay
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

from .udgeom import UDGeom


# -----------------------------------------------------------------------------#
# Utility helpers
# -----------------------------------------------------------------------------#

def _remove_unused(vertices: np.ndarray, faces: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Remove vertices that are not referenced by any face."""
    used = np.unique(faces.flatten())
    new_index = -np.ones(len(vertices), dtype=int)
    new_index[used] = np.arange(len(used))
    new_vertices = vertices[used]
    new_faces = new_index[faces]
    return new_vertices, new_faces


def _divide_faces(mesh: "trimesh.Trimesh", times: int = 1) -> "trimesh.Trimesh":
    """
    Subdivide each triangular face into four by adding midpoints.
    Repeats `times` times (MATLAB divideFaces equivalent).
    """
    if times <= 0:
        return mesh.copy()

    current = mesh
    for _ in range(times):
        verts = current.vertices
        faces = current.faces

        edge_mid_cache = {}

        def midpoint(a: int, b: int) -> int:
            key = tuple(sorted((a, b)))
            if key in edge_mid_cache:
                return edge_mid_cache[key]
            mid = (verts[a] + verts[b]) / 2.0
            edge_mid_cache[key] = len(verts_list)
            verts_list.append(mid)
            return edge_mid_cache[key]

        verts_list = verts.tolist()
        new_faces = []
        for f in faces:
            a, b, c = f
            ab = midpoint(a, b)
            bc = midpoint(b, c)
            ca = midpoint(c, a)
            new_faces.extend(
                [
                    (a, ab, ca),
                    (ab, b, bc),
                    (ca, bc, c),
                    (ab, bc, ca),
                ]
            )

        verts_arr = np.asarray(verts_list, dtype=float)
        faces_arr = np.asarray(new_faces, dtype=int)
        verts_arr, faces_arr = _remove_unused(verts_arr, faces_arr)
        current = trimesh.Trimesh(vertices=verts_arr, faces=faces_arr, process=False)

    return current


def _unit_cube() -> "trimesh.Trimesh":
    """
    Unit cube centered at origin with side length 1 (open bottom to match MATLAB).
    Vertices are in [-0.5, 0.5] for each axis.
    """
    verts = np.array([
        [-0.5, -0.5, -0.5],
        [-0.5,  0.5, -0.5],
        [ 0.5,  0.5, -0.5],
        [ 0.5, -0.5, -0.5],
        [ 0.5, -0.5,  0.5],
        [-0.5, -0.5,  0.5],
        [-0.5,  0.5,  0.5],
        [ 0.5,  0.5,  0.5],
    ], dtype=float)

    faces = np.array([
        [0, 1, 2], [2, 3, 0],  # bottom (open in MATLAB; keep for watertightness)
        [0, 3, 4], [4, 5, 0],  # side
        [0, 5, 6], [6, 1, 0],  # side
        [7, 4, 3], [3, 2, 7],  # side
        [7, 2, 1], [1, 6, 7],  # side
        [7, 6, 5], [5, 4, 7],  # top
    ], dtype=int)

    return trimesh.Trimesh(vertices=verts, faces=faces, process=False)


def _operate_unit_cube(scale: Sequence[float], shifts: np.ndarray, divisions: int) -> "trimesh.Trimesh":
    """
    Scale and shift subdivided unit cube copies, then concatenate.
    """
    base = _divide_faces(_unit_cube(), times=divisions)
    scaled_vertices = base.vertices * np.asarray(scale, dtype=float)
    base_faces = base.faces

    meshes = []
    for shift in shifts:
        verts = scaled_vertices + np.asarray(shift, dtype=float)
        meshes.append(trimesh.Trimesh(vertices=verts, faces=base_faces, process=False))

    return trimesh.util.concatenate(meshes)


def _structured_ground_points(xsize: float, ysize: float, edgelength: float) -> np.ndarray:
    """Create regularly spaced ground points including domain edges."""
    nx = int(round(xsize / edgelength))
    ny = int(round(ysize / edgelength))
    xs = np.linspace(0.0, xsize, nx + 1)
    ys = np.linspace(0.0, ysize, ny + 1)
    pts = np.array([[x, y] for x in xs for y in ys], dtype=float)
    return pts


def _triangulate_ground(points_xy: np.ndarray, xsize: float, ysize: float) -> "trimesh.Trimesh":
    """
    Triangulate 2D points on the ground plane.
    Uses SciPy Delaunay if available, otherwise falls back to a structured grid.
    """
    if SCIPY_AVAILABLE:
        tri = Delaunay(points_xy)
        simplices = tri.simplices
        centroids = points_xy[simplices].mean(axis=1)
        keep = (
            (centroids[:, 0] >= -1e-9)
            & (centroids[:, 0] <= xsize + 1e-9)
            & (centroids[:, 1] >= -1e-9)
            & (centroids[:, 1] <= ysize + 1e-9)
        )
        simplices = simplices[keep]
        vertices_3d = np.column_stack([points_xy, np.zeros(len(points_xy))])
        return trimesh.Trimesh(vertices=vertices_3d, faces=simplices, process=False)

    # Fallback: build structured grid
    nx = int(round(xsize / (points_xy[:, 0].ptp() / max(1, len(np.unique(points_xy[:, 0])) - 1))))
    ny = int(round(ysize / (points_xy[:, 1].ptp() / max(1, len(np.unique(points_xy[:, 1])) - 1))))
    xs = np.linspace(0.0, xsize, nx + 1)
    ys = np.linspace(0.0, ysize, ny + 1)
    vertices = np.array([(x, y, 0.0) for x in xs for y in ys], dtype=float)
    faces: List[Tuple[int, int, int]] = []

    def vid(i: int, j: int) -> int:
        return i * (ny + 1) + j

    for i in range(nx):
        for j in range(ny):
            v0 = vid(i, j)
            v1 = vid(i + 1, j)
            v2 = vid(i + 1, j + 1)
            v3 = vid(i, j + 1)
            faces.append((v0, v1, v2))
            faces.append((v0, v2, v3))

    return trimesh.Trimesh(vertices=vertices, faces=np.array(faces), process=False)


def _generate_ground(existing: Optional["trimesh.Trimesh"], xsize: float, ysize: float, edgelength: float) -> "trimesh.Trimesh":
    """
    Port of MATLAB generateGround: add ground facets and merge with existing mesh.
    """
    ground_points = np.empty((0, 2))
    if existing is not None and len(existing.faces) > 0:
        ground_vert_ids = np.where(np.isclose(existing.vertices[:, 2], 0.0))[0]
        ground_points = existing.vertices[ground_vert_ids][:, :2]

    grid_pts = _structured_ground_points(xsize, ysize, edgelength)
    all_points = np.vstack([ground_points, grid_pts])
    ground_mesh = _triangulate_ground(all_points, xsize, ysize)

    if existing is None or len(existing.faces) == 0:
        return ground_mesh

    combined = trimesh.util.concatenate([existing, ground_mesh])

    # Remove ground faces under buildings (requires watertightness)
    verts = combined.vertices
    faces = combined.faces
    ground_mask = np.all(np.isclose(verts[faces][:, :, 2], 0.0), axis=1)

    if existing.is_watertight:
        ground_centers = verts[faces[ground_mask]].mean(axis=1)
        try:
            inside = existing.contains(ground_centers)
            ground_mask_indices = np.where(ground_mask)[0]
            faces = np.delete(faces, ground_mask_indices[inside], axis=0)
        except Exception:
            pass

    verts, faces = _remove_unused(verts, faces)
    return trimesh.Trimesh(vertices=verts, faces=faces, process=False)


# -----------------------------------------------------------------------------#
# Public API
# -----------------------------------------------------------------------------#

def create_flat_surface(xsize: float, ysize: float, edgelength: float) -> UDGeom:
    """
    Create a flat surface consisting of triangular facets (MATLAB createFlatSurface).
    """
    ground = _generate_ground(None, xsize, ysize, edgelength)
    return UDGeom(stl=ground)


def create_canyons(
    xsize: float,
    ysize: float,
    B: float,
    W: float,
    H: float,
    shift: float,
    edgelength: float,
    rotate90: bool = False,
) -> UDGeom:
    """
    Create one-dimensional street canyons (MATLAB createCanyons).
    """
    L = B + W
    divisions = int(round(B / edgelength))
    if divisions > 2:
        raise ValueError("divisions must be <= 2 (B/edgelength <= 2)")
    if divisions == 2:
        L = L / 2.0

    Nx = xsize / (B + W)
    Ny = ysize / L
    if abs(Nx - round(Nx)) > 1e-9 or abs(Ny - round(Ny)) > 1e-9:
        raise ValueError("The domain size should be a multiple of canyon width/length")
    Nx = int(round(Nx))
    Ny = int(round(Ny))

    # Base canyon vertices (floors, walls, roof)
    points_floor1 = np.array([[0, 0, 0], [W / 2, 0, 0], [W / 2, L, 0], [0, L, 0]], dtype=float)
    points_floor2 = points_floor1 + np.array([W / 2 + B, 0, 0])
    points_wall1 = np.array([[0, 0, 0], [0, 0, H], [0, L, H], [0, L, 0]], dtype=float) + np.array([W / 2, 0, 0])
    points_wall2 = np.flip(points_wall1, axis=0) + np.array([B, 0, 0])
    points_roof = np.array([[0, 0, 0], [B, 0, 0], [B, L, 0], [0, L, 0]], dtype=float) + np.array([0, 0, H]) + np.array([W / 2, 0, 0])
    points_unit = np.vstack([points_floor1, points_wall1, points_roof, points_wall2, points_floor2])

    conn: List[Tuple[int, int, int]] = []
    points: List[Tuple[float, float, float]] = []

    def add_point(pt: np.ndarray) -> int:
        p = tuple(pt.tolist())
        if p in point_index:
            return point_index[p]
        idx = len(points)
        points.append(p)
        point_index[p] = idx
        return idx

    point_index = {}

    for i in range(Nx):
        for j in range(Ny):
            for n in range(5):
                locs = []
                for m in range(4):
                    pt = points_unit[4 * n + m] + np.array([(i) * (B + W), (j) * L, 0])
                    locs.append(add_point(pt))
                conn.append((locs[0], locs[1], locs[2]))
                conn.append((locs[0], locs[2], locs[3]))

    vertices = np.array(points, dtype=float)
    faces = np.array(conn, dtype=int)

    # shift interior points in x (matching MATLAB behavior)
    x_max = vertices[:, 0].max() if len(vertices) else 0.0
    mask_shift = (vertices[:, 0] > 0.0) & (vertices[:, 0] < x_max)
    vertices[mask_shift, 0] += shift

    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)

    if divisions >= 1:
        mesh = _divide_faces(mesh, times=1)

    if rotate90:
        angle_rad = np.pi / 2.0
        R = np.array([[np.cos(angle_rad), -np.sin(angle_rad), 0.0],
                      [np.sin(angle_rad), np.cos(angle_rad), 0.0],
                      [0.0, 0.0, 1.0]])
        rot = np.block([[R, np.zeros((3, 1))], [np.zeros((1, 3)), np.ones((1, 1))]])
        # Match MATLAB: rotate first, then shift by xsize in +x
        mesh.apply_transform(rot)
        mesh.apply_translation([xsize, 0.0, 0.0])

    # Remove existing ground facets (z == 0) before adding merged ground
    #verts = mesh.vertices
    #faces = mesh.faces
    #ground_mask = np.all(np.isclose(verts[faces][:, :, 2], 0.0), axis=1)
    #if np.any(ground_mask):
    #    faces = faces[~ground_mask]
    #    verts, faces = _remove_unused(verts, faces)
    #    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)

    # Ground generation disabled (request): return canyon mesh directly
    return UDGeom(stl=mesh)


def create_cubes(
    xsize: float,
    ysize: float,
    Hx: float,
    Hy: float,
    Hz: float,
    Cx: float,
    Cy: float,
    geom_option: str,
    edgelength: float,
) -> UDGeom:
    """
    Create cubes: single, aligned, or staggered (MATLAB createCubes).

    Parameters
    ----------
    xsize, ysize : float
        Domain size in meters (x, y).
    Hx, Hy, Hz : float
        Cube lengths in x, y and height z.
    Cx, Cy : float
        Spacing between cubes in x and y directions (ignored for single cube).
    geom_option : {'S','AC','SC'}
        'S'  - single cube centred in the domain.
        'AC' - aligned array on a regular grid.
        'SC' - staggered array (alternate rows shifted by half spacing).
    edgelength : float
        Target facet size; controls subdivision level.

    Raises
    ------
    ValueError
        If domain is not an integer multiple of cube+spacing for AC/SC,
        or if geom_option is invalid.
    """
    divisions = int(round(Hx / edgelength))
    tol = 0.0

    opt = geom_option.upper()
    if opt not in {"S", "AC", "SC"}:
        raise ValueError("geom_option must be 'S', 'AC', or 'SC'")

    # Build cube arrays
    if opt == "SC":
        Nx = xsize / (Hx + Cx)
        Ny = ysize / (Hy + Cy)
        if abs(Nx - round(Nx)) > 1e-9 or abs(Ny - round(Ny)) > 1e-9:
            raise ValueError("The domain size should be a multiple of cube width + canyon width")
        Nx = int(round(Nx))
        Ny = int(round(Ny))

        shifts = []
        for i in range(1, Nx + 1):
            for j in range(1, Ny + 1):
                if i % 2 == 1:
                    shifts.append([i * (Cx + Hx) - Hx / 2 - Cx / 2, j * (Cy + Hy) - Hy / 2 - Cy / 2, Hz / 2])
                else:
                    shifts.append([i * (Cx + Hx) - Hx / 2 - Cx / 2, j * (Cy + Hy), Hz / 2])
        array = _operate_unit_cube([Hx, Hy, Hz], np.array(shifts), divisions)
        mesh = array

        # Trim outside domain in y
        verts = mesh.vertices.copy()
        faces = mesh.faces.copy()
        if divisions < 2:
            verts[:, 1] = np.clip(verts[:, 1], 0.0, ysize)

        points_remove = np.where((verts[:, 1] <= 0.0) | (verts[:, 1] >= ysize))[0]
        if len(points_remove) > 0:
            mask = np.all(np.isin(faces, points_remove), axis=1)
            faces = faces[~mask]
        verts, faces = _remove_unused(verts, faces)
        mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)

    elif opt == "AC":
        Nx = xsize / (Hx + Cx)
        Ny = ysize / (Hy + Cy)
        if abs(Nx - round(Nx)) > 1e-9 or abs(Ny - round(Ny)) > 1e-9:
            raise ValueError("The domain size should be a multiple of cube width + canyon width")
        Nx = int(round(Nx))
        Ny = int(round(Ny))

        shifts = []
        for i in range(1, Nx + 1):
            for j in range(1, Ny + 1):
                shifts.append([i * (Cx + Hx) - Hx / 2 - Cx / 2, j * (Cy + Hy) - Hy / 2 - Cy / 2, Hz / 2])
        mesh = _operate_unit_cube([Hx, Hy, Hz], np.array(shifts), divisions)

    else:  # 'S'
        shift = np.array([[xsize / 2.0, ysize / 2.0, 0.0]]) + np.array([0.0, 0.0, Hz / 2.0])
        mesh = _operate_unit_cube([Hx, Hy, Hz], shift, divisions)

    # Apply tiny normal shifts (tol) except -z faces (tol is zero here, kept for parity)
    if tol != 0.0:
        normals = mesh.face_normals
        verts = mesh.vertices.copy()
        faces = mesh.faces
        for idx, nrm in enumerate(normals):
            ids = faces[idx]
            if np.allclose(nrm, [1, 0, 0]):
                verts[ids] += tol * np.array([1, 0, 0])
            elif np.allclose(nrm, [-1, 0, 0]):
                verts[ids] += tol * np.array([-1, 0, 0])
            elif np.allclose(nrm, [0, 1, 0]):
                verts[ids] += tol * np.array([0, 1, 0])
            elif np.allclose(nrm, [0, -1, 0]):
                verts[ids] += tol * np.array([0, -1, 0])
            elif np.allclose(nrm, [0, 0, 1]):
                verts[ids] += tol * np.array([0, 0, 1])
        mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)

    ground = _generate_ground(mesh, xsize, ysize, edgelength)
    return UDGeom(stl=ground)


def create_realistic(
    stlfile: Path,
    xsize: float,
    ysize: float,
    shift: Tuple[float, float, float],
    edgelength: float,
) -> UDGeom:
    """
    Create a realistic urban surface by loading buildings from STL and adding ground.
    Mirrors MATLAB createRealistic.
    """
    building_mesh = trimesh.load_mesh(str(stlfile))
    building_mesh.apply_translation(np.asarray(shift, dtype=float))
    ground = _generate_ground(building_mesh, xsize, ysize, edgelength)
    return UDGeom(stl=ground)


__all__ = [
    "create_flat_surface",
    "create_canyons",
    "create_cubes",
    "create_realistic",
]
