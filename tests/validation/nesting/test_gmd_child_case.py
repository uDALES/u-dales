#!/usr/bin/env python3
"""``make_gmd_child_case`` must build a nested child of the GMD 2024 parent (567).

Unit tests here need no solver and (mostly) no parent run -- they check the
grid, the geometry translation and the namelist transform against numbers
either fixed by construction or read from the parent's own ALREADY-WRITTEN
``namoptions.567``/``prof.inp.567`` (static inputs, not run output). The two
tests that need the parent's own simulation to have produced its
``nesting.out.*`` band files -- the ``&NESTPARENT`` box cross-check (which
only needs the file to exist, not the run) and the full end-to-end build --
are gated by ``UDALES_GMD_PARENT_DIR`` and skip with a reason when it is unset.

No solver, no MPI, for everything except the (optional) end-to-end build,
which calls the repo's IBM preprocessing.
"""

from __future__ import annotations

import os
import re
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import numpy as np

import indoor_geometry as ig
import make_gmd_child_case as gc

#: Synthetic stand-in for the parent's namoptions.567 -- the same sections and
#: keys as the real file (see the module docstring's parent facts), so
#: child_sections is exercised against realistic input without needing the
#: real, machine-specific parent case directory.
_SYNTHETIC_NAMOPTIONS = """
&RUN
iexpnr       = 567
runtime      = 150.
trestart     = 12.5
lwarmstart   = .false.
startfile    = ''
dtmax        = 0.0008
irandom      = 43
randu        = 0.01
ladaptive    = .true.
nprocx       = 16
nprocy       = 16
libm         = .true.
/

&OUTPUT
lfielddump   = .false.
tstatsdump   = 37.5
tsample      = 0.008
/

&DOMAIN
itot         = 1024
jtot         = 512
ktot         = 128
xlen         = 3.42
ylen         = 1.8
/

&PHYSICS
ps           = 101500.00
igrw_damp    = 0
/

&DYNAMICS
ipoiss       = 0
/

&NAMSUBGRID
lvreman      = .true.
/

&BC
BCxm         = 3
BCtopm       = 3
wtsurf       = 0.
wqsurf       = 0.
thls         = 288.
z0           = 0.0009
z0h          = 0.000067
/

&WALLS
nfcts = 352
nsolpts_u = 92016
nbndpts_u = 555572
nfctsecs_u = 522465
iwallmom   = 1
/

&DRIVER
idriver      = 2
driverjobnr  = 567
driverstore  = 187501
lchunkread   = .true.
chunkread_size = 1000
/

&INPS
zsize        = 0.96
lzstretch    = .true.
hlin         = 0.24
dzlin        = 0.0035
lstretchexp  = .true.
stretchconst = 1.
u0           = 4.9524
stl_file     = indoor_object_final.stl
stl_ground   = .true.
diag_neighbs = .true.
/

&NESTPARENT
lnestparent      = .true.
tnestparent      = 0.0012
nestparent_x0    = 0.3607031250
nestparent_y0    = 0.6011718750
nestparent_xsize = 0.7214062500
nestparent_ysize = 0.6011718750
nestparent_nzone = 13
nestparent_linit = .true.
/
"""


def _write_synthetic_parent(tmp: Path) -> Path:
    """A parent case directory good enough for the grid/namelist tests.

    ``prof.inp.567``'s zf is NOT the real (stretched) column -- it is a stand-in,
    monotone and the right length -- since these tests exercise the box/namelist
    machinery, not the parent's actual stretching. The real cross-checks
    (TestBoxMatchesRealParent, TestEndToEnd) read the real file instead, gated
    by UDALES_GMD_PARENT_DIR.
    """
    parent = tmp / "prod567"
    parent.mkdir()
    (parent / f"namoptions.{gc.PARENT_EXPNR}").write_text(_SYNTHETIC_NAMOPTIONS)
    zf = (np.arange(gc.PARENT_KTOT) + 0.5) * (0.96 / gc.PARENT_KTOT)
    with (parent / f"prof.inp.{gc.PARENT_EXPNR}").open("w", encoding="ascii") as f:
        f.write("# synthetic stand-in\n# z thl qt u v e12\n")
        for z in zf:
            f.write(f"{z:.12f} 288.0 0.0 4.9524 0.0 0.0\n")
    (parent / f"lscale.inp.{gc.PARENT_EXPNR}").write_text(
        "# synthetic\n# z ug vg pgx pgy wfls dqtdxls dqtdyls dqtdtls dthlrad\n"
        + "".join(f"{z:.12f}" + " 0.0" * 9 + "\n" for z in zf)
    )
    return parent


def _real_parent_dir() -> "Path | None":
    root = os.environ.get("UDALES_GMD_PARENT_DIR")
    return Path(root) if root else None


_VERTEX_RE = re.compile(r"vertex\s+([\-0-9.eE+]+)\s+([\-0-9.eE+]+)\s+([\-0-9.eE+]+)")


def _ascii_stl_vertices(path: Path) -> np.ndarray:
    text = Path(path).read_text(encoding="ascii")
    return np.asarray(_VERTEX_RE.findall(text), dtype=float)


class TestChildGridNesting(unittest.TestCase):
    """The child's grid must nest in the parent's, at r = 1 and r = 2."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.parent_dir = _write_synthetic_parent(Path(self._tmp.name))
        self.parent = gc.parent_grid(self.parent_dir)

    def tearDown(self):
        self._tmp.cleanup()

    def test_stretched_zh_is_preserved_exactly(self):
        for r in (1, 2):
            with self.subTest(r=r):
                box = gc.ChildBox(refine=r)
                child = gc.child_grid(box, self.parent)
                self.assertTrue(np.array_equal(child.zh, self.parent.zh))

    def test_every_parent_face_in_the_box_is_a_child_face(self):
        for r in (1, 2):
            with self.subTest(r=r):
                box = gc.ChildBox(refine=r)
                child = gc.child_grid(box, self.parent)  # raises if misaligned
                tol = 1.0e-9 * min(np.min(np.diff(self.parent.xh)),
                                   np.min(np.diff(self.parent.yh)))
                for axis, parent_edges, child_edges in (
                    ("x", self.parent.xh, child.xh), ("y", self.parent.yh, child.yh),
                ):
                    inside = parent_edges[(parent_edges >= child_edges[0] - tol)
                                          & (parent_edges <= child_edges[-1] + tol)]
                    for p in inside:
                        nearest = np.min(np.abs(child_edges - p))
                        self.assertLess(nearest, tol,
                                        f"axis {axis}: parent face {p:g} at r={r} is not "
                                        "a child face")

    def test_box_size_scales_with_the_grid_not_the_refinement(self):
        # itot/jtot double at r = 2, but xlen/ylen -- the physical box -- do not.
        b1 = gc.ChildBox(refine=1)
        b2 = gc.ChildBox(refine=2)
        self.assertEqual(b2.itot, 2 * b1.itot)
        self.assertEqual(b2.jtot, 2 * b1.jtot)
        self.assertAlmostEqual(b2.xlen, b1.xlen, places=12)
        self.assertAlmostEqual(b2.ylen, b1.ylen, places=12)
        self.assertAlmostEqual(b2.x0, b1.x0, places=12)
        self.assertAlmostEqual(b2.y0, b1.y0, places=12)


class TestBoxMatchesRealParent(unittest.TestCase):
    """The default box must match the real parent's own &NESTPARENT block."""

    @unittest.skipIf(_real_parent_dir() is None,
                     "set UDALES_GMD_PARENT_DIR to the 567 parent case directory")
    def test_box_origin_and_size_match_nestparent_to_1e9(self):
        parent_dir = _real_parent_dir()
        sections = gc.read_namoptions_sections(parent_dir / "namoptions.567")
        nestparent = sections["NESTPARENT"]
        box = gc.ChildBox(refine=1)
        self.assertAlmostEqual(box.x0, float(nestparent["nestparent_x0"]), places=9)
        self.assertAlmostEqual(box.y0, float(nestparent["nestparent_y0"]), places=9)
        self.assertAlmostEqual(box.xlen, float(nestparent["nestparent_xsize"]), places=9)
        self.assertAlmostEqual(box.ylen, float(nestparent["nestparent_ysize"]), places=9)


class TestStlTranslation(unittest.TestCase):
    """The child's STL must carry the enclosure in the CHILD's own coordinates."""

    def test_translation_puts_the_cavity_at_the_expected_child_coordinates(self):
        box = gc.ChildBox(refine=1)
        t = ig.PAPER_WALL
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "shell.stl"
            ig.write_stl(path, t, ground=False, origin=(box.x0, box.y0))
            v = _ascii_stl_vertices(path)
        (xo0, xo1), (yo0, yo1), (_, zo1) = ig.outer_bounds(t)
        self.assertAlmostEqual(float(v[:, 0].min()), xo0 - box.x0, places=6)
        self.assertAlmostEqual(float(v[:, 0].max()), xo1 - box.x0, places=6)
        self.assertAlmostEqual(float(v[:, 1].min()), yo0 - box.y0, places=6)
        self.assertAlmostEqual(float(v[:, 1].max()), yo1 - box.y0, places=6)
        # z is never re-based: the outer envelope's z-range is unaffected by origin.
        self.assertAlmostEqual(float(v[:, 2].min()), 0.0, places=6)
        self.assertAlmostEqual(float(v[:, 2].max()), zo1, places=6)

    def test_default_origin_reproduces_the_untranslated_shell(self):
        t = ig.PAPER_WALL
        with tempfile.TemporaryDirectory() as tmp:
            bare = Path(tmp) / "bare.stl"
            explicit = Path(tmp) / "explicit.stl"
            ig.write_stl(bare, t, ground=False)
            ig.write_stl(explicit, t, ground=False, origin=(0.0, 0.0))
            self.assertEqual(bare.read_text(), explicit.read_text())


class TestCheckWallAtTheGmdResolution(unittest.TestCase):
    """check_wall, on this case's own dx/dy: 7 mm is refused at r = 1, accepted at r = 2."""

    def test_7mm_is_refused_at_r1_and_accepted_at_r2(self):
        dx1, dy1 = gc.parent_dx(), gc.parent_dy()
        with self.assertRaises(ig.WallTooThin):
            ig.check_wall(0.007, dx1, dy1)
        got = ig.check_wall(0.007, dx1 / 2.0, dy1 / 2.0)
        self.assertAlmostEqual(got, 0.007, places=12)

    def test_the_paper_wall_is_accepted_at_r1(self):
        dx1, dy1 = gc.parent_dx(), gc.parent_dy()
        self.assertAlmostEqual(ig.check_wall(ig.PAPER_WALL, dx1, dy1), ig.PAPER_WALL, places=12)


class TestChildSections(unittest.TestCase):
    """The namelist transform: only the documented keys change."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.parent_dir = _write_synthetic_parent(Path(self._tmp.name))
        self.parent_sections = gc.read_namoptions_sections(
            self.parent_dir / f"namoptions.{gc.PARENT_EXPNR}")

    def tearDown(self):
        self._tmp.cleanup()

    def _sections(self, r: int, t: float = ig.PAPER_WALL, nr: str = "601"):
        box = gc.ChildBox(refine=r)
        dtmax_child = float(self.parent_sections["RUN"]["dtmax"]) / r
        return gc.child_sections(self.parent_sections, nr, box, t, dtmax_child,
                                 runtime=150.0, trestart=12.5, nprocx=8, nprocy=9,
                                 stl_name=f"geom.{nr}.stl"), box, dtmax_child

    def test_bc_and_nesting_switches(self):
        sec, box, dtmax_child = self._sections(r=1)
        self.assertEqual(sec["BC"]["BCxm"], 4)
        self.assertEqual(sec["BC"]["BCym"], 3)
        # BCtopm and the roughness/surface keys are the parent's own, untouched.
        self.assertEqual(sec["BC"]["BCtopm"], 3)
        self.assertAlmostEqual(sec["BC"]["z0"], 0.0009, places=12)
        self.assertIs(sec["NESTING"]["lnesting"], True)
        self.assertIs(sec["NESTING"]["nest_linitfromparent"], True)
        self.assertIs(sec["NESTING"]["nest_lendabort"], True)
        self.assertEqual(sec["NESTING"]["nest_shape"], 1)
        self.assertEqual(sec["NESTING"]["nest_timeinterp"], 2)
        # tau = n_tau dt with n_tau = 2, the campaign's ratio (1.0 s at 0.5 s),
        # not the campaign's 1.0 s itself: here dt is 0.8 ms at r = 1.
        self.assertAlmostEqual(sec["NESTING"]["nest_tau"],
                               gc.N_TAU_STEPS * sec["RUN"]["dtmax"], places=12)
        self.assertAlmostEqual(sec["NESTING"]["nest_tau"], 2 * 0.0008, places=12)

    def test_itot_jtot_dtmax_and_zone_widths(self):
        for r in (1, 2):
            with self.subTest(r=r):
                sec, box, dtmax_child = self._sections(r=r)
                self.assertEqual(sec["DOMAIN"]["itot"], box.itot)
                self.assertEqual(sec["DOMAIN"]["jtot"], box.jtot)
                self.assertEqual(sec["DOMAIN"]["ktot"], gc.PARENT_KTOT)
                self.assertAlmostEqual(sec["RUN"]["dtmax"], dtmax_child, places=15)
                self.assertAlmostEqual(sec["RUN"]["dtmax"],
                                       float(self.parent_sections["RUN"]["dtmax"]) / r,
                                       places=15)
                dx_child = gc.parent_dx() / r
                self.assertAlmostEqual(sec["NESTING"]["nest_guardwidth"],
                                       gc.GUARD_CELLS * dx_child, places=12)
                self.assertAlmostEqual(sec["NESTING"]["nest_zonewidth"],
                                       gc.ZONE_CELLS * dx_child, places=12)

    def test_nest_lparentgeom_only_at_the_paper_wall(self):
        sec_paper, _, _ = self._sections(r=1, t=ig.PAPER_WALL)
        sec_thin, _, _ = self._sections(r=2, t=0.010)
        self.assertIs(sec_paper["NESTING"]["nest_lparentgeom"], True)
        self.assertIs(sec_thin["NESTING"]["nest_lparentgeom"], False)

    def test_driver_is_reduced_to_idriver_off(self):
        sec, _, _ = self._sections(r=1)
        self.assertEqual(list(sec["DRIVER"].items()), [("idriver", 0)])

    def test_nestparent_is_dropped(self):
        sec, _, _ = self._sections(r=1)
        self.assertNotIn("NESTPARENT", sec)

    def test_walls_keeps_iwallmom_but_drops_the_parents_counts(self):
        sec, _, _ = self._sections(r=1)
        self.assertEqual(sec["WALLS"]["iwallmom"], 1)
        for key in sec["WALLS"]:
            self.assertFalse(key.startswith(("nfcts", "nsolpts", "nbndpts", "nfctsecs")),
                             f"{key} should have been dropped; counts come from preprocessing")

    def test_inps_stretching_is_unchanged_only_geometry_keys_move(self):
        sec, box, _ = self._sections(r=1, nr="601")
        for key in ("zsize", "lzstretch", "hlin", "dzlin", "lstretchexp", "stretchconst", "u0"):
            self.assertEqual(sec["INPS"][key], self.parent_sections["INPS"][key], key)
        self.assertEqual(sec["INPS"]["stl_file"], "geom.601.stl")
        self.assertIs(sec["INPS"]["stl_ground"], True)


class TestEndToEnd(unittest.TestCase):
    """The full build (preprocessing + nesting file) against the real, run parent."""

    @unittest.skipIf(_real_parent_dir() is None,
                     "set UDALES_GMD_PARENT_DIR to a 567 parent case directory that has run "
                     "(so nesting.out.*.567.nc band files exist)")
    def test_builds_the_full_case_from_the_real_parent(self):
        parent_dir = _real_parent_dir()
        with tempfile.TemporaryDirectory() as tmp:
            casedir = gc.build(parent_dir, Path(tmp), "901", refine=1)
            self.assertTrue((casedir / "nesting.inp.901.nc").exists())
            self.assertTrue((casedir / "manifest.json").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
