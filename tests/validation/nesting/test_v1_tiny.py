#!/usr/bin/env python3
"""Smoke test: drive the whole V1 pipeline at the ``tiny`` preset.

This is not a physics test -- a 64 x 64 x 32 parent run for 240 s cannot say
anything about whether nesting reproduces a turbulent boundary layer.  It is a
test of the *harness*: that the parent case builds and runs, that the child is
cut out of its dumps with a flux-balanced nesting file and a geometry identical
to the parent sub-region, that the nested run completes with a bounded
divergence, and that the analysis produces finite numbers.  Because the tiny
and production paths differ only by a :class:`config.Preset`, keeping this green
keeps the production campaign runnable.

Run directly, or through ``tests/run_tests.py nesting-validation``:

    module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
    source /rds/general/user/mvr/home/udales/.venv/bin/activate
    python tests/validation/nesting/test_v1_tiny.py

Environment: ``UDALES_BUILD`` (solver binary), ``UDALES_RUNTIME_MODULES``,
``MPIEXEC``, ``UDALES_V1_RUNDIR`` (where to work; kept if set),
``UDALES_V1_KEEP=1`` (keep a temporary run directory too).
"""

from __future__ import annotations

import json
import os
import re
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import numpy as np

import analyse
import caselib
import make_child_case
import make_parent_case
import run_v1
from caselib import run_solver
from config import get_preset


def _solid_mask(path: Path, shape) -> np.ndarray:
    idx = np.loadtxt(path, comments="#", dtype=int, ndmin=2)
    m = np.zeros(shape, dtype=bool)
    m[idx[:, 0] - 1, idx[:, 1] - 1, idx[:, 2] - 1] = True
    return m


class TestV1TinyPipeline(unittest.TestCase):
    """One end-to-end run, with the assertions spread over the stages."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.preset = get_preset("tiny")
        override = os.environ.get("UDALES_V1_RUNDIR")
        if override:
            cls.rundir = Path(override)
            cls.rundir.mkdir(parents=True, exist_ok=True)
            cls._temp = None
        else:
            cls._temp = tempfile.mkdtemp(prefix="udales-v1-tiny-")
            cls.rundir = Path(cls._temp)
        cls.parent_dir = cls.rundir / cls.preset.parent_expnr
        cls.child_dir = cls.rundir / cls.preset.child_expnr
        cls.outdir = cls.rundir / "analysis"

        binary = caselib.solver_binary()
        if not binary.exists():
            raise unittest.SkipTest(f"solver binary not found at {binary}")

        p = cls.preset
        make_parent_case.build(cls.rundir, p)
        run_solver(cls.parent_dir, f"namoptions_spinup.{p.parent_expnr}",
                   p.nprocx * p.nprocy, cls.parent_dir / "spinup.log")
        run_v1._set_startfile(cls.parent_dir / f"namoptions.{p.parent_expnr}",
                              run_v1._restart_file(cls.parent_dir, p.parent_expnr))
        run_solver(cls.parent_dir, f"namoptions.{p.parent_expnr}",
                   p.nprocx * p.nprocy, cls.parent_dir / "production.log")
        make_child_case.build(cls.parent_dir, cls.rundir, p)
        run_solver(cls.child_dir, f"namoptions.{p.child_expnr}",
                   p.child_nprocx * p.child_nprocy, cls.child_dir / "child.log")
        cls.metrics = analyse.run(cls.parent_dir, cls.child_dir, cls.outdir, p,
                                  make_plots=False)
        cls.manifest = json.loads((cls.child_dir / "manifest.json").read_text())
        cls.child_log = (cls.child_dir / "child.log").read_text(errors="replace")

    @classmethod
    def tearDownClass(cls) -> None:
        if cls._temp and not os.environ.get("UDALES_V1_KEEP"):
            shutil.rmtree(cls._temp, ignore_errors=True)
        elif cls._temp:
            print(f"run directory kept at {cls._temp}")

    # -- the case builders -------------------------------------------------- #

    def test_child_geometry_is_the_parent_sub_region(self):
        """The whole experiment is void if the two geometries differ."""
        p = self.preset
        pm = _solid_mask(self.parent_dir / "solid_c.txt", (p.itot, p.jtot, p.ktot))
        cm = _solid_mask(self.child_dir / "solid_c.txt",
                         (p.child_itot, p.child_jtot, p.child_ktot))
        sub = pm[p.child_i0:p.child_i0 + p.child_itot,
                 p.child_j0:p.child_j0 + p.child_jtot, :]
        self.assertTrue(np.array_equal(sub, cm),
                        "the child's solid cells are not the parent's sub-region")
        self.assertGreater(cm.sum(), 0, "the child case has no buildings at all")

    def test_nesting_file_is_flux_balanced(self):
        after = self.manifest["flux_residual_after_correction"]["max_abs_normalised"]
        before = self.manifest["flux_residual_before_correction"]["max_abs_normalised"]
        self.assertLess(after, 1.0e-12, f"corrected file still carries Phi/A = {after:g}")
        # The correction has to be doing real work, or the test is vacuous: the
        # parent dumps are single precision, so the raw residual sits well above
        # the solver's 1e-10 tolerance.
        self.assertGreater(before, after * 100.0)

    def test_initial_condition_is_projected(self):
        ic = self.manifest["initial_condition_divmax"]
        if ic["before_projection"] is None:
            self.skipTest("this preset does not initialise the child from the parent")
        self.assertLess(ic["after_projection"], 1.0e-12)
        self.assertGreater(ic["before_projection"], ic["after_projection"])

    # -- the nested run ----------------------------------------------------- #

    def test_child_run_kept_the_flux_and_the_divergence_bounded(self):
        phi = [float(m) for m in
               re.findall(r"Phi \(norm\) =\s*(\S+)", self.child_log)]
        div = [float(m) for m in
               re.findall(r"divmax, divtot =\s*(\S+)", self.child_log)]
        self.assertTrue(phi, "the child run printed no nesting flux diagnostics")
        self.assertTrue(div, "the child run printed no divergence diagnostics")
        self.assertLess(max(abs(x) for x in phi), 1.0e-9,
                        "the runtime boundary flux residual is not at round-off")
        self.assertLess(max(div), 1.0e-10,
                        "the projection left a divergence far above round-off")

    def test_all_four_faces_are_forced(self):
        for face in ("west", "east", "south", "north"):
            self.assertIn(f"face {face} is forced", self.child_log)

    def test_the_zone_is_building_free(self):
        """The design section 5 rule, as an assertion rather than a comment.

        With ``geometry = 'plaza'`` the child sets ``nest_lparentgeom = .false.``,
        so ``nesting_init`` *aborts* if any solid point has ``W > 0``.  The run
        having completed is therefore already most of the proof; this checks the
        two ends of it -- that the solver did not merely warn, and that the zone
        really is clear in the geometry the case was built from.
        """
        if not self.preset.building_free_zone:
            self.assertIn("solid points inside the relaxation zone", self.child_log)
            self.skipTest("this preset deliberately puts buildings in the zone")
        self.assertNotIn("solid points inside the relaxation zone", self.child_log,
                         "nesting_init found solid points in a zone that is "
                         "supposed to be building-free")
        p = self.preset
        mask = _solid_mask(self.child_dir / "solid_c.txt",
                           (p.child_itot, p.child_jtot, p.child_ktot))
        n = p.zone_cells
        frame = np.ones((p.child_itot, p.child_jtot), dtype=bool)
        frame[n:p.child_itot - n, n:p.child_jtot - n] = False
        self.assertEqual(int(mask[frame].sum()), 0,
                         "solid cells found in the child's guard + ramp band")
        self.assertGreater(int(mask[~frame].sum()), 0,
                           "the child interior has no buildings; V1 would test nothing")

    def test_uniform_layout_reproduces_create_cubes(self):
        """The private-helper coupling in ``caselib.cube_mesh``, made checkable.

        ``cube_mesh`` calls the two helpers ``udgeom.create_cubes`` itself uses,
        because ``create_cubes`` cannot express the plaza layout.  Feeding it the
        *full* layout must therefore give back exactly what ``create_cubes``
        gives, or the two have drifted.
        """
        from udgeom.geometry_generation import create_cubes

        p = self.preset
        mine = caselib.cube_mesh(p, p._full_cube_centres(), p.xlen, p.ylen)
        theirs = create_cubes(p.xlen, p.ylen, p.building_width, p.building_width,
                              p.building_height, p.street_width, p.street_width,
                              "AC", p.edgelength)
        self.assertTrue(np.array_equal(mine.stl.vertices, theirs.stl.vertices))
        self.assertTrue(np.array_equal(mine.stl.faces, theirs.stl.faces))

    # -- the analysis ------------------------------------------------------- #

    def test_analysis_emits_machine_readable_numbers(self):
        for name in ("v1_metrics.json", "profiles.csv",
                     "error_vs_distance_x_umean.csv", "spectrum_z_17m.csv"):
            self.assertTrue((self.outdir / name).exists(), f"{name} was not written")
        m = self.metrics["profile_metrics"]
        for key, value in m.items():
            self.assertTrue(np.isfinite(value), f"{key} is not finite")
        self.assertGreaterEqual(self.metrics["samples"]["parent"], 4)
        self.assertGreaterEqual(self.metrics["samples"]["child"], 4)

    def test_interior_profiles_are_not_wildly_different(self):
        """A loose sanity bound, not a physics claim.

        The tiny run is far too short to say anything about accuracy; this only
        catches a child that has diverged, blown up or been driven by the wrong
        data.  The production preset is where the numbers mean something.
        """
        m = self.metrics["profile_metrics"]
        self.assertLess(m["u_rms_difference_over_ustar"], 2.0)
        self.assertLess(m["tke_rms_difference_over_ustar2"], 5.0)

    def test_error_curves_cover_the_whole_span(self):
        p = self.preset
        c = self.metrics["error_vs_distance"]["x_umean"]
        self.assertEqual(len(c["error"]), p.child_itot - 1)
        self.assertTrue(np.all(np.isfinite(c["error"])))
        self.assertTrue(np.all(np.isfinite(c["noise_floor"])))


if __name__ == "__main__":
    unittest.main(verbosity=2)
