#!/usr/bin/env python3
"""Smoke test: drive the whole V2 sweep at the ``tiny`` preset.

Same contract as ``test_v1_tiny.py`` and for the same reason.  This makes no
physical claim -- a 96 x 96 x 32 parent run for 240 s cannot say whether the
resolved-TKE deficit tracks the zone width or the fetch.  It tests the
*harness*: that a sweep whose points share one parent really can share it, that
the reference child is reused rather than silently re-run, that a child whose
zone contains buildings runs with ``nest_lparentgeom = .true.`` instead of
aborting, that the common comparison block is the same physical region for every
point, and that the summary table comes out with both arms populated and finite.

The tiny and production sweeps differ only by the :class:`config.Sweep` object,
so keeping this green keeps the production campaign runnable -- and the
production sweep's own consistency is checked here too, without running it.

Run directly, or through ``tests/run_tests.py nesting-validation``:

    module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
    source /rds/general/user/mvr/home/udales/.venv/bin/activate
    python tests/validation/nesting/test_v2_tiny.py

Environment: ``UDALES_BUILD`` (solver binary), ``UDALES_RUNTIME_MODULES``,
``MPIEXEC``, ``UDALES_V2_RUNDIR`` (where to work; kept if set),
``UDALES_V2_KEEP=1`` (keep a temporary run directory too).
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

import caselib
import make_child_case
import make_parent_case
import run_v1
import run_v2
from caselib import run_solver
from config import CONVERGED, get_preset, get_sweep


def _solid_mask(path: Path, shape) -> np.ndarray:
    idx = np.loadtxt(path, comments="#", dtype=int, ndmin=2)
    m = np.zeros(shape, dtype=bool)
    m[idx[:, 0] - 1, idx[:, 1] - 1, idx[:, 2] - 1] = True
    return m


class TestSweepConfiguration(unittest.TestCase):
    """Everything that can be checked without running anything.

    These are the checks that protect the reuse of the V1 parent, which is what
    the whole campaign is costed against.  They are cheap, so they run whether
    or not a solver binary is available.
    """

    def test_both_sweeps_are_self_consistent(self):
        for name in ("v2", "v2-tiny"):
            get_sweep(name)  # raises with every problem listed

    def test_every_production_point_regenerates_the_v1_parent(self):
        """The one invariant that makes reusing the parent's dumps legitimate."""
        sweep = get_sweep("v2")
        want = CONVERGED.cube_centres()
        for pt in sweep.points:
            got = pt.preset.cube_centres()
            self.assertEqual(got.shape, want.shape, f"{pt.key}: cube count")
            self.assertTrue(np.array_equal(np.sort(got, axis=0), np.sort(want, axis=0)),
                            f"{pt.key}: cube layout differs from the V1 parent's")
            self.assertEqual(pt.preset.plaza_window, CONVERGED.plaza_window)

    def test_the_reference_point_is_reused_and_shared(self):
        sweep = get_sweep("v2")
        ref = sweep.point("ref")
        self.assertTrue(ref.reuse)
        self.assertEqual(set(ref.arms), {"zone", "size"})
        self.assertIs(ref.preset, CONVERGED)
        self.assertEqual(len(sweep.to_run), len(sweep.points) - 1)

    def test_expnrs_are_unique(self):
        for name in ("v2", "v2-tiny"):
            nrs = [p.preset.child_expnr for p in get_sweep(name).points]
            self.assertEqual(len(nrs), len(set(nrs)), f"{name}: {nrs}")
            self.assertNotIn(get_sweep(name).parent.parent_expnr, nrs)

    def test_nzone_covers_the_zone_at_every_point(self):
        for name in ("v2", "v2-tiny"):
            for pt in get_sweep(name).points:
                p = pt.preset
                self.assertGreaterEqual(p.nzone, p.zone_cells, f"{name}/{pt.key}")

    def test_the_common_block_is_the_same_physical_region(self):
        """Otherwise the 'matched window' column would compare different places."""
        import analyse
        for name in ("v2", "v2-tiny"):
            sweep = get_sweep(name)
            boxes = set()
            for pt in sweep.points:
                p = pt.preset
                ii, jj = analyse.central_indices(p, sweep.common_block_cells)
                boxes.add((round(p.child_origin[0] + ii[0] * p.dx, 6),
                           round(p.child_origin[1] + jj[0] * p.dy, 6),
                           ii.size, jj.size))
            self.assertEqual(len(boxes), 1,
                             f"{name}: the common block lands in {len(boxes)} "
                             f"different places: {boxes}")

    def test_a_point_that_cannot_share_the_parent_is_rejected(self):
        """The guard has to actually fire, or it is decoration."""
        from dataclasses import replace
        from config import Sweep, SweepPoint

        base = get_preset("tiny")
        bad = replace(base, name="tiny-bad", child_expnr="905", itot=64, jtot=64,
                      nprocx=2, nprocy=2, plaza=base.plaza_window)
        sweep = Sweep(name="bad", parent=base, points=(
            SweepPoint("ref", ("zone", "size"), base, reuse=True),
            SweepPoint("bad", ("zone",), bad)))
        with self.assertRaises(ValueError) as ctx:
            sweep.validate()
        self.assertIn("itot", str(ctx.exception))

    def test_every_point_has_a_building_free_zone(self):
        """The whole sweep runs one boundary treatment, so only one thing moves."""
        for name in ("v2", "v2-tiny"):
            for pt in get_sweep(name).points:
                self.assertTrue(pt.preset.building_free_zone,
                                f"{name}/{pt.key}: zone is not building-free")

    def test_the_size_arm_clears_its_zone_rather_than_inheriting_it(self):
        """The size arm only works because child geometry may differ from parent."""
        sweep = get_sweep("v2")
        # The reference child is the one the plaza was carved for: nothing to
        # clear, and it is an exact sub-model of the parent.
        self.assertEqual(sweep.point("ref").preset.n_child_cubes_removed, 0)
        for pt in sweep.arm("zone"):
            self.assertEqual(pt.preset.n_child_cubes_removed, 0,
                             f"{pt.key}: the zone arm should need no clearing")
        # The smaller children do: the parent has cubes where their zones land.
        for key in ("size64", "size96"):
            p = sweep.point(key).preset
            self.assertGreater(p.n_child_cubes_removed, 0, key)
            self.assertLess(p.building_clearance_available, p.zone_clearance,
                            f"{key}: nothing would have needed clearing")
            self.assertTrue(p.clear_child_zone, key)

    def test_no_cleared_cube_reaches_the_region_compared(self):
        """The invariant that keeps the size arm interpretable."""
        for name in ("v2", "v2-tiny"):
            for pt in get_sweep(name).points:
                self.assertEqual(
                    len(pt.preset.removed_cubes_reaching_the_interior()), 0,
                    f"{name}/{pt.key}: a cleared cube reaches the analysis interior")

    def test_the_clearing_invariant_actually_fires(self):
        """A guard that has never rejected anything is not a guard.

        ``TINY``'s 14 m zone is shallower than the 24 m a cube sits in from a
        child face in this array, so clearing a small child's zone there would
        also change its interior -- which is exactly why the tiny sweep carries
        the production zone instead.
        """
        from dataclasses import replace
        base = get_preset("tiny")
        bad = replace(base, name="tiny-clear-bad", child_expnr="905",
                      child_itot=32, child_jtot=32,
                      plaza=base.plaza_window, clear_child_zone=True)
        self.assertGreater(len(bad.removed_cubes_reaching_the_interior()), 0)
        with self.assertRaises(ValueError) as ctx:
            bad.validate()
        self.assertIn("also reaches the analysis interior", str(ctx.exception))

    def test_clearing_off_leaves_the_child_as_the_parents_restriction(self):
        """The V1 behaviour, and the nest_lparentgeom = .true. case, still exist."""
        from dataclasses import replace
        p = replace(get_sweep("v2").point("size64").preset,
                    name="size64-uncleared", clear_child_zone=False)
        self.assertFalse(p.building_free_zone)
        self.assertEqual(p.n_child_cubes_removed, 0)
        self.assertEqual(len(p.child_cube_centres()),
                         len(p.cube_centres_in(*p.child_origin,
                                               p.child_xlen, p.child_ylen)))
        # and the cleared version really does carry fewer cubes than that
        cleared = get_sweep("v2").point("size64").preset
        self.assertLess(len(cleared.child_cube_centres()),
                        len(p.child_cube_centres()))


class TestV2TinySweep(unittest.TestCase):
    """One end-to-end sweep, with the assertions spread over the stages."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.sweep = get_sweep("v2-tiny")
        cls.base = cls.sweep.parent
        override = os.environ.get("UDALES_V2_RUNDIR")
        if override:
            cls.rundir = Path(override)
            cls.rundir.mkdir(parents=True, exist_ok=True)
            cls._temp = None
        else:
            cls._temp = tempfile.mkdtemp(prefix="udales-v2-tiny-")
            cls.rundir = Path(cls._temp)
        cls.parent_dir = cls.rundir / cls.base.parent_expnr

        binary = caselib.solver_binary()
        if not binary.exists():
            raise unittest.SkipTest(f"solver binary not found at {binary}")

        p = cls.base
        # -- the V1 pipeline, which produces the parent and the reference child #
        make_parent_case.build(cls.rundir, p)
        run_solver(cls.parent_dir, f"namoptions_spinup.{p.parent_expnr}",
                   p.nprocx * p.nprocy, cls.parent_dir / "spinup.log")
        run_v1._set_startfile(cls.parent_dir / f"namoptions.{p.parent_expnr}",
                              run_v1._restart_file(cls.parent_dir, p.parent_expnr))
        run_solver(cls.parent_dir, f"namoptions.{p.parent_expnr}",
                   p.nprocx * p.nprocy, cls.parent_dir / "production.log")
        make_child_case.build(cls.parent_dir, cls.rundir, p)
        run_solver(cls.rundir / p.child_expnr, f"namoptions.{p.child_expnr}",
                   p.child_nprocx * p.child_nprocy,
                   cls.rundir / p.child_expnr / "child.log")

        # -- the sweep, through its real driver ------------------------------ #
        # Plots left on deliberately: analyse._plots and sweep_summary._plots are
        # both wrapped in try/except, so a broken figure degrades to a printed
        # warning in a four-hour production job rather than an error.  The only
        # way that gets noticed is a test that asserts the files appear.
        argv = sys.argv
        sys.argv = ["run_v2.py", str(cls.rundir), "--sweep", cls.sweep.name,
                    "--parent-dir", str(cls.parent_dir)]
        try:
            cls.rc = run_v2.main()
        finally:
            sys.argv = argv
        cls.outdir = cls.rundir / "analysis"
        cls.metrics = {
            pt.key: json.loads(
                (cls.outdir / pt.key / run_v2.METRICS_NAME).read_text())
            for pt in cls.sweep.points}
        cls.summary = json.loads((cls.outdir / "sweep_summary.json").read_text())

    @classmethod
    def tearDownClass(cls) -> None:
        if cls._temp and not os.environ.get("UDALES_V2_KEEP"):
            shutil.rmtree(cls._temp, ignore_errors=True)
        elif cls._temp:
            print(f"run directory kept at {cls._temp}")

    def _casedir(self, key: str) -> Path:
        return self.rundir / self.sweep.point(key).preset.child_expnr

    # -- the driver --------------------------------------------------------- #

    def test_the_sweep_ran_every_point(self):
        self.assertEqual(self.rc, 0)
        self.assertEqual(len(self.summary["rows"]), len(self.sweep.points))
        self.assertEqual({r["key"] for r in self.summary["rows"]},
                         {p.key for p in self.sweep.points})

    def test_the_reference_child_was_reused_not_rerun(self):
        """If it had been re-run there would be a second case directory for it."""
        ref = self.sweep.point("ref")
        self.assertTrue(ref.reuse)
        casedir = self.rundir / ref.preset.child_expnr
        # exactly one directory, the one the V1 stage built
        self.assertTrue((casedir / "manifest.json").exists())
        row = next(r for r in self.summary["rows"] if r["key"] == "ref")
        self.assertTrue(row["reused"])

    def test_the_reuse_guard_rejects_the_wrong_run(self):
        """Reuse is where a paired comparison could silently stop being paired."""
        ref = self.sweep.point("ref")
        other = next(p for p in self.sweep.points if not p.reuse)
        with self.assertRaises(RuntimeError):
            run_v2._check_reused(ref, self._casedir(other.key))
        with self.assertRaises(FileNotFoundError):
            run_v2._check_reused(ref, self.rundir / "no-such-directory")

    def test_the_parent_statistics_were_shared_where_the_window_is_shared(self):
        """The zone arm must be measured against one parent accumulation."""
        same_window = [p.key for p in self.sweep.points
                       if (p.preset.child_itot, p.preset.child_jtot)
                       == (self.base.child_itot, self.base.child_jtot)]
        reused = [k for k in same_window if self.metrics[k]["parent_bundle_reused"]]
        self.assertEqual(len(reused), len(same_window) - 1,
                         "the parent sub-region was re-accumulated needlessly")

    # -- the cases the sweep built ------------------------------------------ #

    def test_every_child_matches_the_parent_outside_its_zone(self):
        """The child may differ from the parent, but only inside the band.

        Two assertions, and the pair is the point: **zero** solid cells in the
        child's guard + ramp band, so ``nesting_init``'s ``nest_lparentgeom =
        .false.`` assertion has something real to pass; and solid cells
        **identical** to the parent's everywhere outside it, so the region the
        statistics are taken over is the same geometry in both runs.  Between
        the two lies only the band, where the solution is imposed.
        """
        p0 = self.base
        pm = _solid_mask(self.parent_dir / "solid_c.txt", (p0.itot, p0.jtot, p0.ktot))
        for pt in self.sweep.points:
            p = pt.preset
            cm = _solid_mask(self._casedir(pt.key) / "solid_c.txt",
                             (p.child_itot, p.child_jtot, p.child_ktot))
            sub = pm[p.child_i0:p.child_i0 + p.child_itot,
                     p.child_j0:p.child_j0 + p.child_jtot, :]
            n = p.zone_cells
            band = np.ones((p.child_itot, p.child_jtot), dtype=bool)
            band[n:p.child_itot - n, n:p.child_jtot - n] = False
            self.assertEqual(int(cm[band].sum()), 0,
                             f"{pt.key}: solid cells in the child's guard + ramp band")
            self.assertTrue(np.array_equal(sub[~band], cm[~band]),
                            f"{pt.key}: outside the band the child's solid cells "
                            "differ from the parent's")
            self.assertGreater(int(cm[~band].sum()), 0,
                               f"{pt.key}: no buildings outside the band at all")
            # and the clearing did what the preset said it would
            cleared = int(sub[band].sum()) - int(cm[band].sum())
            self.assertEqual(cleared > 0, p.n_child_cubes_removed > 0, pt.key)

    def test_the_mask_the_statistics_use_excludes_either_runs_buildings(self):
        for pt in self.sweep.points:
            sm = self.metrics[pt.key]["solid_mask"]
            self.assertEqual(sm["solid_in_the_child_only"], 0, pt.key)
            self.assertEqual(sm["solid_in_the_parent_only"] > 0,
                             pt.preset.n_child_cubes_removed > 0, pt.key)
            self.assertGreater(sm["fluid_cells_compared"], 0, pt.key)

    def test_every_nesting_file_is_flux_balanced(self):
        """Balanced in the manifest, and the manifest is what the streamed file holds."""
        from netCDF4 import Dataset
        for pt in self.sweep.points:
            m = json.loads((self._casedir(pt.key) / "manifest.json").read_text())
            after = m["flux_residual_after_correction"]["max_abs_normalised"]
            self.assertLess(after, 1.0e-12, f"{pt.key}: Phi/A = {after:g}")
            # the per-level writer's own tally: every stored level, nothing else
            wd = m["writer_diagnostics"]
            self.assertEqual(wd["ntime"], m["n_parent_levels"], pt.key)
            self.assertEqual(wd["correction"]["residual_max_abs"],
                             m["flux_residual_before_correction"]["max_abs"], pt.key)
            path = self._casedir(pt.key) / f"nesting.inp.{pt.preset.child_expnr}.nc"
            with Dataset(path, "r") as ds:
                self.assertEqual(len(ds.dimensions["time"]), wd["ntime"], pt.key)
                residual = np.asarray(ds.variables["flux_residual"][:], dtype=float)
            self.assertAlmostEqual(float(np.max(np.abs(residual))),
                                   m["flux_residual_after_correction"]["max_abs"],
                                   delta=1.0e-20, msg=pt.key)

    def test_every_child_kept_the_flux_and_divergence_bounded(self):
        for pt in self.sweep.points:
            log = (self._casedir(pt.key) / "child.log").read_text(errors="replace")
            phi = [float(x) for x in re.findall(r"Phi \(norm\) =\s*(\S+)", log)]
            div = [float(x) for x in re.findall(r"divmax, divtot =\s*(\S+)", log)]
            self.assertTrue(phi and div, f"{pt.key}: no nesting diagnostics")
            self.assertLess(max(abs(x) for x in phi), 1.0e-9, pt.key)
            self.assertLess(max(div), 1.0e-10, pt.key)
            for face in ("west", "east", "south", "north"):
                self.assertIn(f"face {face} is forced", log, f"{pt.key}/{face}")

    def test_the_solver_asserts_the_zone_is_clear_at_every_point(self):
        """``nesting_init`` must verify the rule, not warn about it.

        Every point of the sweep runs ``nest_lparentgeom = .false.``, so the
        solver aborts if any solid point has ``W > 0``.  Each child having
        completed is therefore already most of the proof; this checks both ends
        of it -- the namelist really says ``.false.``, and the solver did not
        take the "allowed by nest_lparentgeom" branch.
        """
        marker = "solid points inside the relaxation zone"
        for pt in self.sweep.points:
            p = pt.preset
            self.assertTrue(p.building_free_zone, pt.key)
            log = (self._casedir(pt.key) / "child.log").read_text(errors="replace")
            nml = (self._casedir(pt.key) / f"namoptions.{p.child_expnr}"
                   ).read_text()
            self.assertRegex(nml, r"(?im)^\s*nest_lparentgeom\s*=\s*\.false\.")
            self.assertNotIn(marker, log,
                             f"{pt.key}: solid points in a zone declared clear")

    def test_the_zone_fraction_warning_fires_exactly_where_predicted(self):
        """Expected at the small sizes, and not a fault -- so pin it down."""
        for pt in self.sweep.points:
            log = (self._casedir(pt.key) / "child.log").read_text(errors="replace")
            fired = "zone occupies" in log and "> 15 %" in log
            self.assertEqual(fired, pt.preset.zone_fraction_warns,
                             f"{pt.key}: zone fraction "
                             f"{100 * pt.preset.zone_fraction:.1f} %, warning "
                             f"{'fired' if fired else 'did not fire'}")

    # -- the measurement ---------------------------------------------------- #

    def test_every_point_emits_the_falsification_metrics(self):
        for key, m in self.metrics.items():
            v2 = m["v2"]
            d = v2["tke_deficit"]["above"]
            self.assertIsNotNone(d["mean_relative"], key)
            self.assertTrue(np.isfinite(d["mean_relative"]), key)
            self.assertTrue(np.isfinite(d["mean_spread"]), key)
            self.assertTrue(np.isfinite(
                v2["criterion_a"]["max_interior_umean_error_over_ustar"]), key)
            f = v2["tke_error_vs_fetch"]
            self.assertEqual(len(f["per_face"]), 4, key)
            self.assertTrue(np.isfinite(f["mean_error_at_zone_edge"]), key)
            for band in ("band_16_64m", "band_8_16m"):
                for name, sp in m["spectra"].items():
                    self.assertIn(band, sp["bands"], f"{key}/{name}")
            self.assertTrue(v2["common_block"]["cells"] > 0, key)
            self.assertIsNotNone(
                v2["common_block"]["tke_deficit"]["above"]["mean_relative"], key)

    def test_the_summary_table_has_both_arms_and_the_shared_reference(self):
        text = (self.outdir / "sweep_summary.md").read_text()
        self.assertIn("P1 -- zone width", text)
        self.assertIn("P2 -- child size", text)
        # the reference point is the only row that appears in both arm tables
        self.assertEqual(text.count("| ref "), 2)
        for name in ("sweep_summary.csv", "sweep_summary.json", "sweep_summary.md",
                     "sweep_deficit.png", "sweep_spectra.png"):
            self.assertTrue((self.outdir / name).exists(), name)
        for pt in self.sweep.points:
            for name in ("profiles.png", "tke_deficit.png", "spectra.png",
                         "error_vs_distance.png", "tke_series.png"):
                self.assertTrue((self.outdir / pt.key / name).exists(),
                                f"{pt.key}/{name}")
        v = self.summary["verdict"]
        self.assertGreaterEqual(v["zone_arm"]["n_points"], 2)
        self.assertGreaterEqual(v["size_arm"]["n_points"], 2)
        for arm in ("zone_arm", "size_arm"):
            self.assertTrue(np.isfinite(v[arm]["range_pct"]), arm)
            self.assertEqual(v[arm]["x"], sorted(v[arm]["x"]), arm)

    def test_no_point_is_confounded_and_the_clearing_is_reported(self):
        """Both halves matter: nothing confounded, and the difference disclosed."""
        self.assertEqual(self.summary["verdict"]["confounds"], [])
        self.assertIn("child_geometry_note", self.summary["verdict"])
        cut = {r["key"]: r["cubes_cleared"] for r in self.summary["rows"]}
        for pt in self.sweep.points:
            self.assertEqual(cut[pt.key], pt.preset.n_child_cubes_removed, pt.key)
        self.assertTrue(any(v > 0 for v in cut.values()),
                        "the tiny sweep no longer clears anything, so it stopped "
                        "covering the path the production size arm depends on")

    def test_the_headline_numbers_are_not_wildly_different(self):
        """A loose bound, not a physics claim: it catches a child that diverged."""
        for key, m in self.metrics.items():
            pm = m["profile_metrics"]
            self.assertLess(pm["u_rms_difference_over_ustar"], 2.0, key)
            self.assertLess(pm["tke_rms_difference_over_ustar2"], 5.0, key)


if __name__ == "__main__":
    unittest.main(verbosity=2)
