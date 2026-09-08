#!/usr/bin/env python3
"""Smoke test: drive the V0c matched-control-plus-refinement suite at ``tiny``.

Same contract as ``test_v0_tiny.py``/``test_v0b_tiny.py``: no physical claim --
a tiny, few-minutes run cannot say whether criterion A is genuinely violated at
r > 1, only that the harness builds and runs the r = 1 matched-resolution
control alongside r = 2 and r = 4 through one shared code path, and that the
new pressure diagnostics and the profile-RMS/criterion-A separation the review
asked for (nesting-review-2026-09-08-codex.md, finding 1) actually show up in
the metrics.  What this covers, specific to V0c:

* three points build and run to completion: ``r1`` (the matched control),
  ``r2``, ``r4``, all ``arm = "filtered"``, all with ``prolongation`` left
  unset (the writer's shipped default applies uniformly -- V0c does not sweep
  the reconstruction again, V0b already settled that);
* ``r1`` is a plain slab CUT of the reference's own dumps
  (``manifest["refinement"]["slabs"]`` names ``"cut"``), while ``r2``/``r4``
  are interpolated -- ``RefinedPoint.coarsen == 1`` degenerates
  ``CoarsenedFieldDump`` to the identity, so the writer never reaches the
  prolongation for the control, and the offline divergence identity therefore
  holds for ``r1`` exactly, regardless of which reconstruction is the default;
* the pressure diagnostic (``analyse_v0.runtime_diagnostics``) now reports a
  RESOLVED time series (``nest_statint`` is a real, finite interval on this
  suite, not the inherited ``tstatsdump`` default V0b used) with the startup
  discarded at ``child_spinup``, and both a mean AND a spread -- review finding
  2 -- rather than V0b's two compulsory endpoints;
* the summary table and the per-point metrics carry the mean-flow **profile
  RMS** (``u_rms_difference_over_ustar``) and criterion A (the **maximum
  interior slab RMS**) as two separate columns/keys, not one number doing
  double duty -- review finding 1, `analyse.py:521` and `:692`.

Run directly, or through ``tests/run_tests.py nesting-validation``:

    module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
    source ~/udales/.venv/bin/activate
    python tests/validation/nesting/test_v0c_tiny.py

Environment: ``UDALES_BUILD`` (solver binary), ``UDALES_RUNTIME_MODULES``,
``MPIEXEC``, ``UDALES_V0C_RUNDIR`` (where to work; kept if set),
``UDALES_V0C_KEEP=1`` (keep a temporary run directory too).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import numpy as np

import caselib
from config import CONVERGED, V0C, V0C_FINE, get_suite


# --------------------------------------------------------------------------- #
# Checks that need no run
# --------------------------------------------------------------------------- #


class TestV0cSuiteConfiguration(unittest.TestCase):
    """Everything checkable without a solver, including the production suite."""

    def test_both_suites_are_self_consistent(self):
        for name in ("v0c", "v0c-tiny"):
            get_suite(name)  # raises with every problem listed

    def test_matched_control_present_and_no_coarse_arm(self):
        for name in ("v0c", "v0c-tiny"):
            suite = get_suite(name)
            self.assertFalse(suite.require_paired_arms)
            self.assertTrue(suite.allow_matched_control)
            self.assertTrue(all(p.arm == "filtered" for p in suite.points))
            self.assertEqual(suite.refinements, (1, 2, 4))

    def test_prolongation_is_left_unset_everywhere(self):
        """V0c does not sweep the reconstruction again; V0b already settled it."""
        for name in ("v0c", "v0c-tiny"):
            for pt in get_suite(name).points:
                self.assertIsNone(pt.child.prolongation, pt.key)

    def test_the_control_is_a_degenerate_filtered_point(self):
        """r = 1's 'coarsen' is 1: CoarsenedFieldDump(dump, 1) is the identity."""
        for name in ("v0c", "v0c-tiny"):
            suite = get_suite(name)
            r1 = suite.point("r1")
            self.assertEqual(r1.refine, 1)
            self.assertEqual(r1.coarsen, 1)
            self.assertAlmostEqual(r1.driver.dx, r1.child.dx, places=12)

    def test_r2_and_r4_are_genuinely_coarsened(self):
        for name in ("v0c", "v0c-tiny"):
            suite = get_suite(name)
            for refine in (2, 4):
                pt = suite.point(f"r{refine}")
                self.assertEqual(pt.coarsen, refine)
                self.assertAlmostEqual(pt.driver.dx, refine * pt.child.dx, places=12)

    def test_the_window_is_v1_converged_s_own(self):
        """~5.7x V0b's 1800 s window (roughly halves the sampling floor)."""
        self.assertEqual(V0C_FINE.production, CONVERGED.production)
        self.assertEqual(V0C_FINE.child_spinup, CONVERGED.child_spinup)
        window = V0C_FINE.production - V0C_FINE.child_spinup
        self.assertGreaterEqual(window / 1800.0, 5.0)

    def test_the_final_configuration_is_used(self):
        """The FINAL cadence/interpolant (0.5 s, Catmull-Rom), same as V0b."""
        for name in ("v0c", "v0c-tiny"):
            for pt in get_suite(name).points:
                self.assertEqual(pt.child.cadence, 0.5, pt.key)
                self.assertEqual(pt.child.timeinterp, 2, pt.key)
                self.assertLessEqual(pt.child.c_dump_u0, 2.0, pt.key)

    def test_nest_statint_is_a_real_resolved_interval(self):
        """Review finding 2: not the inherited tstatsdump default V0b used."""
        for name in ("v0c", "v0c-tiny"):
            suite = get_suite(name)
            for pt in suite.points:
                p = pt.child
                self.assertIsNotNone(p.nest_statint, pt.key)
                self.assertGreater(p.nest_statint, 0.0, pt.key)
                n_reports = p.production / p.nest_statint
                n_post_discard = (p.production - p.child_spinup) / p.nest_statint
                self.assertGreater(n_reports, 10, pt.key)
                self.assertGreater(n_post_discard, 5, pt.key)

    def test_expnrs_do_not_collide(self):
        for name in ("v0c", "v0c-tiny"):
            suite = get_suite(name)
            nrs = [p.expnr for p in suite.points] + [suite.reference.parent_expnr]
            self.assertEqual(len(nrs), len(set(nrs)),
                             f"{name}: colliding experiment numbers among {nrs}")
            # the r = 2 / r = 4 driver placeholders are never built, and the
            # r = 1 'driver' reuses the reference's own expnr on purpose.
            self.assertEqual(suite.drivers_to_run, [])

    def test_child_expnrs_match_the_register_comment(self):
        want = {"r1": "983", "r2": "985", "r4": "987"}
        got = {p.key: p.expnr for p in V0C.points}
        self.assertEqual(got, want)

    # -- registration --------------------------------------------------- #

    @staticmethod
    def _groups():
        tests_dir = HERE.parents[1]
        if str(tests_dir) not in sys.path:
            sys.path.insert(0, str(tests_dir))
        import run_tests

        manifest = run_tests._load_manifest()
        return manifest, run_tests

    def test_v0c_is_in_nesting_validation_and_nowhere_else(self):
        manifest, run_tests = self._groups()
        labels = {g: [s["label"] for s in run_tests._expand_groups(manifest, g)]
                  for g in manifest["groups"]}
        v0c = [l for l in labels["nesting-validation"] if "V0c" in l]
        self.assertEqual(len(v0c), 1, f"expected one V0c entry, got {v0c}")
        for group in ("all", "supported", "nesting", "experimental",
                      "python-library", "supported-macos", "lint"):
            self.assertFalse([l for l in labels.get(group, []) if "V0c" in l],
                             f"group '{group}' reaches a V0c entry")


# --------------------------------------------------------------------------- #
# The end-to-end run
# --------------------------------------------------------------------------- #


class TestV0cTinyPipeline(unittest.TestCase):
    """One run of the whole tiny V0c suite, with the assertions spread over it."""

    @classmethod
    def setUpClass(cls) -> None:
        binary = caselib.solver_binary()
        if not binary.exists():
            raise unittest.SkipTest(f"solver binary not found at {binary}")
        cls.suite = get_suite("v0c-tiny")
        override = os.environ.get("UDALES_V0C_RUNDIR")
        if override:
            cls.rundir = Path(override)
            cls.rundir.mkdir(parents=True, exist_ok=True)
            cls._temp = None
        else:
            cls._temp = tempfile.mkdtemp(prefix="udales-v0c-tiny-")
            cls.rundir = Path(cls._temp)
        env = dict(os.environ)
        env.setdefault("MPLCONFIGDIR", str(cls.rundir / ".mpl"))
        completed = subprocess.run(
            [sys.executable, str(HERE / "run_v0.py"), str(cls.rundir),
             "--suite", "v0c-tiny", "--no-plots"],
            cwd=str(HERE), env=env, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
        cls.driver_output = completed.stdout.decode(errors="replace")
        (cls.rundir / "run_v0.log").write_text(cls.driver_output, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError("run_v0.py failed; last lines:\n"
                               + "\n".join(cls.driver_output.splitlines()[-40:]))
        cls.summary = json.loads(
            (cls.rundir / "analysis" / "v0_summary.json").read_text())
        cls.metrics = {}
        cls.manifests = {}
        for pt in cls.suite.points:
            cls.metrics[pt.key] = json.loads(
                (cls.rundir / "analysis" / pt.key / "v0_metrics.json").read_text())
            cls.manifests[pt.key] = json.loads(
                (cls.rundir / pt.expnr / "manifest.json").read_text())

    @classmethod
    def tearDownClass(cls) -> None:
        if cls._temp and not os.environ.get("UDALES_V0C_KEEP"):
            shutil.rmtree(cls._temp, ignore_errors=True)
        elif cls._temp:
            print(f"run directory kept at {cls._temp}")

    # -- all three points built and ran ---------------------------------- #

    def test_all_three_points_built_and_ran(self):
        for pt in self.suite.points:
            self.assertTrue((self.rundir / pt.expnr / "manifest.json").exists(), pt.key)
            self.assertTrue((self.rundir / pt.expnr / "child.log").exists(), pt.key)
            r = self.manifests[pt.key]["refinement"]
            self.assertEqual(r["spatial"], pt.refine, pt.key)

    def test_r1_is_a_slab_cut_r2_r4_are_interpolated(self):
        self.assertIn("cut", self.manifests["r1"]["refinement"]["slabs"])
        self.assertIn("interpolated", self.manifests["r2"]["refinement"]["slabs"])
        self.assertIn("interpolated", self.manifests["r4"]["refinement"]["slabs"])

    def test_prolongation_was_never_requested_but_r1_is_inert_to_it_anyway(self):
        for pt in self.suite.points:
            m = self.manifests[pt.key]
            self.assertIsNone(m["prolongation_requested"], pt.key)

    # -- the control's divergence identity holds exactly, unconditionally - #

    def test_r1_never_reaches_the_prolongation(self):
        """No coarsening, no interpolation: the writer never touches the target.

        ``make_child_case.build`` only computes a separate "parent before
        prolongation" divmax inside ``if driving.interpolates:`` -- at
        ``refine = 1`` the cut block IS the stored initial condition, so
        there is nothing to compare it against and the field is correctly
        left ``None`` rather than a redundant copy of ``before_projection``.
        The identity design section 1.3 states for ``r > 1`` holds here
        trivially, by construction, not by measurement -- which is why r = 2
        and r = 4 (interpolated) are where that identity is actually
        exercised, not the control.
        """
        ic = self.manifests["r1"]["initial_condition_divmax"]
        self.assertIsNone(ic["parent_before_prolongation"])
        self.assertIsNotNone(ic["before_projection"])
        self.assertTrue(np.isfinite(ic["before_projection"]))

    def test_r2_and_r4_do_reach_the_prolongation(self):
        """The identity design section 1.3 states is actually exercised here."""
        for refine in (2, 4):
            ic = self.manifests[f"r{refine}"]["initial_condition_divmax"]
            self.assertIsNotNone(ic["parent_before_prolongation"], refine)
            self.assertIsNotNone(ic["before_projection"], refine)
            self.assertTrue(np.isfinite(ic["before_projection"]), refine)

    def test_the_nesting_file_is_flux_balanced_at_every_ratio(self):
        for pt in self.suite.points:
            after = self.manifests[pt.key]["flux_residual_after_correction"]["max_abs_normalised"]
            self.assertLess(after, 1.0e-12, f"{pt.key}: Phi/A = {after:g}")

    # -- the pressure diagnostic is a resolved, startup-discarded series -- #

    def test_the_pressure_response_is_a_resolved_series(self):
        """Review finding 2: more than V0b's two compulsory endpoints."""
        for pt in self.suite.points:
            rt = self.metrics[pt.key]["v0"]["runtime"]
            self.assertGreater(rt["gradp_zone"]["n"], 2,
                               f"{pt.key}: only {rt['gradp_zone']['n']} reports; "
                               "nest_statint is not giving a resolved series")
            self.assertTrue(np.isfinite(rt["gradp_zone"]["mean"]), pt.key)
            self.assertTrue(np.isfinite(rt["gradp_zone"]["std"]), pt.key)
            self.assertTrue(np.isfinite(rt["gradp_interior"]["mean"]), pt.key)
            self.assertTrue(np.isfinite(rt["gradp_interior"]["std"]), pt.key)

    def test_the_startup_is_discarded_consistently_with_the_stats_window(self):
        for pt in self.suite.points:
            rt = self.metrics[pt.key]["v0"]["runtime"]
            self.assertEqual(rt["gradp_stats_start"], pt.child.child_spinup, pt.key)
            post = rt["gradp_zone_post_startup"]
            self.assertGreater(post["n"], 0, pt.key)
            self.assertLess(post["n"], rt["gradp_zone"]["n"],
                            f"{pt.key}: the startup discard should drop at least one report")
            self.assertEqual(
                post["n"] + rt["gradp_n_discarded_as_startup"], rt["gradp_zone"]["n"], pt.key)
            times = np.asarray(rt["gradp_report_times"])
            self.assertEqual(int((times >= pt.child.child_spinup).sum()), post["n"], pt.key)

    # -- profile RMS and criterion A are separate, both present ----------- #

    def test_profile_rms_and_criterion_a_are_reported_separately(self):
        """Review finding 1: these are different statistics (analyse.py:521, :692).

        ``profile_metrics.u_rms_difference_over_ustar`` is the RMS, over
        height, of the interior- and area-averaged child-minus-truth mean
        profile; ``criterion_a`` is the MAXIMUM, over interior fetch stations
        and faces, of the RMS taken along a single slab
        (``analyse.slab_rms_difference`` / ``decay_length``).  A max-of-slab-
        RMS is not a simple function of an area-averaged RMS (different
        reductions commute differently under spatial cancellation), so this
        only asserts what the review actually asked for: both are present,
        finite, and kept under distinct keys rather than one number silently
        standing in for the other.  Every production point measured so far
        (V1, C0c, V0b r2/r4) has criterion_a several times profile_rms; that
        is not asserted here as a law, only reported for a human to notice if
        a tiny, noisy run ever breaks the pattern.
        """
        for pt in self.suite.points:
            m = self.metrics[pt.key]
            profile_rms = m["profile_metrics"]["u_rms_difference_over_ustar"]
            crit_a = m["v2"]["criterion_a"]["max_interior_umean_error_over_ustar"]
            self.assertTrue(np.isfinite(profile_rms), pt.key)
            self.assertTrue(np.isfinite(crit_a), pt.key)
            print(f"    {pt.key}: profile RMS = {profile_rms:.4f} u*, "
                  f"criterion A = {crit_a:.4f} u*")

    def test_summary_table_carries_both_statistics_as_separate_columns(self):
        for row in self.summary["rows"]:
            self.assertIn("u_rms_difference_over_ustar", row)
            self.assertIn("criterion_a", row)
            self.assertIsNotNone(row["u_rms_difference_over_ustar"], row["key"])
            self.assertIsNotNone(row["criterion_a"], row["key"])
        for name in ("v0_summary.json", "v0_summary.csv", "v0_summary.md"):
            self.assertTrue((self.rundir / "analysis" / name).exists(), name)
        csv_header = (self.rundir / "analysis" / "v0_summary.csv").read_text().splitlines()[0]
        self.assertIn("u_rms_difference_over_ustar", csv_header)
        self.assertIn("criterion_a", csv_header)


# --------------------------------------------------------------------------- #
# Register the tests to run when invoked directly
# --------------------------------------------------------------------------- #


if __name__ == "__main__":
    unittest.main(verbosity=2)
