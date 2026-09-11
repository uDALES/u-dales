#!/usr/bin/env python3
"""Smoke test: drive V0c16 -- V0c's own suite, decomposed for 16 ranks (4x4).

V0c (``submit_cx3_v0c.pbs``, job 4004496, ``ncpus=64:mem=300gb``) sat queued
in ``v1_medium24`` for hours; the 16-core queues turn over far better.  V0c16
is a SEPARATE preset/suite (``config.V0C16``/``V0C16_FINE``, fresh experiment
numbers 992-997) rather than a mutation of V0c's own objects, specifically so
4004496 keeps working if it ever starts.  This test is V0c's own
``test_v0c_tiny.py`` with the suite name and a few numbers changed: the
physics claims (matched control at r = 1, box-filtered r = 2/r = 4, the
pressure diagnostic, profile RMS vs criterion A) are already covered there
and are not re-litigated here.  What THIS test covers, specific to V0c16:

* the suite actually runs 16 ranks (4x4) end to end, on both the fine-truth
  parent (256 x 256 x 64 in production, this test's own tiny 96 x 96 x 32)
  and the children (128 x 128 x 64 in production, tiny 64 x 64 x 32) --
  V0c16_FINE_TINY sets ``nprocx = nprocy = child_nprocx = child_nprocy = 4``
  rather than inheriting TINY's own default 2x2, so this run exercises the
  real decomposition the production job ships, not a smaller one;
* V0c16's experiment numbers (992-997) do not collide with V0c's own
  (982-989) or any other suite's, so the two can share a login node /
  ``$EPHEMERAL`` filesystem (though not the same run directory -- see
  ``submit_cx3_v0c16.pbs``'s header) without stepping on each other;
* the three points build and run to completion under this decomposition, the
  nesting file is flux-balanced at every ratio, and the summary/metrics files
  land where the production job expects them.

Run directly, or through ``tests/run_tests.py nesting-validation``:

    module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
    source ~/udales/.venv/bin/activate
    python tests/validation/nesting/test_v0c16_tiny.py

Environment: ``UDALES_BUILD`` (solver binary), ``UDALES_RUNTIME_MODULES``,
``MPIEXEC``, ``UDALES_V0C16_RUNDIR`` (where to work; kept if set),
``UDALES_V0C16_KEEP=1`` (keep a temporary run directory too).
"""

from __future__ import annotations

import json
import os
import re
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
from config import V0C, V0C16, get_suite


# --------------------------------------------------------------------------- #
# Checks that need no run
# --------------------------------------------------------------------------- #


class TestV0c16SuiteConfiguration(unittest.TestCase):
    """Everything checkable without a solver, including the production suite."""

    def test_both_suites_are_self_consistent(self):
        for name in ("v0c16", "v0c16-tiny"):
            get_suite(name)  # raises with every problem listed

    def test_matched_control_present_and_no_coarse_arm(self):
        for name in ("v0c16", "v0c16-tiny"):
            suite = get_suite(name)
            self.assertFalse(suite.require_paired_arms)
            self.assertTrue(suite.allow_matched_control)
            self.assertTrue(all(p.arm == "filtered" for p in suite.points))
            self.assertEqual(suite.refinements, (1, 2, 4))

    def test_decomposition_is_4x4_everywhere(self):
        """16 ranks: 4x4, not V0c's 8x8 -- both itot and child_itot divide by 4."""
        for name in ("v0c16", "v0c16-tiny"):
            suite = get_suite(name)
            ref = suite.reference
            for attr in ("nprocx", "nprocy", "child_nprocx", "child_nprocy"):
                self.assertEqual(getattr(ref, attr), 4, f"{name}.{attr}")
            self.assertEqual(ref.itot % ref.nprocx, 0, name)
            self.assertEqual(ref.jtot % ref.nprocy, 0, name)
            self.assertEqual(ref.child_itot % ref.child_nprocx, 0, name)
            self.assertEqual(ref.child_jtot % ref.child_nprocy, 0, name)

    def test_expnrs_do_not_collide_with_v0c_or_each_other(self):
        v0c_nrs = {p.expnr for p in V0C.points} | {V0C.reference.parent_expnr}
        for name in ("v0c16", "v0c16-tiny"):
            suite = get_suite(name)
            nrs = [p.expnr for p in suite.points] + [suite.reference.parent_expnr]
            self.assertEqual(len(nrs), len(set(nrs)),
                             f"{name}: colliding experiment numbers among {nrs}")
            self.assertFalse(set(nrs) & v0c_nrs,
                             f"{name} reuses a V0c experiment number: "
                             f"{set(nrs) & v0c_nrs}")
            # the r = 2 / r = 4 driver placeholders are never built, and the
            # r = 1 'driver' reuses the reference's own expnr on purpose.
            self.assertEqual(suite.drivers_to_run, [])

    def test_child_expnrs_match_the_register_comment(self):
        want = {"r1": "993", "r2": "994", "r4": "995"}
        got = {p.key: p.expnr for p in V0C16.points}
        self.assertEqual(got, want)

    def test_prolongation_is_left_unset_everywhere(self):
        for name in ("v0c16", "v0c16-tiny"):
            for pt in get_suite(name).points:
                self.assertIsNone(pt.child.prolongation, pt.key)

    def test_the_control_is_a_degenerate_filtered_point(self):
        for name in ("v0c16", "v0c16-tiny"):
            suite = get_suite(name)
            r1 = suite.point("r1")
            self.assertEqual(r1.refine, 1)
            self.assertEqual(r1.coarsen, 1)
            self.assertAlmostEqual(r1.driver.dx, r1.child.dx, places=12)

    def test_r2_and_r4_are_genuinely_coarsened(self):
        for name in ("v0c16", "v0c16-tiny"):
            suite = get_suite(name)
            for refine in (2, 4):
                pt = suite.point(f"r{refine}")
                self.assertEqual(pt.coarsen, refine)
                self.assertAlmostEqual(pt.driver.dx, refine * pt.child.dx, places=12)

    # -- registration --------------------------------------------------- #

    @staticmethod
    def _groups():
        tests_dir = HERE.parents[1]
        if str(tests_dir) not in sys.path:
            sys.path.insert(0, str(tests_dir))
        import run_tests

        manifest = run_tests._load_manifest()
        return manifest, run_tests

    def test_v0c16_is_in_nesting_validation_and_nowhere_else(self):
        manifest, run_tests = self._groups()
        labels = {g: [s["label"] for s in run_tests._expand_groups(manifest, g)]
                  for g in manifest["groups"]}
        pattern = re.compile(r"V0c16")
        v0c16 = [l for l in labels["nesting-validation"] if pattern.search(l)]
        self.assertEqual(len(v0c16), 1, f"expected one V0c16 entry, got {v0c16}")
        for group in ("all", "supported", "nesting", "experimental",
                      "python-library", "supported-macos", "lint"):
            self.assertFalse([l for l in labels.get(group, []) if pattern.search(l)],
                             f"group '{group}' reaches a V0c16 entry")


# --------------------------------------------------------------------------- #
# The end-to-end run
# --------------------------------------------------------------------------- #


class TestV0c16TinyPipeline(unittest.TestCase):
    """One run of the whole tiny V0c16 suite, at the real 4x4/16-rank decomposition."""

    @classmethod
    def setUpClass(cls) -> None:
        binary = caselib.solver_binary()
        if not binary.exists():
            raise unittest.SkipTest(f"solver binary not found at {binary}")
        cls.suite = get_suite("v0c16-tiny")
        override = os.environ.get("UDALES_V0C16_RUNDIR")
        if override:
            cls.rundir = Path(override)
            cls.rundir.mkdir(parents=True, exist_ok=True)
            cls._temp = None
        else:
            cls._temp = tempfile.mkdtemp(prefix="udales-v0c16-tiny-")
            cls.rundir = Path(cls._temp)
        env = dict(os.environ)
        env.setdefault("MPLCONFIGDIR", str(cls.rundir / ".mpl"))
        completed = subprocess.run(
            # --yes: v0c16-tiny's child runs at 4x4 = 16 ranks, above
            # run_v0.py's nprocs > 8 confirmation threshold (v0c-tiny's own
            # 2x2 = 4 ranks is under it, so test_v0c_tiny.py never needs this).
            [sys.executable, str(HERE / "run_v0.py"), str(cls.rundir),
             "--suite", "v0c16-tiny", "--no-plots", "--yes"],
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
        if cls._temp and not os.environ.get("UDALES_V0C16_KEEP"):
            shutil.rmtree(cls._temp, ignore_errors=True)
        elif cls._temp:
            print(f"run directory kept at {cls._temp}")

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

    def test_the_nesting_file_is_flux_balanced_at_every_ratio(self):
        for pt in self.suite.points:
            after = self.manifests[pt.key]["flux_residual_after_correction"]["max_abs_normalised"]
            self.assertLess(after, 1.0e-12, f"{pt.key}: Phi/A = {after:g}")

    def test_profile_rms_and_criterion_a_are_both_finite(self):
        for pt in self.suite.points:
            m = self.metrics[pt.key]
            profile_rms = m["profile_metrics"]["u_rms_difference_over_ustar"]
            crit_a = m["v2"]["criterion_a"]["max_interior_umean_error_over_ustar"]
            self.assertTrue(np.isfinite(profile_rms), pt.key)
            self.assertTrue(np.isfinite(crit_a), pt.key)

    def test_summary_table_written(self):
        for name in ("v0_summary.json", "v0_summary.csv", "v0_summary.md"):
            self.assertTrue((self.rundir / "analysis" / name).exists(), name)


# --------------------------------------------------------------------------- #
# Register the tests to run when invoked directly
# --------------------------------------------------------------------------- #


if __name__ == "__main__":
    unittest.main(verbosity=2)
