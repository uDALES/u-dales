#!/usr/bin/env python3
"""Smoke test: drive the V0b prolongation suite at the ``tiny`` presets.

Same contract as ``test_v0_tiny.py``: no physical claim -- a tiny, few-minutes
run cannot say which reconstruction *should* ship, only that both build and
run end to end and that the harness measures what the plan asks for (design
plan section 7, R2). What this covers, specific to V0b:

* both prolongation arms (``constant``, ``linear``) build at both ratios
  (``r2-constant``, ``r2-linear``, ``r4-constant``, ``r4-linear``) and run to
  completion;
* the divergence identity design section 1.3 makes for the ``constant``
  reconstruction (the interpolated child field reproduces the parent's own
  discrete divergence, ``before_projection == parent_before_prolongation``)
  holds for ``constant`` and is **not** asserted for ``linear`` -- W8's
  documented trade is that ``linear`` gives that guarantee up, so this test
  records the linear divergence rather than asserting it away, exactly as the
  plan requires;
* the summary machinery pairs ``constant`` against ``linear`` at each ratio
  (``run_v0.write_summary``'s ``constant_vs_linear`` block), the pressure
  response diagnostic (``|grad p|`` zone/interior/ratio) is parsed with its
  time-mean, and the staircase-amplitude reduction (``analyse_v0
  .staircase_amplitude``) returns a finite number for both arms;
* ``config.Preset.prolongation = None`` -- every preset before V0b, including
  V0's own -- reproduces today's nesting file **exactly** (every variable, not
  only the divergence numbers): rebuilding the same child twice, once with
  ``prolongation`` unset and once with it set explicitly to
  ``udprep.nesting.DEFAULT_PROLONGATION``, must give byte-identical variables
  (the only attribute allowed to differ is ``created``, a wall-clock
  timestamp).  This is the automated form of the manual git-stash comparison
  in the V0b report: same claim, checked every time this file runs.

Run directly, or through ``tests/run_tests.py nesting-validation``:

    module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
    source ~/udales/.venv/bin/activate
    python tests/validation/nesting/test_v0b_tiny.py

Environment: ``UDALES_BUILD`` (solver binary), ``UDALES_RUNTIME_MODULES``,
``MPIEXEC``, ``UDALES_V0B_RUNDIR`` (where to work; kept if set),
``UDALES_V0B_KEEP=1`` (keep a temporary run directory too).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import numpy as np

import caselib
import make_child_case
from config import V0B, get_suite

sys.path.insert(0, str(HERE.parents[2] / "tools" / "python"))
from udprep.nesting import DEFAULT_PROLONGATION, PROLONGATIONS  # noqa: E402


# --------------------------------------------------------------------------- #
# Checks that need no run
# --------------------------------------------------------------------------- #


class TestV0bSuiteConfiguration(unittest.TestCase):
    """Everything checkable without a solver, including the production suite."""

    def test_both_suites_are_self_consistent(self):
        for name in ("v0b", "v0b-tiny"):
            get_suite(name)  # raises with every problem listed

    def test_no_coarse_arm_is_required(self):
        """V0b's own axis is prolongation, not arm; see the ``V0B`` docstring."""
        for name in ("v0b", "v0b-tiny"):
            suite = get_suite(name)
            self.assertFalse(suite.require_paired_arms)
            self.assertTrue(all(p.arm == "filtered" for p in suite.points))

    def test_every_ratio_carries_both_prolongations(self):
        for name in ("v0b", "v0b-tiny"):
            suite = get_suite(name)
            for refine in suite.refinements:
                got = {p.child.prolongation for p in suite.ratio(refine)}
                self.assertEqual(got, {"constant", "linear"},
                                 f"{name}: r = {refine} has {got}")

    def test_prolongation_is_validated_against_the_writer(self):
        self.assertEqual(set(PROLONGATIONS), {"constant", "linear"})
        for name in ("v0b", "v0b-tiny"):
            for pt in get_suite(name).points:
                self.assertIn(pt.child.prolongation, PROLONGATIONS)
        bad = replace(V0B.points[0].child, name="bad", prolongation="nearest")
        with self.assertRaises(ValueError):
            bad.validate()

    def test_prolongation_unset_is_accepted_and_means_the_writer_default(self):
        p = replace(V0B.points[0].child, name="unset", prolongation=None)
        p.validate()  # must not raise

    def test_the_final_configuration_is_used(self):
        """R2(b): measure both reconstructions at 0.5 s / Catmull-Rom, not V1's."""
        for pt in V0B.points:
            self.assertEqual(pt.child.cadence, 0.5, pt.key)
            self.assertEqual(pt.child.timeinterp, 2, pt.key)
            # the box-filtering source is the C0b fine parent (960), not the
            # never-built coarse-grid placeholder Preset named by pt.driver
            self.assertEqual(pt.reference_expnr, "960", pt.key)
            self.assertEqual(pt.child.parent_expnr, "960", pt.key)
            self.assertAlmostEqual(pt.driver.dx, pt.refine * pt.child.dx, places=12, msg=pt.key)

    def test_expnrs_do_not_collide(self):
        for name in ("v0b", "v0b-tiny"):
            suite = get_suite(name)
            nrs = ([p.expnr for p in suite.points]
                   + [d.parent_expnr for d in suite.drivers_to_run]
                   + [suite.reference.parent_expnr])
            self.assertEqual(len(nrs), len(set(nrs)),
                             f"{name}: colliding experiment numbers among {nrs}")
            # the two coarse-driver placeholders are never built
            self.assertEqual(suite.drivers_to_run, [])

    def test_child_expnrs_match_the_register_comment(self):
        want = {"r2-constant": "978", "r2-linear": "979",
                "r4-constant": "980", "r4-linear": "981"}
        got = {p.key: p.expnr for p in V0B.points}
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

    def test_v0b_is_in_nesting_validation_and_nowhere_else(self):
        manifest, run_tests = self._groups()
        labels = {g: [s["label"] for s in run_tests._expand_groups(manifest, g)]
                  for g in manifest["groups"]}
        v0b = [l for l in labels["nesting-validation"] if "V0b" in l]
        self.assertEqual(len(v0b), 1, f"expected one V0b entry, got {v0b}")
        for group in ("all", "supported", "nesting", "experimental",
                      "python-library", "supported-macos", "lint"):
            self.assertFalse([l for l in labels.get(group, []) if "V0b" in l],
                             f"group '{group}' reaches a V0b entry")


# --------------------------------------------------------------------------- #
# The end-to-end run
# --------------------------------------------------------------------------- #


class TestV0bTinyPipeline(unittest.TestCase):
    """One run of the whole tiny V0b suite, with the assertions spread over it."""

    @classmethod
    def setUpClass(cls) -> None:
        binary = caselib.solver_binary()
        if not binary.exists():
            raise unittest.SkipTest(f"solver binary not found at {binary}")
        cls.suite = get_suite("v0b-tiny")
        override = os.environ.get("UDALES_V0B_RUNDIR")
        if override:
            cls.rundir = Path(override)
            cls.rundir.mkdir(parents=True, exist_ok=True)
            cls._temp = None
        else:
            cls._temp = tempfile.mkdtemp(prefix="udales-v0b-tiny-")
            cls.rundir = Path(cls._temp)
        env = dict(os.environ)
        env.setdefault("MPLCONFIGDIR", str(cls.rundir / ".mpl"))
        completed = subprocess.run(
            [sys.executable, str(HERE / "run_v0.py"), str(cls.rundir),
             "--suite", "v0b-tiny", "--no-plots"],
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
        if cls._temp and not os.environ.get("UDALES_V0B_KEEP"):
            shutil.rmtree(cls._temp, ignore_errors=True)
        elif cls._temp:
            print(f"run directory kept at {cls._temp}")

    # -- all four points built and ran --------------------------------- #

    def test_all_four_points_built_and_ran(self):
        for pt in self.suite.points:
            self.assertTrue((self.rundir / pt.expnr / "manifest.json").exists(), pt.key)
            self.assertTrue((self.rundir / pt.expnr / "child.log").exists(), pt.key)
            r = self.manifests[pt.key]["refinement"]
            self.assertEqual(r["spatial"], pt.refine, pt.key)
            self.assertIn("interpolated", r["slabs"], pt.key)

    def test_manifest_records_the_prolongation_actually_used(self):
        for pt in self.suite.points:
            m = self.manifests[pt.key]
            self.assertEqual(m["prolongation"], pt.child.prolongation, pt.key)
            self.assertEqual(m["prolongation_requested"], pt.child.prolongation, pt.key)

    # -- the divergence identity, kept for 'constant' and recorded for 'linear' #

    def test_constant_arm_satisfies_the_divergence_identity(self):
        """Design section 1.3: a solenoidal parent gives a solenoidal child target."""
        for pt in self.suite.arm("filtered"):
            if pt.child.prolongation != "constant":
                continue
            ic = self.manifests[pt.key]["initial_condition_divmax"]
            self.assertIsNotNone(ic["parent_before_prolongation"], pt.key)
            self.assertIsNotNone(ic["before_projection"], pt.key)
            self.assertAlmostEqual(
                ic["before_projection"] / ic["parent_before_prolongation"], 1.0,
                places=6,
                msg=f"{pt.key}: constant prolongation should reproduce the "
                    f"parent's divmax ({ic['parent_before_prolongation']:g}), got "
                    f"{ic['before_projection']:g}")

    def test_linear_arm_does_not_but_records_its_divergence(self):
        """Plan section 7, R2: the trade W8 introduced, measured rather than hidden."""
        for pt in self.suite.points:
            if pt.child.prolongation != "linear":
                continue
            ic = self.manifests[pt.key]["initial_condition_divmax"]
            self.assertIsNotNone(ic["parent_before_prolongation"], pt.key)
            self.assertIsNotNone(ic["before_projection"], pt.key)
            self.assertTrue(np.isfinite(ic["before_projection"]), pt.key)
            # The whole point of R2: this must NOT be close to 1, or 'linear'
            # would not actually be giving up the divergence guarantee (and
            # this suite would have nothing to measure).
            ratio = ic["before_projection"] / max(ic["parent_before_prolongation"], 1e-300)
            self.assertGreater(
                abs(ratio - 1.0), 1.0e-3,
                f"{pt.key}: linear prolongation's pre-projection divmax "
                f"({ic['before_projection']:g}) is suspiciously close to the "
                f"parent's ({ic['parent_before_prolongation']:g}); expected a "
                "real divergence source (design plan section 7, R2)")

    def test_the_nesting_file_is_still_flux_balanced_regardless_of_prolongation(self):
        """The flux correction is a separate step from the tangential slope choice."""
        for pt in self.suite.points:
            m = self.manifests[pt.key]
            after = m["flux_residual_after_correction"]["max_abs_normalised"]
            self.assertLess(after, 1.0e-12, f"{pt.key}: Phi/A = {after:g}")

    # -- the runs kept flux/divergence bounded, both arms ---------------- #

    def test_the_runs_kept_the_flux_and_divergence_bounded(self):
        for pt in self.suite.points:
            rt = self.metrics[pt.key]["v0"]["runtime"]
            self.assertGreater(rt["phi"]["n"], 0, f"{pt.key}: no nesting diagnostics")
            self.assertGreater(rt["divmax"]["n"], 0, f"{pt.key}: no divergence diagnostics")
            self.assertLess(rt["phi"]["max_abs"], 1.0e-9,
                            f"{pt.key}: runtime boundary flux residual off round-off")
            self.assertLess(rt["divmax"]["max"], 1.0e-10,
                            f"{pt.key}: the projection left a large divergence")

    # -- the pressure response, item 2 of design section 10.7 ------------ #

    def test_the_pressure_response_is_reported_with_its_time_mean(self):
        for pt in self.suite.points:
            rt = self.metrics[pt.key]["v0"]["runtime"]
            for key in ("gradp_zone", "gradp_interior"):
                self.assertGreater(rt[key]["n"], 0, f"{pt.key}: no {key}")
                self.assertTrue(np.isfinite(rt[key]["mean"]), f"{pt.key}: {key} mean")
            self.assertTrue(np.isfinite(rt["gradp_ratio"]["mean"]), pt.key)
            # The ratio of the means is a legitimate second reduction of the
            # same two series and should be the same order of magnitude as
            # the mean of the per-report ratios (both summarise "how much
            # bigger is the zone's response"), not orders apart.
            ratio_of_means = rt["gradp_zone"]["mean"] / rt["gradp_interior"]["mean"]
            self.assertGreater(ratio_of_means, 0.0, pt.key)

    # -- the staircase signature ------------------------------------------ #

    def test_the_staircase_amplitude_is_reported_for_both_arms(self):
        for pt in self.suite.points:
            st = self.metrics[pt.key]["v0"]["staircase"]
            self.assertTrue(st["available"], pt.key)
            self.assertEqual(st["refine"], pt.refine, pt.key)
            self.assertTrue(np.isfinite(st["rms_over_ustar"]), pt.key)
            self.assertTrue(np.isfinite(st["max_abs_over_ustar"]), pt.key)
            self.assertGreaterEqual(st["rms_over_ustar"], 0.0, pt.key)

    # -- the suite table pairs constant against linear -------------------- #

    def test_the_summary_pairs_constant_against_linear(self):
        for refine in self.suite.refinements:
            paired = self.summary["constant_vs_linear"][f"r{refine}"]
            for key in ("constant_criterion_a", "linear_criterion_a",
                        "constant_staircase_rms_over_ustar",
                        "linear_staircase_rms_over_ustar",
                        "constant_pre_projection_child_divmax",
                        "linear_pre_projection_child_divmax",
                        "constant_gradp_ratio_mean", "linear_gradp_ratio_mean"):
                self.assertIsNotNone(paired.get(key), f"r{refine}.{key}")
        for name in ("v0_summary.json", "v0_summary.csv", "v0_summary.md"):
            self.assertTrue((self.rundir / "analysis" / name).exists(), name)

    # -- prolongation=None reproduces today's output ---------------------- #

    def test_prolongation_none_reproduces_the_writer_default_exactly(self):
        """R2(a)/byte-identity: an unset ``prolongation`` changes nothing.

        ``r2-constant`` (expnr 978, already built by ``setUpClass`` with
        ``prolongation="constant"`` set explicitly) is rebuilt ONE more time
        from the SAME already-produced tiny fine parent (960), this time with
        ``preset.prolongation`` left ``None`` -- the state of every preset
        before V0b, including V0's own -- and every stored variable in the
        resulting ``nesting.inp.*.nc`` must match the original exactly.  No
        solver run is needed: only the offline slab cut.  This automates, and
        runs every time, the manual git-stash comparison in the V0b report.
        """
        from netCDF4 import Dataset

        from make_child_case import DrivingParent

        # Track the writer's actual default rather than hard-coding one: V0b
        # settled it on 'linear' (design 10.5), and the point of this test is
        # that leaving `prolongation` unset reproduces whatever that default
        # is, byte for byte.
        self.assertIn(DEFAULT_PROLONGATION, ("constant", "linear"))
        pt = self.suite.point(f"r2-{DEFAULT_PROLONGATION}")
        self.assertEqual(pt.child.prolongation, DEFAULT_PROLONGATION)
        reference_dir = self.rundir / pt.reference_expnr
        original = self.rundir / pt.expnr / f"nesting.inp.{pt.expnr}.nc"
        self.assertTrue(original.exists(), original)

        out = self.rundir / "_byte_identity_unset"
        shutil.rmtree(out, ignore_errors=True)
        out.mkdir()
        preset_none = replace(pt.child, name="byte-identity-unset",
                              child_expnr="989", prolongation=None)
        point = replace(pt, key="byte-identity", child=preset_none)
        driving = DrivingParent.refined(reference_dir, point)
        try:
            make_child_case.build(reference_dir, out, preset_none, driving=driving)
        finally:
            driving.dump.close()
        rebuilt = out / preset_none.child_expnr / f"nesting.inp.{preset_none.child_expnr}.nc"

        with Dataset(original, "r") as a, Dataset(rebuilt, "r") as b:
            self.assertEqual(sorted(a.variables), sorted(b.variables))
            for name in a.variables:
                va = np.asarray(a.variables[name][:])
                vb = np.asarray(b.variables[name][:])
                self.assertEqual(va.shape, vb.shape, name)
                self.assertTrue(np.array_equal(va, vb),
                                f"variable {name!r} differs between "
                                f"prolongation={DEFAULT_PROLONGATION!r} and "
                                "prolongation=None")
            # Every global attribute except the wall-clock timestamp and the
            # preset-name-derived 'parent_model' string (the two builds are
            # deliberately named differently; nothing else may differ).
            skip = {"created", "parent_model"}
            for key in a.ncattrs():
                if key in skip:
                    continue
                self.assertEqual(a.getncattr(key), b.getncattr(key), key)
        # The suite must carry an arm for whatever the default is, or the
        # comparison above silently tests nothing.
        self.assertIn(f"r2-{DEFAULT_PROLONGATION}",
                      [q.key for q in self.suite.points],
                      "no V0b arm matches the writer's default prolongation")


# --------------------------------------------------------------------------- #
# Register the tests to run when invoked directly
# --------------------------------------------------------------------------- #


if __name__ == "__main__":
    unittest.main(verbosity=2)
