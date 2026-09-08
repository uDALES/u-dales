#!/usr/bin/env python3
"""Smoke test: drive the whole V0 refinement suite at the ``tiny`` presets.

Same contract as ``test_v1_tiny.py`` and ``test_v2_tiny.py``.  It makes **no
physical claim** -- a 24 x 24 x 8 coarse parent run for 240 s cannot say whether
a child recovers the turbulence its parent could not hold, and the tiny suite
deliberately violates one of the design's own conditions (``L_rel >= 2 dx_P``
fails at ``r = 4``), which the configuration reports rather than hides.  What it
tests is the *harness*, and specifically the two things that only exist once the
refinement ratio exceeds one:

* the grids really are related by ``r`` -- the coarse parent's IBM mask is the
  fine reference's block-ANDed down, cell for cell, so the child's buildings are
  the parent's and the Big Brother comparison is not quietly comparing two
  geometries;
* the conservative prolongation of design section 1.3 does what section 1.3
  claims, offline **and** in the running solver: the interpolated child field
  reproduces the parent's discrete divergence, the stored boundary data is flux
  balanced, and the solver's own ``Phi`` and ``divmax`` stay at round-off.

The tiny and production suites differ only by the :class:`config.RefinementSuite`
object, so keeping this green keeps the production campaign runnable -- and the
production suite's own consistency is checked here too, without running it.

Run directly, or through ``tests/run_tests.py nesting-validation``:

    module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
    source /rds/general/user/mvr/home/udales/.venv/bin/activate
    python tests/validation/nesting/test_v0_tiny.py

Environment: ``UDALES_BUILD`` (solver binary), ``UDALES_RUNTIME_MODULES``,
``MPIEXEC``, ``UDALES_V0_RUNDIR`` (where to work; kept if set),
``UDALES_V0_KEEP=1`` (keep a temporary run directory too).
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
from caselib import coarsen_fluid_mask, coarsen_staggered, load_solid_mask
from config import CONVERGED, get_suite

sys.path.insert(0, str(HERE.parents[2] / "tools" / "python"))
from udprep.nesting import (
    DEFAULT_PROLONGATION,  # noqa: E402
    MAX_SPATIAL_REFINEMENT,
    NestGrid,
    conservative_interpolate,
    discrete_divergence,
)


# --------------------------------------------------------------------------- #
# Checks that need no run
# --------------------------------------------------------------------------- #


class TestSuiteConfiguration(unittest.TestCase):
    """Everything checkable without a solver, including the production suite."""

    def test_both_suites_are_self_consistent(self):
        for name in ("v0", "v0-tiny"):
            get_suite(name)  # raises with every problem listed

    def test_every_point_really_refines(self):
        for name in ("v0", "v0-tiny"):
            suite = get_suite(name)
            for pt in suite.points:
                self.assertGreater(pt.refine, 1, f"{name}/{pt.key}")
                self.assertLessEqual(pt.refine, MAX_SPATIAL_REFINEMENT,
                                     f"{name}/{pt.key} exceeds the writer's limit")
                self.assertAlmostEqual(pt.driver.dx, pt.refine * pt.child.dx, places=12)
                for a, b in ((pt.driver.itot * pt.refine, suite.reference.itot),
                             (pt.driver.jtot * pt.refine, suite.reference.jtot),
                             (pt.driver.ktot * pt.refine, suite.reference.ktot)):
                    self.assertEqual(a, b, f"{name}/{pt.key}: coarse grid is not r x coarser")

    def test_every_child_is_the_v1_child(self):
        """The reason the V1 result is the r = 1 row and not a separate experiment."""
        suite = get_suite("v0")
        for pt in suite.points:
            c = pt.child
            self.assertEqual(c.parent_expnr, CONVERGED.parent_expnr)
            for field in ("itot", "jtot", "ktot", "dx", "child_itot", "child_jtot",
                          "guardwidth", "zonewidth", "tau", "nzone", "nwall",
                          "ustar", "u0", "spinup", "production", "dtdump",
                          "child_spinup", "timeinterp", "init_from_parent"):
                self.assertEqual(getattr(c, field), getattr(CONVERGED, field),
                                 f"{pt.key}: child {field} differs from the V1 child")
            self.assertTrue(np.array_equal(np.sort(c.child_cube_centres(), axis=0),
                                           np.sort(CONVERGED.child_cube_centres(), axis=0)),
                            f"{pt.key}: the child's buildings are not the V1 child's")

    def test_the_coarse_parents_carry_the_same_cubes(self):
        for name in ("v0", "v0-tiny"):
            suite = get_suite(name)
            want = suite.reference.cube_centres()
            for pt in suite.points:
                got = pt.driver.cube_centres()
                self.assertTrue(np.array_equal(np.sort(got, axis=0), np.sort(want, axis=0)),
                                f"{name}/{pt.key}: coarse parent layout differs")
                self.assertEqual(pt.driver.plaza_window, suite.reference.plaza_window)

    def test_the_production_ramp_is_resolved_by_both_parents(self):
        """Design section 1.4(c): L_rel >= 2 dx_P, in the *parent's* terms."""
        for pt in get_suite("v0").points:
            self.assertTrue(pt.resolves_the_ramp,
                            f"{pt.key}: L_rel = {pt.child.zonewidth} m < "
                            f"2 dx_P = {2 * pt.driver.dx} m")

    def test_the_tiny_suite_reports_its_own_violation(self):
        """The smoke test's r = 4 point breaks 1.4(c); that must be visible."""
        tiny = get_suite("v0-tiny")
        self.assertFalse(tiny.point("r4-filtered").resolves_the_ramp)
        self.assertFalse(tiny.point("r4-coarse").resolves_the_ramp)

    def test_arms_share_one_coarse_grid_per_ratio(self):
        """Or 'filtered' and 'coarse' would not be comparing like with like."""
        for name in ("v0", "v0-tiny"):
            suite = get_suite(name)
            for refine in suite.refinements:
                grids = {(p.driver.itot, p.driver.jtot, p.driver.ktot, p.driver.dx)
                         for p in suite.ratio(refine)}
                self.assertEqual(len(grids), 1, f"{name}: r = {refine} has {grids}")
                arms = {p.arm for p in suite.ratio(refine)}
                self.assertEqual(arms, {"filtered", "coarse"})

    def test_only_the_coarse_arm_runs_a_parent(self):
        suite = get_suite("v0")
        self.assertEqual(len(suite.drivers_to_run), len(suite.refinements))
        for pt in suite.arm("filtered"):
            self.assertFalse(pt.runs_driver)
            self.assertEqual(pt.coarsen, pt.refine)
        for pt in suite.arm("coarse"):
            self.assertTrue(pt.runs_driver)
            self.assertEqual(pt.coarsen, 1)

    def test_expnrs_do_not_collide(self):
        suite = get_suite("v0")
        nrs = ([p.expnr for p in suite.points]
               + [d.parent_expnr for d in suite.drivers_to_run]
               + [suite.reference.parent_expnr])
        self.assertEqual(len(nrs), len(set(nrs)), f"colliding experiment numbers: {nrs}")

    def test_a_ratio_one_point_is_refused(self):
        from dataclasses import replace
        from config import RefinedPoint
        pt = get_suite("v0").point("r2-filtered")
        bad = RefinedPoint(key="bad", arm="filtered", refine=1,
                           driver=replace(pt.driver, dx=pt.child.dx,
                                          itot=pt.child.itot, jtot=pt.child.jtot,
                                          ktot=pt.child.ktot,
                                          child_itot=pt.child.child_itot,
                                          child_jtot=pt.child.child_jtot),
                           child=pt.child)
        with self.assertRaises(ValueError):
            bad.validate()


class TestCoarseningOperator(unittest.TestCase):
    """The filter that makes the 'filtered' arm meaningful, on its own.

    ``coarsen_staggered`` has to be flux conservative -- otherwise the coarse
    field it produces is not a parent any prolongation can be conservative
    against -- and it has to be the exact left inverse of the prolongation's
    tangential half, or a round trip would move the field and V0 would be
    measuring the round trip rather than the refinement.
    """

    @staticmethod
    def _solenoidal(n: int, dx: float, seed: int = 7):
        rng = np.random.default_rng(seed)
        a = rng.standard_normal((n + 1, n + 1, n + 1, 3))
        u = (np.diff(a[:, :, :, 2], axis=1)[:, :, :-1]
             - np.diff(a[:, :, :, 1], axis=2)[:, :-1, :]) / dx
        v = (np.diff(a[:, :, :, 0], axis=2)[:-1, :, :]
             - np.diff(a[:, :, :, 2], axis=0)[:, :, :-1]) / dx
        w = (np.diff(a[:, :, :, 1], axis=0)[:, :-1, :]
             - np.diff(a[:, :, :, 0], axis=1)[:-1, :, :]) / dx
        return u, v, w

    def test_face_fluxes_are_conserved(self):
        u, v, w = self._solenoidal(16, 2.0)
        ur, vr, wr = u[:-1], v[:, :-1], w[:, :, :-1]
        for f in (2, 4):
            uc, vc, wc = coarsen_staggered(ur, vr, wr, f)
            self.assertEqual(uc.shape, (16 // f,) * 3)
            for i, j, k in ((0, 0, 0), (3, 2, 1), (16 // f - 1,) * 3):
                self.assertAlmostEqual(
                    uc[i, j, k], ur[i * f, j * f:(j + 1) * f, k * f:(k + 1) * f].mean(),
                    places=14, msg=f"u face flux at r = {f}")
                self.assertAlmostEqual(
                    wc[i, j, k], wr[i * f:(i + 1) * f, j * f:(j + 1) * f, k * f].mean(),
                    places=14, msg=f"w face flux at r = {f}")

    def test_a_solenoidal_field_stays_solenoidal_through_the_round_trip(self):
        """Design section 1.3, both directions, on one field."""
        n, dx = 16, 2.0
        u, v, w = self._solenoidal(n, dx)
        fine = NestGrid.uniform(n, n, n, n * dx, n * dx, n * dx)
        self.assertLess(np.abs(discrete_divergence(fine, u, v, w)).max(), 1.0e-12)
        for f in (2, 4):
            m = n // f
            uc = u[::f].reshape(-1, m, f, m, f).mean(axis=(2, 4))
            vc = v[:, ::f].reshape(m, f, -1, m, f).mean(axis=(1, 4))
            wc = w[:, :, ::f].reshape(m, f, m, f, -1).mean(axis=(1, 3))
            coarse = NestGrid.uniform(m, m, m, n * dx, n * dx, n * dx)
            self.assertLess(np.abs(discrete_divergence(coarse, uc, vc, wc)).max(),
                            1.0e-12, f"coarsening broke the divergence at r = {f}")
            # Design 1.3's identity is a property of the divergence-preserving
            # tangential reconstruction, so name it rather than inheriting the
            # writer's default, which V0b settled on 'linear' (see 10.5).
            back = [conservative_interpolate(
                        coarse, pf, comp,
                        *[fine.component_coords(comp, ax) for ax in range(3)],
                        prolongation="constant")
                    for comp, pf in zip("uvw", (uc, vc, wc))]
            self.assertLess(np.abs(discrete_divergence(fine, *back)).max(), 1.0e-12,
                            f"prolongation broke the divergence at r = {f}")
            self.assertTrue(
                np.allclose(back[0][::f].reshape(-1, m, f, m, f).mean(axis=(2, 4)), uc),
                f"coarsen(prolong(x)) != x at r = {f}")

    def test_the_fluid_mask_filter_is_conservative(self):
        fluid = np.ones((8, 8, 8), dtype=bool)
        fluid[0, 0, 0] = False
        out = coarsen_fluid_mask(fluid, 2)
        self.assertEqual(out.shape, (4, 4, 4))
        self.assertFalse(out[0, 0, 0], "a partly solid coarse cell was called fluid")
        self.assertEqual(int((~out).sum()), 1)


# --------------------------------------------------------------------------- #
# The end-to-end run
# --------------------------------------------------------------------------- #


class TestV0TinyPipeline(unittest.TestCase):
    """One run of the whole tiny suite, with the assertions spread over it."""

    @classmethod
    def setUpClass(cls) -> None:
        binary = caselib.solver_binary()
        if not binary.exists():
            raise unittest.SkipTest(f"solver binary not found at {binary}")
        cls.suite = get_suite("v0-tiny")
        override = os.environ.get("UDALES_V0_RUNDIR")
        if override:
            cls.rundir = Path(override)
            cls.rundir.mkdir(parents=True, exist_ok=True)
            cls._temp = None
        else:
            cls._temp = tempfile.mkdtemp(prefix="udales-v0-tiny-")
            cls.rundir = Path(cls._temp)
        env = dict(os.environ)
        env.setdefault("MPLCONFIGDIR", str(cls.rundir / ".mpl"))
        completed = subprocess.run(
            [sys.executable, str(HERE / "run_v0.py"), str(cls.rundir),
             "--suite", "v0-tiny", "--no-plots"],
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
        cls.logs = {}
        for pt in cls.suite.points:
            cls.metrics[pt.key] = json.loads(
                (cls.rundir / "analysis" / pt.key / "v0_metrics.json").read_text())
            cls.manifests[pt.key] = json.loads(
                (cls.rundir / pt.expnr / "manifest.json").read_text())
            cls.logs[pt.key] = (cls.rundir / pt.expnr / "child.log").read_text(
                errors="replace")

    @classmethod
    def tearDownClass(cls) -> None:
        if cls._temp and not os.environ.get("UDALES_V0_KEEP"):
            shutil.rmtree(cls._temp, ignore_errors=True)
        elif cls._temp:
            print(f"run directory kept at {cls._temp}")

    # -- the grids ---------------------------------------------------------- #

    def test_the_coarse_parent_is_the_fine_one_block_averaged(self):
        """The geometric invariant the whole comparison rests on.

        Buildings are aligned to every grid in the suite, so the coarse run's
        IBM mask must be the fine run's block-ANDed down exactly.  If it is not,
        parent and child resolve different obstacles and nothing measured
        downstream means anything.
        """
        ref = self.suite.reference
        fine = load_solid_mask(self.rundir / ref.parent_expnr,
                               (ref.itot, ref.jtot, ref.ktot))
        self.assertGreater(int((~fine).sum()), 0, "the fine reference has no buildings")
        for driver in self.suite.drivers_to_run:
            refine = ref.itot // driver.itot
            coarse = load_solid_mask(self.rundir / driver.parent_expnr,
                                     (driver.itot, driver.jtot, driver.ktot))
            self.assertTrue(np.array_equal(coarse, coarsen_fluid_mask(fine, refine)),
                            f"the {driver.dx:g} m parent's solid cells are not the "
                            "2 m parent's, block-ANDed down")

    def test_the_child_geometry_is_the_fine_reference_sub_region(self):
        ref = self.suite.reference
        fine = load_solid_mask(self.rundir / ref.parent_expnr,
                               (ref.itot, ref.jtot, ref.ktot))
        sub = fine[ref.child_i0:ref.child_i0 + ref.child_itot,
                   ref.child_j0:ref.child_j0 + ref.child_jtot, :ref.child_ktot]
        for pt in self.suite.points:
            cm = load_solid_mask(self.rundir / pt.expnr,
                                 (ref.child_itot, ref.child_jtot, ref.child_ktot))
            self.assertTrue(np.array_equal(cm, sub),
                            f"{pt.key}: the child's solid cells are not the fine "
                            "reference's sub-region")

    # -- the nesting file --------------------------------------------------- #

    def test_the_slabs_were_interpolated_not_cut(self):
        for pt in self.suite.points:
            r = self.manifests[pt.key]["refinement"]
            self.assertEqual(r["spatial"], pt.refine, pt.key)
            self.assertAlmostEqual(r["spatial_from_file"], float(pt.refine), places=9)
            self.assertIn("interpolated", r["slabs"], pt.key)
            self.assertEqual(r["source_expnr"],
                             pt.reference_expnr if pt.arm == "filtered"
                             else pt.driver_expnr)
            self.assertEqual(r["coarsen_factor"], pt.coarsen, pt.key)

    def test_the_nesting_file_is_flux_balanced(self):
        for pt in self.suite.points:
            m = self.manifests[pt.key]
            after = m["flux_residual_after_correction"]["max_abs_normalised"]
            before = m["flux_residual_before_correction"]["max_abs_normalised"]
            self.assertLess(after, 1.0e-12, f"{pt.key}: corrected file has Phi/A = {after:g}")
            self.assertGreater(before, after * 100.0,
                               f"{pt.key}: the correction did no work")

    def test_the_prolongation_reproduces_the_parent_divergence(self):
        """Design section 1.3, end to end on a real field rather than in a unit test.

        The prolongation is piecewise constant tangentially and linear normally,
        so within a parent cell each of du/dx, dv/dy, dw/dz is the parent's --
        the interpolated child field must therefore carry the parent's discrete
        divergence, not merely a conservative average of it.  The two numbers
        come from opposite sides of the interpolation and are compared, not
        merely bounded.
        """
        for pt in self.suite.points:
            m = self.manifests[pt.key]
            ic = m["initial_condition_divmax"]
            self.assertIsNotNone(ic["parent_before_prolongation"], pt.key)
            scheme = m.get("prolongation") or DEFAULT_PROLONGATION
            if scheme == "constant":
                # The exact identity: only the divergence-preserving scheme has it.
                self.assertAlmostEqual(
                    ic["before_projection"] / ic["parent_before_prolongation"], 1.0,
                    places=6,
                    msg=f"{pt.key}: the prolonged field's divmax "
                        f"({ic['before_projection']:g}) is not the parent's "
                        f"({ic['parent_before_prolongation']:g})")
            else:
                # 'linear' trades the identity for accuracy (V0b, section 10.5).
                # What must still hold is that the projection removes what it
                # leaves -- record the size rather than assert it away.
                self.assertGreater(ic["before_projection"],
                                   ic["parent_before_prolongation"], pt.key)
            self.assertLess(ic["after_projection"], 1.0e-12, pt.key)

    # -- the nested run ----------------------------------------------------- #

    def test_the_runs_kept_the_flux_and_the_divergence_bounded(self):
        for pt in self.suite.points:
            rt = self.metrics[pt.key]["v0"]["runtime"]
            self.assertGreater(rt["phi"]["n"], 0, f"{pt.key}: no nesting diagnostics")
            self.assertGreater(rt["divmax"]["n"], 0, f"{pt.key}: no divergence diagnostics")
            self.assertLess(rt["phi"]["max_abs"], 1.0e-9,
                            f"{pt.key}: runtime boundary flux residual off round-off")
            self.assertLess(rt["divmax"]["max"], 1.0e-10,
                            f"{pt.key}: the projection left a large divergence")
            self.assertEqual(set(rt["faces_forced"]),
                             {"west", "east", "south", "north"}, pt.key)

    def test_the_zone_is_building_free_and_asserted(self):
        for pt in self.suite.points:
            self.assertNotIn("solid points inside the relaxation zone", self.logs[pt.key],
                             f"{pt.key}: nesting_init found solid points in the zone")

    # -- the analysis ------------------------------------------------------- #

    def test_the_spectra_are_split_at_the_parent_nyquist(self):
        for pt in self.suite.points:
            split = self.metrics[pt.key]["v0"]["spectra_across_parent_nyquist"]
            self.assertTrue(split, f"{pt.key}: no split spectra")
            for height, sp in split.items():
                self.assertAlmostEqual(sp["parent_nyquist_wavelength_m"],
                                       2.0 * pt.driver.dx, places=9)
                bands = sp["bands"]
                self.assertEqual(set(bands),
                                 {"parent_resolved", "parent_marginal",
                                  "sub_parent_filter"})
                total = sum(b["n_modes"] for b in bands.values())
                self.assertEqual(
                    total, len(self.metrics[pt.key]["spectra"][height]
                               ["wavenumber_rad_per_m"]) - 1,
                    f"{pt.key}/{height}: the three bands do not partition the modes")
                for name, b in bands.items():
                    if b["n_modes"]:
                        self.assertTrue(np.isfinite(b["ratio_of_sums"]),
                                        f"{pt.key}/{height}/{name}")

    def test_the_driving_parent_profile_is_recorded(self):
        """The only way to tell a bad child from a faithful child of a bad parent."""
        for pt in self.suite.points:
            pd = self.metrics[pt.key]["v0"]["parent_deficit"]
            self.assertTrue(pd["available"], pt.key)
            self.assertAlmostEqual(pd["parent_dx_m"], pt.driver.dx, places=9)
            self.assertTrue(np.isfinite(pd["u_rms_difference"]), pt.key)
            self.assertTrue(np.isfinite(pd["tke_rms_difference"]), pt.key)

    def test_the_analysis_emits_finite_numbers(self):
        for pt in self.suite.points:
            m = self.metrics[pt.key]
            for key, value in m["profile_metrics"].items():
                self.assertTrue(np.isfinite(value), f"{pt.key}: {key} is not finite")
            self.assertGreaterEqual(m["samples"]["parent"], 4, pt.key)
            self.assertGreaterEqual(m["samples"]["child"], 4, pt.key)
            outdir = self.rundir / "analysis" / pt.key
            for name in ("v0_metrics.json", "profiles.csv",
                         "spectral_bands_across_parent_nyquist.csv",
                         "driving_parent_vs_truth.csv"):
                self.assertTrue((outdir / name).exists(),
                                f"{pt.key}: {name} was not written")

    def test_the_children_did_not_diverge(self):
        """A loose bound, not a physics claim: the tiny run is far too short."""
        for pt in self.suite.points:
            m = self.metrics[pt.key]["profile_metrics"]
            self.assertLess(m["u_rms_difference_over_ustar"], 5.0, pt.key)
            self.assertLess(m["tke_rms_difference_over_ustar2"], 20.0, pt.key)

    # -- the suite table ---------------------------------------------------- #

    def test_the_summary_has_every_point_and_both_arms_paired(self):
        rows = {r["key"]: r for r in self.summary["rows"]}
        self.assertEqual(set(rows), {p.key for p in self.suite.points})
        self.assertEqual(sorted(self.summary["refinements"]),
                         list(self.suite.refinements))
        for refine in self.suite.refinements:
            paired = self.summary["filtered_vs_coarse"][f"r{refine}"]
            self.assertIsNotNone(paired["cost_of_a_real_parent"],
                                 f"r = {refine}: the two arms were not paired")
        for name in ("v0_summary.json", "v0_summary.csv", "v0_summary.md"):
            self.assertTrue((self.rundir / "analysis" / name).exists(), name)


# --------------------------------------------------------------------------- #
# Suite registration
# --------------------------------------------------------------------------- #


class TestSuiteRegistration(unittest.TestCase):
    """V0 must be reachable by name and unreachable from the gates."""

    @staticmethod
    def _groups():
        tests_dir = HERE.parents[1]
        if str(tests_dir) not in sys.path:
            sys.path.insert(0, str(tests_dir))
        import run_tests

        manifest = run_tests._load_manifest()
        return manifest, run_tests

    def test_v0_is_in_nesting_validation_and_nowhere_else(self):
        manifest, run_tests = self._groups()
        labels = {g: [s["label"] for s in run_tests._expand_groups(manifest, g)]
                  for g in manifest["groups"]}
        # "V0(?![bc])": matches this suite's own entries ("V0 tiny...", "V0
        # refinement...") but not V0b's or V0c's ("V0b tiny prolongation
        # suite...", "V0c tiny matched-control suite..."), which are
        # registered and checked separately in test_v0b_tiny.py/test_v0c_tiny.py.
        v0 = [l for l in labels["nesting-validation"] if re.search(r"V0(?![bc])", l)]
        self.assertEqual(len(v0), 2, f"expected two V0 entries, got {v0}")
        for group in ("all", "supported", "nesting", "experimental",
                      "python-library", "supported-macos", "lint"):
            self.assertFalse([l for l in labels.get(group, []) if re.search(r"V0(?!b)", l)],
                             f"group '{group}' reaches a V0 entry")


if __name__ == "__main__":
    unittest.main(verbosity=2)
