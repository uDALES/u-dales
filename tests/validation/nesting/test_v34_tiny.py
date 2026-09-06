#!/usr/bin/env python3
"""Smoke test: drive V3 and V4 end to end at a size that runs on a login node.

Same contract as ``test_v1_tiny.py`` and ``test_v2_tiny.py``.  It makes **no
physical claim** -- a 90 s window on a 40 s spin-up cannot say anything about an
adjustment length -- and it is not supposed to.  It tests the *harness*: that a
flat parent preprocesses and runs, that a volume-flow-rate controller can be
calibrated from a separate reference run and actually holds that bulk velocity,
that a child whose canopy is generated independently of its parent still asserts
a building-free zone, that the streamwise block reduction lands where the
geometry says it should, and that the V4 comparison refuses to run against a
baseline whose geometry is not the child's.

The configuration half needs no solver at all and checks the **production**
experiments, including the two claims the campaign rests on:

* the V4 child's cubes are the V1 ``converged`` child's cubes, so that V1 really
  is the reference;
* a staggered matched-geometry control is impossible at any size, which is why
  V4 has none.

Run directly, or through ``tests/run_tests.py nesting-validation``:

    module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
    source /rds/general/user/mvr/home/udales/.venv/bin/activate
    python tests/validation/nesting/test_v34_tiny.py

Environment: ``UDALES_BUILD`` (solver binary), ``UDALES_RUNTIME_MODULES``,
``MPIEXEC``, ``UDALES_V34_RUNDIR`` (where to work; kept if set),
``UDALES_V34_KEEP=1`` (keep a temporary run directory too).
"""

from __future__ import annotations

import json
import os
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
import analyse_geometry as ag
import caselib
import config
import presets_geometry as pg
import run_geometry
from config import CONVERGED


def _run_experiment(name: str, extra=()) -> Path:
    """Drive one tiny experiment end to end and return its run directory."""
    env = os.environ.get("UDALES_V34_RUNDIR")
    rundir = Path(env) / name if env else Path(tempfile.mkdtemp(prefix=f"{name}-"))
    argv = sys.argv
    sys.argv = ["run_geometry.py", str(rundir), "--experiment", name] + list(extra)
    try:
        rc = run_geometry.main()
    finally:
        sys.argv = argv
    if rc != 0:
        raise RuntimeError(f"{name} exited {rc}")
    return rundir


class TestGeometryConfiguration(unittest.TestCase):
    """Everything checkable without running anything -- production included."""

    def test_every_experiment_is_self_consistent(self):
        for name in sorted(pg.EXPERIMENTS):
            pg.get_experiment(name)  # raises with every problem listed

    def test_the_v4_child_is_the_v1_converged_child(self):
        """The claim that makes V1 the reference V4 is measured against.

        V4 changes the parent's layout and nothing else, so the child has to be
        V1's child -- not merely the same size, the same cubes in the same
        places.  V1 got its building-free zone by carving a plaza out of the
        parent; V4 gets it by clearing the child.  Those are different
        mechanisms and they have to produce the same layout, so it is checked
        rather than argued.
        """
        got = np.sort(pg.V4_MISMATCH.child_cube_centres(), axis=0)
        want = np.sort(CONVERGED.child_cube_centres(), axis=0)
        self.assertEqual(got.shape, want.shape)
        self.assertTrue(np.allclose(got, want))
        # ... and the plaza and the clearing remove the same 28 cubes.
        self.assertEqual(len(pg.V4_MISMATCH.child_cubes_removed()),
                         CONVERGED.n_cubes_removed)
        for f in ("itot", "jtot", "ktot", "dx", "child_itot", "child_jtot",
                  "guardwidth", "zonewidth", "tau", "nzone", "nwall",
                  "timeinterp", "ustar", "spinup", "production", "dtdump",
                  "child_spinup", "dtmax"):
            self.assertEqual(getattr(pg.V4_MISMATCH, f), getattr(CONVERGED, f), f)

    def test_a_staggered_matched_control_is_refused(self):
        """Why V4 carries no matched-geometry control, as a checked invariant.

        A dropped cube may not occupy 36 m of a 32 m period (26 m of clearance
        plus 8 m of half-width, at each end of the child), and a staggered
        array's two column families are exactly half a period apart, so one of
        them always lands in a blocked window.  If this ever stops firing, the
        argument in ``presets_geometry.V4`` needs revisiting -- not deleting.
        """
        with self.assertRaises(ValueError) as ctx:
            pg.V4_TINY_MATCHED_IMPOSSIBLE.validate()
        self.assertIn("reached the analysis interior", str(ctx.exception))

    def test_every_child_asserts_a_building_free_zone(self):
        for name in sorted(pg.EXPERIMENTS):
            for c in pg.get_experiment(name).children:
                p = c.preset
                self.assertTrue(p.building_free_zone, f"{name}/{c.key}")
                self.assertEqual(len(p.cubes_in_zone()), 0, f"{name}/{c.key}")
                self.assertEqual(len(p.removed_cubes_reaching_the_interior()), 0,
                                 f"{name}/{c.key}")
                self.assertGreater(len(p.cubes_in_analysis_interior()), 0,
                                   f"{name}/{c.key}")

    def test_the_v3_standoffs_are_the_design_tables_plus_a_long_one(self):
        """0/5/15 are section 10.4's; 40 is added so the claim can fail."""
        self.assertEqual(pg.V3_STANDOFFS[:3], (0, 5, 15))
        self.assertGreater(pg.V3_STANDOFFS[-1], 15)
        for c in pg.V3.children:
            p = c.preset
            # The first building face sits standoff cells inside the clear box,
            # which is nest_nwall cells deeper than the ramp.
            self.assertAlmostEqual(
                p.first_row_fetch_m,
                p.nwall * p.dx + p.standoff_cells * p.dx, places=9)
            self.assertGreaterEqual(p.n_rows, 8, c.key)

    def test_the_v3_parent_resolves_nothing_and_the_v4_parent_is_staggered(self):
        self.assertEqual(len(pg.V3_PARENT.cube_centres()), 0)
        aligned = pg.aligned_centres(pg.V4_PARENT.xlen, pg.V4_PARENT.ylen, 16.0, 16.0)
        stag = pg.V4_PARENT.cube_centres()
        # create_cubes' staggered array writes the two clipped half cubes on the
        # spanwise periodic faces as separate entries, so it carries Nx/2 more
        # nominal centres than the aligned array at the same plan area density.
        nx = int(round(pg.V4_PARENT.xlen / pg.V4_PARENT.period))
        self.assertEqual(len(stag), len(aligned) + nx // 2)
        # ... and the layouts really are different.
        self.assertFalse(np.allclose(np.sort(stag[:len(aligned)], axis=0),
                                     np.sort(aligned, axis=0)))

    def test_the_blocks_tile_the_interior_and_are_phase_locked(self):
        for name in sorted(pg.EXPERIMENTS):
            exp = pg.get_experiment(name)
            if exp.kind != "v3":
                continue
            for c in exp.children:
                p = c.preset
                e = ag.block_edges(p)
                lz = p.guardwidth + p.zonewidth
                self.assertGreater(len(e), 1, f"{name}/{c.key}")
                self.assertGreaterEqual(e[0], lz - 1.0e-9, f"{name}/{c.key}")
                self.assertLessEqual(e[-1], p.child_xlen - lz + 1.0e-9,
                                     f"{name}/{c.key}")
                self.assertTrue(np.allclose(np.diff(e), p.period))
                # phase-locked: some edge is the first building face
                self.assertTrue(np.any(np.isclose(e, p.canopy_x_range[0])),
                                f"{name}/{c.key}")

    def test_the_settling_rule_needs_every_later_point(self):
        f = np.array([0.0, 1.0, 2.0, 3.0])
        self.assertEqual(ag._first_settled(f, np.array([1.0, 0.0, 0.0, 0.0]), 0.5), 1.0)
        # a single late excursion pushes it back
        self.assertEqual(ag._first_settled(f, np.array([0.0, 1.0, 0.0, 0.0]), 0.5), 2.0)
        self.assertIsNone(ag._first_settled(f, np.array([0.0, 0.0, 0.0, 1.0]), 0.5))

    def test_the_presets_register_themselves_with_config(self):
        for name in ("v3-parent", "v3-reference", "v4-mismatch"):
            self.assertIn(name, config.PRESETS)
            self.assertIs(config.get_preset(name), pg.GEO_PRESETS[name])

    def test_the_new_suite_entries_stay_off_the_gates(self):
        sys.path.insert(0, str(HERE.parents[1]))
        import run_tests
        manifest = run_tests._load_manifest()
        for group in ("all", "supported", "nesting"):
            labels = [s["label"] for s in run_tests._expand_groups(manifest, group)]
            self.assertFalse([l for l in labels if l.startswith("nesting-validation")],
                             f"group '{group}' reaches the validation campaign")
        labels = [s["label"]
                  for s in run_tests._expand_groups(manifest, "nesting-validation")]
        self.assertTrue([l for l in labels if "V3" in l], labels)
        self.assertTrue([l for l in labels if "V4" in l], labels)


class _TinyRun(unittest.TestCase):
    experiment = ""
    rundir: Path = None
    _temp = False

    @classmethod
    def setUpClass(cls) -> None:
        if not cls.experiment:
            raise unittest.SkipTest("base class")
        binary = caselib.solver_binary()
        if not binary.exists():
            raise unittest.SkipTest(f"no solver binary at {binary}")
        cls._temp = not os.environ.get("UDALES_V34_RUNDIR")
        cls.rundir = _run_experiment(cls.experiment)
        cls.exp = pg.get_experiment(cls.experiment)

    @classmethod
    def tearDownClass(cls) -> None:
        if cls.rundir and cls._temp and not os.environ.get("UDALES_V34_KEEP"):
            shutil.rmtree(cls.rundir, ignore_errors=True)

    def periodic_stats(self, key):
        return json.loads(
            (self.rundir / "analysis" / f"periodic_{key}.json").read_text())

    def child_dir(self, key):
        return self.rundir / self.exp.child(key).preset.child_expnr

    def child_metrics(self, key, name):
        return json.loads((self.rundir / "analysis" / key / name).read_text())

    # -- shared assertions -------------------------------------------------- #

    def test_every_child_ran_with_the_zone_asserted_clear(self):
        for c in self.exp.default_children:
            p = c.preset
            nml = (self.child_dir(c.key) / f"namoptions.{p.child_expnr}").read_text()
            self.assertRegex(nml, r"(?im)^\s*nest_lparentgeom\s*=\s*\.false\.")
            log = (self.child_dir(c.key) / "child.log").read_text(errors="replace")
            self.assertNotIn("solid points inside the relaxation zone", log, c.key)
            self.assertIn("face west is forced", log, c.key)
            self.assertIn("face north is forced", log, c.key)

    def test_every_child_kept_the_flux_and_divergence_bounded(self):
        for c in self.exp.default_children:
            log = (self.child_dir(c.key) / "child.log").read_text(errors="replace")
            phis = [abs(float(l.split("=")[1])) for l in log.splitlines()
                    if "Phi (norm)" in l]
            divs = [abs(float(l.split("=")[1].split()[0])) for l in log.splitlines()
                    if "divmax, divtot" in l]
            self.assertTrue(phis, c.key)
            self.assertLess(max(phis), 1.0e-10, f"{c.key}: Phi {max(phis):.2e}")
            self.assertLess(max(divs), 1.0e-10, f"{c.key}: divmax {max(divs):.2e}")

    def test_the_childs_solid_mask_is_clear_in_the_zone(self):
        """The other end of ``nest_lparentgeom = .false.``, checked from the data."""
        for c in self.exp.default_children:
            p = c.preset
            shape = (p.child_itot - 1, p.child_jtot - 1, p.child_ktot - 1)
            fluid = caselib.load_solid_mask(self.child_dir(c.key), shape)
            solid = ~fluid
            n = p.zone_cells
            for sl in (np.s_[:n, :, :], np.s_[-n:, :, :],
                       np.s_[:, :n, :], np.s_[:, -n:, :]):
                self.assertEqual(int(solid[sl].sum()), 0,
                                 f"{c.key}: solid cells inside the zone")
            ii, jj = analyse.interior_indices(p)
            self.assertGreater(int(solid[np.ix_(ii, jj)].sum()), 0,
                               f"{c.key}: no buildings in the interior")

    def test_the_summary_exists_and_carries_a_verdict(self):
        out = self.rundir / "analysis" / f"{self.exp.kind}_summary.json"
        self.assertTrue(out.exists())
        summary = json.loads(out.read_text())
        self.assertTrue(summary["rows"])
        self.assertIn("verdict", summary)
        self.assertTrue((self.rundir / "analysis"
                         / f"{self.exp.kind}_summary.md").exists())


class TestV3Tiny(_TinyRun):
    experiment = "v3-tiny"

    def test_the_flat_parent_holds_the_reference_bulk_velocity(self):
        """The calibration the whole V3 design rests on.

        The flat parent is asked to run at the canopy's own bulk velocity so
        that the child's interior is in global momentum balance and the only
        thing out of equilibrium is the near-surface profile.  If the controller
        did not hold it, V3 would be measuring a bulk deceleration instead of an
        internal boundary layer.
        """
        ref = self.periodic_stats("reference")
        par = self.periodic_stats("parent")
        self.assertEqual(par["forcing"]["source"], "measured from 'reference'")
        self.assertAlmostEqual(par["forcing"]["uflowrate"],
                               ref["profiles"]["bulk_u"], places=9)
        rel = abs(par["profiles"]["bulk_u"] / ref["profiles"]["bulk_u"] - 1.0)
        self.assertLess(rel, 0.05, f"the parent's bulk is {100 * rel:.1f} % off "
                                   "the target it was given")

    def test_the_parent_really_has_no_buildings(self):
        par = self.periodic_stats("parent")
        self.assertEqual(par["parent_layout"], "none")
        self.assertEqual(par["samples"]["n_solid_cells"], 0)
        ref = self.periodic_stats("reference")
        self.assertGreater(ref["samples"]["n_solid_cells"], 0)

    def test_the_imposed_profile_is_out_of_equilibrium_with_the_canopy(self):
        """There is something for the child to adjust to -- the premise of V3."""
        m = self.child_metrics("standoff0", "v3_metrics.json")
        self.assertIsNotNone(m["mismatch"]["canopy_u_relative"])
        self.assertGreater(abs(m["mismatch"]["canopy_u_relative"]), 0.2)

    def test_the_blocks_cover_the_standoff_and_the_canopy(self):
        for c in self.exp.default_children:
            m = self.child_metrics(c.key, "v3_metrics.json")
            blocks = m["blocks"]["core"]
            self.assertTrue(blocks, c.key)
            n_free = sum(1 for b in blocks if not b["has_cube"])
            n_rows = sum(1 for b in blocks if b["has_cube"])
            self.assertEqual(n_rows, c.preset.n_rows, c.key)
            # A standoff long enough to hold a whole period must show up as at
            # least one building-free block ahead of the canopy; a zero standoff
            # must not.
            expected_free = int(c.preset.standoff_cells * c.preset.dx
                                >= c.preset.period)
            self.assertGreaterEqual(n_free, expected_free, c.key)
            if c.preset.standoff_cells == 0:
                self.assertEqual(n_free, 0, c.key)
            fetch = [b["fetch_from_zone_m"] for b in blocks]
            self.assertEqual(fetch, sorted(fetch), c.key)

    def test_the_core_scope_is_narrower_than_the_full_canopy(self):
        m = self.child_metrics("standoff0", "v3_metrics.json")
        core = m["blocks"]["core"][0]["n_fluid_cells"]
        full = m["blocks"]["canopy"][0]["n_fluid_cells"]
        self.assertLess(core, full)

    def test_the_adjustment_block_is_populated(self):
        for c in self.exp.default_children:
            m = self.child_metrics(c.key, "v3_metrics.json")
            adj = m["adjustment"]
            for q in ("u_canopy", "uw_at_roof"):
                self.assertIn(q, adj, f"{c.key}/{q}")
                self.assertIn("vs_equilibrium", adj[q])
                self.assertIn("vs_last_row", adj[q])
            self.assertEqual(len(m["ibl"]["delta_i_m"]), len(m["blocks"]["core"]))

    def test_the_child_versus_parent_block_is_written_but_labelled(self):
        """It is context here, not a criterion -- so it must not be v1_metrics."""
        d = self.rundir / "analysis" / "standoff0"
        self.assertTrue((d / run_geometry.LEGACY_NAME).exists())
        self.assertFalse((d / "v1_metrics.json").exists())


class TestV4Tiny(_TinyRun):
    experiment = "v4-tiny"

    def test_the_two_parents_carry_different_layouts(self):
        a = self.periodic_stats("aligned")
        s = self.periodic_stats("staggered")
        self.assertEqual(a["parent_layout"], "aligned")
        self.assertEqual(s["parent_layout"], "staggered")
        self.assertGreater(a["samples"]["n_solid_cells"], 0)
        self.assertGreater(s["samples"]["n_solid_cells"], 0)
        # Same plan area density, so the same number of solid cells -- counted
        # from solid_c.txt, not from the reduced arrays the statistics use.
        # ``cell_centred`` drops the last cell in each direction, and the
        # staggered array's clipped half cubes put solid cells in exactly that
        # last spanwise row, so the reduced counts differ by ~1 % between the two
        # layouts while the geometries carry identical solid volume.
        counts = []
        for key in ("aligned", "staggered"):
            nr = self.exp.periodic_run(key).preset.parent_expnr
            text = (self.rundir / nr / "solid_c.txt").read_text()
            counts.append(sum(1 for l in text.splitlines()
                              if l.strip() and not l.startswith("#")))
        self.assertEqual(counts[0], counts[1])

    def test_the_two_children_are_geometrically_identical(self):
        """What makes the difference between them attributable to the parent."""
        shapes = []
        for key in ("baseline", "mismatch"):
            p = self.exp.child(key).preset
            shape = (p.child_itot - 1, p.child_jtot - 1, p.child_ktot - 1)
            shapes.append(caselib.load_solid_mask(self.child_dir(key), shape))
        self.assertTrue(np.array_equal(shapes[0], shapes[1]))
        self.assertGreater(int((~shapes[0]).sum()), 0)

    def test_the_comparison_is_against_the_baseline_child(self):
        m = self.child_metrics("mismatch", "v4_metrics.json")
        self.assertIn("comparison", m)
        c = m["comparison"]
        self.assertEqual(c["labels"]["baseline"], "baseline")
        self.assertTrue(np.isfinite(c["umean_difference"]["rms_over_ustar"]))
        self.assertIsNotNone(
            c["criterion_a_prime"]["max_interior_umean_error_over_ustar"])
        self.assertTrue(c["spectra_child_over_baseline"])

    def test_a_baseline_with_the_wrong_geometry_is_refused(self):
        """The guard that stops V4 comparing two different canopies.

        Pointed at a *staggered* child (here, the aligned parent's own case
        directory relabelled) the comparison must abort rather than average two
        geometries together.  Simulated by handing it the parent's grid, which
        has a different mask shape.
        """
        c = self.exp.child("mismatch")
        with self.assertRaises(Exception):
            run_geometry.analyse_child_v4(
                self.rundir, self.exp, c, self.child_dir("mismatch"),
                self.rundir / self.exp.periodic_run("staggered").preset.parent_expnr,
                self.rundir / "analysis" / "_guard", False,
                (self.child_dir("baseline"),
                 pg.get_experiment("v3-tiny").child("standoff0").preset), {})


if __name__ == "__main__":
    unittest.main(verbosity=2)
