#!/usr/bin/env python3
"""Smoke test: drive both C0 cadence sweeps at the ``tiny`` preset.

Same contract as ``test_v2_tiny.py``: no physical claim -- a 96 x 96 x 32 parent
run for 240 s cannot say whether the boundary cadence causes the V1 TKE deficit
-- but every piece of harness C0 relies on is exercised end to end, and the
production sweeps' configuration is checked without running them.

What is exercised:

* **C0a** (``c0-tiny``): the reference child reused, three children cut from the
  same 3 s dumps at 6 s and 9 s (linear) and 3 s (Catmull-Rom).  The
  subsampling has to reduce what is *read*, not merely what is stored, so the
  manifest's ``stride``, the nesting file's time axis and the child's runtime
  are all checked against the parent's dump count.  The Catmull-Rom point is
  also the first running check that ``nest_timeinterp = 2`` is flux safe now
  that the limiter is gone: its ``Phi`` must be at round-off like the linear
  points'.
* **C0b** (``c0b-tiny``): a second parent **warm-started** from the tiny
  spin-up's restart files, dumping at 0.5 s, and three children sliced from
  those dumps at 0.5, 1.5 and 3 s -- built and run by ``run_v2.py`` itself
  through ``--parent-restart-dir``.  The parent's clock, dump interval and
  symlinked restart set are checked, and so is that the children dump at 3 s
  while the parent dumps at 0.5 s and the analysis still samples both at the
  same interval.

Run directly, or through ``tests/run_tests.py nesting-validation``:

    module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
    source /rds/general/user/mvr/home/udales/.venv/bin/activate
    python tests/validation/nesting/test_c0_tiny.py

Environment: ``UDALES_BUILD`` (solver binary), ``UDALES_RUNTIME_MODULES``,
``MPIEXEC``, ``UDALES_C0_RUNDIR`` (where to work; kept if set),
``UDALES_C0_KEEP=1`` (keep a temporary run directory too).
"""

from __future__ import annotations

import json
import os
import re
import shutil
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
import make_parent_case
import run_v1
import run_v2
from caselib import FieldDump, run_solver
from config import C0, C0B, C0_FINE, CONVERGED, TINY, get_preset, get_sweep

C0_EXPNRS = {str(n) for n in range(960, 970)}


def _nc_times(path: Path) -> np.ndarray:
    from netCDF4 import Dataset
    with Dataset(path, "r") as ds:
        return np.asarray(ds.variables["time"][:], dtype=float)


class TestC0Configuration(unittest.TestCase):
    """Everything about the production sweeps that needs no run."""

    def test_all_four_sweeps_are_self_consistent(self):
        for name in ("c0", "c0-tiny", "c0b", "c0b-tiny"):
            get_sweep(name)  # raises with every problem listed

    def test_c0a_is_the_v1_child_at_other_cadences_off_the_same_parent(self):
        sweep = get_sweep("c0")
        self.assertEqual(sweep.arms, ("cadence",))
        ref = sweep.point("ref")
        self.assertTrue(ref.reuse)
        self.assertIs(ref.preset, CONVERGED)
        self.assertEqual(len(sweep.to_run), 3)
        for pt in sweep.points:
            p = pt.preset
            self.assertEqual(p.parent_expnr, "903", pt.key)
            for field in ("child_itot", "child_jtot", "zonewidth", "guardwidth",
                          "nzone", "tau", "child_spinup", "production", "dtdump",
                          "init_from_parent", "clear_child_zone", "child_dtdump"):
                self.assertEqual(getattr(p, field), getattr(CONVERGED, field),
                                 f"{pt.key}: {field}")
            self.assertTrue(np.array_equal(p.cube_centres(), CONVERGED.cube_centres()))
            self.assertTrue(np.array_equal(p.child_cube_centres(),
                                           CONVERGED.child_cube_centres()))
        self.assertEqual(sweep.point("cad6").preset.cadence_stride, 2)
        self.assertEqual(sweep.point("cad9").preset.cadence_stride, 3)
        cr = sweep.point("cr3").preset
        self.assertEqual((cr.cadence, cr.cadence_stride, cr.timeinterp), (3.0, 1, 2))
        for key in ("cad6", "cad9"):
            self.assertEqual(sweep.point(key).preset.timeinterp, 1)

    def test_c0b_parent_is_the_converged_parent_continued_at_half_second_dumps(self):
        differs = {"name", "parent_expnr", "child_expnr", "production", "dtdump",
                   "cadence", "child_dtdump", "child_spinup", "stride"}
        for field in CONVERGED.__dataclass_fields__:
            if field in differs:
                continue
            self.assertEqual(getattr(C0_FINE, field), getattr(CONVERGED, field), field)
        # the restart the plan names was written at the end of the 10800 s spin-up
        self.assertEqual(C0_FINE.t_start, 10800.0)
        self.assertEqual(C0_FINE.t_end, 13200.0)
        self.assertEqual((C0_FINE.dtdump, C0_FINE.child_dtdump), (0.5, 3.0))
        self.assertEqual(C0_FINE.analysis_parent_stride, 6)
        self.assertEqual(C0_FINE.parent_expnr, "960")
        sweep = get_sweep("c0b")
        self.assertEqual(len(sweep.to_run), len(sweep.points))
        self.assertEqual([p.preset.cadence for p in sweep.points],
                         [0.5, 1.0, 1.5, 3.0, 6.0, 9.0])
        self.assertEqual([p.preset.cadence_stride for p in sweep.points],
                         [1, 2, 3, 6, 12, 18])
        self.assertIs(sweep.point("cad0.5").preset, C0_FINE)
        # 1800 s of statistics after the discard, as the plan sizes it
        self.assertAlmostEqual(C0_FINE.production - C0_FINE.child_spinup, 1800.0)

    def test_c0_experiment_numbers_collide_with_nothing(self):
        from config import SUITES, SWEEPS
        others = set()
        for name, s in SWEEPS.items():
            if name.startswith("c0"):
                continue
            others |= {p.preset.child_expnr for p in s.points} | {s.parent.parent_expnr}
        for s in SUITES.values():
            others |= {p.expnr for p in s.points} | {p.driver_expnr for p in s.points}
            others.add(s.reference.parent_expnr)
        geometry = set(re.findall(r'"(9\d\d)"', (HERE / "presets_geometry.py").read_text()))
        used = set()
        for name in ("c0", "c0b", "c0-tiny", "c0b-tiny"):
            s = get_sweep(name)
            used |= {p.preset.child_expnr for p in s.points if not p.reuse}
            used.add(s.parent.parent_expnr)
        used -= {"903", "904"}  # C0a's parent and reused reference are V1's
        # V0b's reference is deliberately the same directory as C0's own fine
        # parent (960, $EPHEMERAL/nesting-c0b/960) -- no new parent run, by
        # design (nesting-plan-2026-09-06.md section 7).  Same kind of
        # intentional reuse as 903/904 above, not a collision to catch.
        used -= {"960"}
        self.assertTrue(used <= C0_EXPNRS, used)
        self.assertFalse(used & others, used & others)
        self.assertFalse(used & geometry, used & geometry)

    def test_a_cadence_that_is_not_a_multiple_of_dtdump_is_rejected(self):
        for bad in ({"cadence": 4.0}, {"cadence": 0.0}, {"cadence": -3.0},
                    {"child_dtdump": 1.0}, {"cadence": 1.5}):
            with self.assertRaises(ValueError, msg=str(bad)) as ctx:
                replace(TINY, name="tiny-bad", **bad).validate()
            self.assertRegex(str(ctx.exception), "cadence|child_dtdump")
        # and a legitimate multiple passes, with the right stride
        ok = replace(TINY, name="tiny-ok", cadence=9.0)
        ok.validate()
        self.assertEqual(ok.cadence_stride, 3)

    def test_cadence_defaults_to_the_dump_interval(self):
        self.assertEqual(TINY.cadence, TINY.dtdump)
        self.assertEqual(CONVERGED.cadence_stride, 1)
        self.assertEqual(CONVERGED.child_dtdump, CONVERGED.dtdump)
        self.assertEqual(CONVERGED.analysis_parent_stride, CONVERGED.stride)

    def test_the_dump_courant_number(self):
        # C_dump = u cadence / dx: 3 m/s * 3 s / 2 m
        self.assertAlmostEqual(CONVERGED.c_dump_u0, 4.5)
        self.assertAlmostEqual(CONVERGED.dump_courant(3.58), 5.37)
        self.assertAlmostEqual(C0_FINE.c_dump_u0, 0.75)
        text = CONVERGED.describe()
        self.assertIn("boundary cadence     3 s", text)
        self.assertIn("C_dump = u0 * cadence / dx = 4.50 at u0 = 3 m/s", text)
        self.assertIn("Catmull-Rom", get_sweep("c0").point("cr3").preset.describe())

    def test_registration(self):
        """C0 must be reachable by name and unreachable from the gates."""
        tests_dir = HERE.parents[1]
        if str(tests_dir) not in sys.path:
            sys.path.insert(0, str(tests_dir))
        import run_tests
        manifest = run_tests._load_manifest()
        labels = {g: [s["label"] for s in run_tests._expand_groups(manifest, g)]
                  for g in manifest["groups"]}
        c0 = [l for l in labels["nesting-validation"] if "C0" in l]
        self.assertEqual(len(c0), 3, f"expected three C0 entries, got {c0}")
        for group in ("all", "supported", "nesting", "experimental",
                      "python-library", "supported-macos", "lint"):
            self.assertFalse([l for l in labels.get(group, []) if "C0" in l],
                             f"group '{group}' reaches a C0 entry")


class TestC0Tiny(unittest.TestCase):
    """Both tiny sweeps end to end, through ``run_v2.py``."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.c0a = get_sweep("c0-tiny")
        cls.c0b = get_sweep("c0b-tiny")
        cls.base = cls.c0a.parent            # TINY
        cls.fine = cls.c0b.parent            # C0_FINE_TINY
        override = os.environ.get("UDALES_C0_RUNDIR")
        if override:
            cls.rundir = Path(override)
            cls.rundir.mkdir(parents=True, exist_ok=True)
            cls._temp = None
        else:
            cls._temp = tempfile.mkdtemp(prefix="udales-c0-tiny-")
            cls.rundir = Path(cls._temp)
        cls.parent_dir = cls.rundir / cls.base.parent_expnr
        cls.fine_dir = cls.rundir / cls.fine.parent_expnr

        binary = caselib.solver_binary()
        if not binary.exists():
            raise unittest.SkipTest(f"solver binary not found at {binary}")

        p = cls.base
        # -- the V1 pipeline: tiny parent (spin-up + 3 s dumps), reference child #
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

        # -- C0a: off the same dumps, reference reused ----------------------- #
        cls.rc_a = cls._run_v2(str(cls.rundir), "--sweep", cls.c0a.name,
                               "--parent-dir", str(cls.parent_dir))
        # -- C0b: parent warm-started from the spin-up restart, then sliced --- #
        cls.rc_b = cls._run_v2(str(cls.rundir), "--sweep", cls.c0b.name,
                               "--parent-restart-dir", str(cls.parent_dir))
        cls.outdir = cls.rundir / "analysis"
        cls.metrics = {
            pt.key: json.loads((cls.outdir / pt.key / run_v2.METRICS_NAME).read_text())
            for sweep in (cls.c0a, cls.c0b) for pt in sweep.points}
        # run_v2 writes one summary per run directory; the second run overwrote
        # the first's, so the C0a table is rebuilt from the per-point metrics.
        import sweep_summary
        cls.summary_b = json.loads((cls.outdir / "sweep_summary.json").read_text())
        cls.summary_a = sweep_summary.build(cls.rundir, cls.c0a,
                                            metrics_name=run_v2.METRICS_NAME)
        cls.outdir_a = cls.rundir / "analysis-c0a"
        sweep_summary.write(cls.outdir_a, cls.c0a, cls.summary_a, make_plots=True)

    @staticmethod
    def _run_v2(*argv: str) -> int:
        saved = sys.argv
        sys.argv = ["run_v2.py", *argv]
        try:
            return run_v2.main()
        finally:
            sys.argv = saved

    @classmethod
    def tearDownClass(cls) -> None:
        if cls._temp and not os.environ.get("UDALES_C0_KEEP"):
            shutil.rmtree(cls._temp, ignore_errors=True)
        elif cls._temp:
            print(f"run directory kept at {cls._temp}")

    def _casedir(self, sweep, key: str) -> Path:
        return self.rundir / sweep.point(key).preset.child_expnr

    def _manifest(self, sweep, key: str) -> dict:
        return json.loads((self._casedir(sweep, key) / "manifest.json").read_text())

    # -- C0a ---------------------------------------------------------------- #

    def test_both_sweeps_ran_every_point(self):
        self.assertEqual((self.rc_a, self.rc_b), (0, 0))
        self.assertEqual({r["key"] for r in self.summary_a["rows"]},
                         {p.key for p in self.c0a.points})
        self.assertEqual({r["key"] for r in self.summary_b["rows"]},
                         {p.key for p in self.c0b.points})

    def test_the_subsampling_reads_every_nth_level_and_nothing_else(self):
        """Manifest stride, nesting-file time axis and runtime all agree."""
        n_dumped = FieldDump(self.parent_dir, self.base.parent_expnr, self.base.dx).ntime
        ref = self._manifest(self.c0a, "ref")
        self.assertEqual(ref["cadence"]["stride"], 1)
        self.assertEqual(ref["n_parent_levels"], n_dumped)
        for key, stride in (("cad6", 2), ("cad9", 3), ("cr3", 1)):
            m = self._manifest(self.c0a, key)
            p = self.c0a.point(key).preset
            cad = m["cadence"]
            self.assertEqual(cad["stride"], stride, key)
            self.assertEqual(cad["seconds"], p.cadence, key)
            self.assertEqual(cad["n_levels_dumped"], n_dumped, key)
            want = len(range(0, n_dumped, stride))
            self.assertEqual(cad["n_levels_used"], want, key)
            self.assertEqual(m["n_parent_levels"], want, key)
            # the streaming writer appended exactly the subsampled levels: it
            # never saw the ones the cadence skipped
            self.assertEqual(m["writer_diagnostics"]["ntime"], want, key)
            times = _nc_times(self._casedir(self.c0a, key) / f"nesting.inp.{p.child_expnr}.nc")
            self.assertEqual(times.size, want, f"{key}: nesting file levels")
            self.assertEqual(times[0], 0.0, key)
            self.assertTrue(np.all(np.diff(times) > 0), f"{key}: times not increasing")
            # the stored cadence is the requested one to within a time step
            self.assertLess(abs(np.median(np.diff(times)) - p.cadence), p.dtmax, key)
            self.assertAlmostEqual(m["parent_dt_median"], cad["parent_dt"], places=12)
            self.assertAlmostEqual(cad["parent_dt"], float(np.median(np.diff(times))), places=9)
            # runtime stops two levels short of the subsampled record
            self.assertAlmostEqual(m["runtime"], float(times[-3]), places=6, msg=key)
            # and the child ran to that runtime, dumping at its own interval
            ct = FieldDump(self._casedir(self.c0a, key), p.child_expnr, p.dx).times
            self.assertLess(abs(np.median(np.diff(ct)) - p.child_dtdump), p.dtmax, key)
        self.assertEqual(self._manifest(self.c0a, "cad6")["n_parent_levels"],
                         (n_dumped + 1) // 2)

    def test_the_catmull_rom_child_is_flux_safe(self):
        """Finding N1 is history: the unlimited interpolant keeps Phi at round-off."""
        nml = (self._casedir(self.c0a, "cr3") / "namoptions.963").read_text()
        self.assertRegex(nml, r"(?im)^\s*nest_timeinterp\s*=\s*2\b")
        for sweep in (self.c0a, self.c0b):
            for pt in sweep.points:
                log = (self._casedir(sweep, pt.key) / "child.log").read_text(errors="replace")
                phi = [float(x) for x in re.findall(r"Phi \(norm\) =\s*(\S+)", log)]
                div = [float(x) for x in re.findall(r"divmax, divtot =\s*(\S+)", log)]
                self.assertTrue(phi and div, f"{pt.key}: no nesting diagnostics")
                self.assertLess(max(abs(x) for x in phi), 1.0e-9, pt.key)
                self.assertLess(max(div), 1.0e-10, pt.key)

    def test_the_reference_was_reused_and_the_guard_sees_the_cadence(self):
        row = next(r for r in self.summary_a["rows"] if r["key"] == "ref")
        self.assertTrue(row["reused"])
        ref = self.c0a.point("ref")
        with self.assertRaises(RuntimeError) as ctx:
            run_v2._check_reused(ref, self._casedir(self.c0a, "cad6"))
        self.assertIn("cadence", str(ctx.exception))

    def test_the_c0a_table_is_produced(self):
        text = (self.outdir_a / "sweep_summary.md").read_text()
        self.assertIn("C0 -- boundary-data cadence", text)
        self.assertIn("Band ratios (child/parent) at every sampled height", text)
        for key in ("ref", "cad6", "cad9", "cr3"):
            self.assertEqual(text.count(f"| {key} "), 2, key)
        for name in ("sweep_summary.csv", "sweep_summary.json", "sweep_summary.md",
                     "sweep_deficit.png", "sweep_spectra.png", "sweep_cadence_bands.png"):
            self.assertTrue((self.outdir_a / name).exists(), name)
        v = self.summary_a["verdict"]["cadence_arm"]
        self.assertEqual(v["n_points"], 4)
        self.assertEqual(v["x"], sorted(v["x"]))
        self.assertEqual(v["keys"], ["ref", "cr3", "cad6", "cad9"])
        for band in ("band_8_16_at_z2h", "band_16_64_at_z2h"):
            self.assertTrue(all(np.isfinite(x) for x in v[band]), band)
        # the per-height table has every sampled height for every point
        heights = {n: hb["z_over_h"] for r in self.summary_a["rows"]
                   for n, hb in r["bands_by_height"].items()}
        self.assertEqual(len(heights), len(self.base.spectra_heights))
        csv = (self.outdir_a / "sweep_summary.csv").read_text().splitlines()[0]
        self.assertIn("band_8_16m@", csv)
        self.assertIn("cadence_s", csv)

    # -- C0b ---------------------------------------------------------------- #

    def test_the_fine_parent_was_warm_started_from_the_spinup_restart(self):
        nr, src = self.fine.parent_expnr, self.base.parent_expnr
        nml = (self.fine_dir / f"namoptions.{nr}").read_text()
        self.assertRegex(nml, r"(?im)^\s*lwarmstart\s*=\s*\.true\.")
        m = re.search(r"(?im)^\s*startfile\s*=\s*'(initd\d{8}_000_000\.%s)'" % nr, nml)
        self.assertIsNotNone(m, "startfile not set to the linked restart")
        self.assertRegex(nml, r"(?im)^\s*tfielddump\s*=\s*0\.5\b")
        self.assertRegex(nml, r"(?im)^\s*runtime\s*=\s*120\b")
        self.assertFalse((self.fine_dir / f"namoptions_spinup.{nr}").exists())
        links = sorted(self.fine_dir.glob(f"initd????????_???_???.{nr}"))
        self.assertEqual(len(links), self.fine.nprocx * self.fine.nprocy)
        for link in links:
            self.assertTrue(link.is_symlink(), link.name)
            target = link.resolve()
            self.assertEqual(target.parent, self.parent_dir.resolve())
            self.assertEqual(target.name, link.name[:-len(nr)] + src)
        info = json.loads((self.fine_dir / "preset.json").read_text())
        self.assertEqual(info["startfile"], m.group(1))
        self.assertEqual(Path(info["warmstart_from"]), self.parent_dir.resolve())
        # the clock continued from the restart: dumps start at t_start, every 0.5 s
        times = FieldDump(self.fine_dir, nr, self.fine.dx).times
        self.assertGreaterEqual(times[0], self.fine.t_start)
        self.assertLess(times[0] - self.fine.t_start, 2 * self.fine.dtdump)
        self.assertLess(abs(np.median(np.diff(times)) - self.fine.dtdump), self.fine.dtmax)
        self.assertGreaterEqual(times.size, 0.8 * self.fine.production / self.fine.dtdump)
        # and it was not re-run by the second stage of run_v2
        self.assertEqual(len(list(self.fine_dir.glob("production.log"))), 1)

    def test_the_c0b_children_are_sliced_from_the_fine_dumps(self):
        n_dumped = FieldDump(self.fine_dir, self.fine.parent_expnr, self.fine.dx).ntime
        for key, stride in (("cad0.5", 1), ("cad1.5", 3), ("cad3", 6)):
            m = self._manifest(self.c0b, key)
            p = self.c0b.point(key).preset
            self.assertEqual(m["parent_expnr"], self.fine.parent_expnr)
            self.assertEqual(m["cadence"]["stride"], stride, key)
            self.assertEqual(m["n_parent_levels"], len(range(0, n_dumped, stride)), key)
            self.assertEqual(m["writer_diagnostics"]["ntime"], m["n_parent_levels"], key)
            times = _nc_times(self._casedir(self.c0b, key) / f"nesting.inp.{p.child_expnr}.nc")
            self.assertEqual(times.size, m["n_parent_levels"], key)
            self.assertEqual(times[0], 0.0)
            self.assertTrue(np.all(np.diff(times) > 0), key)
            self.assertLess(abs(np.median(np.diff(times)) - p.cadence), p.dtmax, key)
            # the child dumps at 3 s whatever the parent did
            ct = FieldDump(self._casedir(self.c0b, key), p.child_expnr, p.dx).times
            self.assertLess(abs(np.median(np.diff(ct)) - 3.0), p.dtmax, key)
            self.assertEqual(m["cadence"]["child_dtdump"], 3.0)
            self.assertEqual(m["cadence"]["dump_dt"], 0.5)

    def test_the_analysis_samples_parent_and_child_at_the_same_interval(self):
        for pt in self.c0b.points:
            s = self.metrics[pt.key]["samples"]
            # parent every 6th 0.5 s level, child every 3 s level: same count to
            # within the jitter of adaptive stepping at the window's edges
            self.assertLessEqual(abs(s["parent"] - s["child"]), 2, f"{pt.key}: {s}")
            self.assertGreaterEqual(s["child"], 4, pt.key)
            cfg = self.metrics[pt.key]["v2"]["configuration"]
            self.assertEqual(cfg["parent_dtdump_s"], 0.5)
            self.assertEqual(cfg["child_dtdump_s"], 3.0)
            self.assertEqual(cfg["cadence_s"], pt.preset.cadence)

    def test_the_c0b_table_is_produced(self):
        text = (self.outdir / "sweep_summary.md").read_text()
        self.assertIn("'c0b-tiny'", text)
        self.assertIn("Band ratios (child/parent) at every sampled height", text)
        v = self.summary_b["verdict"]["cadence_arm"]
        self.assertEqual(v["n_points"], 3)
        self.assertEqual(v["x"], [0.5, 1.5, 3.0])
        self.assertTrue((self.outdir / "sweep_cadence_bands.png").exists())
        self.assertEqual(self.summary_b["verdict"]["confounds"], [])

    def test_every_point_emits_finite_headline_numbers(self):
        for key, m in self.metrics.items():
            d = m["v2"]["tke_deficit"]["above"]
            self.assertTrue(np.isfinite(d["mean_relative"]), key)
            for name, sp in m["spectra"].items():
                for band in ("band_16_64m", "band_8_16m"):
                    self.assertTrue(np.isfinite(sp["bands"][band]["mean_of_ratios"]),
                                    f"{key}/{name}/{band}")
            pm = m["profile_metrics"]
            self.assertLess(pm["u_rms_difference_over_ustar"], 2.0, key)
            self.assertLess(pm["tke_rms_difference_over_ustar2"], 5.0, key)


if __name__ == "__main__":
    unittest.main(verbosity=2)
