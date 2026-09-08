#!/usr/bin/env python3
"""Smoke test: drive the V6 (long-run mass-drift) pipeline at tiny scale.

Two things are checked, and they are different kinds of claim:

* :class:`TestV6TinyPipeline` runs the real harness -- ``config.V6``'s parent
  and child, built and run exactly as ``run_v6.py`` would, but with a runtime
  only a little past the parent's short record instead of ~1e5 steps.  This
  cannot say anything about whether mass drifts over a long run (that is what
  the production submission is for); it can and does say that the harness
  builds a case that runs to completion with ``nest_lendabort = .false.``,
  that the documented freeze (``src/modnesting.f90``'s ``check_record_end``)
  fires exactly once and does not abort the run, and that
  ``analyse_v6`` parses the resulting log into finite, sane numbers.
* :class:`TestLinfitAndVerdict` is a pure-Python unit test of
  ``analyse_v6.summarize_series``'s pass/fail logic against synthetic series
  with a known slope (or none), independent of the solver entirely -- the
  tiny run above is far too short and noisy to exercise both branches of the
  verdict on its own.

Run directly, or through ``tests/run_tests.py nesting-validation``:

    module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
    source /rds/general/user/mvr/home/udales/.venv/bin/activate
    python tests/validation/nesting/test_v6_tiny.py

Environment: ``UDALES_BUILD`` (solver binary), ``UDALES_RUNTIME_MODULES``,
``MPIEXEC``, ``UDALES_V6_RUNDIR`` (where to work; kept if set),
``UDALES_V6_KEEP=1`` (keep a temporary run directory too).
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

import analyse_v6
import caselib
import make_child_case
import make_parent_case
from caselib import run_solver, set_namoption
from config import get_preset
from run_v1 import _restart_file, _set_startfile

#: Well past the parent's ~8-level, ~15-24 s record, so the freeze is
#: actually exercised -- but still a login-node-in-minutes number, not the
#: production ~1e5-step target.
_TINY_RUNTIME_S = 60.0
#: Small enough to get several nesting_stats/chkdiv samples inside 60 s.
_TINY_STAT_INTERVAL_S = 4.0


class TestV6TinyPipeline(unittest.TestCase):
    """One end-to-end run at tiny scale, past the record, well short of 1e5 steps."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.preset = get_preset("v6")
        override = os.environ.get("UDALES_V6_RUNDIR")
        if override:
            cls.rundir = Path(override)
            cls.rundir.mkdir(parents=True, exist_ok=True)
            cls._temp = None
        else:
            cls._temp = tempfile.mkdtemp(prefix="udales-v6-tiny-")
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
        _set_startfile(cls.parent_dir / f"namoptions.{p.parent_expnr}",
                       _restart_file(cls.parent_dir, p.parent_expnr))
        run_solver(cls.parent_dir, f"namoptions.{p.parent_expnr}",
                   p.nprocx * p.nprocy, cls.parent_dir / "production.log")

        make_child_case.build(cls.parent_dir, cls.rundir, p)
        cls.manifest = json.loads((cls.child_dir / "manifest.json").read_text())
        namoptions = cls.child_dir / f"namoptions.{p.child_expnr}"
        set_namoption(namoptions, "runtime", _TINY_RUNTIME_S)
        set_namoption(namoptions, "nest_statint", _TINY_STAT_INTERVAL_S)
        set_namoption(namoptions, "tcheck", _TINY_STAT_INTERVAL_S)

        run_solver(cls.child_dir, f"namoptions.{p.child_expnr}",
                   p.child_nprocx * p.child_nprocy, cls.child_dir / "child.log")
        cls.child_log = (cls.child_dir / "child.log").read_text(errors="replace")
        cls.result = analyse_v6.run(cls.child_dir, cls.outdir, make_plots=False,
                                    expected_runtime_s=_TINY_RUNTIME_S)

    @classmethod
    def tearDownClass(cls) -> None:
        if cls._temp and not os.environ.get("UDALES_V6_KEEP"):
            shutil.rmtree(cls._temp, ignore_errors=True)
        elif cls._temp:
            print(f"run directory kept at {cls._temp}")

    # -- the construction ---------------------------------------------------- #

    def test_preset_sets_lendabort_false(self):
        self.assertFalse(self.preset.nest_lendabort)

    def test_runtime_is_past_the_stored_record(self):
        """The whole point: the run has to outlast what make_child_case cut."""
        self.assertGreater(_TINY_RUNTIME_S, self.manifest["runtime"],
                           "the tiny runtime does not exceed the parent record; "
                           "the freeze would never be exercised")

    def test_freeze_warning_fires_exactly_once(self):
        """nendwarn in modnesting.f90: warn once, then go quiet and hold."""
        self.assertEqual(
            self.result["n_freeze_warnings"], 1,
            f"expected exactly one freeze warning, log carries "
            f"{self.result['n_freeze_warnings']}:\n{self.child_log}")

    def test_run_did_not_abort(self):
        """nest_lendabort = .false. means this must never fire."""
        self.assertFalse(self.result["aborted_past_record"])
        self.assertNotIn("ERROR", self.child_log)

    def test_all_four_faces_are_forced(self):
        for face in ("west", "east", "south", "north"):
            self.assertIn(f"face {face} is forced", self.child_log)

    # -- the diagnostics ------------------------------------------------------ #

    def test_diagnostics_were_sampled(self):
        self.assertGreater(self.result["n_checksim_samples"], 1,
                           "no divmax/divtot samples parsed -- tcheck too coarse?")
        self.assertGreater(self.result["n_nesting_samples"], 1,
                           "no nesting_stats samples parsed -- nest_statint too coarse?")

    def test_divergence_is_bounded_and_finite(self):
        for key in ("divmax", "divtot"):
            r = self.result["series"][key]
            self.assertGreater(r["n"], 0)
            self.assertTrue(np.isfinite(r["max_abs"]), f"{key} is not finite")
            self.assertLess(r["max_abs"], 1.0e-6,
                            f"{key} = {r['max_abs']:.3e} is far above round-off")

    def test_flux_residual_is_bounded_and_finite(self):
        r = self.result["series"]["phi"]
        self.assertGreater(r["n"], 0)
        self.assertTrue(np.isfinite(r["max_abs"]))
        self.assertLess(r["max_abs"], 1.0e-6,
                        f"Phi = {r['max_abs']:.3e} is far above round-off")

    def test_every_series_has_a_verdict(self):
        for key in analyse_v6.SERIES:
            r = self.result["series"][key]
            self.assertIn("verdict", r)
            self.assertIn(r["n"], range(0, 10**6))  # just: present and an int

    def test_summary_text_is_produced(self):
        text = analyse_v6.summary(self.result)
        self.assertIn("V6", text)
        self.assertIn("Overall:", text)


class TestLinfitAndVerdict(unittest.TestCase):
    """Pure-Python check of the trend-fit and pass/fail logic, no solver involved.

    The tiny run above is 60 s of a noisy, transient-laden LES -- it cannot be
    relied on to produce both a clean "flat" series and a clean "drifting" one
    on demand.  These synthetic series can.
    """

    def test_flat_noisy_series_passes(self):
        rng = np.random.default_rng(0)
        t = np.linspace(0.0, 1000.0, 200)
        y = 1.0e-9 + 1.0e-11 * rng.standard_normal(t.size)
        r = analyse_v6.summarize_series("flat", t, y)
        self.assertTrue(str(r["verdict"]).startswith("PASS"), r)

    def test_monotone_drift_fails(self):
        rng = np.random.default_rng(1)
        t = np.linspace(0.0, 1000.0, 200)
        # slope large enough to move the series by several multiples of its
        # own noise floor over the run -- an unambiguous drift.
        y = 1.0e-9 + 5.0e-11 * t + 1.0e-11 * rng.standard_normal(t.size)
        r = analyse_v6.summarize_series("drifting", t, y)
        self.assertTrue(str(r["verdict"]).startswith("FAIL"), r)
        self.assertGreater(abs(r["t_stat"]), 3.0)

    def test_too_few_samples_is_inconclusive_not_a_false_pass(self):
        r = analyse_v6.summarize_series("sparse", [0.0, 1.0], [1.0e-9, 1.0e-9])
        self.assertEqual(r["verdict"], "INCONCLUSIVE")

    def test_empty_series_reports_no_data(self):
        r = analyse_v6.summarize_series("empty", [], [])
        self.assertEqual(r["verdict"], "NO DATA")

    def test_log_parser_on_a_synthetic_fragment(self):
        text = (
            " Time of Day: 120000.000    Time of Simulation:         3.00000    dt:  0.500000000\n"
            "divmax, divtot =   1.00E-09 2.00E-09\n"
            " modnesting: t          =        3.000\n"
            " modnesting: Phi (norm) = 1.5000E-10  (largest |Phi| since the previous report)\n"
            " modnesting: Phi lid    = 1.0000E-10  closed faces = 2.0000E-10  (largest since the previous report)\n"
            " modnesting: zone misfit rms [m/s] = 3.0000E-05\n"
            " modnesting: |grad p| zone = 1.2000E+00  interior = 1.1000E+00  ratio =    1.091\n"
            " modnesting: WARNING t =  3.00000E+00 is past the last parent time level at "
            "t =  2.40000E+01; the boundary now freezes on that level (nest_lendabort = .false.)\n"
        )
        p = HERE / "_v6_log_fragment_test.tmp"
        p.write_text(text, encoding="ascii")
        try:
            parsed = analyse_v6.parse_log(p)
        finally:
            p.unlink(missing_ok=True)
        self.assertEqual(parsed["checksim_t"], [3.0])
        self.assertEqual(parsed["divmax"], [1.00e-09])
        self.assertEqual(parsed["divtot"], [2.00e-09])
        self.assertEqual(len(parsed["nest_records"]), 1)
        rec = parsed["nest_records"][0]
        self.assertEqual(rec["t"], 3.0)
        self.assertAlmostEqual(rec["phi"], 1.5e-10)
        self.assertAlmostEqual(rec["phi_lid"], 1.0e-10)
        self.assertAlmostEqual(rec["phi_closed"], 2.0e-10)
        self.assertAlmostEqual(rec["misfit_rms"], 3.0e-05)
        self.assertAlmostEqual(rec["gradp_zone"], 1.2)
        self.assertAlmostEqual(rec["gradp_interior"], 1.1)
        self.assertAlmostEqual(rec["gradp_ratio"], 1.091)
        self.assertEqual(parsed["n_freeze_warnings"], 1)
        self.assertFalse(parsed["aborted_past_record"])


class TestOverallVerdict(unittest.TestCase):
    """Unit tests for ``analyse_v6.analyse()``'s aggregate PASS/FAIL/NO-DATA/
    INCONCLUSIVE verdict -- no solver involved, and independent of
    :class:`TestLinfitAndVerdict` above, which only checks one series' verdict
    in isolation.

    These guard the bug fixed alongside them: the aggregate used to read
    "PASS" off of "no verdict starts with FAIL", so an empty log (no verdicts
    at all) or an incomplete one (all INCONCLUSIVE) passed vacuously --
    ``analyse_v6.analyse(analyse_v6.parse_log(Path('/dev/null')))
    ['overall_verdict']`` was ``'PASS'``. PASS now requires positive evidence
    (complete, finite series; enough samples; the intended duration reached;
    the freeze warning seen exactly once; no abort past the record), and
    every one of the cases below must NOT come back PASS.
    """

    @staticmethod
    def _flat_parsed(n: int = 20, duration: float = 1000.0,
                     n_freeze: int = 1, aborted: bool = False):
        """A hand-built ``parse_log``-shaped dict: flat, noise-free series
        (so per-series verdicts are unambiguous), letting each test vary
        exactly one thing the aggregate is supposed to gate on.
        """
        t = list(np.linspace(0.0, duration, n)) if n > 0 else []
        zeros = [0.0] * n
        dt = [duration / n] * n if n > 0 else []
        nest_records = [
            {"t": ti, "phi": 1e-12, "phi_lid": 1e-12, "phi_closed": 1e-12,
             "misfit_rms": 1e-9, "gradp_zone": 1.0, "gradp_interior": 1.0,
             "gradp_ratio": 1.0}
            for ti in t
        ]
        return {
            "checksim_t": t, "checksim_dt": dt,
            "divmax": list(zeros), "divtot": list(zeros),
            "nest_records": nest_records,
            "n_freeze_warnings": n_freeze, "aborted_past_record": aborted,
        }

    def test_empty_log_is_not_pass(self):
        parsed = analyse_v6.parse_log(Path("/dev/null"))
        result = analyse_v6.analyse(parsed)
        self.assertNotEqual(result["overall_verdict"], "PASS")
        self.assertEqual(result["overall_verdict"], "NO-DATA")

    def test_sparse_log_is_not_pass(self):
        parsed = self._flat_parsed(n=2, duration=1.0)
        result = analyse_v6.analyse(parsed)
        self.assertNotEqual(result["overall_verdict"], "PASS")
        self.assertEqual(result["overall_verdict"], "INCONCLUSIVE")

    def test_all_inconclusive_series_is_not_pass(self):
        # A single sample: every series is below summarize_series' own n >= 3
        # floor, so every one of them comes back INCONCLUSIVE.
        parsed = self._flat_parsed(n=1, duration=0.0)
        result = analyse_v6.analyse(parsed)
        self.assertNotEqual(result["overall_verdict"], "PASS")
        self.assertIn(result["overall_verdict"], ("NO-DATA", "INCONCLUSIVE"))

    def test_missing_freeze_warning_is_not_pass(self):
        parsed = self._flat_parsed(n=20, duration=1000.0, n_freeze=0)
        result = analyse_v6.analyse(parsed, expected_runtime_s=1000.0)
        self.assertNotEqual(result["overall_verdict"], "PASS")
        self.assertEqual(result["overall_verdict"], "INCONCLUSIVE")

    def test_repeated_freeze_warning_is_not_pass(self):
        parsed = self._flat_parsed(n=20, duration=1000.0, n_freeze=2)
        result = analyse_v6.analyse(parsed, expected_runtime_s=1000.0)
        self.assertNotEqual(result["overall_verdict"], "PASS")
        self.assertEqual(result["overall_verdict"], "INCONCLUSIVE")

    def test_aborted_run_is_not_pass(self):
        parsed = self._flat_parsed(n=20, duration=1000.0, aborted=True)
        result = analyse_v6.analyse(parsed, expected_runtime_s=1000.0)
        self.assertNotEqual(result["overall_verdict"], "PASS")
        self.assertEqual(result["overall_verdict"], "FAIL")

    def test_short_of_intended_duration_is_not_pass(self):
        parsed = self._flat_parsed(n=20, duration=100.0)
        result = analyse_v6.analyse(parsed, expected_runtime_s=100000.0)
        self.assertNotEqual(result["overall_verdict"], "PASS")
        self.assertEqual(result["overall_verdict"], "INCONCLUSIVE")

    def test_too_few_samples_below_explicit_minimum_is_not_pass(self):
        parsed = self._flat_parsed(n=5, duration=1000.0)
        result = analyse_v6.analyse(parsed, expected_runtime_s=1000.0, min_samples=10)
        self.assertNotEqual(result["overall_verdict"], "PASS")
        self.assertEqual(result["overall_verdict"], "INCONCLUSIVE")

    def test_well_formed_run_passes(self):
        parsed = self._flat_parsed(n=20, duration=1000.0)
        result = analyse_v6.analyse(parsed, expected_runtime_s=1000.0)
        self.assertEqual(result["overall_verdict"], "PASS", result["overall_reasons"])

    def test_well_formed_run_with_no_expected_runtime_still_passes(self):
        # expected_runtime_s is optional -- a caller exploring a log with no
        # known intended duration should not be blocked on that check alone.
        parsed = self._flat_parsed(n=20, duration=1000.0)
        result = analyse_v6.analyse(parsed)
        self.assertEqual(result["overall_verdict"], "PASS", result["overall_reasons"])

    def test_drifting_series_fails_regardless_of_other_criteria(self):
        parsed = self._flat_parsed(n=20, duration=1000.0)
        parsed["divmax"] = [1e-12 + 5e-9 * ti for ti in parsed["checksim_t"]]
        result = analyse_v6.analyse(parsed, expected_runtime_s=1000.0)
        self.assertEqual(result["overall_verdict"], "FAIL")

    def test_steps_estimate_is_labelled_approximate(self):
        parsed = self._flat_parsed(n=20, duration=1000.0)
        result = analyse_v6.analyse(parsed)
        self.assertTrue(result["steps_estimate_is_approximate"])
        self.assertIn("approx.", analyse_v6.summary(result))


if __name__ == "__main__":
    unittest.main(verbosity=2)
