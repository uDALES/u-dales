#!/usr/bin/env python3
"""Smoke test of the parent-side zone dump (``&NESTPARENT``, plan item D1).

Runs the ``tiny`` parent once with **both** ``lfielddump`` and ``lnestparent``
at the same cadence, then builds the child's ``nesting.inp`` from each source
and requires the two files to be bit-identical: the band files carry exactly
what the slab cut takes from the full dumps, promoted the same way.  It also
reports the storage ratio between the two outputs, checks that the parent
printed its I/O cost, and builds a refined child (``dx_child = dx_parent / 2``,
the V0 tiny ``r2-coarse`` point) from a coarse driver that wrote only the band
sized by :meth:`config.Preset.nestparent_nzone` -- with the interior ``NaN``, a
band too thin for the prolongation stencil would fail that build.

Run directly, or through ``tests/run_tests.py nesting-validation``:

    module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
    source /rds/general/user/mvr/home/udales/.venv/bin/activate
    python tests/validation/nesting/test_nestparent_tiny.py

Environment: ``UDALES_BUILD`` (solver binary), ``UDALES_RUNTIME_MODULES``,
``MPIEXEC``, ``UDALES_NESTPARENT_RUNDIR`` (where to work; kept if set),
``UDALES_NESTPARENT_KEEP=1`` (keep a temporary run directory too).
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
from caselib import NestParent, run_solver
from config import get_preset, get_suite
from make_child_case import DrivingParent


def _run_parent(rundir: Path, preset, **build_kwargs) -> Path:
    """Build, spin up and run the production phase of one parent."""
    parent_dir = rundir / preset.parent_expnr
    make_parent_case.build(rundir, preset, **build_kwargs)
    nr = preset.parent_expnr
    nproc = preset.nprocx * preset.nprocy
    run_solver(parent_dir, f"namoptions_spinup.{nr}", nproc, parent_dir / "spinup.log")
    run_v1._set_startfile(parent_dir / f"namoptions.{nr}",
                          run_v1._restart_file(parent_dir, nr))
    run_solver(parent_dir, f"namoptions.{nr}", nproc, parent_dir / "production.log")
    return parent_dir


def _size(paths) -> int:
    return sum(p.stat().st_size for p in paths)



def _set_namelist(path: Path, **values) -> None:
    """Set existing scalar keys in a namelist in place."""
    text = path.read_text(encoding="ascii")
    for key, value in values.items():
        text, n = re.subn(rf"(?m)^(\s*{key}\s*=\s*).*$",
                          lambda m: m.group(1) + str(value), text)
        if n != 1:
            raise RuntimeError(f"{key}: {n} substitutions in {path}")
    path.write_text(text, encoding="ascii")


def _band_times(casedir: Path, expnr: str):
    import netCDF4
    path = sorted(casedir.glob(f"nesting.out.???.???.{expnr}.nc"))[0]
    with netCDF4.Dataset(path) as ds:
        return [float(t) for t in ds.variables["time"][:]]


def _init_time(casedir: Path, expnr: str) -> float:
    import netCDF4
    path = sorted(casedir.glob(f"nesting.out.init.???.???.{expnr}.nc"))[0]
    with netCDF4.Dataset(path) as ds:
        return float(ds.variables["time"][:])


class TestNestParentTiny(unittest.TestCase):
    """One parent run with both outputs, two child builds, one refined build."""

    @classmethod
    def setUpClass(cls) -> None:
        binary = caselib.solver_binary()
        if not binary.exists():
            raise unittest.SkipTest(f"solver binary not found at {binary}")
        cls.preset = replace(get_preset("tiny"), name="tiny-nestparent", parent_output="both")
        cls.preset.validate()
        override = os.environ.get("UDALES_NESTPARENT_RUNDIR")
        if override:
            cls.rundir = Path(override)
            cls.rundir.mkdir(parents=True, exist_ok=True)
            cls._temp = None
        else:
            cls._temp = tempfile.mkdtemp(prefix="udales-nestparent-tiny-")
            cls.rundir = Path(cls._temp)
        p = cls.preset
        cls.parent_dir = _run_parent(cls.rundir, p)
        cls.production_log = (cls.parent_dir / "production.log").read_text(errors="replace")

        cls.child_dirs = {}
        cls.manifests = {}
        for source in ("fielddump", "nestparent"):
            outdir = cls.rundir / f"from_{source}"
            driving = DrivingParent.matched(cls.parent_dir, p, source=source)
            try:
                casedir = make_child_case.build(cls.parent_dir, outdir, p, driving=driving)
            finally:
                driving.dump.close()
            cls.child_dirs[source] = casedir
            cls.manifests[source] = json.loads((casedir / "manifest.json").read_text())

    @classmethod
    def tearDownClass(cls) -> None:
        if cls._temp and not os.environ.get("UDALES_NESTPARENT_KEEP"):
            shutil.rmtree(cls._temp, ignore_errors=True)
        elif cls._temp:
            print(f"run directory kept at {cls._temp}")

    # -- (1) bit identity ---------------------------------------------------- #

    def test_the_two_sources_give_bit_identical_nesting_files(self):
        from netCDF4 import Dataset

        nr = self.preset.child_expnr
        a = self.child_dirs["fielddump"] / f"nesting.inp.{nr}.nc"
        b = self.child_dirs["nestparent"] / f"nesting.inp.{nr}.nc"
        names = ["time", "net_volume_flux", "flux_residual",
                 "u_init", "v_init", "w_init"]
        names += [f"{c}_{f}" for f in ("west", "east", "south", "north") for c in "uvw"]
        with Dataset(a, "r") as da, Dataset(b, "r") as db:
            self.assertEqual(len(da.dimensions["time"]), len(db.dimensions["time"]))
            self.assertGreaterEqual(len(da.dimensions["time"]), 4)
            for name in names:
                x = np.asarray(da.variables[name][:])
                y = np.asarray(db.variables[name][:])
                self.assertEqual(x.shape, y.shape, name)
                self.assertTrue(np.array_equal(x, y),
                                f"{name} differs between the fielddump and nestparent builds "
                                f"(max |diff| = {np.max(np.abs(x - y)):.3e})")
            for key in ("itot", "jtot", "ktot", "nzone", "has_initial_condition",
                        "fluid_lateral_area"):
                self.assertEqual(da.getncattr(key), db.getncattr(key), key)

    def test_the_manifests_agree_except_for_the_source(self):
        mf, mn = self.manifests["fielddump"], self.manifests["nestparent"]
        self.assertEqual(mf["driving_source"], "fielddump")
        self.assertEqual(mn["driving_source"], "nestparent")
        self.assertEqual(set(mf), set(mn), "the two manifests carry different keys")
        for key in ("n_parent_levels", "t_offset", "runtime", "parent_dt_median",
                    "flux_residual_before_correction", "flux_residual_after_correction",
                    "initial_condition_divmax", "cadence"):
            self.assertEqual(mf[key], mn[key], key)
        self.assertLess(mn["flux_residual_after_correction"]["max_abs_normalised"], 1.0e-12)

    # -- (2) storage --------------------------------------------------------- #

    def test_the_band_files_are_much_smaller_than_the_field_dumps(self):
        nr = self.preset.parent_expnr
        field = _size(self.parent_dir.glob(f"fielddump.???.???.{nr}.nc"))
        band = _size(self.parent_dir.glob(f"nesting.out.???.???.{nr}.nc"))
        init = _size(self.parent_dir.glob(f"nesting.out.init.???.???.{nr}.nc"))
        self.assertGreater(band, 0)
        self.assertGreater(init, 0)
        ratio = field / (band + init)
        p = self.preset
        with NestParent(self.parent_dir, nr, p.dx) as d:
            nb, ni, nj = d.nzone, d.ni, d.nj
        ideal = (p.itot * p.jtot) / (ni * nj - (ni - 2 * nb) * (nj - 2 * nb))
        print(f"\n[nestparent] fielddump {field / 1e6:.2f} MB, nestparent band {band / 1e6:.2f} MB "
              f"+ init {init / 1e6:.2f} MB: ratio {ratio:.2f}x "
              f"(cell-count ideal for this geometry {ideal:.2f}x; production case ~15x)")
        # The tiny geometry (96^2 parent, 64^2 box, 8-cell band) can only give
        # ~5x by cell count; corner overlaps and the one-off init block take a
        # little more.  The production case's ~15x is a matter of geometry.
        self.assertGreater(ratio, 3.0, f"nestparent only {ratio:.2f}x smaller")

    # -- the parent's own accounting ----------------------------------------- #

    def test_the_parent_printed_its_io_cost(self):
        first = re.search(r"nestparent: dump 1 at t = \S+ s: (\S+) MB written \(all ranks\) "
                          r"in (\S+) s", self.production_log)
        total = re.search(r"nestparent: (\d+) dumps, (\S+) MB in total \(all ranks\), (\S+) s",
                          self.production_log)
        self.assertIsNotNone(first, "no per-dump cost line in the production log")
        self.assertIsNotNone(total, "no run-total cost line in the production log")
        print(f"\n[nestparent] parent: first dump {first.group(1)} MB in {first.group(2)} s; "
              f"{total.group(1)} dumps, {total.group(2)} MB, {total.group(3)} s in writes")
        self.assertEqual(int(total.group(1)), self.manifests["nestparent"]["cadence"]["n_levels_dumped"])
        # Exactly one report per 100 dumps plus the first: the tiny run has < 100.
        self.assertEqual(len(re.findall(r"nestparent: dump \d+ at t", self.production_log)), 1)

    def test_the_band_is_the_zone_plus_one_cell(self):
        p = self.preset
        with NestParent(self.parent_dir, p.parent_expnr, p.dx) as d:
            self.assertEqual(d.nzone, p.nestparent_nzone())
            self.assertEqual(d.nzone, p.zone_cells + 1)
            self.assertEqual((d.i0, d.j0, d.ni, d.nj),
                             (p.child_i0, p.child_j0, p.child_itot, p.child_jtot))
            u, v, w = d.read_level(0)
            # NaN in the interior, finite in the band
            nb = d.nzone
            self.assertTrue(np.isnan(u[nb + 1:-nb - 1, nb + 1:-nb - 1]).all())
            self.assertTrue(np.isfinite(u[:nb + 1]).all() and np.isfinite(u[-nb - 1:]).all())
            self.assertTrue(np.isfinite(v[:, :nb + 1]).all() and np.isfinite(v[:, -nb - 1:]).all())
            self.assertTrue(np.isfinite(w[:nb]).all() and np.isfinite(w[:, -nb:]).all())
            self.assertEqual(w.shape[2], p.ktot + 1)
            self.assertTrue((w[:nb, :, -1] == 0.0).all(), "w at the lid is not zero")

    # -- (3) a refined child from a band-only driver -------------------------- #

    def test_a_refined_child_builds_from_a_band_only_driver(self):
        """``dx_child = dx_parent / 2``: the V0 tiny ``r2-coarse`` point, driver
        writing **only** the band sized for the fine child.  With the interior
        ``NaN``, the build raises if any slab or prolongation stencil reached
        past the band -- so completing it is the index-range assertion."""
        suite = get_suite("v0-tiny")
        point = next(pt for pt in suite.points if pt.key == "r2-coarse")
        driver = replace(point.driver, parent_output="nestparent")
        point = replace(point, driver=driver)
        child = point.child
        self.assertEqual(driver.dx, 2 * child.dx)
        self.assertEqual(driver.nestparent_nzone(child),
                         int(np.ceil(child.nzone * child.dx / driver.dx)) + 1)
        rundir = self.rundir / "refined"
        driver_dir = _run_parent(rundir, driver, nestparent_child=child)
        self.assertFalse(list(driver_dir.glob(f"fielddump.*.{driver.parent_expnr}.nc")),
                         "the band-only driver wrote field dumps")
        driving = DrivingParent.refined(driver_dir, point, source="nestparent")
        try:
            self.assertTrue(driving.interpolates)
            casedir = make_child_case.build(driver_dir, rundir, child, driving=driving)
        finally:
            driving.dump.close()
        manifest = json.loads((casedir / "manifest.json").read_text())
        self.assertEqual(manifest["driving_source"], "nestparent")
        self.assertEqual(manifest["refinement"]["spatial"], 2)
        self.assertLess(manifest["flux_residual_after_correction"]["max_abs_normalised"],
                        1.0e-12)
        from udprep.nesting import validate_nesting_file

        attrs = validate_nesting_file(casedir / f"nesting.inp.{child.child_expnr}.nc")
        self.assertEqual(int(attrs["has_initial_condition"]), 1)
        # A band one cell thinner is NOT enough: the outermost zone cell's
        # slope needs its inner neighbour.  Shown by masking the band down.
        with NestParent(driver_dir, driver.parent_expnr, driver.dx) as d:
            u, v, w = d.read_level(0)
            nb = d.nzone
            for arr in (u, v, w):
                arr[nb - 1:-(nb - 1), nb - 1:-(nb - 1), :] = np.nan
            from udprep.nesting import slabs_from_parent
            slabs = slabs_from_parent(driving.grid(child), u, v, w,
                                      child=make_child_case.child_grid(child),
                                      nzone=child.nzone)
            with self.assertRaises(ValueError):
                caselib.check_finite_slabs(slabs, nb - 1, driver.dx)



class TestNestParentRestart(unittest.TestCase):
    """A parent continued across two jobs must leave a usable input pair.

    ``open_band_file`` reopens an existing band file and appends, keeping its
    original time origin, but ``write_init_file`` used to recreate
    ``nesting.out.init.*.nc`` with ``NF90_CLOBBER`` on every startup, because
    ``linitdone`` was reset from ``nestparent_linit`` alone and never consulted
    whether this run was a continuation.  A normal continuation therefore
    replaced the original full-domain snapshot with one taken at the restart
    time while keeping the first segment's band times -- and
    :mod:`make_child_case` requires the snapshot to match the FIRST band
    level, so the pair became unusable *and* the original snapshot was gone.
    The append cursor also used ``modstat_nc``'s ``>=`` convention, which
    overwrote the boundary sample taken at exactly the restart time.

    A production parent is the case that needs this, and until now could not
    be built: ``make_parent_case`` gave the production phase ``trestart = 1e9``
    unconditionally, so it wrote no restart of its own and no campaign parent
    could be continued -- which is how the defect above reached production.
    ``Preset.production_trestart`` now overrides that, and this suite sets it.
    """

    SEG1_RUNTIME = 6.0
    SEG2_RUNTIME = 12.0
    CADENCE = 1.5

    @classmethod
    def setUpClass(cls) -> None:
        binary = caselib.solver_binary()
        if not binary.exists():
            raise unittest.SkipTest(f"solver binary not found at {binary}")
        # production_trestart is what makes a continuable parent possible at
        # all: without it the production phase gets trestart = 1e9 and writes no
        # restart, which is why nothing in the campaign exercised this path.
        cls.preset = replace(get_preset("tiny"), name="tiny-nestparent-restart",
                             parent_output="nestparent",
                             production_trestart=cls.SEG1_RUNTIME)
        cls.preset.validate()
        cls._temp = tempfile.mkdtemp(prefix="udales-nestparent-restart-")
        cls.rundir = Path(cls._temp)
        p, nr = cls.preset, cls.preset.parent_expnr
        nproc = p.nprocx * p.nprocy

        cls.parent_dir = make_parent_case.build(cls.rundir, p)
        nml = cls.parent_dir / f"namoptions.{nr}"
        run_solver(cls.parent_dir, f"namoptions_spinup.{nr}", nproc,
                   cls.parent_dir / "spinup.log")

        # segment 1 -- a short production phase that does write a restart
        run_v1._set_startfile(nml, run_v1._restart_file(cls.parent_dir, nr))
        # runtime and cadence still need patching (they are the preset's own
        # production window), but trestart now comes from the preset.
        _set_namelist(nml, runtime=cls.SEG1_RUNTIME, tnestparent=cls.CADENCE)
        # captured here: segment 2 patches trestart back to 1e9 below, so the
        # file cannot be read for this after setUpClass has finished.
        cls.seg1_trestart = caselib.read_namoption(nml, "trestart")
        cls.seg1_restart = run_v1._restart_file(cls.parent_dir, nr)
        run_solver(cls.parent_dir, f"namoptions.{nr}", nproc, cls.parent_dir / "seg1.log")
        cls.band1 = _band_times(cls.parent_dir, nr)
        cls.init1 = _init_time(cls.parent_dir, nr)

        # segment 2 -- continue in the SAME directory, as a queued job would
        run_v1._set_startfile(nml, run_v1._restart_file(cls.parent_dir, nr))
        _set_namelist(nml, runtime=cls.SEG2_RUNTIME, trestart=1.0e9)
        run_solver(cls.parent_dir, f"namoptions.{nr}", nproc, cls.parent_dir / "seg2.log")
        cls.band2 = _band_times(cls.parent_dir, nr)
        cls.init2 = _init_time(cls.parent_dir, nr)
        cls.seg2_log = (cls.parent_dir / "seg2.log").read_text(errors="replace")

    @classmethod
    def tearDownClass(cls) -> None:
        if cls._temp and not os.environ.get("UDALES_NESTPARENT_KEEP"):
            shutil.rmtree(cls._temp, ignore_errors=True)
        elif cls._temp:
            print(f"run directory kept at {cls._temp}")

    def test_the_original_initial_block_survives_the_restart(self):
        self.assertAlmostEqual(self.init2, self.init1, places=4,
                               msg=f"the snapshot moved {self.init1} -> {self.init2}")
        self.assertAlmostEqual(self.init2, self.band2[0], places=4,
                               msg="the snapshot must match the first band level")

    def test_the_preset_made_a_continuable_parent(self):
        """The restart the second segment starts from must exist because the
        PRESET asked for it, not because the test patched the namelist."""
        self.assertEqual(self.preset.production_trestart, self.SEG1_RUNTIME)
        self.assertIsNotNone(self.seg1_trestart,
                             "no trestart line in the production namelist")
        self.assertAlmostEqual(float(self.seg1_trestart), self.SEG1_RUNTIME, places=6)
        self.assertTrue(self.seg1_restart,
                        "the production phase wrote no restart to continue from")

    def test_the_continuation_says_so(self):
        self.assertIn("continuing an existing output series", self.seg2_log)
        self.assertIn("kept initial block matches the first band level", self.seg2_log)

    def test_the_first_segments_levels_are_untouched(self):
        self.assertEqual([round(t, 4) for t in self.band2[:len(self.band1)]],
                         [round(t, 4) for t in self.band1])

    def test_the_sample_at_the_restart_time_is_kept(self):
        """The last level of segment 1 is at the restart time; it must remain."""
        self.assertIn(round(self.band1[-1], 4), [round(t, 4) for t in self.band2])

    def test_the_record_is_contiguous_and_grew(self):
        self.assertGreater(len(self.band2), len(self.band1))
        self.assertEqual(sorted(self.band2), self.band2)
        self.assertEqual(len(set(round(t, 4) for t in self.band2)), len(self.band2))

    def test_a_child_builds_from_the_continued_pair(self):
        """The end-to-end consequence: make_child_case pairs init with band[0]."""
        p = self.preset
        outdir = self.rundir / "child_from_continued"
        driving = DrivingParent.matched(self.parent_dir, p, source="nestparent")
        try:
            casedir = make_child_case.build(self.parent_dir, outdir, p, driving=driving)
        finally:
            driving.dump.close()
        manifest = json.loads((casedir / "manifest.json").read_text())
        self.assertGreater(manifest["n_parent_levels"], 0)
        self.assertTrue((casedir / f"nesting.inp.{p.child_expnr}.nc").exists()
                        or any(casedir.glob("nesting.inp.*.nc")))


if __name__ == "__main__":
    unittest.main(verbosity=2)
