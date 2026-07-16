"""Pipeline-level robustness tests for udprep (review findings P6, P8).

P6: run_ibm must reload the freshly written solid masks and facet sections into
    the UDBase instance, so a first-ever preprocessing run does not feed the
    radiation step stale (empty) Sc/facsec.
P8: write_changed_params must only swallow the benign missing-namoptions case
    (as a visible warning); every other save_param failure must propagate.
"""
import sys
import types
import unittest
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

import numpy as np

from _common import PYTHON_DIR  # noqa: E402

from udbase import UDBase  # noqa: E402
from udprep.udprep_ibm import IBMSection  # noqa: E402
from udprep._section import Section  # noqa: E402

_DOMAIN_NAMOPTIONS = "\n".join(
    [
        "&DOMAIN",
        " itot = 4",
        " jtot = 3",
        " ktot = 2",
        " xlen = 40.0",
        " ylen = 30.0",
        " zsize = 20.0",
        "/",
    ]
) + "\n"

class TestRunIbmReloadsOutputs(unittest.TestCase):
    """P6: run_ibm reloads solid masks and facet sections into sim."""

    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)
        (self.workdir / "namoptions.001").write_text(_DOMAIN_NAMOPTIONS, encoding="ascii")

    def _fresh_sim(self):
        """A UDBase with no IBM output files -> empty facsec, Sc is None."""
        return UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

    def _make_section(self, sim):
        values = {
            "libm": True,
            "gen_geom": True,
            "geom_path": "",
            "stl_ground": True,
            "diag_neighbs": True,
            "nompthreads": 1,
            "ibmtol": 5e-4,
            "ray_dir_u": [0.0, 0.0, 1.0],
            "ray_dir_v": [0.0, 0.0, 1.0],
            "ray_dir_w": [0.0, 0.0, 1.0],
            "ray_dir_c": [0.0, 0.0, 1.0],
        }
        section = IBMSection("ibm", values, sim=sim, defaults={})
        # run_ibm's _require_sim wants a geometry when libm is on; the heavy
        # computation is mocked, so a placeholder object is enough.
        sim.geom = object()
        return section

    def _write_plausible_outputs(self):
        """Mimic the IBM backend writing small solid/facet-section files."""
        # solid_*.txt: header + one 1-based (i, j, k) index inside the domain.
        for grid in ("u", "v", "w", "c"):
            (self.workdir / f"solid_{grid}.txt").write_text(
                "# i j k\n1 1 1\n", encoding="ascii"
            )
            # fluid_boundary_*.txt: header + one 1-based (i, j, k) fluid cell.
            (self.workdir / f"fluid_boundary_{grid}.txt").write_text(
                "# i j k\n1 1 1\n", encoding="ascii"
            )
            # facet_sections_*.txt: header + (facid, area, fluid-loc, distance).
            (self.workdir / f"facet_sections_{grid}.txt").write_text(
                "# facid area loc distance\n1 0.5 1 0.1\n", encoding="ascii"
            )

    def test_run_ibm_populates_facsec_and_solid_masks(self):
        sim = self._fresh_sim()
        # Precondition: a fresh case has empty facet sections and no solid mask.
        self.assertEqual(sim.facsec, {})
        self.assertIsNone(sim.Sc)

        section = self._make_section(sim)

        # Mock one level deeper than run_ibm: the f2py backend writes the output
        # files and returns the 13-element count vector, so run_ibm's own
        # bookkeeping and reload code still execute.
        counts = np.ones(13, dtype=int)

        def fake_f2py():
            self._write_plausible_outputs()
            return counts

        with mock.patch.object(section, "_run_ibm_via_f2py", side_effect=fake_f2py):
            section.run_ibm(backend="f2py")

        # After run_ibm the reloaded state must be populated.
        self.assertIn("c", sim.facsec)
        self.assertGreater(len(sim.facsec), 0)
        self.assertIsNotNone(sim.Sc)
        self.assertTrue(sim.Sc.any())

    def test_run_ibm_reload_runs_for_legacy_backend(self):
        sim = self._fresh_sim()
        section = self._make_section(sim)

        counts = np.ones(13, dtype=int)

        def fake_legacy():
            self._write_plausible_outputs()

        with mock.patch.object(section, "_write_ibm_input_files"), \
                mock.patch.object(section, "_run_ibm_via_legacy", side_effect=fake_legacy), \
                mock.patch.object(section, "_update_counts_from_info_fort"):
            # _update_counts_from_info_fort is mocked out (no info_fort.txt), so
            # only the reload is under test for the legacy path.
            section.run_ibm(backend="legacy")

        self.assertIn("c", sim.facsec)
        self.assertIsNotNone(sim.Sc)
        self.assertTrue(sim.Sc.any())

class TestWriteChangedParamsErrorHandling(unittest.TestCase):
    """P8: write_changed_params only swallows the benign missing-file case."""

    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)

    def _section(self):
        sim = types.SimpleNamespace(path=self.workdir, expnr="001", alpha=2)
        return Section("dummy", {"alpha": 2}, sim=sim, defaults={"alpha": 1})

    def test_missing_namoptions_warns_and_does_not_raise(self):
        # No namoptions.001 on disk -> save_param raises FileNotFoundError, which
        # is the one benign case: it must surface as a visible warning, not crash.
        section = self._section()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            section.write_changed_params()
        self.assertTrue(
            any("skipping writeback" in str(w.message) for w in caught),
            msg=f"expected a skipping-writeback warning, got {[str(w.message) for w in caught]}",
        )

    def test_genuine_bug_propagates(self):
        # A save_param that raises for a non-file reason (here a TypeError bug)
        # must propagate rather than be silently swallowed.
        section = self._section()
        with mock.patch.object(section, "save_param", side_effect=TypeError("boom")):
            with self.assertRaises(TypeError):
                section.write_changed_params()

    def test_permission_error_propagates(self):
        # A write that fails for permissions is a real problem and must not be
        # hidden; only FileNotFoundError is treated as benign.
        section = self._section()
        with mock.patch.object(section, "save_param", side_effect=PermissionError("denied")):
            with self.assertRaises(PermissionError):
                section.write_changed_params()

if __name__ == "__main__":
    unittest.main()
