import sys
import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

import numpy as np

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udprep.udprep import Section, SectionSpec, SKIP, UDPrep  # noqa: E402


class DummySection(Section):
    def ping(self):
        self.called = True

    def accept_named(self, force=False):
        self.force_seen = force

    def accept_kwargs(self, **kwargs):
        self.kwargs_seen = kwargs


class DummySim:
    def __init__(self, expnr="999", path=None):
        self.expnr = expnr
        self.path = Path(path or ".")
        self.saved = []
        self.foo = 10
        self.bar = 5

    def save_param(self, key, value):
        self.saved.append((key, value))


class FakeSection:
    def __init__(self, name, values, sim=None, defaults=None):
        self._name = name
        self.sim = sim
        self._defaults = defaults or {}
        self.values = dict(values)
        for key, val in values.items():
            setattr(self, key, val)
        self.run_all_calls = []
        self.write_changed_params_calls = 0
        self.show_changed_params_calls = 0

    def run_all(self, **kwargs):
        self.run_all_calls.append(kwargs)

    def write_changed_params(self):
        self.write_changed_params_calls += 1

    def show_changed_params(self):
        self.show_changed_params_calls += 1

    def __repr__(self):
        return f"{self._name}:"


class TestSectionCore(unittest.TestCase):
    def test_resolve_default_supports_fraction_reference(self):
        ctx = types.SimpleNamespace(foo=12.0, bar=3.0)
        value = Section.resolve_default(
            "grid",
            "ratio",
            ctx,
            {"grid": {"ratio": "foo / bar"}},
            {},
        )
        self.assertEqual(value, 4.0)

    def test_resolve_default_falls_back_to_callable(self):
        ctx = types.SimpleNamespace(foo=8)
        value = Section.resolve_default(
            "grid",
            "bar",
            ctx,
            {},
            {"bar": lambda obj: obj.foo + 1},
        )
        self.assertEqual(value, 9)

    def test_resolve_default_returns_skip_for_missing(self):
        value = Section.resolve_default("grid", "missing", object(), {}, {})
        self.assertIs(value, SKIP)

    def test_run_steps_filters_kwargs_by_signature(self):
        section = DummySection("dummy", {})
        section.run_steps(
            "dummy",
            [
                ("accept_named", section.accept_named),
                ("accept_kwargs", section.accept_kwargs),
            ],
            force=True,
            ignored=123,
        )
        self.assertTrue(section.force_seen)
        self.assertEqual(section.kwargs_seen, {"force": True, "ignored": 123})

    def test_write_changed_params_skips_arrays_and_nested_collections(self):
        sim = DummySim()
        section = Section(
            "dummy",
            {
                "alpha": 2,
                "arr": np.array([1, 2]),
                "nested": [1, np.array([2, 3])],
            },
            sim=sim,
            defaults={"alpha": 1, "arr": np.array([0]), "nested": []},
        )
        section.write_changed_params()
        self.assertEqual(sim.saved, [("alpha", 2)])

    def test_changed_params_detects_numpy_arrays(self):
        section = Section(
            "dummy",
            {"arr": np.array([1, 2]), "same": np.array([3])},
            defaults={"arr": np.array([1, 3]), "same": np.array([3])},
        )
        changed = section._changed_params()
        self.assertEqual(len(changed), 1)
        self.assertEqual(changed[0][0], "arr")


class TestUDPrepCore(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)

    def _fake_udbase_module(self):
        module = types.ModuleType("udbase")

        class FakeUDBase(DummySim):
            def __init__(self, expnr, path=None, load_geometry=True):
                super().__init__(expnr=expnr, path=path)
                self.load_geometry = load_geometry
                self.existing = 4

        module.UDBase = FakeUDBase
        return module

    def test_init_populates_sections_from_defaults_and_sim(self):
        specs = [
            SectionSpec(
                name="alpha",
                fields=["foo", "bar"],
                defaults={"foo": lambda ctx: ctx.existing + 1, "bar": 7},
                section_cls=FakeSection,
            ),
        ]

        fake_module = self._fake_udbase_module()
        with mock.patch.dict(sys.modules, {"udbase": fake_module}):
            with mock.patch.object(UDPrep, "SECTION_SPECS", specs):
                prep = UDPrep("123", path=self.workdir, load_geometry=False)
        prep.SECTION_SPECS = specs
        self.assertEqual(prep.alpha.foo, 10)
        self.assertEqual(prep.alpha.bar, 5)

    def test_addvar_updates_specific_section(self):
        specs = [
            SectionSpec(
                name="alpha",
                fields=["foo", "bar"],
                defaults={"foo": 1, "bar": 2},
                section_cls=FakeSection,
            ),
        ]
        fake_module = self._fake_udbase_module()
        with mock.patch.dict(sys.modules, {"udbase": fake_module}):
            with mock.patch.object(UDPrep, "SECTION_SPECS", specs):
                prep = UDPrep("123", path=self.workdir)
        prep.SECTION_SPECS = specs
        prep.addvar("bar", 9, section="alpha")
        self.assertEqual(prep.alpha.bar, 9)
        self.assertEqual(prep.alpha.values["foo"], 10)

    def test_addvar_rejects_unknown_section_and_field(self):
        specs = [
            SectionSpec(name="alpha", fields=["foo"], defaults={"foo": 1}, section_cls=FakeSection),
        ]
        fake_module = self._fake_udbase_module()
        with mock.patch.dict(sys.modules, {"udbase": fake_module}):
            with mock.patch.object(UDPrep, "SECTION_SPECS", specs):
                prep = UDPrep("123", path=self.workdir)
        with self.assertRaises(ValueError):
            prep.addvar("foo", 1, section="missing")
        with self.assertRaises(ValueError):
            prep.addvar("bar", 1, section="alpha")

    def test_run_all_respects_section_gates(self):
        specs = [
            SectionSpec("grid", ["x"], {"x": 1}, FakeSection),
            SectionSpec("seb", ["x"], {"x": 1}, FakeSection),
            SectionSpec("ibm", ["libm", "gen_geom"], {"libm": True, "gen_geom": True}, FakeSection),
            SectionSpec("ic", ["x"], {"x": 1}, FakeSection),
            SectionSpec(
                "vegetation",
                ["ltrees", "ltreesfile"],
                {"ltrees": True, "ltreesfile": False},
                FakeSection,
            ),
            SectionSpec(
                "scalars",
                ["nsv", "lscasrc", "lscasrcl"],
                {"nsv": 0, "lscasrc": False, "lscasrcl": False},
                FakeSection,
            ),
            SectionSpec("radiation", ["x"], {"x": 1}, FakeSection),
        ]
        fake_module = self._fake_udbase_module()
        with mock.patch.dict(sys.modules, {"udbase": fake_module}):
            with mock.patch.object(UDPrep, "SECTION_SPECS", specs):
                prep = UDPrep("123", path=self.workdir)
        prep.SECTION_SPECS = specs

        prep.sim.lEB = True
        prep.ibm.generate_lscale = mock.Mock()
        prep.ibm.write_lscale = mock.Mock()
        prep.ibm.generate_factypes = mock.Mock()
        prep.ibm.write_factypes = mock.Mock()
        prep.ibm.run_ibm_fortran = mock.Mock()
        prep.ibm.copy_geom_outputs = mock.Mock()
        prep.ibm.write_facets = mock.Mock()
        prep.ibm.write_facetarea = mock.Mock()
        prep.vegetation.save = mock.Mock()
        prep.radiation.run_all = mock.Mock()
        prep.seb.run_all = mock.Mock()
        prep.vegetation.ltrees = True
        prep.vegetation.ltreesfile = False

        prep.run_all(force=True)

        self.assertEqual(len(prep.grid.run_all_calls), 1)
        self.assertEqual(len(prep.ic.run_all_calls), 1)
        self.assertEqual(len(prep.vegetation.run_all_calls), 1)
        self.assertEqual(len(prep.scalars.run_all_calls), 0)
        prep.ibm.generate_lscale.assert_called_once()
        prep.ibm.write_lscale.assert_called_once()
        prep.ibm.generate_factypes.assert_called_once()
        prep.ibm.write_factypes.assert_called_once()
        prep.ibm.run_ibm_fortran.assert_called_once()
        prep.ibm.write_facets.assert_called_once()
        prep.ibm.write_facetarea.assert_called_once()
        prep.vegetation.save.assert_called_once()
        prep.radiation.run_all.assert_called_once_with(force=True)
        prep.seb.run_all.assert_called_once()

    def test_write_changed_params_and_show_changed_params_visit_all_sections(self):
        specs = [
            SectionSpec("alpha", ["foo"], {"foo": 1}, FakeSection),
            SectionSpec("beta", ["bar"], {"bar": 2}, FakeSection),
        ]
        fake_module = self._fake_udbase_module()
        with mock.patch.dict(sys.modules, {"udbase": fake_module}):
            with mock.patch.object(UDPrep, "SECTION_SPECS", specs):
                prep = UDPrep("123", path=self.workdir)
        prep.SECTION_SPECS = specs
        prep.write_changed_params()
        prep.show_changed_params()
        self.assertEqual(prep.alpha.write_changed_params_calls, 1)
        self.assertEqual(prep.beta.write_changed_params_calls, 1)
        self.assertEqual(prep.alpha.show_changed_params_calls, 1)
        self.assertEqual(prep.beta.show_changed_params_calls, 1)


if __name__ == "__main__":
    unittest.main()
