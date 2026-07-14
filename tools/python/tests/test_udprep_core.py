import contextlib
import sys
import types
import unittest
import warnings
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
from udprep.udprep_bcs import SPEC as BCS_SPEC  # noqa: E402
from udprep.udprep_ibm import IBMSection  # noqa: E402
from udprep.udprep_radiation import RadiationSection  # noqa: E402


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
        with mock.patch.object(section, "save_param") as mock_sp:
            section.write_changed_params()
        mock_sp.assert_called_once_with("alpha", 2)

    def test_changed_params_detects_numpy_arrays(self):
        section = Section(
            "dummy",
            {"arr": np.array([1, 2]), "same": np.array([3])},
            defaults={"arr": np.array([1, 3]), "same": np.array([3])},
        )
        changed = section._changed_params()
        self.assertEqual(len(changed), 1)
        self.assertEqual(changed[0][0], "arr")


class TestIBMSection(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)

    def _make_section(self, **overrides):
        sim = types.SimpleNamespace(
            path=self.workdir,
            expnr="321",
            geom=object(),
        )
        values = {
            "libm": True,
            "gen_geom": True,
            "geom_path": "",
            "iwallmom": 2,
            "ltempeq": False,
            "lmoist": False,
            "lwritefac": False,
            "stl_ground": True,
            "diag_neighbs": True,
            "nompthreads": 8,
            "ibmtol": 5e-4,
            "ray_dir_u": [0.0, 0.0, 1.0],
            "ray_dir_v": [0.0, 0.0, 1.0],
            "ray_dir_w": [0.0, 0.0, 1.0],
            "ray_dir_c": [0.0, 0.0, 1.0],
        }
        values.update(overrides)
        return IBMSection("ibm", values, sim=sim, defaults={})

    def _recording_section(self, **overrides):
        section = self._make_section(**overrides)
        calls = []
        section.run_ibm = mock.Mock(side_effect=lambda backend="f2py": calls.append(("run_ibm", backend)))
        section.write_facets = mock.Mock(side_effect=lambda: calls.append(("write_facets", None)))
        section.write_facets_unused = mock.Mock(side_effect=lambda: calls.append(("write_facets_unused", None)))
        section.write_facetarea = mock.Mock(side_effect=lambda: calls.append(("write_facetarea", None)))
        section.generate_factypes = mock.Mock(side_effect=lambda: calls.append(("generate_factypes", None)))
        section.write_factypes = mock.Mock(side_effect=lambda: calls.append(("write_factypes", None)))
        section.copy_geom_outputs = mock.Mock(side_effect=lambda: calls.append(("copy_geom_outputs", None)))
        return section, calls

    def test_run_all_orchestrates_ibm_preprocessing_modes(self):
        with self.subTest("generate geometry and default factypes"):
            section, calls = self._recording_section()
            section.run_all(backend="legacy")

            self.assertEqual(
                calls,
                [
                    ("run_ibm", "legacy"),
                    ("write_facets", None),
                    ("write_facetarea", None),
                    ("generate_factypes", None),
                    ("write_factypes", None),
                ],
            )

        with self.subTest("write unused facets when c facet sections are enabled"):
            section, calls = self._recording_section(ltempeq=True)
            section.run_all()

            self.assertEqual(
                calls,
                [
                    ("run_ibm", "f2py"),
                    ("write_facets", None),
                    ("write_facetarea", None),
                    ("write_facets_unused", None),
                    ("generate_factypes", None),
                    ("write_factypes", None),
                ],
            )

        with self.subTest("preserve existing factypes"):
            (self.workdir / "factypes.inp.321").write_text("existing\n", encoding="ascii")
            section, calls = self._recording_section()
            with self.assertWarns(UserWarning):
                section.run_all()

            self.assertEqual(
                calls,
                [
                    ("run_ibm", "f2py"),
                    ("write_facets", None),
                    ("write_facetarea", None),
                ],
            )

        with self.subTest("copy existing geometry outputs"):
            section, calls = self._recording_section(gen_geom=False)
            section.run_all()

            self.assertEqual(calls, [("copy_geom_outputs", None)])

        with self.subTest("copy existing geometry outputs with c sections"):
            section, calls = self._recording_section(gen_geom=False, ltempeq=True)
            section.run_all()

            self.assertEqual(calls, [("copy_geom_outputs", None)])

        with self.subTest("skip when IBM disabled"):
            section, calls = self._recording_section(libm=False)
            section.run_all()

            self.assertEqual(calls, [])

    def test_write_facets_unused_records_facets_without_c_sections(self):
        section = self._make_section(ltempeq=True)
        section.sim.nfcts = 5
        section.sim.calculate_facet_sections_c = True
        (self.workdir / "facet_sections_c.txt").write_text(
            "# facet area flux point distance\n"
            "1 0.5 10 0.1\n"
            "3 0.7 11 0.2\n"
            "3 0.2 12 0.3\n",
            encoding="ascii",
        )

        section.write_facets_unused()

        self.assertEqual((self.workdir / "facets_unused.321").read_text(encoding="ascii"), "2\n4\n5\n")

    def test_write_facets_unused_writes_empty_file_when_all_facets_are_used(self):
        section = self._make_section(ltempeq=True)
        section.sim.nfcts = 3
        section.sim.calculate_facet_sections_c = True
        (self.workdir / "facet_sections_c.txt").write_text(
            "# facet area flux point distance\n"
            "1 0.5 10 0.1\n"
            "2 0.7 11 0.2\n"
            "3 0.2 12 0.3\n",
            encoding="ascii",
        )

        section.write_facets_unused()

        self.assertEqual((self.workdir / "facets_unused.321").read_text(encoding="ascii"), "")

    def test_copy_geom_outputs_copies_existing_facets_unused(self):
        geom_path = self.workdir / "geom"
        geom_path.mkdir()
        for name in ("solid_c.txt", "fluid_boundary_c.txt", "facet_sections_c.txt"):
            (geom_path / name).write_text("# header\n1 2 3\n", encoding="ascii")
        for name in ("facetarea.inp.999", "facets.inp.999", "factypes.inp.999"):
            (geom_path / name).write_text("# header\n1\n", encoding="ascii")
        (geom_path / "facets_unused.999").write_text("2,4\n", encoding="ascii")

        section = self._make_section(
            gen_geom=False,
            geom_path=str(geom_path),
            iwallmom=1,
            ltempeq=True,
        )
        section.sim.geom = types.SimpleNamespace(stl=types.SimpleNamespace(faces=[0, 1, 2, 3]))
        section.sim.calculate_facet_sections_uvw = False
        section.sim.calculate_facet_sections_c = True
        section.save_param = mock.Mock()

        section.copy_geom_outputs()

        self.assertEqual((self.workdir / "facets_unused.321").read_text(encoding="ascii"), "2,4\n")

    def test_write_facets_reads_case_relative_headerless_type_file(self):
        section = self._make_section(read_types=True, types_path="geom_facet_types.txt")
        section.sim.geom = types.SimpleNamespace(
            stl=types.SimpleNamespace(
                faces=np.array([[0, 1, 2], [2, 3, 0], [0, 3, 4]], dtype=int),
                face_normals=np.array(
                    [
                        [0.0, 0.0, 1.0],
                        [0.0, 1.0, 0.0],
                        [1.0, 0.0, 0.0],
                    ],
                    dtype=float,
                ),
            )
        )
        (self.workdir / "geom_facet_types.txt").write_text("1\n2\n4\n", encoding="ascii")

        section.write_facets()

        self.assertEqual(section.sim.facet_types.tolist(), [1, 2, 4])
        lines = (self.workdir / "facets.inp.321").read_text(encoding="ascii").splitlines()
        self.assertEqual(len(lines), 4)
        self.assertTrue(lines[1].startswith("   1"))
        self.assertTrue(lines[3].startswith("   4"))

    def test_write_facets_does_not_overwrite_existing_file(self):
        # Regression: an existing facets.inp is authored input. A re-run must not
        # clobber hand-assigned materials (previously it overwrote unconditionally,
        # destroying e.g. green-roof/asphalt-floor types with all-concrete).
        section = self._make_section()
        section.sim.geom = types.SimpleNamespace(
            stl=types.SimpleNamespace(
                faces=np.array([[0, 1, 2], [2, 3, 0], [0, 3, 4]], dtype=int),
                face_normals=np.array(
                    [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=float
                ),
            )
        )
        existing = (
            "# type, normal\n"
            "  -1 0.0000 0.0000 1.0000\n"
            "   1 0.0000 1.0000 0.0000\n"
            "  12 1.0000 0.0000 0.0000\n"
        )
        facets_path = self.workdir / "facets.inp.321"
        facets_path.write_text(existing, encoding="ascii")

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            section.write_facets()

        # File is left untouched, and its types are loaded into state.
        self.assertEqual(facets_path.read_text(encoding="ascii"), existing)
        self.assertEqual(section.sim.facet_types.tolist(), [-1, 1, 12])
        self.assertEqual(section.sim.nfcts, 3)
        self.assertTrue(any("NOT overwriting" in str(w.message) for w in caught))

    def test_load_facet_types_accepts_one_line_header(self):
        types_path = self.workdir / "types.txt"
        types_path.write_text("type\n1\n2\n4\n", encoding="ascii")

        values = IBMSection._load_facet_types(types_path, 3)

        self.assertEqual(values.tolist(), [1, 2, 4])

    def test_load_facet_types_rejects_wrong_count(self):
        types_path = self.workdir / "types.txt"
        types_path.write_text("1\n2\n", encoding="ascii")

        with self.assertRaisesRegex(ValueError, "expected 3"):
            IBMSection._load_facet_types(types_path, 3)

    def test_ray_direction_parses_lists_and_namoptions_strings(self):
        section = self._make_section(
            ray_dir_u=[1.0, 0.0, 0.0],
            ray_dir_v="0.0, 1.0, 0.0",
            ray_dir_w="0.0 0.0 -1.0",
            ray_dir_c=np.array([1.0, 1.0, 1.0]),
        )

        np.testing.assert_allclose(section._ray_direction("ray_dir_u"), [1.0, 0.0, 0.0])
        np.testing.assert_allclose(section._ray_direction("ray_dir_v"), [0.0, 1.0, 0.0])
        np.testing.assert_allclose(section._ray_direction("ray_dir_w"), [0.0, 0.0, -1.0])
        np.testing.assert_allclose(section._ray_direction("ray_dir_c"), [1.0, 1.0, 1.0])

    def test_ray_direction_rejects_invalid_values(self):
        section = self._make_section(ray_dir_u="1.0, 2.0")
        with self.assertRaisesRegex(ValueError, "exactly three"):
            section._ray_direction("ray_dir_u")

        section = self._make_section(ray_dir_u="0.0, 0.0, 0.0")
        with self.assertRaisesRegex(ValueError, "non-zero"):
            section._ray_direction("ray_dir_u")

        section = self._make_section(ray_dir_u="1.0, nope, 0.0")
        with self.assertRaisesRegex(ValueError, "numeric"):
            section._ray_direction("ray_dir_u")

    def test_write_ibm_input_files_uses_configured_ray_directions(self):
        section = self._make_section(
            ibmtol=1e-4,
            nompthreads=3,
            diag_neighbs=False,
            ray_dir_u="1.0, 0.0, 0.0",
            ray_dir_v=[0.0, 1.0, 0.0],
            ray_dir_w=[0.0, 0.0, -1.0],
            ray_dir_c=[1.0, 1.0, 1.0],
        )
        stl = types.SimpleNamespace(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2]], dtype=int),
        )
        section.sim.geom = types.SimpleNamespace(
            stl=stl,
            face_incenters=np.array([[0.25, 0.25, 0.0]]),
            face_normals=np.array([[0.0, 0.0, 1.0]]),
        )
        section.sim.dx = 2.0
        section.sim.dy = 3.0
        section.sim.itot = 4
        section.sim.jtot = 5
        section.sim.ktot = 2
        section.sim.zt = np.array([0.5, 1.5])
        section.sim.zm = np.array([0.0, 1.0])
        section.sim.BCxm = 1
        section.sim.BCym = 2

        section._write_ibm_input_files()

        lines = (self.workdir / "inmypoly_inp_info.txt").read_text(encoding="ascii").splitlines()
        rays = [[float(part) for part in lines[idx].split()] for idx in range(3, 7)]
        self.assertEqual(
            rays,
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, -1.0],
                [1.0, 1.0, 1.0],
            ],
        )


class TestUDPrepCore(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        # validate_expnr requires: dir name == namoptions suffix == iexpnr value.
        self.workdir = Path(self.temp_dir.name) / "123"
        self.workdir.mkdir()
        (self.workdir / "namoptions.123").write_text("&RUN\niexpnr = 123\n/\n")

    def _fake_udbase_module(self):
        module = types.ModuleType("udbase")

        class FakeUDBase(DummySim):
            def __init__(self, expnr, path=None, load_geometry=True, suppress_load_warnings=False):
                super().__init__(expnr=expnr, path=path)
                self.load_geometry = load_geometry
                self.existing = 4

        module.UDBase = FakeUDBase
        return module

    def _make_run_all_prep(self, *, libm=True, radiation_lEB=True, ltrees=False):
        specs = [
            SectionSpec("grid", ["x"], {"x": 1}, FakeSection),
            SectionSpec("forcing", ["x"], {"x": 1}, FakeSection),
            SectionSpec("seb", ["x"], {"x": 1}, FakeSection),
            SectionSpec("ibm", ["libm", "gen_geom"], {"libm": libm, "gen_geom": True}, FakeSection),
            SectionSpec(
                "vegetation",
                ["ltrees", "ltreesfile"],
                {"ltrees": ltrees, "ltreesfile": False},
                FakeSection,
            ),
            SectionSpec(
                "scalars",
                ["nsv", "lscasrc", "lscasrcl"],
                {"nsv": 0, "lscasrc": False, "lscasrcl": False},
                FakeSection,
            ),
            SectionSpec("radiation", ["lEB"], {"lEB": radiation_lEB}, FakeSection),
        ]
        fake_module = self._fake_udbase_module()
        with mock.patch.dict(sys.modules, {"udbase": fake_module}):
            with mock.patch.object(UDPrep, "SECTION_SPECS", specs):
                prep = UDPrep("123", path=self.workdir)
        prep.SECTION_SPECS = specs
        prep.ibm.libm = libm
        prep.radiation.lEB = radiation_lEB
        prep.vegetation.ltrees = ltrees
        prep.vegetation.ltreesfile = False
        prep.scalars.nsv = 0
        return prep

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

    def test_bcs_section_uses_defaults_and_preserves_namoptions_overrides(self):
        fake_module = self._fake_udbase_module()
        with mock.patch.dict(sys.modules, {"udbase": fake_module}):
            with mock.patch.object(UDPrep, "SECTION_SPECS", [BCS_SPEC]):
                prep = UDPrep("123", path=self.workdir, load_geometry=False)

        self.assertIsInstance(prep.bcs, Section)
        self.assertEqual(prep.bcs.BCxm, 1)
        self.assertEqual(prep.bcs.BCym, 1)
        self.assertEqual(prep.sim.BCxm, 1)
        self.assertEqual(prep.sim.BCym, 1)

        class OverrideUDBase(DummySim):
            def __init__(self, expnr, path=None, load_geometry=True, suppress_load_warnings=False):
                super().__init__(expnr=expnr, path=path)
                self.BCxm = 2
                self.BCym = 3

        override_module = types.ModuleType("udbase")
        override_module.UDBase = OverrideUDBase
        with mock.patch.dict(sys.modules, {"udbase": override_module}):
            with mock.patch.object(UDPrep, "SECTION_SPECS", [BCS_SPEC]):
                prep = UDPrep("123", path=self.workdir, load_geometry=False)

        self.assertEqual(prep.bcs.BCxm, 2)
        self.assertEqual(prep.bcs.BCym, 3)
        self.assertEqual(prep.sim.BCxm, 2)
        self.assertEqual(prep.sim.BCym, 3)

    def test_run_all_respects_section_gates(self):
        prep = self._make_run_all_prep(libm=True, radiation_lEB=True, ltrees=True)
        prep.sim.lEB = True
        prep.vegetation.save = mock.Mock()
        prep.radiation.run_all = mock.Mock()
        prep.seb.run_all = mock.Mock()

        prep.run_all(force=True, ibm_backend="legacy")

        self.assertEqual(len(prep.grid.run_all_calls), 1)
        self.assertEqual(len(prep.forcing.run_all_calls), 1)
        self.assertEqual(len(prep.vegetation.run_all_calls), 1)
        self.assertEqual(len(prep.scalars.run_all_calls), 0)
        self.assertEqual(prep.ibm.run_all_calls, [{"backend": "legacy"}])
        prep.radiation.run_all.assert_called_once_with(force=True)
        prep.seb.run_all.assert_called_once()

    def test_run_all_uses_radiation_lEB_for_current_radiation_gate(self):
        prep = self._make_run_all_prep(libm=True, radiation_lEB=True)
        prep.sim.lEB = False
        prep.radiation.run_all = mock.Mock()
        prep.seb.run_all = mock.Mock()

        prep.run_all(force=True, ibm_backend="legacy")

        self.assertEqual(prep.ibm.run_all_calls, [{"backend": "legacy"}])
        prep.radiation.run_all.assert_called_once_with(force=True)
        prep.seb.run_all.assert_called_once()

    def test_run_all_skips_radiation_and_seb_when_radiation_lEB_is_false(self):
        prep = self._make_run_all_prep(libm=True, radiation_lEB=False)
        prep.sim.lEB = True
        prep.radiation.run_all = mock.Mock()
        prep.seb.run_all = mock.Mock()

        prep.run_all(force=True, ibm_backend="legacy")

        self.assertEqual(prep.ibm.run_all_calls, [{"backend": "legacy"}])
        prep.radiation.run_all.assert_not_called()
        prep.seb.run_all.assert_not_called()

    def test_run_all_skips_radiation_and_seb_when_ibm_is_disabled(self):
        prep = self._make_run_all_prep(libm=False, radiation_lEB=True)
        prep.radiation.run_all = mock.Mock()
        prep.seb.run_all = mock.Mock()

        prep.run_all(force=True, ibm_backend="legacy")

        self.assertEqual(prep.ibm.run_all_calls, [])
        prep.radiation.run_all.assert_not_called()
        prep.seb.run_all.assert_not_called()

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


class TestRadiationSection(unittest.TestCase):
    def test_run_short_wave_skip_removes_full_vf_text_intermediate(self):
        with TemporaryDirectory() as temp_dir:
            case_dir = Path(temp_dir)
            sim = DummySim(expnr="321", path=case_dir)
            sim.lvfsparse = False
            section = RadiationSection("radiation", {"view3d_out": 0}, sim=sim, defaults={})

            (case_dir / "Sdir.txt").write_text("1.0\n", encoding="ascii")
            (case_dir / "netsw.inp.321").write_text("1.0\n", encoding="ascii")
            (case_dir / "vf.nc.inp.321").write_bytes(b"netcdf placeholder")
            vf_path = case_dir / "vf.txt"
            vf_path.write_text("stale view3d text output\n", encoding="ascii")

            section.run_short_wave()

            self.assertFalse(vf_path.exists())

    def test_write_netsw_only_writes_sveg_when_nonempty(self):
        with TemporaryDirectory() as temp_dir:
            case_dir = Path(temp_dir)
            sim = DummySim(expnr="321", path=case_dir)
            section = RadiationSection("radiation", {}, sim=sim, defaults={})

            section.write_netsw(np.array([1.0]), s_veg=np.array([]))
            self.assertFalse((case_dir / "sveg.inp.321").exists())

            section.write_netsw(np.array([1.0]), s_veg=np.array([2.0]))
            self.assertTrue((case_dir / "sveg.inp.321").exists())

    def test_scanline_knet_uses_sdir_file_precision(self):
        section = RadiationSection("radiation", {}, sim=DummySim())
        full_sdir = np.array([1.234, 5.675, 9.999], dtype=float)

        def fake_calc_direct_sw(*_args, **_kwargs):
            # Mirror the REAL solver contract: scanline returns an empty
            # veg-absorption array (np.zeros(0)), never None.
            return full_sdir.copy(), np.zeros(0, dtype=float), {}

        seen = {}

        def fake_calc_reflections_sw(sdir, *_args, **_kwargs):
            seen["sdir"] = sdir.copy()
            return sdir + 1.0

        section.calc_direct_sw = fake_calc_direct_sw
        section.calc_reflections_sw = fake_calc_reflections_sw

        sdir, knet, s_veg = section._compute_knet(
            np.array([0.0, 0.0, 1.0]),
            800.0,
            250.0,
            "scanline",
            0.1,
            True,
            np.zeros(3),
            object(),
            np.ones(3),
            None,
        )

        expected = np.round(full_sdir, 2)
        np.testing.assert_allclose(seen["sdir"], expected)
        np.testing.assert_allclose(sdir, expected)
        np.testing.assert_allclose(knet, expected + 1.0)
        self.assertIsInstance(s_veg, np.ndarray)
        self.assertEqual(s_veg.size, 0)

    # ------------------------------------------------------------------
    # Item 1 (P1): run_short_wave_timedep vegetation array-shape contract
    # ------------------------------------------------------------------
    def _run_timedep(self, tmp, nfcts, ltrees, s_veg_value):
        """Drive run_short_wave_timedep with the solver mocked to honour the
        REAL contract: _compute_knet returns an s_veg array (never None)."""
        sim = types.SimpleNamespace(
            path=Path(tmp),
            expnr="001",
            ltrees=ltrees,
            geom=types.SimpleNamespace(
                stl=types.SimpleNamespace(face_normals=np.zeros((nfcts, 3)))
            ),
        )
        sim.assign_prop_to_fac = lambda _prop: np.ones(nfcts)

        section = RadiationSection("radiation", {}, sim=sim, defaults={})
        section.ltimedepsw = True
        section.lEB = False
        section.isolar = 2
        section.ishortwave = 2
        section.directsw_method = "moller"
        section.runtime = 20.0
        section.dtSP = 10.0
        section.year = 2020
        section.month = 6
        section.day = 21
        section.hour = 12
        section.minute = 0
        section.second = 0

        section._solar_state_time = mock.Mock(
            return_value=(np.array([0.0, 0.0, 1.0]), 30.0, 0.0, 800.0, 100.0)
        )
        section._compute_knet = mock.Mock(
            return_value=(np.ones(nfcts), np.ones(nfcts), s_veg_value)
        )
        section._write_sdir_nc = mock.Mock()
        section.write_timedepsw = mock.Mock()
        section.write_timedepsveg = mock.Mock()

        section.run_short_wave_timedep()
        return section

    def test_timedep_no_veg_completes_without_valueerror(self):
        with TemporaryDirectory() as tmp:
            section = self._run_timedep(
                tmp, nfcts=4, ltrees=False, s_veg_value=np.zeros(0, dtype=float)
            )
        section.write_timedepsw.assert_called_once()
        _tSP, knet = section.write_timedepsw.call_args.args
        self.assertEqual(knet.shape, (4, 3))
        section._write_sdir_nc.assert_called_once()
        section.write_timedepsveg.assert_not_called()

    def test_timedep_veg_stores_nveg_rows(self):
        nfcts, nveg = 4, 7
        with TemporaryDirectory() as tmp:
            section = self._run_timedep(
                tmp, nfcts=nfcts, ltrees=True, s_veg_value=np.arange(nveg, dtype=float)
            )
        section.write_timedepsveg.assert_called_once()
        _tSP, sveg = section.write_timedepsveg.call_args.args
        self.assertEqual(sveg.shape, (nveg, 3))
        np.testing.assert_allclose(sveg[:, 0], np.arange(nveg))

    # ------------------------------------------------------------------
    # Item 2 (P3): calc_view_factors honours the per-call maxD argument
    # ------------------------------------------------------------------
    @staticmethod
    def _vf_section(tmp):
        sim = types.SimpleNamespace(
            path=Path(tmp),
            expnr="001",
            stl_file="geom.stl",
            geom=types.SimpleNamespace(
                stl=types.SimpleNamespace(faces=np.zeros((3, 3), dtype=int))
            ),
        )
        section = RadiationSection("radiation", {}, sim=sim, defaults={})
        section.view3d_out = 0
        section.lvfsparse = False
        section.maxD = 1000.0
        return section

    @staticmethod
    def _enter_view3d_mocks(stack):
        fake_vf = mock.MagicMock()
        fake_vf.nnz = 5
        fake_vf.toarray.return_value = np.zeros((3, 3))
        mocks = stack.enter_context(
            mock.patch.multiple(
                "udprep.udprep_radiation",
                stl_to_view3d=mock.DEFAULT,
                resolve_view3d_exe=mock.DEFAULT,
                run_view3d=mock.DEFAULT,
                read_view3d_output=mock.DEFAULT,
                compute_svf=mock.DEFAULT,
                write_svf=mock.DEFAULT,
                write_vf=mock.DEFAULT,
                write_vfsparse=mock.DEFAULT,
            )
        )
        mocks["read_view3d_output"].return_value = fake_vf
        mocks["compute_svf"].return_value = np.zeros(3)
        return mocks

    def test_calc_view_factors_uses_per_call_maxD(self):
        with TemporaryDirectory() as tmp, contextlib.ExitStack() as stack:
            mocks = self._enter_view3d_mocks(stack)
            section = self._vf_section(tmp)
            section.calc_view_factors(maxD=123.0)
            self.assertEqual(mocks["stl_to_view3d"].call_args.kwargs["maxD"], 123.0)

    def test_calc_view_factors_maxD_none_uses_self_maxD(self):
        with TemporaryDirectory() as tmp, contextlib.ExitStack() as stack:
            mocks = self._enter_view3d_mocks(stack)
            section = self._vf_section(tmp)
            section.calc_view_factors(maxD=None)
            self.assertEqual(mocks["stl_to_view3d"].call_args.kwargs["maxD"], 1000.0)

    def test_calc_view_factors_maxD_participates_in_cache_key(self):
        with TemporaryDirectory() as tmp, contextlib.ExitStack() as stack:
            mocks = self._enter_view3d_mocks(stack)
            section = self._vf_section(tmp)
            section.calc_view_factors(maxD=123.0)
            section.calc_view_factors(maxD=123.0)  # same maxD -> memory cache hit
            self.assertEqual(mocks["stl_to_view3d"].call_count, 1)
            section.calc_view_factors(maxD=456.0)  # different maxD -> recompute
            self.assertEqual(mocks["stl_to_view3d"].call_count, 2)


if __name__ == "__main__":
    unittest.main()
