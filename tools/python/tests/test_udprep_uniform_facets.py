"""Tests for the udprep hook that uniformises STL facet sizes (lfacetuniform).

Uses a small case directory with a box STL whose walls are far larger than the
target edge length, so the hook has to write a new STL, a new facet-type file,
the face map and the sidecar, point stl_file/types_path at them (in memory and
in namoptions) and reload the geometry.
"""
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import trimesh

from _common import PYTHON_DIR  # noqa: F401

from udbase import UDBase
from udgeom import FacetCountError
from udprep.udprep import UDPrep


def _box_stl(path: Path) -> int:
    v = np.array(
        [
            [0.0, 0.0, 0.0], [30.0, 0.0, 0.0], [30.0, 30.0, 0.0], [0.0, 30.0, 0.0],
            [0.0, 0.0, 30.0], [30.0, 0.0, 30.0], [30.0, 30.0, 30.0], [0.0, 30.0, 30.0],
        ]
    )

    def quad(a, b, c, d):
        return [[a, b, c], [a, c, d]]

    f = quad(0, 1, 5, 4) + quad(1, 2, 6, 5) + quad(2, 3, 7, 6) + quad(3, 0, 4, 7) + quad(4, 5, 6, 7)
    mesh = trimesh.Trimesh(vertices=v, faces=np.asarray(f), process=False)
    mesh.export(str(path))
    return len(mesh.faces)


def _namoptions(**inps) -> str:
    lines = [
        "&RUN", " lwarmstart = .false.", "/",
        "&DOMAIN", " itot = 8", " jtot = 8", " ktot = 8", " xlen = 80.0", " ylen = 80.0", "/",
        "&WALLS", " nfcts = 0", "/",
        "&INPS", " zsize = 80.0", " stl_file = 'box.stl'", " read_types = .true.", " types_path = 'box_types.txt'",
    ]
    lines += [f" {k} = {v}" for k, v in inps.items()]
    lines += ["/", ""]
    return "\n".join(lines)


class TestUniformizeFacetsHook(unittest.TestCase):
    def setUp(self):
        self._tmp = TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.case = Path(self._tmp.name)
        self.n_in = _box_stl(self.case / "box.stl")
        (self.case / "box_types.txt").write_text(
            "# types\n" + "\n".join(["1"] * 8 + ["2"] * 2) + "\n", encoding="ascii"
        )

    def _prep(self, **inps):
        (self.case / "namoptions.001").write_text(_namoptions(**inps), encoding="ascii")
        sim = UDBase("001", self.case, load_geometry=True, suppress_load_warnings=True)
        return UDPrep(sim)

    def test_defaults_are_registered(self):
        prep = self._prep()
        self.assertFalse(prep.ibm.lfacetuniform)
        self.assertEqual(prep.ibm.facet_target_edge, 0.0)
        self.assertEqual(prep.ibm.facet_size_band, 4.0)
        self.assertTrue(prep.ibm.facet_allow_split)
        self.assertTrue(prep.ibm.facet_allow_merge)
        self.assertEqual(prep.ibm.facet_count_max, 0)

    def test_hook_writes_outputs_and_switches_the_case_to_them(self):
        prep = self._prep(lfacetuniform=".true.", facet_target_edge="10.0", facet_count_max="1000")
        prep.ibm.uniformize_facets()

        out_stl = self.case / "box_uniform.stl"
        out_types = self.case / "box_types_uniform.txt"
        sidecar = self.case / "box_uniform.json"
        facemap = self.case / "box_uniform_facemap.txt"
        for path in (out_stl, out_types, sidecar, facemap):
            self.assertTrue(path.is_file(), path)

        # Original inputs untouched.
        self.assertEqual(len(trimesh.load_mesh(str(self.case / "box.stl")).faces), self.n_in)

        new_mesh = trimesh.load_mesh(str(out_stl))
        self.assertGreater(len(new_mesh.faces), self.n_in)
        types = np.loadtxt(out_types, dtype=int)
        self.assertEqual(len(types), len(new_mesh.faces))
        roof = np.isclose(new_mesh.triangles_center[:, 2], 30.0)
        self.assertTrue(np.all(types[roof] == 2))
        self.assertTrue(np.all(types[~roof] == 1))
        fmap = np.loadtxt(facemap, dtype=int)
        self.assertEqual(len(fmap), len(new_mesh.faces))
        self.assertEqual(fmap.min(), 1)
        self.assertEqual(fmap.max(), self.n_in)

        meta = json.loads(sidecar.read_text())
        self.assertEqual(meta["source"]["stl"], "box.stl")
        self.assertEqual(meta["settings"]["facet_target_edge"], 10.0)
        self.assertEqual(meta["report"]["n_faces_out"], len(new_mesh.faces))

        # The case now works on the uniform geometry, in memory and in namoptions.
        self.assertEqual(prep.sim.stl_file, "box_uniform.stl")
        self.assertEqual(prep.sim.types_path, "box_types_uniform.txt")
        self.assertEqual(prep.sim.geom.n_faces, len(new_mesh.faces))
        text = (self.case / "namoptions.001").read_text()
        self.assertIn("stl_file = 'box_uniform.stl'", text)
        self.assertIn("types_path = 'box_types_uniform.txt'", text)
        self.assertNotIn("&INP\n", text)  # updated in place inside &INPS, no stray block
        self.assertEqual(text.count("stl_file"), 1)

        # Rerun from the rewritten namoptions: stl_file is already a product -> no-op.
        sim2 = UDBase("001", self.case, load_geometry=True, suppress_load_warnings=True)
        prep2 = UDPrep(sim2)
        prep2.ibm.uniformize_facets()
        self.assertEqual(prep2.sim.stl_file, "box_uniform.stl")
        self.assertFalse((self.case / "box_uniform_uniform.stl").exists())

        # Same product with other settings is refused with guidance.
        (self.case / "namoptions.001").write_text(
            (self.case / "namoptions.001").read_text().replace("facet_target_edge = 10.0", "facet_target_edge = 5.0")
        )
        prep3 = UDPrep(UDBase("001", self.case, load_geometry=True, suppress_load_warnings=True))
        with self.assertRaises(ValueError) as ctx:
            prep3.ibm.uniformize_facets()
        self.assertIn("Set stl_file back", str(ctx.exception))

    def test_existing_output_is_reused_only_when_sidecar_matches(self):
        prep = self._prep(lfacetuniform=".true.", facet_target_edge="10.0", facet_count_max="1000")
        prep.ibm.uniformize_facets()
        out_stl = self.case / "box_uniform.stl"
        mtime = out_stl.stat().st_mtime_ns

        # Fresh case pointing at the source again, same settings -> reuse, no rewrite.
        (self.case / "namoptions.001").write_text(
            _namoptions(lfacetuniform=".true.", facet_target_edge="10.0", facet_count_max="1000"), encoding="ascii"
        )
        prep2 = UDPrep(UDBase("001", self.case, load_geometry=True, suppress_load_warnings=True))
        prep2.ibm.uniformize_facets()
        self.assertEqual(out_stl.stat().st_mtime_ns, mtime)
        self.assertEqual(prep2.sim.stl_file, "box_uniform.stl")

        # Different settings -> the existing file is not overwritten.
        (self.case / "namoptions.001").write_text(
            _namoptions(lfacetuniform=".true.", facet_target_edge="5.0", facet_count_max="1000"), encoding="ascii"
        )
        prep3 = UDPrep(UDBase("001", self.case, load_geometry=True, suppress_load_warnings=True))
        with self.assertRaises(FileExistsError):
            prep3.ibm.uniformize_facets()
        self.assertEqual(out_stl.stat().st_mtime_ns, mtime)

    def test_count_limit_is_enforced(self):
        prep = self._prep(lfacetuniform=".true.", facet_target_edge="2.0", facet_count_max="50")
        with self.assertRaises(FacetCountError) as ctx:
            prep.ibm.uniformize_facets()
        self.assertIn("facet_count_max", str(ctx.exception))
        self.assertFalse((self.case / "box_uniform.stl").exists())

    def test_run_all_schedules_the_step_only_when_switched_on(self):
        prep = self._prep()
        prep.ibm.run_steps = lambda label, steps, **kw: setattr(prep.ibm, "_seen", [name for name, _ in steps])
        prep.ibm.run_all()
        self.assertNotIn("uniformize_facets", prep.ibm._seen)
        prep = self._prep(lfacetuniform=".true.")
        prep.ibm.run_steps = lambda label, steps, **kw: setattr(prep.ibm, "_seen", [name for name, _ in steps])
        prep.ibm.run_all()
        self.assertEqual(prep.ibm._seen[0], "uniformize_facets")
        self.assertEqual(prep.ibm._seen[1], "run_ibm")

    def test_stale_facets_inp_is_refused(self):
        prep = self._prep(lfacetuniform=".true.", facet_target_edge="10.0", facet_count_max="1000")
        # facets.inp written for the 10-face source geometry.
        with (self.case / "facets.inp.001").open("w") as f:
            f.write("# type, normal\n")
            for _ in range(self.n_in):
                f.write("1 0.0 0.0 1.0\n")
        prep.ibm.uniformize_facets()
        with self.assertRaises(ValueError) as ctx:
            prep.ibm.write_facets()
        self.assertIn("different STL", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
