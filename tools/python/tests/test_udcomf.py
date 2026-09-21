"""Step 5 package and UDBase-facade tests."""

import importlib
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from udbase import UDBase
from udcomf import UDComf


class TestUDComfPackage(unittest.TestCase):
    def test_modules_import(self):
        for name in (
            "udcomf.udcomf_radiation",
            "udcomf.udcomf_export",
            "udcomf.udcomf_io",
            "udcomf.checkpoints",
            "udcomf.thermalcomfort",
        ):
            with self.subTest(module=name):
                self.assertIsNotNone(importlib.import_module(name).__doc__)

    def test_udbase_attaches_facade_to_same_case(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "namoptions.001").write_text(
                "&DOMAIN\n itot = 2\n jtot = 2\n ktot = 2\n"
                " xlen = 20.0\n ylen = 20.0\n zsize = 20.0\n/\n",
                encoding="ascii",
            )
            sim = UDBase(1, path, load_geometry=False, suppress_load_warnings=True)

            self.assertIsInstance(sim.comf, UDComf)
            self.assertIs(sim.comf.sim, sim)
            self.assertIsNot(sim.comf, sim.vis)
