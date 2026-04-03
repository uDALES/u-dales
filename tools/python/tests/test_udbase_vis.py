import sys
import unittest
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udbase import UDBase  # noqa: E402
from udvis import UDVis  # noqa: E402


class RecordingVis:
    def __init__(self):
        self.calls = []

    def plot_fac(self, **kwargs):
        self.calls.append(("plot_fac", kwargs))
        return "plot_fac_result"

    def plot_2dmap(self, **kwargs):
        self.calls.append(("plot_2dmap", kwargs))
        return "plot_2dmap_result"


class AliasTrackingVis(UDVis):
    def __init__(self):
        super().__init__(sim=object())
        self.called = None

    def plot_veg(self, veg=None, show=False):
        self.called = {"veg": veg, "show": show}
        return "plot_veg_result"


class TestUDBaseVisualizationCompatibility(unittest.TestCase):
    def test_udbase_plot_fac_forwards_to_vis_facade(self):
        sim = UDBase.__new__(UDBase)
        sim.vis = RecordingVis()

        result = sim.plot_fac(var="dummy_var", building_ids=[1, 2], show=False)

        self.assertEqual(result, "plot_fac_result")
        self.assertEqual(
            sim.vis.calls,
            [("plot_fac", {"var": "dummy_var", "building_ids": [1, 2], "show": False})],
        )

    def test_udbase_plot_2dmap_forwards_to_vis_facade(self):
        sim = UDBase.__new__(UDBase)
        sim.vis = RecordingVis()

        result = sim.plot_2dmap(val=1.0, labels=["A"], show=False)

        self.assertEqual(result, "plot_2dmap_result")
        self.assertEqual(
            sim.vis.calls,
            [("plot_2dmap", {"val": 1.0, "labels": ["A"], "show": False})],
        )

    def test_plot_trees_remains_alias_for_plot_veg(self):
        vis = AliasTrackingVis()

        result = vis.plot_trees(show=True)

        self.assertEqual(result, "plot_veg_result")
        self.assertEqual(vis.called, {"veg": None, "show": True})


if __name__ == "__main__":
    unittest.main()
