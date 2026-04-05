"""
Visualization facade for Python uDALES tooling.

`UDBase` remains the canonical stateful simulation object. `UDVis` provides a
stable plotting entry point on top of that state via `sim.vis`, while keeping
backend choice as an internal implementation detail that can change over time.
"""

from .udbase_vis import UDVis

__all__ = ["UDVis"]
