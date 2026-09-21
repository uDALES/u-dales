"""Facade for thermal-comfort postprocessing of a UDBase case."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from udbase import UDBase


class UDComf:
    """Hold the existing case state for future comfort calculations."""

    def __init__(self, sim: UDBase):
        self.sim = sim
