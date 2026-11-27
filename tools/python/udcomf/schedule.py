"""Statistics averaging schedule from the solver's OUTPUT namelist.

Every record in a time-averaged statistics file is the mean over the
``tstatsdump`` seconds that precede its saved timestamp. The solver starts the
accumulators from zero after each dump and after a restart, and writes a
record only once ``tstatsdump`` seconds have accumulated, so a saved record is
always a full-length mean; a restart shifts the dump instants but never
shortens a window. ``tstatstart`` and ``tstatsgap`` place the dumps on the
clock and are used here only to report records that fall off that schedule.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import f90nml
import numpy as np


# Solver defaults from modglobal.f90, used when a key is absent in namoptions.
_DEFAULTS = {"tstatstart": 0.0, "tstatsdump": 10000.0, "tstatsgap": 0.0, "tsample": 5.0}


@dataclass(frozen=True)
class StatsSchedule:
    tstatstart: float
    tstatsdump: float
    tstatsgap: float
    tsample: float

    @property
    def period(self) -> float:
        """Time between consecutive dumps of an uninterrupted run."""
        return self.tstatsdump + self.tstatsgap

    def window_bounds(self, times: np.ndarray) -> np.ndarray:
        """``[t - tstatsdump, t)`` for each saved record time ``t``."""
        times = np.asarray(times, dtype=float)
        return np.column_stack((times - self.tstatsdump, times))

    def window_midpoints(self, times: np.ndarray) -> np.ndarray:
        return np.asarray(times, dtype=float) - 0.5 * self.tstatsdump

    def off_schedule(self, times: np.ndarray, *, tolerance: float | None = None) -> np.ndarray:
        """True where a dump instant is not ``tstatstart + n*period`` for ``n >= 1``.

        Such records are legitimate after a warm start inside a window: the
        solver then dumps ``tstatsdump`` seconds after the restart instead of
        at the scheduled instant. The record is still a full-length mean.
        """
        times = np.asarray(times, dtype=float)
        tolerance = self.tsample if tolerance is None else float(tolerance)
        offset = times - self.tstatstart
        n = np.round(offset / self.period)
        return (n < 1) | (np.abs(offset - n * self.period) > tolerance)


def read_stats_schedule(case_dir: Path) -> StatsSchedule:
    """Read the statistics schedule of a case, with the solver defaults."""
    case = Path(case_dir)
    output = f90nml.read(case / f"namoptions.{case.name}").get("output", {})
    values = {}
    for key, default in _DEFAULTS.items():
        try:
            values[key] = float(output.get(key, default))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"OUTPUT.{key} must be a number") from exc
    if not np.isfinite(list(values.values())).all():
        raise ValueError("Statistics schedule settings must be finite")
    if values["tstatsdump"] <= 0:
        raise ValueError("tstatsdump must be positive")
    if values["tstatsgap"] < 0 or values["tstatstart"] < 0:
        raise ValueError("tstatsgap and tstatstart must be nonnegative")
    if values["tsample"] <= 0:
        raise ValueError("tsample must be positive")
    return StatsSchedule(**values)


def validate_statistics_times(times: np.ndarray, schedule: StatsSchedule, *,
                              tolerance: float | None = None) -> np.ndarray:
    """Return saved dump instants after checking that their windows do not overlap."""
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or times.size == 0 or not np.isfinite(times).all():
        raise ValueError("Statistics times must be a nonempty finite vector")
    if np.any(np.diff(times) <= 0):
        raise ValueError("Statistics times must increase strictly")
    tolerance = schedule.tsample if tolerance is None else float(tolerance)
    if np.any(np.diff(times) < schedule.tstatsdump - tolerance):
        raise ValueError(
            "Statistics records are closer than tstatsdump; their averaging windows overlap"
        )
    return times
