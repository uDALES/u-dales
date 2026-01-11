from __future__ import annotations

import time
from datetime import datetime, timedelta
from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt

# Add the uDALES python path
udbase_path = Path("C:/Users/mvr/OneDrive - Imperial College London/codes/uDALES/u-dales").resolve()
tools_path = (udbase_path / "tools" / "python").resolve()
if tools_path not in sys.path:
    sys.path.insert(0, str(tools_path))

from udprep import solar  # noqa: E402


def _angle_diff_deg(a: float, b: float) -> float:
    diff = (a - b + 180.0) % 360.0 - 180.0
    return abs(diff)


start = time.perf_counter()

longitude = -0.13
latitude = 51.5
timezone = 0.0
elevation = 0.0
xazimuth = 90.0

start_time = datetime(2011, 9, 30, 6, 0, 0)
steps = 25
dt = timedelta(minutes=30)
times = [start_time + n * dt for n in range(steps)]

zen_py = []
zen_pv = []
az_py = []
az_pv = []
I_py = []
I_pv = []
Dsky_py = []
Dsky_pv = []

for t in times:
    _, zen, az, I, Dsky = solar.solar_state(
        t,
        longitude,
        latitude,
        timezone,
        elevation,
        xazimuth=xazimuth,
        backend="python",
    )
    _, zen2, az2, I2, Dsky2 = solar.solar_state(
        t,
        longitude,
        latitude,
        timezone,
        elevation,
        xazimuth=xazimuth,
        backend="pvlib",
    )
    zen_py.append(float(zen))
    zen_pv.append(float(zen2))
    az_py.append(float(az))
    az_pv.append(float(az2))
    I_py.append(float(I))
    I_pv.append(float(I2))
    Dsky_py.append(float(Dsky))
    Dsky_pv.append(float(Dsky2))

zen_py = np.asarray(zen_py)
zen_pv = np.asarray(zen_pv)
az_py = np.asarray(az_py)
az_pv = np.asarray(az_pv)
I_py = np.asarray(I_py)
I_pv = np.asarray(I_pv)
Dsky_py = np.asarray(Dsky_py)
Dsky_pv = np.asarray(Dsky_pv)

zen_diff = np.array([_angle_diff_deg(a, b) for a, b in zip(zen_py, zen_pv)])
az_diff = np.array([_angle_diff_deg(a, b) for a, b in zip(az_py, az_pv)])

elapsed = time.perf_counter() - start
print(f"solar_state comparison runtime: {elapsed:.3f} s")
print(f"zenith diff: max {zen_diff.max():.3f} deg, mean {zen_diff.mean():.3f} deg")
print(f"azimuth diff: max {az_diff.max():.3f} deg, mean {az_diff.mean():.3f} deg")
print(f"I diff: max {np.max(np.abs(I_py - I_pv)):.2f} W/m2")
print(f"Dsky diff: max {np.max(np.abs(Dsky_py - Dsky_pv)):.2f} W/m2")

# Plot angles
fig_angles, ax = plt.subplots()
zen_color = "tab:blue"
az_color = "tab:orange"
ax.plot(times, zen_py, "o", color=zen_color, label="zenith (python)")
ax.plot(times, zen_pv, color=zen_color, label="zenith (pvlib)")
ax.plot(times, az_py, "o", color=az_color, label="azimuth (python)")
ax.plot(times, az_pv, color=az_color, label="azimuth (pvlib)")
ax.set_title("Solar angles comparison")
ax.set_xlabel("time")
ax.set_ylabel("angle [deg]")
ax.legend()

# Plot irradiance
fig_irr, ax = plt.subplots()
I_color = "tab:green"
Dsky_color = "tab:red"
ax.plot(times, I_py, "o", color=I_color, label="I (python)")
ax.plot(times, I_pv, color=I_color, label="I (pvlib)")
ax.plot(times, Dsky_py, "o", color=Dsky_color, label="Dsky (python)")
ax.plot(times, Dsky_pv, color=Dsky_color, label="Dsky (pvlib)")
ax.set_title("Irradiance comparison")
ax.set_xlabel("time")
ax.set_ylabel("irradiance [W/m2]")
ax.legend()

plt.show()
