from __future__ import annotations

"""
Solar position + irradiance helpers (ported from tools/SEB/SPA + ASHRAE.m).

Backends:
- python: direct port of MATLAB SPA implementation.
- pvlib: use pvlib for angles, ASHRAE for irradiance.
"""

from datetime import datetime, timezone as dt_timezone, timedelta
import math
from typing import Any, Dict, Tuple

import numpy as np


# -----------------------------------------------------------------------------
# ASHRAE coefficients (tools/SEB/ASHRAE.m)
# -----------------------------------------------------------------------------

_ASHRAE = {
    1: (1230.0, 0.142, 0.058),
    2: (1215.0, 0.144, 0.060),
    3: (1186.0, 0.156, 0.071),
    4: (1136.0, 0.180, 0.097),
    5: (1104.0, 0.196, 0.121),
    6: (1088.0, 0.205, 0.134),
    7: (1085.0, 0.207, 0.136),
    8: (1107.0, 0.201, 0.122),
    9: (1151.0, 0.177, 0.092),
    10: (1192.0, 0.160, 0.073),
    11: (1221.0, 0.149, 0.063),
    12: (1234.0, 0.142, 0.057),
}


def solar_strength_ashrae(month: int, zenith_deg: float) -> Tuple[float, float]:
    """
    ASHRAE clear-sky model.

    Returns
    -------
    I : float
        Direct normal irradiance [W/m^2].
    Dsky : float
        Diffuse sky irradiance [W/m^2].
    """
    if month not in _ASHRAE:
        raise ValueError("month must be in 1..12")
    a, b, c = _ASHRAE[month]
    cos_zen = math.cos(math.radians(zenith_deg))
    if cos_zen <= 0.0:
        return 0.0, 0.0
    direct = a * math.exp(-b / cos_zen)
    return direct, c * direct


# -----------------------------------------------------------------------------
# SPA constants (tools/SEB/SPA/spa_const.m)
# -----------------------------------------------------------------------------

SPA_ZA = 0
SPA_ZA_INC = 1
SPA_ZA_RTS = 2
SPA_ALL = 3

SUN_RADIUS = 0.26667
L_COUNT = 6
B_COUNT = 2
R_COUNT = 5
Y_COUNT = 63
L_MAX_SUBCOUNT = 64
B_MAX_SUBCOUNT = 5
R_MAX_SUBCOUNT = 40
TERM_A = 0
TERM_B = 1
TERM_C = 2
TERM_X0 = 0
TERM_X1 = 1
TERM_X2 = 2
TERM_X3 = 3
TERM_X4 = 4
TERM_COUNT = 3
TERM_PSI_A = 0
TERM_PSI_B = 1
TERM_EPS_C = 2
TERM_EPS_D = 3
TERM_PE_COUNT = 4
JD_MINUS = 0
JD_ZERO = 1
JD_PLUS = 2
JD_COUNT = 3
SUN_TRANSIT = 0
SUN_RISE = 1
SUN_SET = 2
SUN_COUNT = 3
TERM_X_COUNT = 5
TERM_Y_COUNT = 5

l_subcount = [64, 34, 20, 7, 3, 1]
b_subcount = [5, 2]
r_subcount = [40, 10, 6, 2, 1]

L_TERMS = [
    np.array(
        [
            [175347046.0, 0.0, 0.0],
            [3341656.0, 4.6692568, 6283.07585],
            [34894.0, 4.6261, 12566.1517],
            [3497.0, 2.7441, 5753.3849],
            [3418.0, 2.8289, 3.5231],
            [3136.0, 3.6277, 77713.7715],
            [2676.0, 4.4181, 7860.4194],
            [2343.0, 6.1352, 3930.2097],
            [1324.0, 0.7425, 11506.7698],
            [1273.0, 2.0371, 529.691],
            [1199.0, 1.1096, 1577.3435],
            [990.0, 5.233, 5884.927],
            [902.0, 2.045, 26.298],
            [857.0, 3.508, 398.149],
            [780.0, 1.179, 5223.694],
            [753.0, 2.533, 5507.553],
            [505.0, 4.583, 18849.228],
            [492.0, 4.205, 775.523],
            [357.0, 2.92, 0.067],
            [317.0, 5.849, 11790.629],
            [284.0, 1.899, 796.298],
            [271.0, 0.315, 10977.079],
            [243.0, 0.345, 5486.778],
            [206.0, 4.806, 2544.314],
            [205.0, 1.869, 5573.143],
            [202.0, 2.458, 6069.777],
            [156.0, 0.833, 213.299],
            [132.0, 3.411, 2942.463],
            [126.0, 1.083, 20.775],
            [115.0, 0.645, 0.98],
            [103.0, 0.636, 4694.003],
            [102.0, 0.976, 15720.839],
            [102.0, 4.267, 7.114],
            [99.0, 6.21, 2146.17],
            [98.0, 0.68, 155.42],
            [86.0, 5.98, 161000.69],
            [85.0, 1.3, 6275.96],
            [85.0, 3.67, 71430.7],
            [80.0, 1.81, 17260.15],
            [79.0, 3.04, 12036.46],
            [75.0, 1.76, 5088.63],
            [74.0, 3.5, 3154.69],
            [74.0, 4.68, 801.82],
            [70.0, 0.83, 9437.76],
            [62.0, 3.98, 8827.39],
            [61.0, 1.82, 7084.9],
            [57.0, 2.78, 6286.6],
            [56.0, 4.39, 14143.5],
            [56.0, 3.47, 6279.55],
            [52.0, 0.19, 12139.55],
            [52.0, 1.33, 1748.02],
            [51.0, 0.28, 5856.48],
            [49.0, 0.49, 1194.45],
            [41.0, 5.37, 8429.24],
            [41.0, 2.4, 19651.05],
            [39.0, 6.17, 10447.39],
            [37.0, 6.04, 10213.29],
            [37.0, 2.57, 1059.38],
            [36.0, 1.71, 2352.87],
            [36.0, 1.78, 6812.77],
            [33.0, 0.59, 17789.85],
            [30.0, 0.44, 83996.85],
            [30.0, 2.74, 1349.87],
            [25.0, 3.16, 4690.48],
        ],
        dtype=float,
    ),
    np.array(
        [
            [628331966747.0, 0.0, 0.0],
            [206059.0, 2.678235, 6283.07585],
            [4303.0, 2.6351, 12566.1517],
            [425.0, 1.59, 3.523],
            [119.0, 5.796, 26.298],
            [109.0, 2.966, 1577.344],
            [93.0, 2.59, 18849.23],
            [72.0, 1.14, 529.69],
            [68.0, 1.87, 398.15],
            [67.0, 4.41, 5507.55],
            [59.0, 2.89, 5223.69],
            [56.0, 2.17, 155.42],
            [45.0, 0.4, 796.3],
            [36.0, 0.47, 775.52],
            [29.0, 2.65, 7.11],
            [21.0, 5.34, 0.98],
            [19.0, 1.85, 5486.78],
            [19.0, 4.97, 213.3],
            [17.0, 2.99, 6275.96],
            [16.0, 0.03, 2544.31],
            [16.0, 1.43, 2146.17],
            [15.0, 1.21, 10977.08],
            [12.0, 2.83, 1748.02],
            [12.0, 3.26, 5088.63],
            [12.0, 5.27, 1194.45],
            [12.0, 2.08, 4694.0],
            [11.0, 0.77, 553.57],
            [10.0, 1.3, 6286.6],
            [10.0, 4.24, 1349.87],
            [9.0, 2.7, 242.73],
            [9.0, 5.64, 951.72],
            [8.0, 5.3, 2352.87],
            [6.0, 2.65, 9437.76],
            [6.0, 4.67, 4690.48],
        ],
        dtype=float,
    ),
    np.array(
        [
            [52919.0, 0.0, 0.0],
            [8720.0, 1.0721, 6283.0758],
            [309.0, 0.867, 12566.152],
            [27.0, 0.05, 3.52],
            [16.0, 5.19, 26.3],
            [16.0, 3.68, 155.42],
            [10.0, 0.76, 18849.23],
            [9.0, 2.06, 77713.77],
            [7.0, 0.83, 775.52],
            [5.0, 4.66, 1577.34],
            [4.0, 1.03, 7.11],
            [4.0, 3.44, 5573.14],
            [3.0, 5.14, 796.3],
            [3.0, 6.05, 5507.55],
            [3.0, 1.19, 242.73],
            [3.0, 6.12, 529.69],
            [3.0, 0.31, 398.15],
            [3.0, 2.28, 553.57],
            [2.0, 4.38, 5223.69],
            [2.0, 3.75, 0.98],
        ],
        dtype=float,
    ),
    np.array(
        [
            [289.0, 5.844, 6283.076],
            [35.0, 0.0, 0.0],
            [17.0, 5.49, 12566.15],
            [3.0, 5.2, 155.42],
            [1.0, 4.72, 3.52],
            [1.0, 5.3, 18849.23],
            [1.0, 5.97, 242.73],
        ],
        dtype=float,
    ),
    np.array(
        [
            [114.0, 3.142, 0.0],
            [8.0, 4.13, 6283.08],
            [1.0, 3.84, 12566.15],
        ],
        dtype=float,
    ),
    np.array([[1.0, 3.14, 0.0]], dtype=float),
]

B_TERMS = [
    np.array(
        [
            [280.0, 3.199, 84334.662],
            [102.0, 5.422, 5507.553],
            [80.0, 3.88, 5223.69],
            [44.0, 3.7, 2352.87],
            [32.0, 4.0, 1577.34],
        ],
        dtype=float,
    ),
    np.array([[9.0, 3.9, 5507.55], [6.0, 1.73, 5223.69]], dtype=float),
]

R_TERMS = [
    np.array(
        [
            [100013989.0, 0.0, 0.0],
            [1670700.0, 3.0984635, 6283.07585],
            [13956.0, 3.05525, 12566.1517],
            [3084.0, 5.1985, 77713.7715],
            [1628.0, 1.1739, 5753.3849],
            [1576.0, 2.8469, 7860.4194],
            [925.0, 5.453, 11506.77],
            [542.0, 4.564, 3930.21],
            [472.0, 3.661, 5884.927],
            [346.0, 0.964, 5507.553],
            [329.0, 5.9, 5223.694],
            [307.0, 0.299, 5573.143],
            [243.0, 4.273, 11790.629],
            [212.0, 5.847, 1577.344],
            [186.0, 5.022, 10977.079],
            [175.0, 3.012, 18849.228],
            [110.0, 5.055, 5486.778],
            [98.0, 0.89, 6069.78],
            [86.0, 5.69, 15720.84],
            [86.0, 1.27, 161000.69],
            [65.0, 0.27, 17260.15],
            [63.0, 0.92, 529.69],
            [57.0, 2.01, 83996.85],
            [56.0, 5.24, 71430.7],
            [49.0, 3.25, 2544.31],
            [47.0, 2.58, 775.52],
            [45.0, 5.54, 9437.76],
            [43.0, 6.01, 6275.96],
            [39.0, 5.36, 4694.0],
            [38.0, 2.39, 8827.39],
            [37.0, 0.83, 19651.05],
            [37.0, 4.9, 12139.55],
            [36.0, 1.67, 12036.46],
            [35.0, 1.84, 2942.46],
            [33.0, 0.24, 7084.9],
            [32.0, 0.18, 5088.63],
            [32.0, 1.78, 398.15],
            [28.0, 1.21, 6286.6],
            [28.0, 1.9, 6279.55],
            [26.0, 4.59, 10447.39],
        ],
        dtype=float,
    ),
    np.array(
        [
            [103019.0, 1.10749, 6283.07585],
            [1721.0, 1.0644, 12566.1517],
            [702.0, 3.142, 0.0],
            [32.0, 1.02, 18849.23],
            [31.0, 2.84, 5507.55],
            [25.0, 1.32, 5223.69],
            [18.0, 1.42, 1577.34],
            [10.0, 5.91, 10977.08],
            [9.0, 1.42, 6275.96],
            [9.0, 0.27, 5486.78],
        ],
        dtype=float,
    ),
    np.array(
        [
            [4359.0, 5.7846, 6283.0758],
            [124.0, 5.579, 12566.152],
            [12.0, 3.14, 0.0],
            [9.0, 3.63, 77713.77],
            [6.0, 1.87, 5573.14],
            [3.0, 5.47, 18849.23],
        ],
        dtype=float,
    ),
    np.array([[145.0, 4.273, 6283.076], [7.0, 3.92, 12566.15]], dtype=float),
    np.array([[4.0, 2.56, 6283.08]], dtype=float),
]

Y_TERMS = np.array(
    [
        [0, 0, 0, 0, 1],
        [-2, 0, 0, 2, 2],
        [0, 0, 0, 2, 2],
        [0, 0, 0, 0, 2],
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [-2, 1, 0, 2, 2],
        [0, 0, 0, 2, 1],
        [0, 0, 1, 2, 2],
        [-2, -1, 0, 2, 2],
        [-2, 0, 1, 0, 0],
        [-2, 0, 0, 2, 1],
        [0, 0, -1, 2, 2],
        [2, 0, 0, 0, 0],
        [0, 0, 1, 0, 1],
        [2, 0, -1, 2, 2],
        [0, 0, -1, 0, 1],
        [0, 0, 1, 2, 1],
        [-2, 0, 2, 0, 0],
        [0, 0, -2, 2, 1],
        [2, 0, 0, 2, 2],
        [0, 0, 2, 2, 2],
        [0, 0, 2, 0, 0],
        [-2, 0, 1, 2, 2],
        [0, 0, 0, 2, 0],
        [-2, 0, 0, 2, 0],
        [0, 0, -1, 2, 1],
        [0, 2, 0, 0, 0],
        [2, 0, -1, 0, 1],
        [-2, 2, 0, 2, 2],
        [0, 1, 0, 0, 1],
        [-2, 0, 1, 0, 1],
        [0, -1, 0, 0, 1],
        [0, 0, 2, -2, 0],
        [2, 0, -1, 2, 1],
        [2, 0, 1, 2, 2],
        [0, 1, 0, 2, 2],
        [-2, 1, 1, 0, 0],
        [0, -1, 0, 2, 2],
        [2, 0, 0, 2, 1],
        [2, 0, 1, 0, 0],
        [-2, 0, 2, 2, 2],
        [-2, 0, 1, 2, 1],
        [2, 0, -2, 0, 1],
        [2, 0, 0, 0, 1],
        [0, -1, 1, 0, 0],
        [-2, -1, 0, 2, 1],
        [-2, 0, 0, 0, 1],
        [0, 0, 2, 2, 1],
        [-2, 0, 2, 0, 1],
        [-2, 1, 0, 2, 1],
        [0, 0, 1, -2, 0],
        [-1, 0, 1, 0, 0],
        [-2, 1, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 1, 2, 0],
        [0, 0, -2, 2, 2],
        [-1, -1, 1, 0, 0],
        [0, 1, 1, 0, 0],
        [0, -1, 1, 2, 2],
        [2, -1, -1, 2, 2],
        [0, 0, 3, 2, 2],
        [2, -1, 0, 2, 2],
    ],
    dtype=float,
)

PE_TERMS = np.array(
    [
        [-171996.0, -174.2, 92025.0, 8.9],
        [-13187.0, -1.6, 5736.0, -3.1],
        [-2274.0, -0.2, 977.0, -0.5],
        [2062.0, 0.2, -895.0, 0.5],
        [1426.0, -3.4, 54.0, -0.1],
        [712.0, 0.1, -7.0, 0.0],
        [-517.0, 1.2, 224.0, -0.6],
        [-386.0, -0.4, 200.0, 0.0],
        [-301.0, 0.0, 129.0, -0.1],
        [217.0, -0.5, -95.0, 0.3],
        [-158.0, 0.0, 0.0, 0.0],
        [129.0, 0.1, -70.0, 0.0],
        [123.0, 0.0, -53.0, 0.0],
        [63.0, 0.0, 0.0, 0.0],
        [63.0, 0.1, -33.0, 0.0],
        [-59.0, 0.0, 26.0, 0.0],
        [-58.0, -0.1, 32.0, 0.0],
        [-51.0, 0.0, 27.0, 0.0],
        [48.0, 0.0, 0.0, 0.0],
        [46.0, 0.0, -24.0, 0.0],
        [-38.0, 0.0, 16.0, 0.0],
        [-31.0, 0.0, 13.0, 0.0],
        [29.0, 0.0, 0.0, 0.0],
        [29.0, 0.0, -12.0, 0.0],
        [26.0, 0.0, 0.0, 0.0],
        [-22.0, 0.0, 0.0, 0.0],
        [21.0, 0.0, -10.0, 0.0],
        [17.0, -0.1, 0.0, 0.0],
        [16.0, 0.0, -8.0, 0.0],
        [-16.0, 0.1, 7.0, 0.0],
        [-15.0, 0.0, 9.0, 0.0],
        [-13.0, 0.0, 7.0, 0.0],
        [-12.0, 0.0, 6.0, 0.0],
        [11.0, 0.0, 0.0, 0.0],
        [-10.0, 0.0, 5.0, 0.0],
        [-8.0, 0.0, 3.0, 0.0],
        [7.0, 0.0, -3.0, 0.0],
        [-7.0, 0.0, 0.0, 0.0],
        [-7.0, 0.0, 3.0, 0.0],
        [-7.0, 0.0, 3.0, 0.0],
        [6.0, 0.0, 0.0, 0.0],
        [6.0, 0.0, -3.0, 0.0],
        [6.0, 0.0, -3.0, 0.0],
        [-6.0, 0.0, 3.0, 0.0],
        [-6.0, 0.0, 3.0, 0.0],
        [5.0, 0.0, 0.0, 0.0],
        [-5.0, 0.0, 3.0, 0.0],
        [-5.0, 0.0, 3.0, 0.0],
        [-5.0, 0.0, 3.0, 0.0],
        [4.0, 0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0, 0.0],
        [-4.0, 0.0, 0.0, 0.0],
        [-4.0, 0.0, 0.0, 0.0],
        [-4.0, 0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0, 0.0],
        [-3.0, 0.0, 0.0, 0.0],
        [-3.0, 0.0, 0.0, 0.0],
        [-3.0, 0.0, 0.0, 0.0],
        [-3.0, 0.0, 0.0, 0.0],
        [-3.0, 0.0, 0.0, 0.0],
        [-3.0, 0.0, 0.0, 0.0],
        [-3.0, 0.0, 0.0, 0.0],
    ],
    dtype=float,
)


# -----------------------------------------------------------------------------
# SPA helpers (ported from solarPosition.m)
# -----------------------------------------------------------------------------


def _limit_degrees(degrees: float) -> float:
    degrees = degrees / 360.0
    limited = 360.0 * (degrees - math.floor(degrees))
    if limited < 0.0:
        limited += 360.0
    return limited


def _limit_degrees180pm(degrees: float) -> float:
    degrees = degrees / 360.0
    limited = 360.0 * (degrees - math.floor(degrees))
    if limited < -180.0:
        limited += 360.0
    elif limited > 180.0:
        limited -= 360.0
    return limited


def _limit_degrees180(degrees: float) -> float:
    degrees = degrees / 180.0
    limited = 180.0 * (degrees - math.floor(degrees))
    if limited < 0.0:
        limited += 180.0
    return limited


def _limit_zero2one(value: float) -> float:
    limited = value - math.floor(value)
    if limited < 0.0:
        limited += 1.0
    return limited


def _limit_minutes(minutes: float) -> float:
    limited = minutes
    if limited < -20.0:
        limited += 1440.0
    elif limited > 20.0:
        limited -= 1440.0
    return limited


def _dayfrac_to_local_hr(dayfrac: float, timezone: float) -> float:
    return 24.0 * _limit_zero2one(dayfrac + timezone / 24.0)


def _third_order_polynomial(a: float, b: float, c: float, d: float, x: float) -> float:
    return ((a * x + b) * x + c) * x + d


def _validate_inputs(spa: Dict[str, Any]) -> int:
    if spa["year"] < -2000 or spa["year"] > 6000:
        return 1
    if spa["month"] < 1 or spa["month"] > 12:
        return 2
    if spa["day"] < 1 or spa["day"] > 31:
        return 3
    if spa["hour"] < 0 or spa["hour"] > 24:
        return 4
    if spa["minute"] < 0 or spa["minute"] > 59:
        return 5
    if spa["second"] < 0 or spa["second"] >= 60:
        return 6
    if spa["pressure"] < 0 or spa["pressure"] > 5000:
        return 12
    if spa["temperature"] <= -273 or spa["temperature"] > 6000:
        return 13
    if spa["delta_ut1"] <= -1 or spa["delta_ut1"] >= 1:
        return 17
    if spa["hour"] == 24 and spa["minute"] > 0:
        return 5
    if spa["hour"] == 24 and spa["second"] > 0:
        return 6
    if abs(spa["delta_t"]) > 8000:
        return 7
    if abs(spa["timezone"]) > 18:
        return 8
    if abs(spa["longitude"]) > 180:
        return 9
    if abs(spa["latitude"]) > 90:
        return 10
    if abs(spa["atmos_refract"]) > 5:
        return 16
    if spa["elevation"] < -6500000:
        return 11
    if spa["function"] in (SPA_ZA_INC, SPA_ALL):
        if abs(spa["slope"]) > 360:
            return 14
        if abs(spa["azm_rotation"]) > 360:
            return 15
    return 0


def _julian_day(
    year: int,
    month: int,
    day: int,
    hour: int,
    minute: int,
    second: float,
    dut1: float,
    tz: float,
) -> float:
    day_decimal = day + (hour - tz + (minute + (second + dut1) / 60.0) / 60.0) / 24.0
    if month < 3:
        month += 12
        year -= 1
    julian_day = (
        math.floor(365.25 * (year + 4716))
        + math.floor(30.6001 * (month + 1))
        + day_decimal
        - 1524.5
    )
    if julian_day > 2299160:
        a = math.floor(year / 100)
        julian_day = julian_day + (2 - a + math.floor(a / 4))
    return julian_day


def _julian_century(jd: float) -> float:
    return (jd - 2451545.0) / 36525.0


def _julian_ephemeris_day(jd: float, delta_t: float) -> float:
    return jd + delta_t / 86400.0


def _julian_ephemeris_century(jde: float) -> float:
    return (jde - 2451545.0) / 36525.0


def _julian_ephemeris_millennium(jce: float) -> float:
    return jce / 10.0


def _earth_periodic_term_summation(terms: np.ndarray, count: int, jme: float) -> float:
    total = 0.0
    for i in range(count):
        total += terms[i, TERM_A] * math.cos(terms[i, TERM_B] + terms[i, TERM_C] * jme)
    return total


def _earth_values(term_sum: np.ndarray, count: int, jme: float) -> float:
    total = 0.0
    for i in range(count):
        total += term_sum[i] * (jme ** i)
    total /= 1e8
    return total


def _earth_heliocentric_longitude(jme: float) -> float:
    sums = np.zeros(L_COUNT, dtype=float)
    for i in range(L_COUNT):
        sums[i] = _earth_periodic_term_summation(L_TERMS[i], l_subcount[i], jme)
    return _limit_degrees(math.degrees(_earth_values(sums, L_COUNT, jme)))


def _earth_heliocentric_latitude(jme: float) -> float:
    sums = np.zeros(B_COUNT, dtype=float)
    for i in range(B_COUNT):
        sums[i] = _earth_periodic_term_summation(B_TERMS[i], b_subcount[i], jme)
    return math.degrees(_earth_values(sums, B_COUNT, jme))


def _earth_radius_vector(jme: float) -> float:
    sums = np.zeros(R_COUNT, dtype=float)
    for i in range(R_COUNT):
        sums[i] = _earth_periodic_term_summation(R_TERMS[i], r_subcount[i], jme)
    return _earth_values(sums, R_COUNT, jme)


def _geocentric_longitude(l: float) -> float:
    theta = l + 180.0
    if theta >= 360.0:
        theta -= 360.0
    return theta


def _geocentric_latitude(b: float) -> float:
    return -b


def _mean_elongation_moon_sun(jce: float) -> float:
    return _third_order_polynomial(1 / 189474, -0.0019142, 445267.11148, 297.85036, jce)


def _mean_anomaly_sun(jce: float) -> float:
    return _third_order_polynomial(-1 / 300000, -0.0001603, 35999.05034, 357.52772, jce)


def _mean_anomaly_moon(jce: float) -> float:
    return _third_order_polynomial(1 / 56250, 0.0086972, 477198.867398, 134.96298, jce)


def _argument_latitude_moon(jce: float) -> float:
    return _third_order_polynomial(1 / 327270, -0.0036825, 483202.017538, 93.27191, jce)


def _ascending_longitude_moon(jce: float) -> float:
    return _third_order_polynomial(1 / 450000, 0.0020708, -1934.136261, 125.04452, jce)


def _xy_term_summation(i: int, x: np.ndarray) -> float:
    total = 0.0
    for j in range(TERM_Y_COUNT):
        total += x[j] * Y_TERMS[i, j]
    return total


def _nutation_longitude_and_obliquity(jce: float, x: np.ndarray) -> Tuple[float, float]:
    sum_psi = 0.0
    sum_epsilon = 0.0
    for i in range(Y_COUNT):
        xy_term_sum = math.radians(_xy_term_summation(i, x))
        sum_psi += (PE_TERMS[i, TERM_PSI_A] + jce * PE_TERMS[i, TERM_PSI_B]) * math.sin(xy_term_sum)
        sum_epsilon += (PE_TERMS[i, TERM_EPS_C] + jce * PE_TERMS[i, TERM_EPS_D]) * math.cos(xy_term_sum)
    del_psi = sum_psi / 36000000.0
    del_eps = sum_epsilon / 36000000.0
    return del_psi, del_eps


def _ecliptic_mean_obliquity(jme: float) -> float:
    u = jme / 10.0
    return (
        84381.448
        + u
        * (
            -4680.93
            + u
            * (
                -1.55
                + u
                * (
                    1999.25
                    + u * (-51.38 + u * (-249.67 + u * (-39.05 + u * (7.12 + u * (27.87 + u * (5.79 + u * 2.45))))))
                )
            )
        )
    )


def _ecliptic_true_obliquity(delta_epsilon: float, epsilon0: float) -> float:
    return delta_epsilon + epsilon0 / 3600.0


def _aberration_correction(r: float) -> float:
    return -20.4898 / (3600.0 * r)


def _apparent_sun_longitude(theta: float, delta_psi: float, delta_tau: float) -> float:
    return theta + delta_psi + delta_tau


def _greenwich_mean_sidereal_time(jd: float, jc: float) -> float:
    return _limit_degrees(280.46061837 + 360.98564736629 * (jd - 2451545.0) + jc * jc * (0.000387933 - jc / 38710000.0))


def _greenwich_sidereal_time(nu0: float, delta_psi: float, epsilon: float) -> float:
    return nu0 + delta_psi * math.cos(math.radians(epsilon))


def _geocentric_right_ascension(lamda: float, epsilon: float, beta: float) -> float:
    lamda_rad = math.radians(lamda)
    epsilon_rad = math.radians(epsilon)
    return _limit_degrees(
        math.degrees(
            math.atan2(
                math.sin(lamda_rad) * math.cos(epsilon_rad) - math.tan(math.radians(beta)) * math.sin(epsilon_rad),
                math.cos(lamda_rad),
            )
        )
    )


def _geocentric_declination(beta: float, epsilon: float, lamda: float) -> float:
    beta_rad = math.radians(beta)
    epsilon_rad = math.radians(epsilon)
    return math.degrees(
        math.asin(math.sin(beta_rad) * math.cos(epsilon_rad) + math.cos(beta_rad) * math.sin(epsilon_rad) * math.sin(math.radians(lamda)))
    )


def _observer_hour_angle(nu: float, longitude: float, alpha_deg: float) -> float:
    return _limit_degrees(nu + longitude - alpha_deg)


def _sun_equatorial_horizontal_parallax(r: float) -> float:
    return 8.794 / (3600.0 * r)


def _right_ascension_parallax_and_topocentric_dec(
    latitude: float, elevation: float, xi: float, h: float, delta: float
) -> Tuple[float, float]:
    lat_rad = math.radians(latitude)
    xi_rad = math.radians(xi)
    h_rad = math.radians(h)
    delta_rad = math.radians(delta)
    u = math.atan(0.99664719 * math.tan(lat_rad))
    y = 0.99664719 * math.sin(u) + elevation * math.sin(lat_rad) / 6378140.0
    x = math.cos(u) + elevation * math.cos(lat_rad) / 6378140.0

    delta_alpha_rad = math.atan2(-x * math.sin(xi_rad) * math.sin(h_rad), math.cos(delta_rad) - x * math.sin(xi_rad) * math.cos(h_rad))
    delta_prime = math.degrees(
        math.atan2(
            (math.sin(delta_rad) - y * math.sin(xi_rad)) * math.cos(delta_alpha_rad),
            math.cos(delta_rad) - x * math.sin(xi_rad) * math.cos(h_rad),
        )
    )
    delta_alpha = math.degrees(delta_alpha_rad)
    return delta_alpha, delta_prime


def _topocentric_right_ascension(alpha_deg: float, delta_alpha: float) -> float:
    return alpha_deg + delta_alpha


def _topocentric_local_hour_angle(h: float, delta_alpha: float) -> float:
    return h - delta_alpha


def _topocentric_elevation_angle(latitude: float, delta_prime: float, h_prime: float) -> float:
    lat_rad = math.radians(latitude)
    delta_prime_rad = math.radians(delta_prime)
    return math.degrees(
        math.asin(math.sin(lat_rad) * math.sin(delta_prime_rad) + math.cos(lat_rad) * math.cos(delta_prime_rad) * math.cos(math.radians(h_prime)))
    )


def _atmospheric_refraction_correction(pressure: float, temperature: float, atmos_refract: float, e0: float) -> float:
    del_e = 0.0
    if e0 >= -1.0 * (SUN_RADIUS + atmos_refract):
        del_e = (pressure / 1010.0) * (283.0 / (273.0 + temperature)) * 1.02 / (60.0 * math.tan(math.radians(e0 + 10.3 / (e0 + 5.11))))
    return del_e


def _topocentric_elevation_angle_corrected(e0: float, delta_e: float) -> float:
    return e0 + delta_e


def _topocentric_zenith_angle(e: float) -> float:
    return 90.0 - e


def _topocentric_azimuth_angle_astro(h_prime: float, latitude: float, delta_prime: float) -> float:
    h_prime_rad = math.radians(h_prime)
    lat_rad = math.radians(latitude)
    return _limit_degrees(
        math.degrees(
            math.atan2(
                math.sin(h_prime_rad),
                math.cos(h_prime_rad) * math.sin(lat_rad) - math.tan(math.radians(delta_prime)) * math.cos(lat_rad),
            )
        )
    )


def _topocentric_azimuth_angle(azimuth_astro: float) -> float:
    return _limit_degrees(azimuth_astro + 180.0)


def _surface_incidence_angle(zenith: float, azimuth_astro: float, azm_rotation: float, slope: float) -> float:
    zenith_rad = math.radians(zenith)
    slope_rad = math.radians(slope)
    return math.degrees(
        math.acos(
            math.cos(zenith_rad) * math.cos(slope_rad)
            + math.sin(slope_rad) * math.sin(zenith_rad) * math.cos(math.radians(azimuth_astro - azm_rotation))
        )
    )


def _sun_mean_longitude(jme: float) -> float:
    return _limit_degrees(
        280.4664567
        + jme * (360007.6982779 + jme * (0.03032028 + jme * (1 / 49931 + jme * (-1 / 15300 + jme * (-1 / 2000000)))))
    )


def _eot(m: float, alpha: float, del_psi: float, epsilon: float) -> float:
    return _limit_minutes(4.0 * (m - 0.0057183 - alpha + del_psi * math.cos(math.radians(epsilon))))


def _approx_sun_transit_time(alpha_zero: float, longitude: float, nu: float) -> float:
    return (alpha_zero - longitude - nu) / 360.0


def _sun_hour_angle_at_rise_set(latitude: float, delta_zero: float, h0_prime: float) -> float:
    latitude_rad = math.radians(latitude)
    delta_zero_rad = math.radians(delta_zero)
    argument = (math.sin(math.radians(h0_prime)) - math.sin(latitude_rad) * math.sin(delta_zero_rad)) / (
        math.cos(latitude_rad) * math.cos(delta_zero_rad)
    )
    if abs(argument) <= 1.0:
        return _limit_degrees180(math.degrees(math.acos(argument)))
    return -99999.0


def _approx_sun_rise_and_set(m_rts: np.ndarray, h0: float) -> np.ndarray:
    h0_dfrac = h0 / 360.0
    m_rts[SUN_RISE] = _limit_zero2one(m_rts[SUN_TRANSIT] - h0_dfrac)
    m_rts[SUN_SET] = _limit_zero2one(m_rts[SUN_TRANSIT] + h0_dfrac)
    m_rts[SUN_TRANSIT] = _limit_zero2one(m_rts[SUN_TRANSIT])
    return m_rts


def _rts_alpha_delta_prime(ad: np.ndarray, n: float) -> float:
    a = ad[JD_ZERO] - ad[JD_MINUS]
    b = ad[JD_PLUS] - ad[JD_ZERO]
    if abs(a) >= 2:
        a = _limit_zero2one(a)
    if abs(b) >= 2:
        b = _limit_zero2one(b)
    return ad[JD_ZERO] + n * (a + b + (b - a) * n) / 2.0


def _rts_sun_altitude(latitude: float, delta_prime: float, h_prime: float) -> float:
    latitude_rad = math.radians(latitude)
    delta_prime_rad = math.radians(delta_prime)
    return math.degrees(
        math.asin(math.sin(latitude_rad) * math.sin(delta_prime_rad) + math.cos(latitude_rad) * math.cos(delta_prime_rad) * math.cos(math.radians(h_prime)))
    )


def _sun_rise_and_set(m_rts: np.ndarray, h_rts: np.ndarray, delta_prime: np.ndarray, latitude: float, h_prime: np.ndarray, h0_prime: float, sun: int) -> float:
    return m_rts[sun] + (h_rts[sun] - h0_prime) / (360.0 * math.cos(math.radians(delta_prime[sun])) * math.cos(math.radians(latitude)) * math.sin(math.radians(h_prime[sun])))


def _calculate_geocentric_sun_right_ascension_and_declination(spa: Dict[str, Any]) -> Dict[str, Any]:
    x = np.zeros(TERM_X_COUNT, dtype=float)

    spa["jc"] = _julian_century(spa["jd"])
    spa["jde"] = _julian_ephemeris_day(spa["jd"], spa["delta_t"])
    spa["jce"] = _julian_ephemeris_century(spa["jde"])
    spa["jme"] = _julian_ephemeris_millennium(spa["jce"])

    spa["l"] = _earth_heliocentric_longitude(spa["jme"])
    spa["b"] = _earth_heliocentric_latitude(spa["jme"])
    spa["r"] = _earth_radius_vector(spa["jme"])

    spa["theta"] = _geocentric_longitude(spa["l"])
    spa["beta"] = _geocentric_latitude(spa["b"])

    spa["x0"] = _mean_elongation_moon_sun(spa["jce"])
    x[TERM_X0] = spa["x0"]
    spa["x1"] = _mean_anomaly_sun(spa["jce"])
    x[TERM_X1] = spa["x1"]
    spa["x2"] = _mean_anomaly_moon(spa["jce"])
    x[TERM_X2] = spa["x2"]
    spa["x3"] = _argument_latitude_moon(spa["jce"])
    x[TERM_X3] = spa["x3"]
    spa["x4"] = _ascending_longitude_moon(spa["jce"])
    x[TERM_X4] = spa["x4"]

    spa["del_psi"], spa["del_epsilon"] = _nutation_longitude_and_obliquity(spa["jce"], x)

    spa["epsilon0"] = _ecliptic_mean_obliquity(spa["jme"])
    spa["epsilon"] = _ecliptic_true_obliquity(spa["del_epsilon"], spa["epsilon0"])

    spa["del_tau"] = _aberration_correction(spa["r"])
    spa["lamda"] = _apparent_sun_longitude(spa["theta"], spa["del_psi"], spa["del_tau"])
    spa["nu0"] = _greenwich_mean_sidereal_time(spa["jd"], spa["jc"])
    spa["nu"] = _greenwich_sidereal_time(spa["nu0"], spa["del_psi"], spa["epsilon"])

    spa["alpha"] = _geocentric_right_ascension(spa["lamda"], spa["epsilon"], spa["beta"])
    spa["delta"] = _geocentric_declination(spa["beta"], spa["epsilon"], spa["lamda"])
    return spa


def _calculate_eot_and_sun_rise_transit_set(spa: Dict[str, Any]) -> Dict[str, Any]:
    alpha = np.zeros(JD_COUNT, dtype=float)
    delta = np.zeros(JD_COUNT, dtype=float)
    m_rts = np.zeros(SUN_COUNT, dtype=float)
    nu_rts = np.zeros(SUN_COUNT, dtype=float)
    h_rts = np.zeros(SUN_COUNT, dtype=float)
    alpha_prime = np.zeros(SUN_COUNT, dtype=float)
    delta_prime = np.zeros(SUN_COUNT, dtype=float)
    h_prime = np.zeros(SUN_COUNT, dtype=float)
    h0_prime = -1.0 * (SUN_RADIUS + spa["atmos_refract"])

    sun_rts = dict(spa)
    m = _sun_mean_longitude(spa["jme"])
    spa["eot"] = _eot(m, spa["alpha"], spa["del_psi"], spa["epsilon"])
    sun_rts["hour"] = 0
    sun_rts["minute"] = 0
    sun_rts["second"] = 0
    sun_rts["delta_ut1"] = 0
    sun_rts["timezone"] = 0

    sun_rts["jd"] = _julian_day(
        sun_rts["year"],
        sun_rts["month"],
        sun_rts["day"],
        sun_rts["hour"],
        sun_rts["minute"],
        sun_rts["second"],
        sun_rts["delta_ut1"],
        sun_rts["timezone"],
    )

    sun_rts = _calculate_geocentric_sun_right_ascension_and_declination(sun_rts)
    nu = sun_rts["nu"]

    sun_rts["delta_t"] = 0
    sun_rts["jd"] = sun_rts["jd"] - 1
    for i in range(JD_COUNT):
        sun_rts = _calculate_geocentric_sun_right_ascension_and_declination(sun_rts)
        alpha[i] = sun_rts["alpha"]
        delta[i] = sun_rts["delta"]
        sun_rts["jd"] = sun_rts["jd"] + 1

    m_rts[SUN_TRANSIT] = _approx_sun_transit_time(alpha[JD_ZERO], spa["longitude"], nu)
    h0 = _sun_hour_angle_at_rise_set(spa["latitude"], delta[JD_ZERO], h0_prime)

    if h0 >= 0:
        m_rts = _approx_sun_rise_and_set(m_rts, h0)
        for i in range(SUN_COUNT):
            nu_rts[i] = nu + 360.985647 * m_rts[i]
            n = m_rts[i] + spa["delta_t"] / 86400.0
            alpha_prime[i] = _rts_alpha_delta_prime(alpha, n)
            delta_prime[i] = _rts_alpha_delta_prime(delta, n)
            h_prime[i] = _limit_degrees180pm(nu_rts[i] + spa["longitude"] - alpha_prime[i])
            h_rts[i] = _rts_sun_altitude(spa["latitude"], delta_prime[i], h_prime[i])

        spa["srha"] = h_prime[SUN_RISE]
        spa["ssha"] = h_prime[SUN_SET]
        spa["sta"] = h_rts[SUN_TRANSIT]

        spa["suntransit"] = _dayfrac_to_local_hr(m_rts[SUN_TRANSIT] - h_prime[SUN_TRANSIT] / 360.0, spa["timezone"])
        spa["sunrise"] = _dayfrac_to_local_hr(
            _sun_rise_and_set(m_rts, h_rts, delta_prime, spa["latitude"], h_prime, h0_prime, SUN_RISE), spa["timezone"]
        )
        spa["sunset"] = _dayfrac_to_local_hr(
            _sun_rise_and_set(m_rts, h_rts, delta_prime, spa["latitude"], h_prime, h0_prime, SUN_SET), spa["timezone"]
        )
    else:
        spa["srha"] = -99999.0
        spa["ssha"] = -99999.0
        spa["sta"] = -99999.0
        spa["suntransit"] = -99999.0
        spa["sunrise"] = -99999.0
        spa["sunset"] = -99999.0

    return spa


def _spa_calculate(spa: Dict[str, Any]) -> Tuple[int, Dict[str, Any]]:
    result = _validate_inputs(spa)
    if result != 0:
        return result, spa

    spa["jd"] = _julian_day(
        spa["year"],
        spa["month"],
        spa["day"],
        spa["hour"],
        spa["minute"],
        spa["second"],
        spa["delta_ut1"],
        spa["timezone"],
    )

    spa = _calculate_geocentric_sun_right_ascension_and_declination(spa)

    spa["h"] = _observer_hour_angle(spa["nu"], spa["longitude"], spa["alpha"])
    spa["xi"] = _sun_equatorial_horizontal_parallax(spa["r"])

    spa["del_alpha"], spa["delta_prime"] = _right_ascension_parallax_and_topocentric_dec(
        spa["latitude"], spa["elevation"], spa["xi"], spa["h"], spa["delta"]
    )
    spa["alpha_prime"] = _topocentric_right_ascension(spa["alpha"], spa["del_alpha"])
    spa["h_prime"] = _topocentric_local_hour_angle(spa["h"], spa["del_alpha"])

    spa["e0"] = _topocentric_elevation_angle(spa["latitude"], spa["delta_prime"], spa["h_prime"])
    spa["del_e"] = _atmospheric_refraction_correction(spa["pressure"], spa["temperature"], spa["atmos_refract"], spa["e0"])
    spa["e"] = _topocentric_elevation_angle_corrected(spa["e0"], spa["del_e"])

    spa["zenith"] = _topocentric_zenith_angle(spa["e"])
    spa["azimuth_astro"] = _topocentric_azimuth_angle_astro(spa["h_prime"], spa["latitude"], spa["delta_prime"])
    spa["azimuth"] = _topocentric_azimuth_angle(spa["azimuth_astro"])

    if spa["function"] in (SPA_ZA_INC, SPA_ALL):
        spa["incidence"] = _surface_incidence_angle(spa["zenith"], spa["azimuth_astro"], spa["azm_rotation"], spa["slope"])
    if spa["function"] in (SPA_ZA_RTS, SPA_ALL):
        spa = _calculate_eot_and_sun_rise_transit_set(spa)

    return 0, spa


def solar_position_python(time_of_day: datetime, longitude: float, latitude: float, timezone: float, elevation: float) -> Dict[str, Any]:
    """
    Port of tools/SEB/SPA/solarPosition.m. Returns the SPA struct as a dict.
    """
    spa: Dict[str, Any] = {}
    spa["year"] = time_of_day.year
    spa["month"] = time_of_day.month
    spa["day"] = time_of_day.day
    spa["hour"] = time_of_day.hour
    spa["minute"] = time_of_day.minute
    spa["second"] = time_of_day.second + time_of_day.microsecond / 1.0e6
    spa["timezone"] = timezone
    spa["delta_ut1"] = 0.0
    spa["delta_t"] = 0.0
    spa["longitude"] = longitude
    spa["latitude"] = latitude
    spa["elevation"] = elevation
    spa["pressure"] = 1010.0
    spa["temperature"] = 10.0
    spa["slope"] = 0.0
    spa["azm_rotation"] = 0.0
    spa["atmos_refract"] = 0.5667
    spa["function"] = SPA_ALL

    result, spa = _spa_calculate(spa)
    if result != 0:
        raise ValueError(f"SPA error code: {result}")
    return spa


def nsun_from_angles(zenith_deg: float, azimuth_deg: float) -> np.ndarray:
    zenith_rad = math.radians(zenith_deg)
    azimuth_rad = math.radians(azimuth_deg)
    nsun = np.array(
        [
            math.sin(zenith_rad) * math.cos(azimuth_rad),
            math.sin(zenith_rad) * -math.sin(azimuth_rad),
            math.cos(zenith_rad),
        ],
        dtype=float,
    )
    return nsun


def solar_state_python(
    time_of_day: datetime,
    longitude: float,
    latitude: float,
    timezone: float,
    elevation: float,
    *,
    xazimuth: float = 0.0,
) -> Tuple[np.ndarray, float, float, float, float]:
    """
    Solar state using the Python SPA port + ASHRAE clear-sky.
    """
    spa = solar_position_python(time_of_day, longitude, latitude, timezone, elevation)
    zenith = float(spa["zenith"])
    azimuth = float(spa["azimuth"])
    azimuth_local = azimuth - float(xazimuth)
    nsun = nsun_from_angles(zenith, azimuth_local)
    I, Dsky = solar_strength_ashrae(time_of_day.month, zenith)
    return nsun, zenith, azimuth_local, I, Dsky


def solar_state_pvlib(
    time_of_day: datetime,
    longitude: float,
    latitude: float,
    timezone: float,
    elevation: float,
    *,
    xazimuth: float = 0.0,
) -> Tuple[np.ndarray, float, float, float, float]:
    """
    Solar state using pvlib for angles + ASHRAE for irradiance.
    """
    try:
        import pandas as pd
        import pvlib
    except ImportError as exc:
        raise ImportError("pvlib backend requested but pvlib is not installed") from exc

    tzinfo = dt_timezone(timedelta(hours=timezone))
    if time_of_day.tzinfo is None:
        time_of_day = time_of_day.replace(tzinfo=tzinfo)
    else:
        time_of_day = time_of_day.astimezone(tzinfo)

    times = pd.DatetimeIndex([time_of_day])
    sp = pvlib.solarposition.get_solarposition(times, latitude, longitude, altitude=elevation)
    zenith = float(sp["zenith"].iloc[0])
    azimuth = float(sp["azimuth"].iloc[0])
    azimuth_local = azimuth - float(xazimuth)
    nsun = nsun_from_angles(zenith, azimuth_local)
    I, Dsky = solar_strength_ashrae(time_of_day.month, zenith)
    return nsun, zenith, azimuth_local, I, Dsky


def solar_state(
    time_of_day: datetime,
    longitude: float,
    latitude: float,
    timezone: float,
    elevation: float,
    *,
    xazimuth: float = 0.0,
    backend: str = "python",
) -> Tuple[np.ndarray, float, float, float, float]:
    """
    Wrapper to select backend ("python" or "pvlib").
    """
    backend_key = backend.strip().lower()
    if backend_key == "python":
        return solar_state_python(time_of_day, longitude, latitude, timezone, elevation, xazimuth=xazimuth)
    if backend_key == "pvlib":
        return solar_state_pvlib(time_of_day, longitude, latitude, timezone, elevation, xazimuth=xazimuth)
    raise ValueError(f"Unknown solar backend: {backend}")
