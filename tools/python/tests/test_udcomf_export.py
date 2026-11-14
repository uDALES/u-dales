"""Pre-run comfort-output checks on small synthetic namelists."""

from datetime import date
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import f90nml

from udcomf.udcomf_export import check_comfort_preflight


def _case(path: Path, *, first_level: float = 0.75, gap: int = 2700) -> Path:
    case = path / "123"
    case.mkdir()
    (case / "namoptions.123").write_text(
        "&RUN\n"
        " runtime=108060.\n"
        "/\n"
        "&OUTPUT\n"
        " receptor_height=1.1\n"
        " ltdump=.true.\n"
        " tstatstart=21600.\n"
        " tstatsdump=900.\n"
        f" tstatsgap={gap}.\n"
        " tsample=5.\n"
        "/\n"
        "&DOMAIN\n"
        " itot=10\n jtot=8\n ktot=4\n"
        "/\n"
        "&PHYSICS\n"
        " ltempeq=.true.\n lmoist=.true.\n"
        " ltimedepsw=.true.\n ltimedeplw=.true.\n"
        "/\n"
        "&ENERGYBALANCE\n lEB=.true.\n lwriteEBfiles=.true.\n/\n"
        "&SOLAR\n"
        " year=2023\n month=8\n day=20\n hour=18\n"
        " minute=0\n second=0\n timezone=0.\n"
        "/\n",
        encoding="ascii",
    )
    (case / "prof.inp.123").write_text(
        "# profile\n# z thl\n"
        + "".join(f"{z} 300\n" for z in (first_level, 2.25, 9.75, 11.25)),
        encoding="ascii",
    )
    return case


class TestComfortPreflight(unittest.TestCase):
    def test_15_minute_hourly_schedule_and_bracket_are_ready(self):
        with TemporaryDirectory() as tmp:
            case = _case(Path(tmp))
            result = check_comfort_preflight(case, date(2023, 8, 21))
            self.assertTrue(result.ready, result.blockers)
            self.assertEqual(result.analysis_times_s[0], 25200.0)
            self.assertEqual(result.analysis_times_s[-1], 108000.0)
            self.assertEqual(result.estimated_stats_bytes, 10 * 8 * 4 * 28 * 4 * 24)

    def test_bad_cadence_and_unbracketed_receptor_fail(self):
        with TemporaryDirectory() as tmp:
            case = _case(Path(tmp), first_level=1.25, gap=13500)
            result = check_comfort_preflight(case, date(2023, 8, 21))
            self.assertFalse(result.ready)
            self.assertTrue(any("ws_local" in issue for issue in result.blockers))
            self.assertTrue(any("hourly" in issue for issue in result.blockers))

    def test_multi_height_kslice_only_requires_each_bracket(self):
        with TemporaryDirectory() as tmp:
            case = _case(Path(tmp))
            path = case / "namoptions.123"
            namelist = f90nml.read(path)
            namelist["output"].update(
                ltdump=False, ltkslicedump=True, nkslice=2, kslice=[1, 2],
                nreceptor_heights=2, receptor_heights=[1.1, 1.5],
            )
            namelist.write(path, force=True)
            report = check_comfort_preflight(case, date(2023, 8, 21))
            self.assertTrue(report.ready, report.blockers)
            self.assertEqual(report.receptor_heights_m, (1.1, 1.5))
            self.assertLess(report.estimated_stats_bytes, 10 * 8 * 4 * 28 * 4 * 24)

            namelist["output"]["nkslice"] = 1
            namelist["output"]["kslice"] = [2]
            namelist.write(path, force=True)
            report = check_comfort_preflight(case, date(2023, 8, 21))
            self.assertFalse(report.ready)
            self.assertTrue(any("kslice is missing" in issue for issue in report.blockers))
