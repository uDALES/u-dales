"""Exercise air-temperature output in the moist NetCDF output paths."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import f90nml
import numpy as np
from netCDF4 import Dataset


ROOT = Path(__file__).resolve().parents[3]
CASE = 994
FIXTURE = ROOT / "tests" / "system" / "experiments" / str(CASE)
EXECUTABLE = Path(os.environ.get("UDALES_EXE", ROOT / "bin" / "u-dales")).resolve()


class AirTemperatureOutputTest(unittest.TestCase):
    def _run_case(
        self, case_dir: Path, fieldvars: str, slicevars: str, probevars: str, *, ltempeq: bool = True
    ) -> None:
        if not EXECUTABLE.is_file():
            self.skipTest(f"uDALES executable not found: {EXECUTABLE}")
        mpiexec = shutil.which(os.environ.get("MPIEXEC", "mpiexec"))
        if mpiexec is None:
            self.skipTest("MPI launcher not found")

        shutil.copytree(FIXTURE, case_dir, dirs_exist_ok=True)

        namelist = f90nml.read(case_dir / f"namoptions.{CASE}")
        run = namelist["run"]
        run.update(runtime=3.0, dtmax=1.0)
        namelist["physics"]["ltempeq"] = ltempeq
        output = namelist["output"]
        output.pop("tfielddump", None)
        output.update(
            fieldvars=fieldvars,
            slicevars=slicevars,
            probevars=probevars,
            lislicedump=True,
            ljslicedump=True,
            lkslicedump=True,
            ltislicedump=True,
            ltjslicedump=True,
            ltkslicedump=True,
            lprobedump=True,
            lxydump=True,
            lytdump=True,
            lydump=True,
            tinstantdump=1.0,
            tstatsdump=1.0,
            tsample=1.0,
            nislice=1,
            njslice=1,
            nkslice=1,
            islice=[32],
            jslice=[16],
            kslice=[12],
            nprobe=1,
        )
        namelist.write(case_dir / f"namoptions.{CASE}", force=True)
        (case_dir / f"probe.inp.{CASE}").write_text("# i j k\n32 16 12\n", encoding="ascii")

        profile_path = case_dir / f"prof.inp.{CASE}"
        profile = np.loadtxt(profile_path, comments="#")
        profile[:, 2] = 0.012
        np.savetxt(profile_path, profile, header="SDBL flow\nz thl qt u v tke", comments="# ")

        version = subprocess.run([mpiexec, "--version"], capture_output=True, text=True, check=False).stdout
        command = [mpiexec]
        if re.search(r"Open MPI|OpenRTE", version, flags=re.IGNORECASE):
            command.append("--oversubscribe")
        command += ["-n", "8", str(EXECUTABLE), f"namoptions.{CASE}"]
        result = subprocess.run(command, cwd=case_dir, capture_output=True, text=True, timeout=120, check=False)
        if ltempeq:
            self.assertEqual(result.returncode, 0, result.stdout[-4000:] + result.stderr[-4000:])
        else:
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("'ta' output requires ltempeq=.true.", result.stdout + result.stderr)

    def test_moist_output_families(self) -> None:
        with tempfile.TemporaryDirectory(prefix="udales_air_temp_") as tmp:
            case_dir = Path(tmp)
            self._run_case(case_dir, "th,ta,qt,ql", "th,ta,qt", "th,ta,qt")

            families = (
                "ins_field", "ins_islice", "ins_jslice", "ins_kslice", "ins_probe",
                "stats_t", "stats_islice", "stats_jslice", "stats_kslice",
                "stats_xyt", "stats_xy", "stats_yt", "stats_y",
            )
            for family in families:
                with self.subTest(family=family):
                    files = list(case_dir.glob(f"{family}.*.{CASE}.nc"))
                    files += list(case_dir.glob(f"{family}.{CASE}.nc"))
                    self.assertTrue(files, f"No {family} files found")
                    for path in files:
                        with Dataset(path) as ds:
                            self.assertIn("thl", ds.variables)
                            self.assertIn("tha", ds.variables)
                            self.assertEqual(ds["tha"].getncattr("units"), "K")
                            values = np.asarray(ds["tha"][:])
                            self.assertTrue(np.isfinite(values).all())
                            self.assertGreater(values.min(), 150.0)
                            self.assertLess(values.max(), 350.0)

            field_file = next(case_dir.glob(f"ins_field.*.{CASE}.nc"))
            with Dataset(field_file) as ds:
                thl = np.asarray(ds["thl"][:])
                ql = np.asarray(ds["ql"][:])
                temp = np.asarray(ds["tha"][:])
                self.assertGreater(ds["qt"][:].max(), 0.0)
                self.assertGreaterEqual(ql.min(), 0.0)
                exner = (temp - (2.26e6 / 1004.0) * ql) / thl
                self.assertTrue(np.isfinite(exner).all())
                self.assertTrue(((exner > 0.9) & (exner < 1.1)).all())
                self.assertLess(np.max(np.abs(exner - exner.mean(axis=(-2, -1), keepdims=True))), 2e-6)

            with (
                Dataset(case_dir / f"ins_field.001.{CASE}.nc") as field,
                Dataset(case_dir / f"ins_islice.001.{CASE}.nc") as islice,
                Dataset(case_dir / f"ins_jslice.000.{CASE}.nc") as jslice,
                Dataset(case_dir / f"ins_kslice.001.{CASE}.nc") as kslice,
                Dataset(case_dir / f"ins_probe.{CASE}.nc") as probe,
            ):
                expected = float(field["tha"][0, 11, 15, 15])
                self.assertAlmostEqual(float(islice["tha"][0, 11, 15, 0]), expected, places=4)
                self.assertAlmostEqual(float(jslice["tha"][0, 11, 0, 31]), expected, places=4)
                self.assertAlmostEqual(float(kslice["tha"][0, 0, 15, 15]), expected, places=4)
                self.assertAlmostEqual(float(probe["tha"][0, 0]), expected, places=4)

            with (
                Dataset(case_dir / f"stats_t.001.{CASE}.nc") as full,
                Dataset(case_dir / f"stats_islice.001.{CASE}.nc") as islice,
                Dataset(case_dir / f"stats_jslice.000.{CASE}.nc") as jslice,
                Dataset(case_dir / f"stats_kslice.001.{CASE}.nc") as kslice,
            ):
                expected = float(full["tha"][0, 11, 15, 15])
                self.assertAlmostEqual(float(islice["tha"][0, 11, 15, 0]), expected, places=4)
                self.assertAlmostEqual(float(jslice["tha"][0, 11, 0, 31]), expected, places=4)
                self.assertAlmostEqual(float(kslice["tha"][0, 0, 15, 15]), expected, places=4)

    def test_instantaneous_ta_selection_is_independent(self) -> None:
        baseline_thl = {}
        selections = (
            ("th", ("thl",), ("tha",)),
            ("th,ta", ("thl", "tha"), ()),
            ("ta", ("tha",), ("thl",)),
        )
        for selected, present, absent in selections:
            with self.subTest(selected=selected), tempfile.TemporaryDirectory(prefix="udales_air_temp_") as tmp:
                case_dir = Path(tmp)
                self._run_case(case_dir, selected, selected, selected)
                for family in ("ins_field", "ins_islice", "ins_jslice", "ins_kslice", "ins_probe"):
                    files = list(case_dir.glob(f"{family}.*.{CASE}.nc"))
                    files += list(case_dir.glob(f"{family}.{CASE}.nc"))
                    self.assertTrue(files)
                    for path in files:
                        with Dataset(path) as ds:
                            for variable in present:
                                self.assertIn(variable, ds.variables)
                            for variable in absent:
                                self.assertNotIn(variable, ds.variables)
                            if selected == "th":
                                baseline_thl[path.name] = np.asarray(ds["thl"][:])
                            elif selected == "th,ta":
                                np.testing.assert_array_equal(ds["thl"][:], baseline_thl[path.name])

    def test_ta_requires_temperature_equation(self) -> None:
        with tempfile.TemporaryDirectory(prefix="udales_air_temp_") as tmp:
            self._run_case(Path(tmp), "ta", "ta", "ta", ltempeq=False)


if __name__ == "__main__":
    unittest.main()
