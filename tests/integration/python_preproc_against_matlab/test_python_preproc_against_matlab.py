import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
CLEANSIM = REPO_ROOT / "tools" / "cleansim.sh"
WRITE_INPUTS = REPO_ROOT / "tools" / "write_inputs.py"
CASE_SPECS = [
    {
        "case_dir": REPO_ROOT / "tests" / "cases" / "100",
        "outputs": [
            "lscale.inp.100",
            "prof.inp.100",
            "factypes.inp.100",
            "facets.inp.100",
            "facetarea.inp.100",
            "solid_c.txt",
            "solid_u.txt",
            "solid_v.txt",
            "solid_w.txt",
            "fluid_boundary_c.txt",
            "fluid_boundary_u.txt",
            "fluid_boundary_v.txt",
            "fluid_boundary_w.txt",
            "facet_sections_c.txt",
            "facet_sections_u.txt",
            "facet_sections_v.txt",
            "facet_sections_w.txt",
        ],
    },
    {
        "case_dir": REPO_ROOT / "tests" / "cases" / "064",
        "outputs": [
            "lscale.inp.064",
            "prof.inp.064",
            "factypes.inp.064",
            "facets.inp.064",
            "facetarea.inp.064",
            "solid_c.txt",
            "solid_u.txt",
            "solid_v.txt",
            "solid_w.txt",
            "fluid_boundary_c.txt",
            "fluid_boundary_u.txt",
            "fluid_boundary_v.txt",
            "fluid_boundary_w.txt",
            "facet_sections_c.txt",
            "facet_sections_u.txt",
            "facet_sections_v.txt",
            "facet_sections_w.txt",
            "svf.inp.064",
            "vfsparse.inp.064",
            "netsw.inp.064",
            "Tfacinit.inp.064",
        ],
    },
]


def _atol_for_output(relpath: str) -> float:
    if relpath.startswith("facet_sections_"):
        # These text files are emitted with 4 decimal places for area and
        # distance, so cross-platform runs can legitimately differ in the last
        # printed digit after floating-point rounding.
        return 1.1e-4
    if relpath.startswith("netsw.inp."):
        return 2.0e-3
    return 1.0e-10


class TestPythonPreprocAgainstMatlab(unittest.TestCase):
    """Compare Python preprocessing against committed MATLAB reference outputs."""

    @staticmethod
    def _load_numeric_table(path: Path) -> np.ndarray:
        data = np.loadtxt(path, comments="#")
        if np.ndim(data) == 0:
            return np.asarray([[data]])
        return np.atleast_2d(data)

    def test_reference_cases_match_committed_matlab_outputs(self) -> None:
        for spec in CASE_SPECS:
            case_source = spec["case_dir"]
            outputs = spec["outputs"]
            with self.subTest(case=case_source.name):
                with tempfile.TemporaryDirectory(prefix="udales-preproc-test-", dir="/tmp") as temp_root:
                    temp_root = Path(temp_root)
                    temp_case = temp_root / case_source.name
                    shutil.copytree(case_source, temp_case)

                    clean = subprocess.run(
                        ["bash", str(CLEANSIM), str(temp_case), "--preprocess-only", "--delete"],
                        cwd=REPO_ROOT,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                    )
                    if clean.returncode != 0:
                        self.fail(f"cleansim failed for case {case_source.name}:\n{clean.stdout}")

                    result = subprocess.run(
                        [sys.executable, str(WRITE_INPUTS), str(temp_case)],
                        cwd=REPO_ROOT,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                    )
                    if result.returncode != 0:
                        self.fail(f"Python preprocessing failed for case {case_source.name}:\n{result.stdout}")

                    mismatches = []
                    for relpath in outputs:
                        reference = case_source / relpath
                        regenerated = temp_case / relpath
                        if not reference.exists():
                            mismatches.append(f"{relpath}: missing committed reference output")
                            continue
                        if not regenerated.exists():
                            mismatches.append(f"{relpath}: missing regenerated output")
                            continue

                        ref_data = self._load_numeric_table(reference)
                        gen_data = self._load_numeric_table(regenerated)
                        if ref_data.shape != gen_data.shape:
                            mismatches.append(
                                f"{relpath}: shape mismatch reference {ref_data.shape} vs regenerated {gen_data.shape}"
                            )
                            continue
                        try:
                            np.testing.assert_allclose(
                                gen_data,
                                ref_data,
                                rtol=0.0,
                                atol=_atol_for_output(relpath),
                            )
                        except AssertionError as exc:
                            mismatches.append(f"{relpath}: {exc}")

                    if mismatches:
                        self.fail(
                            f"Preprocessing comparison failed for case {case_source.name}:\n"
                            + "\n".join(mismatches)
                        )


if __name__ == "__main__":
    unittest.main()
