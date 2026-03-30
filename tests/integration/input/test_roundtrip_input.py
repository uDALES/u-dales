"""Standalone namelist round-trip validation test.

This test:
1. prepares an input case
2. runs u-dales in TEST_ROUNDTRIP mode so the last rank writes `namoptions_roundtrip`
3. converts both the input namelist and generated namelist via `ud_nam2json`
4. compares generated values against explicit input values or schema defaults
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace


def _normalize_value(v):
    return v


def _compare_generated(generated, input_cfg, schema):
    input_mismatches = []
    default_mismatches = []
    missing_defaults = []

    schema_props = schema.get("properties", {})

    def find_schema_section(name):
        return schema_props.get(name)

    def find_prop_entry(sec_props, pname):
        if not isinstance(sec_props, dict):
            return None
        return sec_props.get(pname)

    def find_input_value(input_section, pname):
        if not isinstance(input_section, dict):
            return None
        return input_section.get(pname)

    for section, params in generated.items():
        if not isinstance(params, dict):
            continue
        input_section = input_cfg.get(section, {})
        sec_schema = find_schema_section(section) or {}
        sec_props = sec_schema.get("properties", {}) if isinstance(sec_schema, dict) else {}

        for pname, pval in params.items():
            gen_val = _normalize_value(pval)

            inp_raw = find_input_value(input_section, pname)
            if inp_raw is not None:
                inp_val = _normalize_value(inp_raw)
                try:
                    if isinstance(inp_val, (int, float)) and isinstance(gen_val, (int, float)):
                        denom = max(abs(float(inp_val)), 1e-12)
                        relerr = abs(float(gen_val) - float(inp_val)) / denom
                        if relerr > 1e-3:
                            input_mismatches.append((section, pname, "input", inp_val, gen_val))
                    else:
                        if inp_val != gen_val:
                            input_mismatches.append((section, pname, "input", inp_val, gen_val))
                except Exception:
                    if inp_val != gen_val:
                        input_mismatches.append((section, pname, "input", inp_val, gen_val))
            else:
                p_schema = find_prop_entry(sec_props, pname) or {}
                if section.upper() == "RUN" and pname.lower() in ("nprocx", "nprocy"):
                    continue
                if "default" in p_schema:
                    default = p_schema["default"]
                    default_norm = _normalize_value(default)
                    try:
                        if isinstance(default_norm, (int, float)) and isinstance(gen_val, (int, float)):
                            denom = max(abs(float(default_norm)), 1e-12)
                            relerr = abs(float(gen_val) - float(default_norm)) / denom
                            if relerr > 1e-3:
                                default_mismatches.append((section, pname, "default", default_norm, gen_val))
                        else:
                            if default_norm != gen_val:
                                default_mismatches.append((section, pname, "default", default_norm, gen_val))
                    except Exception:
                        if default_norm != gen_val:
                            default_mismatches.append((section, pname, "default", default_norm, gen_val))
                else:
                    missing_defaults.append((section, pname, gen_val))

    return input_mismatches + default_mismatches, missing_defaults


def _set_runmode_and_fill_missing_groups(text):
    required_groups = [
        "PHYSICS",
        "DYNAMICS",
        "BC",
        "INLET",
        "DRIVER",
        "WALLS",
        "ENERGYBALANCE",
        "SCALARS",
        "CHEMISTRY",
        "TREES",
        "PURIFS",
        "HEATPUMP",
        "OUTPUT",
        "NAMSUBGRID",
    ]

    if "runmode" not in text.lower():
        match = re.search(r"(?ims)^&RUN\s*(.*?)^\s*/", text)
        if match is not None:
            body = match.group(1)
            replacement = "&RUN\n" + body.rstrip() + "\nrunmode      = 1001\n/"
            text = text[:match.start()] + replacement + text[match.end():]

    for group in required_groups:
        if f"&{group}" not in text:
            text = text.rstrip() + f"\n\n&{group}\n/\n"

    return text


def _copy_case_directory(src_case_dir, workdir):
    for item in src_case_dir.iterdir():
        dst = workdir / item.name
        if item.is_dir():
            if item.name == "warmstart_files":
                continue
            shutil.copytree(item, dst)
        else:
            shutil.copy2(item, dst)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run u-DALES namelist round-trip test")
    parser.add_argument("--comparator-only", action="store_true", help="Skip running u-DALES and compare an existing generated file")
    parser.add_argument("--generated", help="Path to generated file to compare (overrides defaults)")
    parser.add_argument("--mode", choices=["default", "random"], default="default", help="Test mode. 'random' is currently not implemented for namelist input")
    args = parser.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[3]

    udales_path = shutil.which("u-dales")
    ud_nam2json_path = shutil.which("ud_nam2json")
    ud_topdir = os.environ.get("UD_TOPDIR")

    if udales_path:
        udales_exe = Path(udales_path)
    elif ud_topdir:
        udales_exe = Path(ud_topdir) / "bin" / "u-dales"
    else:
        udales_exe = repo_root / "bin" / "u-dales"

    if ud_nam2json_path:
        ud_nam2json = Path(ud_nam2json_path)
    elif ud_topdir:
        ud_nam2json = Path(ud_topdir) / "bin" / "ud_nam2json"
    else:
        ud_nam2json = repo_root / "bin" / "ud_nam2json"

    if ud_topdir:
        schema_path = Path(ud_topdir) / "docs" / "schemas" / "udales_input_schema.json"
        input_dir = Path(ud_topdir) / "tests" / "integration" / "input"
    else:
        schema_path = repo_root / "docs" / "schemas" / "udales_input_schema.json"
        input_dir = repo_root / "tests" / "integration" / "input"

    if not schema_path.exists():
        print(f"SKIP: schema not found: {schema_path}")
        return 0

    if args.mode == "random":
        print("SKIP: random namelist generation is not implemented on this branch")
        return 0

    comparator_only = bool(args.comparator_only)

    if not comparator_only:
        if shutil.which("mpirun") is None:
            print("SKIP: mpirun not available")
            return 0
        if not udales_exe.exists():
            print(f"SKIP: u-DALES executable not found: {udales_exe}")
            return 0
        if not ud_nam2json.exists():
            print(f"SKIP: ud_nam2json not found: {ud_nam2json}")
            return 0

    src_case_dir = repo_root / "examples" / "102"
    src_namelist = input_dir / "namoptions.default"
    if not src_namelist.exists():
        src_namelist = src_case_dir / "namoptions.102"
    if not src_namelist.exists():
        print(f"SKIP: namelist input not found: {src_namelist}")
        return 0

    tempdir = tempfile.TemporaryDirectory(prefix="udales-tmp2-")
    workdir = Path(tempdir.name)

    with open(schema_path, "r") as f:
        schema = json.load(f)

    overall_ret = 0
    tag = "default"

    if comparator_only:
        if args.generated:
            generated_json_path = Path(args.generated)
        else:
            generated_json_path = input_dir / "parameters_roundtrip"
        if not generated_json_path.exists():
            print(f"ERROR: generated round-trip file not found for comparator-only: {generated_json_path}")
            return 2
        input_json_path = input_dir / "parameters_input"
        if not input_json_path.exists():
            print(f"ERROR: comparator-only input JSON not found: {input_json_path}")
            return 2
        with open(input_json_path, "r") as f:
            input_config = json.load(f)
    else:
        _copy_case_directory(src_case_dir, workdir)

        source_text = src_namelist.read_text()
        work_input = workdir / "namoptions.102"
        work_input.write_text(_set_runmode_and_fill_missing_groups(source_text))

        proc_in = subprocess.run(
            [str(ud_nam2json), str(work_input), "parameters_input"],
            cwd=str(workdir),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=30,
        )
        if proc_in.returncode != 0:
            print("ud_nam2json failed on input namelist:")
            print(proc_in.stdout)
            print(proc_in.stderr)
            return 2

        input_json_path = workdir / "parameters_input"
        with open(input_json_path, "r") as f:
            input_config = json.load(f)

        nprocx = int(input_config.get("RUN", {}).get("nprocx", 1))
        nprocy = int(input_config.get("RUN", {}).get("nprocy", 1))
        total_procs = nprocx * nprocy

        try:
            cmd1 = f"mpirun -np {total_procs} {udales_exe} {work_input.name}"
            print("Running:", cmd1)
            proc = subprocess.run(
                cmd1,
                cwd=str(workdir),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                timeout=300,
            )
        except Exception as e:
            print("ERROR: running u-DALES failed to start:", e)
            return 2

        if proc.returncode != 0:
            print("u-DALES failed:")
            print(proc.stdout)
            print(proc.stderr)
            return 2

        try:
            cmd2 = f"{ud_nam2json} namoptions_roundtrip parameters_roundtrip"
            print("Running:", cmd2, "in", str(workdir))
            proc2 = subprocess.run(
                cmd2,
                cwd=str(workdir),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                timeout=30,
            )
        except Exception as e:
            print("ERROR: running ud_nam2json failed to start:", e)
            return 2

        if proc2.returncode != 0:
            print("ud_nam2json failed:")
            print(proc2.stdout)
            print(proc2.stderr)
            return 2

        generated_json_path = workdir / "parameters_roundtrip"
        if not generated_json_path.exists():
            print(f"ERROR: expected json file {generated_json_path} not found")
            return 2

    mismatches = []
    missing_defaults = []
    proc_cmp = SimpleNamespace(returncode=0, stdout="", stderr="")
    try:
        gen = json.loads(generated_json_path.read_text())
        mismatches, missing_defaults = _compare_generated(gen, input_config, schema)
        proc_cmp.stdout = "Generated in-memory diff report"
        if mismatches or missing_defaults:
            proc_cmp.returncode = 1
    except Exception as e:
        proc_cmp.returncode = 2
        proc_cmp.stderr = str(e)

    if proc_cmp.stdout:
        print(proc_cmp.stdout)
    if proc_cmp.stderr:
        print("Comparator stderr:", proc_cmp.stderr)
    if proc_cmp.returncode != 0:
        print("Comparator found differences or failed for", tag, "; return code:", proc_cmp.returncode)

    out_md_repo = input_dir / f"test_roundtrip_{tag}.md"
    try:
        rows = []
        for sec, name, source, expected, found in mismatches:
            key = f"{sec.upper()}.{name.lower()}"
            if source == "input":
                rows.append((key, found, expected, ""))
            else:
                rows.append((key, found, "", expected))
        for sec, name, found in missing_defaults:
            key = f"{sec.upper()}.{name.lower()}"
            rows.append((key, found, "", ""))

        with out_md_repo.open("w") as f:
            f.write("| variable | input_value | namoptions_value | schema_value |\n")
            f.write("|---|---:|---:|---:|\n")
            for var, found, input_val, schema_val in rows:
                f.write(f"| {var} | {input_val} | {found} | {schema_val} |\n")
    except Exception as e:
        print("Warning: could not write Markdown:", e)

    if proc_cmp.returncode != 0:
        overall_ret = proc_cmp.returncode

    if overall_ret != 0:
        return overall_ret

    print("PASS: generated namelist matches input values and schema defaults")
    return 0


if __name__ == "__main__":
    sys.exit(main())
