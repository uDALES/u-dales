"""Standalone namelist round-trip validation test.

This test:
1. prepares an input case
2. runs u-dales in TEST_ROUNDTRIP mode so the last rank writes `namoptions_roundtrip`
3. parses both the input namelist and generated namelist via `f90nml`
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

try:
    import f90nml
except Exception:
    f90nml = None


def _normalize_value(v):
    return v


def _normalize_namelist(nml):
    normalized = {}
    for section, params in nml.items():
        if not isinstance(params, dict):
            continue
        sec_key = str(section).upper()
        normalized.setdefault(sec_key, {})
        for key, val in params.items():
            normalized[sec_key][str(key).lower()] = val
    return normalized


def _normalize_generated_value(expected, actual):
    if isinstance(expected, str) and isinstance(actual, list):
        if all(isinstance(item, (str, int, float)) for item in actual):
            return ''.join(str(item) for item in actual)
    return actual


def _compare_generated(generated, input_cfg, schema, compare_defaults=True):
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
        if section.upper() == "INFO":
            continue
        input_section = input_cfg.get(section, {})
        sec_schema = find_schema_section(section) or {}
        sec_props = sec_schema.get("properties", {}) if isinstance(sec_schema, dict) else {}

        for pname, pval in params.items():
            gen_val = _normalize_value(pval)

            inp_raw = find_input_value(input_section, pname)
            if inp_raw is not None:
                inp_val = _normalize_value(inp_raw)
                gen_val = _normalize_generated_value(inp_val, gen_val)
                try:
                    if isinstance(inp_val, list) and isinstance(gen_val, list):
                        if gen_val[:len(inp_val)] != inp_val:
                            input_mismatches.append((section, pname, "input", inp_val, gen_val))
                    elif isinstance(inp_val, (int, float)) and isinstance(gen_val, (int, float)):
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
                if not compare_defaults:
                    continue
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
    parser.add_argument("--mode", choices=["default", "random"], default="default", help="Test mode: default uses namoptions.default, random generates and uses namoptions.random")
    parser.add_argument("--compare-defaults", action="store_true", help="Also compare generated values against schema defaults")
    args = parser.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[3]

    udales_path = shutil.which("u-dales")
    if f90nml is None:
        print("ERROR: f90nml is required for this test but is not installed")
        return 2
    ud_topdir = os.environ.get("UD_TOPDIR")

    if udales_path:
        udales_exe = Path(udales_path)
    elif ud_topdir:
        udales_exe = Path(ud_topdir) / "bin" / "u-dales"
    else:
        udales_exe = repo_root / "bin" / "u-dales"

    if ud_topdir:
        schema_path = Path(ud_topdir) / "docs" / "schemas" / "udales_input_schema.json"
        input_dir = Path(ud_topdir) / "tests" / "integration" / "input"
    else:
        schema_path = repo_root / "docs" / "schemas" / "udales_input_schema.json"
        input_dir = repo_root / "tests" / "integration" / "input"

    if not schema_path.exists():
        print(f"SKIP: schema not found: {schema_path}")
        return 0

    default_input_nml = input_dir / "namoptions.default"
    random_input_nml = input_dir / "namoptions.random"
    if not default_input_nml.exists():
        print(f"SKIP: default input namelist not found: {default_input_nml}")
        return 0

    tag = args.mode
    selected_input_nml = default_input_nml
    if args.mode == "random":
        gen_script = input_dir / "generate_random_config.py"
        if not gen_script.exists():
            print(f"SKIP: random config generator not found: {gen_script}")
            return 0
        print("Generating namoptions.random using generate_random_config.py")
        proc_gen = subprocess.run(
            [sys.executable, str(gen_script)],
            cwd=str(input_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=60,
        )
        if proc_gen.returncode != 0:
            print("generate_random_config.py failed:")
            print(proc_gen.stdout)
            print(proc_gen.stderr)
            return 2
        selected_input_nml = random_input_nml
        if not selected_input_nml.exists():
            print(f"ERROR: generated random input namelist not found: {selected_input_nml}")
            return 2

    comparator_only = bool(args.comparator_only)

    if not comparator_only:
        if shutil.which("mpirun") is None:
            print("ERROR: mpirun not available")
            return 2
        if not udales_exe.exists():
            print(f"SKIP: u-DALES executable not found: {udales_exe}")
            return 0
    src_case_dir = repo_root / "examples" / "102"
    if not src_case_dir.exists():
        print(f"SKIP: example case directory not found: {src_case_dir}")
        return 0

    tempdir = tempfile.TemporaryDirectory(prefix="udales-tmp2-")
    workdir = Path(tempdir.name)

    with open(schema_path, "r") as f:
        schema = json.load(f)

    overall_ret = 0

    if comparator_only:
        if args.generated:
            generated_namelist = Path(args.generated)
        else:
            generated_namelist = input_dir / "namoptions_roundtrip"
        if not generated_namelist.exists():
            print(f"ERROR: generated round-trip file not found for comparator-only: {generated_namelist}")
            return 2
        input_nml_path = selected_input_nml
        if not input_nml_path.exists():
            print(f"ERROR: comparator-only selected input namelist not found: {input_nml_path}")
            return 2

        work_input = workdir / input_nml_path.name
        shutil.copy2(input_nml_path, work_input)

        source_text = work_input.read_text()
        work_input.write_text(_set_runmode_and_fill_missing_groups(source_text))

        try:
            input_config = _normalize_namelist(f90nml.read(str(work_input)))
        except Exception as e:
            print("ERROR: reading input namelist failed:", e)
            return 2
        try:
            gen = _normalize_namelist(f90nml.read(str(generated_namelist)))
        except Exception as e:
            print("ERROR: reading generated namelist failed:", e)
            return 2
    else:
        _copy_case_directory(src_case_dir, workdir)

        work_input = workdir / selected_input_nml.name
        shutil.copy2(selected_input_nml, work_input)

        source_text = work_input.read_text()
        work_input.write_text(_set_runmode_and_fill_missing_groups(source_text))

        try:
            input_config = _normalize_namelist(f90nml.read(str(work_input)))
        except Exception as e:
            print("ERROR: reading input namelist failed:", e)
            return 2

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

        generated_namelist = workdir / "namoptions_roundtrip"
        if not generated_namelist.exists():
            alt_generated_namelist = workdir / "namoptions_json"
            if alt_generated_namelist.exists():
                generated_namelist = alt_generated_namelist
            else:
                print("ERROR: expected generated namelist not found (namoptions_roundtrip or namoptions_json)")
                return 2

        try:
            gen = _normalize_namelist(f90nml.read(str(generated_namelist)))
        except Exception as e:
            print("ERROR: reading generated namelist failed:", e)
            return 2

    mismatches = []
    missing_defaults = []
    proc_cmp = SimpleNamespace(returncode=0, stdout="", stderr="")
    try:
        mismatches, missing_defaults = _compare_generated(
            gen, input_config, schema, compare_defaults=args.compare_defaults
        )
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
