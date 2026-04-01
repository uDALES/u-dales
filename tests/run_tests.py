#!/usr/bin/env python3

"""Dispatch curated test selections from a manifest."""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional


TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent
MANIFEST_PATH = TESTS_DIR / "test_suites.yml"
DEFAULT_VENV_PYTHON = REPO_ROOT.parent / ".venv" / "bin" / "python"
MPLCONFIGDIR = Path("/tmp") / "udales-matplotlib"


def _load_manifest() -> Dict[str, Any]:
    manifest: Dict[str, Any] = {"groups": {}}
    current_group: Optional[Dict[str, Any]] = None
    current_suite: Optional[Dict[str, Any]] = None
    current_list_key: Optional[str] = None

    with MANIFEST_PATH.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip()
            stripped = line.strip()

            if not stripped or stripped.startswith("#"):
                continue

            indent = len(line) - len(line.lstrip(" "))

            if indent == 0 and stripped.startswith("description:"):
                manifest["description"] = stripped.split(":", 1)[1].strip()
                continue

            if indent == 0 and stripped == "groups:":
                continue

            if indent == 2 and stripped.endswith(":"):
                group_name = stripped[:-1]
                current_group = {"suites": [], "includes": []}
                manifest["groups"][group_name] = current_group
                current_suite = None
                current_list_key = None
                continue

            if current_group is None:
                raise RuntimeError(f"Invalid test manifest structure in {MANIFEST_PATH}: {line}")

            if indent == 4 and stripped.endswith(":"):
                current_list_key = stripped[:-1]
                if current_list_key not in current_group:
                    current_group[current_list_key] = []
                current_suite = None
                continue

            if indent == 6 and stripped.startswith("- "):
                item = stripped[2:]
                if current_list_key == "includes":
                    current_group["includes"].append(item)
                    continue

                if current_list_key == "suites":
                    if ":" not in item:
                        raise RuntimeError(f"Invalid suite entry in {MANIFEST_PATH}: {line}")
                    key, value = item.split(":", 1)
                    current_suite = {key.strip(): value.strip().strip('"'), "command": []}
                    current_group["suites"].append(current_suite)
                    continue

            if indent == 8 and current_suite is not None:
                if stripped == "command:":
                    continue
                if ":" in stripped:
                    key, value = stripped.split(":", 1)
                    current_suite[key.strip()] = value.strip().strip('"')
                    continue

            if indent == 10 and stripped.startswith("- ") and current_suite is not None:
                current_suite["command"].append(stripped[2:].strip().strip('"'))
                continue

            raise RuntimeError(f"Unsupported YAML structure in {MANIFEST_PATH}: {line}")

    if "groups" not in manifest or not isinstance(manifest["groups"], dict):
        raise RuntimeError(f"Invalid test manifest: {MANIFEST_PATH}")
    return manifest


def _expand_groups(
    manifest: Dict[str, Any],
    selection: str,
    seen: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
    if seen is None:
        seen = []
    if selection in seen:
        cycle = " -> ".join(seen + [selection])
        raise RuntimeError(f"Cyclic test group definition in {MANIFEST_PATH}: {cycle}")

    groups = manifest["groups"]
    if selection not in groups:
        raise RuntimeError(f"Unknown test group '{selection}' in {MANIFEST_PATH}")

    group = groups[selection]
    suites: List[Dict[str, Any]] = []

    for include in group.get("includes", []):
        suites.extend(_expand_groups(manifest, include, seen + [selection]))

    for suite in group.get("suites", []):
        suites.append(suite)

    return suites


def _format_command(command: List[str], variables: Dict[str, str]) -> List[str]:
    return [part.format(**variables) for part in command]


def _run_command(
    label: str,
    suite_class: str,
    purpose: str,
    command: List[str],
    env: Dict[str, str],
) -> int:
    print(f"\n==> {label}")
    print(f"class: {suite_class}")
    print(f"purpose: {purpose}")
    print(" ".join(command))
    completed = subprocess.run(command, cwd=REPO_ROOT, env=env)
    status = "PASS" if completed.returncode == 0 else "FAIL"
    print(f"result: {status}")
    if completed.returncode != 0:
        print(f"{label} failed with exit code {completed.returncode}", file=sys.stderr)
    return completed.returncode


def _select_python_interpreter() -> str:
    if DEFAULT_VENV_PYTHON.exists():
        return str(DEFAULT_VENV_PYTHON)
    return sys.executable


def _build_child_env() -> Dict[str, str]:
    MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(MPLCONFIGDIR)
    return env


def main() -> int:
    manifest = _load_manifest()
    groups = manifest["groups"]

    parser = argparse.ArgumentParser(
        description=manifest.get(
            "description",
            "Run curated test selections by support level from tests/test_suites.yml.",
        )
    )
    parser.add_argument(
        "selection",
        choices=sorted(groups.keys()),
        help="Which curated test selection to run.",
    )
    parser.add_argument(
        "--branch-a",
        default="master",
        help="Reference branch for supported regression tests (default: master).",
    )
    parser.add_argument(
        "--branch-b",
        default="HEAD",
        help="Comparison branch or revision for supported regression tests (default: HEAD).",
    )
    parser.add_argument(
        "--build-type",
        choices=["Debug", "Release"],
        default="Release",
        help="Build type for supported regression tests (default: Release).",
    )
    args = parser.parse_args()

    variables = {
        "python": _select_python_interpreter(),
        "repo_root": str(REPO_ROOT),
        "tests_dir": str(TESTS_DIR),
        "branch_a": args.branch_a,
        "branch_b": args.branch_b,
        "build_type": args.build_type,
    }

    suites = _expand_groups(manifest, args.selection)
    child_env = _build_child_env()
    exit_codes = []
    suite_results = []

    for suite in suites:
        label = suite["label"]
        suite_class = suite.get("class", "unspecified")
        purpose = suite.get("purpose", "unspecified")
        command = _format_command(suite["command"], variables)
        code = _run_command(label, suite_class, purpose, command, child_env)
        exit_codes.append(code)
        suite_results.append((label, suite_class, purpose, code))

    print("\nSummary")
    for label, suite_class, purpose, code in suite_results:
        status = "PASS" if code == 0 else "FAIL"
        print(f"- {label} [{suite_class}, {purpose}]: {status}")

    overall = "PASS" if all(code == 0 for code in exit_codes) else "FAIL"
    print(f"overall: {overall}")

    return 0 if overall == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
