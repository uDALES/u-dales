#!/usr/bin/env python3

"""Gate compiler warnings from a build log against a recorded baseline.

Why this exists
---------------
CI builds uDALES with gfortran; most development on CX3 happens with Intel.
The two compilers disagree about what is an error (``use mpi`` gives gfortran an
implicit interface per MPI routine per file, so mixing scalar and rank-1
``MPI_ALLREDUCE`` arguments in one file is a hard error there and silent under
ifort) and about what is a warning. "It builds clean here" is therefore a
statement about one compiler, and the only place that used to be checked was a
push to CI.

This check moves that feedback to ``tests/run_tests.py`` so it is seen before
pushing. It reuses the CI parser -- ``.github/scripts/summarise_warnings.sh
--list`` -- rather than re-implementing it, so the local gate and the CI report
can never drift apart.

Cost
----
Parsing a log is instant; producing one is a full build (minutes). So this does
NOT build. It reads a log a build already wrote, and fails loudly when there is
no usable log rather than passing vacuously. ``tools/build_executable.sh`` and
the CI build step both tee one to ``<build dir>/build.log``.

Baseline, not a ratchet
-----------------------
``tests/lint/build_warnings_baseline.txt`` records the warnings already present
in files nobody in front of this check wrote. A bare "zero warnings" gate would
fail on day one for reasons unrelated to the change under test. Only warnings
ABOVE the baseline for a given (file, class) fail; a count below the baseline is
reported as a suggestion to refresh it, never as a failure. Refresh with
``--update``.

Local gate, CI report
---------------------
This GATES locally and only REPORTS under GitHub Actions. That is deliberate and
matches the policy already written into .github/scripts/summarise_warnings.sh:
CI does not pin its compilers (`apt install gfortran`, `brew install gcc` on
rolling runner images), so a runner image bump legitimately changes the warning
set and would turn CI red on an unrelated PR. Locally the compiler IS pinned --
by the module stack in .github/skills/udales-exec/references/clusters.md -- so a
baseline means something. CI keeps its own reporting step for the same log.
Force report-only anywhere with UDALES_WARNINGS_REPORT_ONLY=1.

Usage
-----
    python tests/lint/check_build_warnings.py [--log PATH] [--update]

Exit status: 0 pass, 1 new warnings, 2 no usable build log.
"""

import argparse
import glob
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

TESTS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = TESTS_DIR.parent
PARSER = REPO_ROOT / ".github" / "scripts" / "summarise_warnings.sh"
BASELINE = Path(__file__).resolve().parent / "build_warnings_baseline.txt"

# Warning flags are only switched on in the Debug configurations (see the
# per-compiler blocks in CMakeLists.txt), so a Release log is silent whatever
# the source says and must never be accepted as evidence.
REQUIRED_BUILD_TYPE = "debug"

# The compiler CI uses, and therefore the one whose warnings gate a merge.
PREFERRED_COMPILER = "GNU"

RECIPE = """
Produce one with a Debug build. To reproduce the CI compiler on CX3:

    module purge && module load tools/prod
    module load foss/2023a netCDF-Fortran/4.6.1-gompi-2023a \\
                FFTW/3.3.10-GCC-12.3.0 CMake/3.26.3-GCCcore-12.3.0
    mkdir -p build/gnu && cd build/gnu
    FC=mpif90 cmake ../.. -DCMAKE_BUILD_TYPE=Debug \\
      -DNETCDF_DIR=$EBROOTNETCDF -DNETCDF_FORTRAN_DIR=$EBROOTNETCDFMINFORTRAN
    make -j8 2>&1 | tee build.log

or, for the Intel stack, `./tools/build_executable.sh icl debug`, which tees
build/debug/build.log for you. Point this check at a specific log with
--log PATH or the UDALES_BUILD_LOG environment variable.
""".strip()


class Meta(object):
    """What a build directory says about itself."""

    def __init__(self, log: Path, compiler: Optional[str], build_type: Optional[str]):
        self.log = log
        self.compiler = compiler
        self.build_type = build_type

    def __repr__(self) -> str:
        return "{} (compiler={}, build type={})".format(
            self.log, self.compiler or "unknown", self.build_type or "unknown"
        )


def _build_dir_meta(log: Path) -> Meta:
    """Read the compiler id and build type CMake recorded next to a log.

    Both come from CMake's own files rather than from the log text: the log of a
    non-verbose `make` does not name the compiler, and guessing it from the
    warning format would be circular (an Intel log has no [-W...] tags, which is
    indistinguishable from a clean gfortran build).
    """
    build_dir = log.parent
    compiler = None
    for path in sorted(build_dir.glob("CMakeFiles/*/CMakeFortranCompiler.cmake")):
        match = re.search(r'set\(CMAKE_Fortran_COMPILER_ID\s+"([^"]+)"', path.read_text())
        if match:
            compiler = match.group(1)
            break

    build_type = None
    cache = build_dir / "CMakeCache.txt"
    if cache.is_file():
        match = re.search(r"^CMAKE_BUILD_TYPE:\w+=(.*)$", cache.read_text(), re.M)
        if match:
            build_type = match.group(1).strip() or None

    return Meta(log, compiler, build_type)


def _candidate_logs(explicit: Optional[str]) -> List[Path]:
    if explicit:
        return [Path(explicit)]
    env = os.environ.get("UDALES_BUILD_LOG")
    if env:
        return [Path(env)]
    found = [Path(p) for p in glob.glob(str(REPO_ROOT / "build" / "*" / "build.log"))]
    return sorted(found)


def _pick_log(explicit: Optional[str], baselined: List[str]) -> Tuple[Optional[Meta], List[str]]:
    """Return the log to check, plus the reasons every rejected candidate lost."""
    notes = []
    usable = []
    for log in _candidate_logs(explicit):
        if not log.is_file():
            notes.append("{}: no such file".format(log))
            continue
        meta = _build_dir_meta(log)
        if meta.compiler is None:
            notes.append("{}: cannot tell which compiler wrote it "
                         "(no CMakeFiles/*/CMakeFortranCompiler.cmake)".format(log))
            continue
        if (meta.build_type or "").lower() != REQUIRED_BUILD_TYPE:
            notes.append("{}: build type is {}, and warning flags are only on in "
                         "Debug".format(log, meta.build_type or "unset"))
            continue
        if meta.compiler not in baselined:
            notes.append("{}: no baseline section for compiler {} in {}".format(
                log, meta.compiler, BASELINE.name))
            continue
        usable.append(meta)

    if not usable:
        return None, notes

    # Prefer the compiler CI gates on; among equals, the most recent build.
    usable.sort(key=lambda m: (m.compiler != PREFERRED_COMPILER,
                               -m.log.stat().st_mtime))
    for meta in usable[1:]:
        notes.append("{}: not chosen (using {} instead)".format(meta.log, usable[0].log))
    return usable[0], notes


def _newest_source_mtime() -> Tuple[float, Optional[Path]]:
    newest, which = 0.0, None
    paths = list((REPO_ROOT / "src").rglob("*.f90")) + [REPO_ROOT / "CMakeLists.txt"]
    for path in paths:
        if not path.is_file():
            continue
        mtime = path.stat().st_mtime
        if mtime > newest:
            newest, which = mtime, path
    return newest, which


def _log_is_complete(log: Path) -> bool:
    """True when the log covers every src/*.f90, i.e. it is a full build.

    An incremental build writes a log naming only the files it recompiled. That
    is still worth checking -- those are the files just edited -- but it is not
    evidence about the rest of the tree, so the "your baseline is stale" hint is
    only offered for a complete log.
    """
    text = log.read_text(errors="replace")
    compiled = set(re.findall(r"([A-Za-z0-9_./-]+\.f90)\.o", text))
    compiled = set(Path(p).name for p in compiled)
    sources = set(p.name for p in (REPO_ROOT / "src").rglob("*.f90"))
    return sources.issubset(compiled)


def _parse(log: Path) -> List[Tuple[str, str, str, str]]:
    """(file, line, class, message) for every warning, via the CI parser."""
    completed = subprocess.run(
        ["bash", str(PARSER), "--list", str(log)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
    )
    if completed.returncode != 0:
        raise RuntimeError("{} --list failed: {}".format(PARSER, completed.stderr.strip()))
    records = []
    for line in completed.stdout.splitlines():
        parts = line.split("|", 3)
        if len(parts) == 4:
            records.append((parts[0], parts[1], parts[2], parts[3]))
    return records


def _read_baseline() -> Tuple[List[str], Dict[Tuple[str, str, str], int]]:
    compilers, counts = [], {}
    if not BASELINE.is_file():
        raise RuntimeError("missing baseline file {}".format(BASELINE))
    for raw in BASELINE.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        fields = line.split()
        if fields[0] == "compiler" and len(fields) == 2:
            compilers.append(fields[1])
            continue
        if len(fields) != 4:
            raise RuntimeError("malformed baseline line in {}: {!r}".format(BASELINE, raw))
        compiler, source, warning_class, count = fields
        counts[(compiler, source, warning_class)] = int(count)
    if not compilers:
        raise RuntimeError("{} declares no `compiler <ID>` lines".format(BASELINE))
    return compilers, counts


def _write_baseline(compilers: List[str], counts: Dict[Tuple[str, str, str], int]) -> None:
    lines = [
        "# Compiler warnings already present in the tree, per compiler.",
        "#",
        "# Consumed by tests/lint/check_build_warnings.py. Only counts ABOVE these",
        "# fail; this is a baseline, not a ratchet, so an unrelated change is never",
        "# blocked by warnings it did not introduce.",
        "#",
        "# Format:",
        "#   compiler <CMake Fortran compiler id>   -- a compiler this file covers;",
        "#                                            a Debug log from any other",
        "#                                            compiler is rejected, not",
        "#                                            silently passed",
        "#   <compiler> <source file> <-Wclass> <count>",
        "#",
        "# To refresh after legitimately adding or removing a warning:",
        "#   1. full Debug build with that compiler, teeing to <build dir>/build.log",
        "#   2. python tests/lint/check_build_warnings.py --update --log <that log>",
        "#   3. commit the diff, and say in the message why the new entries are",
        "#      acceptable -- an entry added here is a warning nobody will be told",
        "#      about again.",
        "#",
        "# Counts are per (file, class) rather than per line so that editing a file",
        "# does not invalidate the baseline for every warning below the edit.",
        "#",
        "# Measured against the pinned local module stacks documented in",
        "# .github/skills/udales-exec/references/clusters.md (GNU: foss/2023a,",
        "# gfortran 12.3; Intel: intel/2021a via tools/build_executable.sh icl).",
        "# GitHub Actions does NOT pin its compilers, so the gate is report-only",
        "# there -- see the header of check_build_warnings.py.",
        "",
    ]
    for compiler in sorted(compilers):
        lines.append("compiler {}".format(compiler))
    lines.append("")
    for (compiler, source, warning_class), count in sorted(counts.items()):
        lines.append("{} {} {} {}".format(compiler, source, warning_class, count))
    BASELINE.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--log", default=None,
                        help="Build log to check (default: UDALES_BUILD_LOG, else "
                             "the newest usable build/*/build.log).")
    parser.add_argument("--update", action="store_true",
                        help="Rewrite the baseline from this log instead of checking "
                             "against it. Requires a complete (full-build) log.")
    parser.add_argument("--report-only", action="store_true",
                        help="Print the comparison but always exit 0. Implied under "
                             "GitHub Actions and by UDALES_WARNINGS_REPORT_ONLY=1.")
    args = parser.parse_args()

    report_only = (args.report_only
                   or os.environ.get("UDALES_WARNINGS_REPORT_ONLY") == "1"
                   or os.environ.get("GITHUB_ACTIONS") == "true")

    print("==> compiler warning gate (tests/lint/check_build_warnings.py)")
    if report_only:
        print("  mode:       REPORT ONLY (exit 0 regardless) -- see this file's header "
              "for why CI does not gate on warnings")

    baselined, baseline = _read_baseline()
    meta, notes = _pick_log(args.log, baselined)

    # Report-only mode still prints everything; only the verdict is softened, so
    # the wording must not claim a failure that did not happen.
    tag = "NOTE (report only)" if report_only else "FAIL"

    if meta is None:
        sys.stdout.flush()
        print("{}: no usable build log.".format(tag), file=sys.stderr)
        print("", file=sys.stderr)
        for note in notes:
            print("  rejected {}".format(note), file=sys.stderr)
        if not notes:
            print("  no build/*/build.log found under {}".format(REPO_ROOT), file=sys.stderr)
        print("", file=sys.stderr)
        print(RECIPE, file=sys.stderr)
        print("", file=sys.stderr)
        if not report_only:
            print("This check deliberately fails rather than passing when it did "
                  "not run: a check that is silent when it saw nothing is worse "
                  "than no check.", file=sys.stderr)
        return 0 if report_only else 2

    for note in notes:
        print("  note: rejected {}".format(note))

    newest_src, which_src = _newest_source_mtime()
    log_mtime = meta.log.stat().st_mtime
    complete = _log_is_complete(meta.log)

    print("  log:        {}".format(meta.log))
    print("  compiler:   {}{}".format(
        meta.compiler,
        "" if meta.compiler == PREFERRED_COMPILER
        else "  (CI gates on {}; this log cannot see its warnings)".format(PREFERRED_COMPILER)))
    print("  build type: {}".format(meta.build_type))
    print("  coverage:   {}".format(
        "full build" if complete else "partial (incremental build log)"))

    if log_mtime < newest_src:
        sys.stdout.flush()
        print("{}: the build log is older than the sources it should "
              "describe.".format(tag), file=sys.stderr)
        print("  {} was modified after {} was written.".format(which_src, meta.log),
              file=sys.stderr)
        print("  Rebuild before running this check.", file=sys.stderr)
        return 0 if report_only else 2

    records = _parse(meta.log)
    counts = Counter((meta.compiler, rec[0], rec[2]) for rec in records)

    if args.update:
        if not complete:
            sys.stdout.flush()
            print("FAIL: --update needs a complete (full-build) log; {} covers only "
                  "the files it recompiled.".format(meta.log), file=sys.stderr)
            return 2
        merged = dict((k, v) for k, v in baseline.items() if k[0] != meta.compiler)
        merged.update(counts)
        compilers = sorted(set(baselined) | set([meta.compiler]))
        _write_baseline(compilers, merged)
        print("  updated {} for compiler {} ({} entries)".format(
            BASELINE, meta.compiler, len(counts)))
        return 0

    regressions = []
    for key, count in sorted(counts.items()):
        allowed = baseline.get(key, 0)
        if count > allowed:
            regressions.append((key, count, allowed))

    total = sum(counts.values())
    allowed_total = sum(v for k, v in baseline.items() if k[0] == meta.compiler)
    print("  warnings:   {} in this log, baseline allows {} for {}".format(
        total, allowed_total, meta.compiler))

    if regressions:
        print("")
        sys.stdout.flush()
        print("{}: warnings above the baseline.".format(tag), file=sys.stderr)
        print("", file=sys.stderr)
        for (compiler, source, warning_class), count, allowed in regressions:
            print("  {}: {} x {} (baseline {})".format(
                source, count, warning_class, allowed), file=sys.stderr)
            for rec in records:
                if rec[0] == source and rec[2] == warning_class:
                    print("      {}:{}: {}".format(rec[0], rec[1], rec[3]), file=sys.stderr)
        print("", file=sys.stderr)
        print("Fix them, or -- if they are correct as they stand -- record them:", file=sys.stderr)
        print("  python tests/lint/check_build_warnings.py --update --log {}".format(
            meta.log), file=sys.stderr)
        print("and say in the commit message why. See the header of {}.".format(
            BASELINE.name), file=sys.stderr)
        return 0 if report_only else 1

    if complete:
        stale = [(k, v) for k, v in sorted(baseline.items())
                 if k[0] == meta.compiler and counts.get(k, 0) < v]
        if stale:
            print("")
            print("  NOTE: {} baseline entries are now over-generous (warnings were "
                  "fixed).".format(len(stale)))
            for (compiler, source, warning_class), allowed in stale:
                print("        {} {}: {} now, {} allowed".format(
                    source, warning_class, counts.get((compiler, source, warning_class), 0),
                    allowed))
            print("        Refresh with --update when convenient. Not a failure.")

    print("PASS: no warnings above the baseline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
