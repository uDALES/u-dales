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

Per compiler major version
--------------------------
A warning set is a property of a compiler *version*: gfortran 12 and 13 do not
report the same things, and a runner image bump would turn CI red on an
unrelated PR if one baseline were applied to every version. So the baseline is
recorded per ``(compiler id, major version)``, and this check

* GATES -- fails on a warning above the baseline -- when the log's compiler
  major has a recorded section, whether locally or under GitHub Actions;
* is REPORT-ONLY, and says so, when it has not: the comparison is printed
  against nothing and the exit status is 0.

The compiler id and version are read from the ``CMakeFortranCompiler.cmake``
CMake writes next to the log. Recorded as of 2026-09: gfortran 12 (CX3,
``foss/2023a``), gfortran 13 (``ubuntu-latest``, 13.2.0) and gfortran 16
(``macos-latest``, Homebrew GCC 16.2.0), so both CI legs gate until their
image moves to a major this file has not seen -- at which point the report
names the missing major and someone records it with ``--update``.

Force report-only anywhere with ``UDALES_WARNINGS_REPORT_ONLY=1``.

Under GitHub Actions a Release leg has no Debug log to check (the warning
flags only exist in Debug); that is reported as not applicable and exits 0.
Any other reason for having no usable log is still a failure there too.

Usage
-----
    python tests/lint/check_build_warnings.py [--log PATH] [--update]
        [--compiler ID:VERSION]

``--compiler`` names the compiler when the log did not come from a local build
directory -- a log harvested from a CI job with ``gh run view <id> --log``,
say -- and is only meaningful with a Debug log.

Exit status: 0 pass (or report-only), 1 new warnings, 2 no usable build log.
"""

import argparse
import glob
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

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

#: Compilers whose warnings summarise_warnings.sh can actually parse. Its awk
#: rules key on a trailing "[-Wclass]" tag, which is a gfortran/clang format;
#: ifort emits "warning #NNNN" and parses to nothing. A log from anything else
#: is not evidence of zero warnings, it is an absence of evidence, so treat it
#: as unusable rather than passing on it.
PARSEABLE_COMPILERS = ("GNU",)

#: (compiler id, major) -- one section of the baseline.
Section = Tuple[str, str]

RECIPE = """
Produce one with a Debug build. To reproduce the CI compiler on CX3:

    module purge && module load tools/prod
    module load foss/2023a netCDF-Fortran/4.6.1-gompi-2023a \\
                FFTW/3.3.10-GCC-12.3.0 CMake/3.26.3-GCCcore-12.3.0
    mkdir -p build/gnu && cd build/gnu
    FC=mpif90 cmake ../.. -DCMAKE_BUILD_TYPE=Debug \\
      -DNETCDF_DIR=$EBROOTNETCDF -DNETCDF_FORTRAN_DIR=$EBROOTNETCDFMINFORTRAN
    make -j8 2>&1 | tee build.log

Point this check at a specific log with --log PATH or the UDALES_BUILD_LOG
environment variable.

An Intel log will not do, even though `./tools/build_executable.sh icl debug`
tees one to build/debug/build.log. The parser keys on gfortran's "[-Wclass]"
tag; ifort's "warning #NNNN" parses to nothing, which would make an Intel log
look like a clean build rather than an unreadable one.
""".strip()


class Meta(object):
    """What a build directory says about itself."""

    def __init__(self, log: Path, compiler: Optional[str], version: Optional[str],
                 build_type: Optional[str]):
        self.log = log
        self.compiler = compiler
        self.version = version
        self.build_type = build_type

    @property
    def major(self) -> Optional[str]:
        if not self.version:
            return None
        return self.version.split(".")[0]

    @property
    def section(self) -> Section:
        return (self.compiler or "unknown", self.major or "unknown")

    def __repr__(self) -> str:
        return "{} (compiler={} {}, build type={})".format(
            self.log, self.compiler or "unknown", self.version or "?",
            self.build_type or "unknown")


def _build_dir_meta(log: Path, override: Optional[Tuple[str, str]]) -> Meta:
    """Read the compiler id, its version and the build type CMake recorded next to a log.

    All three come from CMake's own files rather than from the log text: the
    log of a non-verbose `make` does not name the compiler, and guessing it
    from the warning format would be circular (an Intel log has no [-W...]
    tags, which is indistinguishable from a clean gfortran build).

    ``override`` -- from ``--compiler ID:VERSION`` -- stands in for a log that
    has no build directory, and asserts it is a Debug log.
    """
    build_dir = log.parent
    compiler, version = None, None
    for path in sorted(build_dir.glob("CMakeFiles/*/CMakeFortranCompiler.cmake")):
        text = path.read_text()
        match = re.search(r'set\(CMAKE_Fortran_COMPILER_ID\s+"([^"]+)"', text)
        if match:
            compiler = match.group(1)
        match = re.search(r'set\(CMAKE_Fortran_COMPILER_VERSION\s+"([^"]+)"', text)
        if match:
            version = match.group(1)
        if compiler:
            break

    build_type = None
    cache = build_dir / "CMakeCache.txt"
    if cache.is_file():
        match = re.search(r"^CMAKE_BUILD_TYPE:\w+=(.*)$", cache.read_text(), re.M)
        if match:
            build_type = match.group(1).strip() or None

    if override is not None and compiler is None:
        compiler, version = override
        build_type = build_type or "Debug"

    return Meta(log, compiler, version, build_type)


def _candidate_logs(explicit: Optional[str]) -> List[Path]:
    if explicit:
        return [Path(explicit)]
    env = os.environ.get("UDALES_BUILD_LOG")
    if env:
        return [Path(env)]
    found = [Path(p) for p in glob.glob(str(REPO_ROOT / "build" / "*" / "build.log"))]
    return sorted(found)


def _pick_log(explicit: Optional[str], override: Optional[Tuple[str, str]]
              ) -> Tuple[Optional[Meta], List[str], bool]:
    """The log to check, why every rejected candidate lost, and whether the
    only reason for having none is that they were all Release logs."""
    notes = []
    usable = []
    only_release = True
    for log in _candidate_logs(explicit):
        if not log.is_file():
            notes.append("{}: no such file".format(log))
            only_release = False
            continue
        meta = _build_dir_meta(log, override)
        if meta.compiler is None:
            notes.append("{}: cannot tell which compiler wrote it "
                         "(no CMakeFiles/*/CMakeFortranCompiler.cmake; --compiler "
                         "ID:VERSION says so for a log without a build directory)".format(log))
            only_release = False
            continue
        if (meta.build_type or "").lower() != REQUIRED_BUILD_TYPE:
            notes.append("{}: build type is {}, and warning flags are only on in "
                         "Debug".format(log, meta.build_type or "unset"))
            continue
        if meta.compiler not in PARSEABLE_COMPILERS:
            notes.append("{}: {} warnings are not in a format this check can parse, "
                         "so the log cannot show absence of warnings".format(
                             log, meta.compiler))
            only_release = False
            continue
        if meta.version is None:
            notes.append("{}: CMake recorded no compiler version, so the baseline "
                         "section cannot be chosen".format(log))
            only_release = False
            continue
        usable.append(meta)

    if not usable:
        return None, notes, only_release and bool(notes)

    # Prefer the compiler CI gates on; among equals, the most recent build.
    usable.sort(key=lambda m: (m.compiler != PREFERRED_COMPILER,
                               -m.log.stat().st_mtime))
    for meta in usable[1:]:
        notes.append("{}: not chosen (using {} instead)".format(meta.log, usable[0].log))
    return usable[0], notes, False


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


def _read_baseline() -> Tuple[Set[Section], Dict[Tuple[str, str, str, str], int]]:
    """The recorded (compiler, major) sections and the counts under them."""
    sections: Set[Section] = set()
    counts: Dict[Tuple[str, str, str, str], int] = {}
    if not BASELINE.is_file():
        raise RuntimeError("missing baseline file {}".format(BASELINE))
    for raw in BASELINE.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        fields = line.split()
        if fields[0] == "compiler" and len(fields) == 3:
            sections.add((fields[1], fields[2]))
            continue
        if len(fields) != 5:
            raise RuntimeError("malformed baseline line in {}: {!r} (expected "
                               "`compiler <ID> <major>` or `<ID> <major> <file> "
                               "<-Wclass> <count>`)".format(BASELINE, raw))
        compiler, major, source, warning_class, count = fields
        if (compiler, major) not in sections:
            raise RuntimeError("{}: entry {!r} precedes its `compiler {} {}` line".format(
                BASELINE, raw, compiler, major))
        counts[(compiler, major, source, warning_class)] = int(count)
    if not sections:
        raise RuntimeError("{} declares no `compiler <ID> <major>` lines".format(BASELINE))
    return sections, counts


def _write_baseline(sections: Set[Section],
                    counts: Dict[Tuple[str, str, str, str], int]) -> None:
    lines = [
        "# Compiler warnings already present in the tree, per compiler major version.",
        "#",
        "# Consumed by tests/lint/check_build_warnings.py. Only counts ABOVE these",
        "# fail; this is a baseline, not a ratchet, so an unrelated change is never",
        "# blocked by warnings it did not introduce.",
        "#",
        "# Format:",
        "#   compiler <CMake Fortran compiler id> <major version>",
        "#       -- a compiler version this file covers. A Debug log from a compiler",
        "#          major NOT listed here is checked in report-only mode (the check",
        "#          says so and exits 0); a log from a listed one is gated.",
        "#   <compiler> <major> <source file> <-Wclass> <count>",
        "#",
        "# To refresh after legitimately adding or removing a warning:",
        "#   1. full Debug build with that compiler, teeing to <build dir>/build.log",
        "#      (or harvest a CI job's log with `gh run view <id> --log`)",
        "#   2. python tests/lint/check_build_warnings.py --update --log <that log>",
        "#      (add --compiler GNU:<version> for a log with no build directory)",
        "#   3. commit the diff, and say in the message why the new entries are",
        "#      acceptable -- an entry added here is a warning nobody will be told",
        "#      about again.",
        "#",
        "# Counts are per (file, class) rather than per line so that editing a file",
        "# does not invalidate the baseline for every warning below the edit.",
        "#",
        "# Where each section was measured (2026-09):",
        "#   GNU 12 -- CX3, foss/2023a, gfortran 12.3.0 (the local recipe in",
        "#             .github/skills/udales-exec/references/clusters.md)",
        "#   GNU 13 -- GitHub ubuntu-latest, gfortran 13.2.0 (4:13.2.0-7ubuntu1)",
        "#   GNU 16 -- GitHub macos-latest, Homebrew GCC 16.2.0",
        "# CI does not pin its compilers, so when a runner image moves to a major",
        "# not listed here the gate turns report-only there and names the gap.",
        "",
    ]
    for compiler, major in sorted(sections, key=lambda s: (s[0], int(s[1]) if s[1].isdigit() else 0)):
        lines.append("compiler {} {}".format(compiler, major))
    lines.append("")
    for (compiler, major, source, warning_class), count in sorted(
            counts.items(), key=lambda kv: (kv[0][0], int(kv[0][1]) if kv[0][1].isdigit() else 0,
                                            kv[0][2], kv[0][3])):
        lines.append("{} {} {} {} {}".format(compiler, major, source, warning_class, count))
    BASELINE.write_text("\n".join(lines) + "\n")


def _parse_compiler_override(value: Optional[str]) -> Optional[Tuple[str, str]]:
    if not value:
        return None
    if ":" not in value:
        raise SystemExit("--compiler expects ID:VERSION, e.g. GNU:13.2.0")
    compiler, version = value.split(":", 1)
    return compiler.strip(), version.strip()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--log", default=None,
                        help="Build log to check (default: UDALES_BUILD_LOG, else "
                             "the newest usable build/*/build.log).")
    parser.add_argument("--update", action="store_true",
                        help="Rewrite this compiler major's baseline section from the "
                             "log instead of checking against it. Requires a complete "
                             "(full-build) log.")
    parser.add_argument("--report-only", action="store_true",
                        help="Print the comparison but always exit 0. Also implied by "
                             "UDALES_WARNINGS_REPORT_ONLY=1, and by a log from a "
                             "compiler major the baseline does not record.")
    parser.add_argument("--compiler", default=None, metavar="ID:VERSION",
                        help="Compiler id and version of a log that has no CMake build "
                             "directory next to it (e.g. one harvested from CI); such a "
                             "log is taken to be a Debug log.")
    args = parser.parse_args()

    override = _parse_compiler_override(args.compiler)
    forced = args.report_only or os.environ.get("UDALES_WARNINGS_REPORT_ONLY") == "1"
    in_ci = os.environ.get("GITHUB_ACTIONS") == "true"

    print("==> compiler warning gate (tests/lint/check_build_warnings.py)")

    sections, baseline = _read_baseline()
    recorded = ", ".join("{} {}".format(c, m) for c, m in sorted(sections))
    meta, notes, only_release = _pick_log(args.log, override)

    if meta is None:
        sys.stdout.flush()
        if in_ci and only_release:
            print("  NOT APPLICABLE: only a Release log is available and the warning "
                  "flags exist only in Debug; the Debug leg of this matrix gates.")
            for note in notes:
                print("  rejected {}".format(note))
            return 0
        tag = "NOTE (report only)" if forced else "FAIL"
        print("{}: no usable build log.".format(tag), file=sys.stderr)
        print("", file=sys.stderr)
        for note in notes:
            print("  rejected {}".format(note), file=sys.stderr)
        if not notes:
            print("  no build/*/build.log found under {}".format(REPO_ROOT), file=sys.stderr)
        print("", file=sys.stderr)
        print(RECIPE, file=sys.stderr)
        print("", file=sys.stderr)
        if not forced:
            print("This check deliberately fails rather than passing when it did "
                  "not run: a check that is silent when it saw nothing is worse "
                  "than no check.", file=sys.stderr)
        return 0 if forced else 2

    known = meta.section in sections
    report_only = forced or not known
    if forced:
        mode = "REPORT ONLY (forced by --report-only / UDALES_WARNINGS_REPORT_ONLY=1)"
    elif not known:
        mode = ("REPORT ONLY: no baseline recorded for {} {} (recorded: {}). Record one "
                "with --update if this compiler should gate.".format(
                    meta.compiler, meta.major, recorded))
    else:
        mode = "GATE: baseline recorded for {} {}{}".format(
            meta.compiler, meta.major, " (GitHub Actions)" if in_ci else "")
    print("  mode:       {}".format(mode))

    # Report-only mode still prints everything; only the verdict is softened, so
    # the wording must not claim a failure that did not happen.
    tag = "NOTE (report only)" if report_only else "FAIL"

    for note in notes:
        print("  note: rejected {}".format(note))

    newest_src, which_src = _newest_source_mtime()
    log_mtime = meta.log.stat().st_mtime
    complete = _log_is_complete(meta.log)

    print("  log:        {}".format(meta.log))
    print("  compiler:   {} {} (major {})".format(meta.compiler, meta.version, meta.major))
    print("  build type: {}".format(meta.build_type))
    print("  coverage:   {}".format(
        "full build" if complete else "partial (incremental build log)"))

    if log_mtime < newest_src and override is None:
        sys.stdout.flush()
        print("{}: the build log is older than the sources it should "
              "describe.".format(tag), file=sys.stderr)
        print("  {} was modified after {} was written.".format(which_src, meta.log),
              file=sys.stderr)
        print("  Rebuild before running this check.", file=sys.stderr)
        return 0 if report_only else 2

    records = _parse(meta.log)
    counts = Counter((meta.compiler, meta.major, rec[0], rec[2]) for rec in records)

    if args.update:
        if not complete:
            sys.stdout.flush()
            print("FAIL: --update needs a complete (full-build) log; {} covers only "
                  "the files it recompiled.".format(meta.log), file=sys.stderr)
            return 2
        merged = dict((k, v) for k, v in baseline.items() if k[:2] != meta.section)
        merged.update(counts)
        _write_baseline(sections | {meta.section}, merged)
        print("  updated {} for {} {} ({} entries)".format(
            BASELINE, meta.compiler, meta.major, len(counts)))
        return 0

    regressions = []
    for key, count in sorted(counts.items()):
        allowed = baseline.get(key, 0)
        if count > allowed:
            regressions.append((key, count, allowed))

    total = sum(counts.values())
    allowed_total = sum(v for k, v in baseline.items() if k[:2] == meta.section)
    print("  warnings:   {} in this log, baseline allows {} for {} {}".format(
        total, allowed_total, meta.compiler, meta.major))

    if regressions:
        print("")
        sys.stdout.flush()
        print("{}: warnings above the baseline.".format(tag), file=sys.stderr)
        print("", file=sys.stderr)
        for (compiler, major, source, warning_class), count, allowed in regressions:
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

    if complete and known:
        stale = [(k, v) for k, v in sorted(baseline.items())
                 if k[:2] == meta.section and counts.get(k, 0) < v]
        if stale:
            print("")
            print("  NOTE: {} baseline entries are now over-generous (warnings were "
                  "fixed).".format(len(stale)))
            for (compiler, major, source, warning_class), allowed in stale:
                print("        {} {}: {} now, {} allowed".format(
                    source, warning_class,
                    counts.get((compiler, major, source, warning_class), 0), allowed))
            print("        Refresh with --update when convenient. Not a failure.")

    print("PASS: no warnings above the baseline." if known else
          "PASS (report only): nothing was gated, see the mode line above.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
