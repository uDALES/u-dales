"""Build the pre-branch ("baseline") solver the I1 no-op test compares against.

I1 asks whether ``lnesting = .false.`` leaves the solver bit-identical to the
code before the branch.  That needs a second binary, built from the base
branch with the **same compiler, build type and library paths** as the binary
under test -- a Release baseline against a Debug branch build measures the
optimisation level, not the branch.  Until now the baseline was an uncommitted
``build/u-dales.baseline`` somebody had to produce by hand, which is why I1
only ever ran on one machine.

This module builds it:

1. ``UDALES_BASELINE_REF`` (default ``origin/master``) is resolved to a commit.
   ``git fetch origin`` first; a stale local ``master`` has burned this project
   before, and CI only has ``origin/master``.
2. A detached ``git worktree`` of that commit is created under ``build/``
   (with the ``2decomp-fft`` submodule), so nothing is checked out in place --
   unlike ``tests/regression/david_tests``, which switches the working tree.
3. It is configured with the compiler, build type, NetCDF/FFTW locations and
   MPI settings read back from the ``CMakeCache.txt`` next to ``UDALES_BUILD``,
   plus ``UDALES_CMAKE_ARGS`` (which ``.github/scripts/setup_mpi_env.sh``
   exports in CI), and built.
4. A stamp records the commit and the configuration, so a second call with
   nothing changed is free.

``UDALES_BASELINE`` still short-circuits all of this: point it at a binary and
that is what I1 uses.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shlex
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import _launch as launch

REPO_ROOT = launch.REPO_ROOT
DEFAULT_REF = "origin/master"

#: CMake cache entries copied from the build under test to the baseline build.
#: Anything that selects a compiler, a library or a build type; nothing that
#: CMake derives on its own.
CACHE_KEYS = (
    "CMAKE_BUILD_TYPE",
    "CMAKE_Fortran_COMPILER",
    "CMAKE_C_COMPILER",
    "MPI_Fortran_COMPILER",
    "MPI_C_COMPILER",
    "MPIEXEC_EXECUTABLE",
    "NETCDF_DIR",
    "NETCDF_FORTRAN_DIR",
    "FFTW_ROOT",
    "FFTW_INCLUDE_DIRS",
    "FFTW_DOUBLE_LIB",
    "FFTW_FLOAT_LIB",
    "CMAKE_POLICY_VERSION_MINIMUM",
)


class BaselineError(RuntimeError):
    pass


def _git(*args: str, cwd: Path = REPO_ROOT) -> str:
    done = subprocess.run(["git", *args], cwd=cwd, check=False,
                          capture_output=True, text=True)
    if done.returncode != 0:
        raise BaselineError(f"git {' '.join(args)} failed in {cwd}:\n{done.stderr.strip()}")
    return done.stdout.strip()


def resolve_ref(ref: str) -> str:
    """The commit ``ref`` names, after a fetch; ``origin/<ref>`` as a fallback.

    CI's checkout has no local ``master``, only ``origin/master``; a developer
    may have a local ``master`` that is weeks stale.  Fetching first and
    preferring what was asked for, then the remote-tracking name, covers both.
    """
    subprocess.run(["git", "fetch", "origin"], cwd=REPO_ROOT, check=False,
                   capture_output=True, text=True)
    for candidate in (ref, f"origin/{ref}"):
        done = subprocess.run(["git", "rev-parse", "--verify", "--quiet", f"{candidate}^{{commit}}"],
                              cwd=REPO_ROOT, check=False, capture_output=True, text=True)
        if done.returncode == 0:
            return done.stdout.strip()
    raise BaselineError(f"cannot resolve baseline ref {ref!r} (nor origin/{ref})")


def cache_settings(build_exe: Path) -> Dict[str, str]:
    """The compiler/library choices recorded next to the binary under test."""
    cache = build_exe.parent / "CMakeCache.txt"
    if not cache.is_file():
        raise BaselineError(f"no CMakeCache.txt next to {build_exe}; the baseline must be "
                            "configured like the build under test, and that is where "
                            "the configuration is read from")
    found: Dict[str, str] = {}
    for line in cache.read_text(errors="replace").splitlines():
        m = re.match(r"^([A-Za-z0-9_]+):[A-Z]+=(.*)$", line)
        if m and m.group(1) in CACHE_KEYS and m.group(2).strip():
            found[m.group(1)] = m.group(2).strip()
    if "CMAKE_Fortran_COMPILER" not in found or "CMAKE_BUILD_TYPE" not in found:
        raise BaselineError(f"{cache} does not record a Fortran compiler and a build type")
    return found


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def ensure_baseline(build_exe: Path, ref: Optional[str] = None,
                    jobs: Optional[int] = None) -> Path:
    """Return a baseline executable configured like ``build_exe``, building it if needed."""
    explicit = os.environ.get("UDALES_BASELINE", "").strip()
    if explicit:
        path = Path(explicit)
        if not path.is_file():
            raise BaselineError(f"UDALES_BASELINE={explicit} is not a file")
        print(f"[baseline] using UDALES_BASELINE={path}", flush=True)
        return path

    ref = ref or os.environ.get("UDALES_BASELINE_REF", DEFAULT_REF)
    commit = resolve_ref(ref)
    settings = cache_settings(build_exe)
    build_type = settings["CMAKE_BUILD_TYPE"]
    # One build directory per (ref, compiler, build type): a CMake tree
    # configured for one compiler cannot be re-configured for another in
    # place, and an Intel Debug baseline next to a gfortran Debug one is the
    # normal case on CX3.
    compiler = Path(settings["CMAKE_Fortran_COMPILER"]).name
    tag = f"{_safe(ref)}-{_safe(compiler)}-{build_type.lower()}"
    src_dir = REPO_ROOT / "build" / f"baseline-src-{_safe(ref)}"
    build_dir = REPO_ROOT / "build" / f"baseline-{tag}"
    exe = build_dir / "u-dales"
    stamp = build_dir / "baseline.stamp"

    extra = os.environ.get("UDALES_CMAKE_ARGS", "").strip()
    config = {"commit": commit, "settings": settings, "extra": extra,
              "modules": launch.RUNTIME_MODULES}
    digest = hashlib.sha1(json.dumps(config, sort_keys=True).encode()).hexdigest()
    if exe.is_file() and stamp.is_file() and stamp.read_text().strip() == digest:
        print(f"[baseline] reusing {exe} ({ref} @ {commit[:10]}, {build_type})", flush=True)
        return exe

    _checkout(src_dir, commit)
    _configure_and_build(src_dir, build_dir, settings, extra, jobs)
    if not exe.is_file():
        raise BaselineError(f"the baseline build finished but {exe} does not exist")
    stamp.write_text(digest + "\n")
    print(f"[baseline] built {exe} from {ref} @ {commit[:10]} ({build_type})", flush=True)
    return exe


def _checkout(src_dir: Path, commit: str) -> None:
    if src_dir.is_dir():
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=src_dir, check=False,
                              capture_output=True, text=True).stdout.strip()
        if head != commit:
            _git("worktree", "remove", "--force", str(src_dir))
    if not src_dir.is_dir():
        # A previous worktree at this path may have been deleted by hand.
        subprocess.run(["git", "worktree", "prune"], cwd=REPO_ROOT, check=False,
                       capture_output=True, text=True)
        _git("worktree", "add", "--detach", str(src_dir), commit)
    # The main checkout's copy of the submodule is a valid --reference, which
    # saves the clone on a machine without network access to GitHub.
    reference = REPO_ROOT / "2decomp-fft" / ".git"
    args = ["submodule", "update", "--init", "--recursive"]
    if reference.exists():
        args += ["--reference", str(REPO_ROOT / "2decomp-fft")]
    _git(*args, "--", "2decomp-fft", cwd=src_dir)
    if not (src_dir / "2decomp-fft" / "CMakeLists.txt").is_file():
        raise BaselineError(f"2decomp-fft did not initialise in {src_dir}")


def _configure_and_build(src_dir: Path, build_dir: Path, settings: Dict[str, str],
                         extra: str, jobs: Optional[int]) -> None:
    build_dir.mkdir(parents=True, exist_ok=True)
    defines: List[str] = [f"-D{key}={value}" for key, value in sorted(settings.items())]
    if extra:
        defines += shlex.split(extra)
    jobs = jobs or min(8, os.cpu_count() or 2)
    log = build_dir / "build.log"
    # Through the same shell prefix the solver runs under, so the module stack
    # that provides the compiler is loaded; on a runner it is a no-op.
    command = (
        f"{launch.shell_prefix()}cd {shlex.quote(str(build_dir))} && "
        f"cmake -S {shlex.quote(str(src_dir))} -B . {' '.join(shlex.quote(d) for d in defines)} "
        f"&& cmake --build . -- -j {jobs}"
    )
    print(f"[baseline] building {settings['CMAKE_BUILD_TYPE']} baseline in {build_dir} "
          f"(log: {log})", flush=True)
    done = launch.run_in_shell(command)
    log.write_text((done.stdout or "") + (done.stderr or ""))
    if done.returncode != 0:
        raise BaselineError(
            f"baseline build failed (exit {done.returncode}); see {log}\n"
            + launch.tail("output", (done.stdout or "") + (done.stderr or ""), 40))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("build_exe", type=Path,
                        help="the branch binary whose configuration the baseline copies")
    parser.add_argument("--ref", default=None, help=f"git ref (default {DEFAULT_REF})")
    args = parser.parse_args()
    print(ensure_baseline(args.build_exe.resolve(), args.ref))
