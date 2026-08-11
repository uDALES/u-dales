"""Lightweight benchmark timing for radiation preprocessing."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import socket
import subprocess
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

try:
    import resource
except ImportError:  # pragma: no cover - unavailable on Windows
    resource = None


TIMESTAMP_FIELDS = (
    "index",
    "total",
    "model_time_seconds",
    "simulation_datetime",
    "mode",
    "ghi_w_m2",
    "dni_w_m2",
    "dsky_w_m2",
    "zenith_degrees",
    "azimuth_local_degrees",
    "direct_wall_seconds",
    "direct_cpu_seconds",
    "net_wall_seconds",
    "net_cpu_seconds",
    "step_wall_seconds",
    "step_cpu_seconds",
    "cumulative_wall_seconds",
    "peak_rss_mib",
    "completed_at_utc",
)


@dataclass(frozen=True)
class TimingSample:
    wall_seconds: float
    cpu_seconds: float
    peak_rss_mib: float | None


class Stopwatch:
    """Measure elapsed monotonic wall time and process CPU time."""

    def __init__(self) -> None:
        self._wall_start = time.perf_counter()
        self._cpu_start = time.process_time()

    def sample(self) -> TimingSample:
        return TimingSample(
            wall_seconds=time.perf_counter() - self._wall_start,
            cpu_seconds=time.process_time() - self._cpu_start,
            peak_rss_mib=peak_rss_mib(),
        )


def peak_rss_mib() -> float | None:
    """Return peak resident memory for this process in MiB."""
    if resource is None:
        return None
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform == "darwin":
        return value / (1024.0 * 1024.0)
    return value / 1024.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_metadata(repo_root: Path) -> dict[str, Any]:
    metadata: dict[str, Any] = {}
    try:
        commit = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(repo_root), "status", "--porcelain=v1"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.splitlines()
    except (OSError, subprocess.CalledProcessError):
        return metadata
    metadata["git_commit"] = commit
    metadata["git_dirty"] = bool(status)
    metadata["git_status"] = status
    return metadata


def benchmark_metadata(repo_root: Path, case_dir: Path) -> dict[str, Any]:
    """Capture enough runtime provenance to reproduce a benchmark."""
    metadata: dict[str, Any] = {
        "case_dir": str(case_dir),
        "command": sys.argv,
        "hostname": socket.gethostname(),
        "pid": os.getpid(),
        "platform": platform.platform(),
        "python_executable": sys.executable,
        "python_version": platform.python_version(),
    }
    if hasattr(os, "sched_getaffinity"):
        metadata["cpu_affinity"] = sorted(os.sched_getaffinity(0))
    metadata.update(_git_metadata(repo_root))

    extensions = sorted(
        (repo_root / "tools/python/udprep").glob("directshortwave_f2py*.so")
    )
    if extensions:
        extension = extensions[0]
        metadata["directshortwave_extension"] = str(extension)
        metadata["directshortwave_extension_sha256"] = _sha256(extension)
    return metadata


class ShortwaveTimingRecorder:
    """Write durable per-timestamp progress and an atomic JSON summary."""

    def __init__(
        self,
        prefix: Path,
        *,
        metadata: dict[str, Any],
        overwrite: bool,
    ) -> None:
        prefix = Path(prefix).expanduser().resolve()
        self.csv_path = Path(f"{prefix}.timestamps.csv")
        self.summary_path = Path(f"{prefix}.summary.json")
        existing = [path for path in (self.csv_path, self.summary_path) if path.exists()]
        if existing and not overwrite:
            raise FileExistsError(
                "Timing output already exists. Pass --overwrite to replace: "
                + ", ".join(str(path) for path in existing)
            )

        self.csv_path.parent.mkdir(parents=True, exist_ok=True)
        self.summary_path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = self.csv_path.open("w", encoding="ascii", newline="")
        self._writer = csv.DictWriter(self._handle, fieldnames=TIMESTAMP_FIELDS)
        self._writer.writeheader()
        self._handle.flush()

        self._run_timer = Stopwatch()
        self._summary: dict[str, Any] = {
            "schema_version": 1,
            "status": "running",
            "started_at_utc": datetime.now(timezone.utc).isoformat(),
            "timestamps_completed": 0,
            "stages": [],
            "metadata": metadata,
        }
        self._write_summary()

    def elapsed_wall_seconds(self) -> float:
        return self._run_timer.sample().wall_seconds

    def update_metadata(self, **values: Any) -> None:
        self._summary["metadata"].update(values)
        self._write_summary()

    def record_stage(self, name: str, sample: TimingSample) -> None:
        self._summary["stages"].append(
            {
                "name": name,
                "wall_seconds": sample.wall_seconds,
                "cpu_seconds": sample.cpu_seconds,
                "peak_rss_mib": sample.peak_rss_mib,
                "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            }
        )
        self._write_summary()

    def record_timestamp(self, values: dict[str, Any]) -> None:
        row = {field: values.get(field, "") for field in TIMESTAMP_FIELDS}
        row["cumulative_wall_seconds"] = self.elapsed_wall_seconds()
        row["peak_rss_mib"] = peak_rss_mib()
        row["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
        self._writer.writerow(row)
        self._handle.flush()
        self._summary["timestamps_completed"] += 1
        self._write_summary()

    def finish(
        self,
        *,
        status: str,
        total: TimingSample,
        outputs: dict[str, str | None] | None = None,
        error: str | None = None,
    ) -> None:
        self._summary.update(
            {
                "status": status,
                "completed_at_utc": datetime.now(timezone.utc).isoformat(),
                "total_wall_seconds": total.wall_seconds,
                "total_cpu_seconds": total.cpu_seconds,
                "peak_rss_mib": total.peak_rss_mib,
            }
        )
        if outputs is not None:
            self._summary["outputs"] = outputs
        if error is not None:
            self._summary["error"] = error
        self._write_summary()
        self._handle.close()

    def _write_summary(self) -> None:
        temp_path = self.summary_path.with_name(
            f".{self.summary_path.name}.{os.getpid()}.tmp"
        )
        with temp_path.open("w", encoding="ascii") as handle:
            json.dump(self._summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
        temp_path.replace(self.summary_path)


@contextmanager
def timed_stage(
    recorder: ShortwaveTimingRecorder | None, name: str
) -> Iterator[None]:
    if recorder is None:
        yield
        return
    timer = Stopwatch()
    try:
        yield
    finally:
        recorder.record_stage(name, timer.sample())
