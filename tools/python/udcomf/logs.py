"""One log for the whole thermal-comfort chain.

Library modules log through ``stage_logger(<stage>)`` and never configure
handlers, so they stay silent inside other programs. The command-line
wrappers call :func:`configure` once; it writes every message to the case log
file and to the console in the form

    12:01:03 INFO    [radiation ] 1.1 m shortwave: 12 tiles, 4 resumed
    12:03:40 INFO    [radiation ]   tile 5/12 done (41 %)

Stage headers sit at the left margin; the details of a stage are indented
under it. Use :func:`header` for the former and the stage logger for the
latter.
"""

from __future__ import annotations

import logging
from pathlib import Path
import sys

ROOT = "udcomf"
STAGES = ("pipeline", "preflight", "merge", "atmosphere", "radiation", "exchange", "indices")
_WIDTH = max(len(stage) for stage in STAGES)


class StageFormatter(logging.Formatter):
    """``time level [stage] message`` with stage details indented."""

    def format(self, record: logging.LogRecord) -> str:
        parts = record.name.split(".")
        stage = parts[1] if len(parts) > 1 else "pipeline"
        indent = "" if getattr(record, "header", False) else "  "
        text = f"{self.formatTime(record, '%H:%M:%S')} {record.levelname:<7} [{stage:<{_WIDTH}}] {indent}{record.getMessage()}"
        if record.exc_info:
            text += "\n" + self.formatException(record.exc_info)
        return text


def stage_logger(stage: str) -> logging.Logger:
    """Logger of one stage; its messages are indented under the stage header."""
    if stage not in STAGES:
        raise ValueError(f"Unknown log stage {stage!r}; choose from {list(STAGES)}")
    return logging.getLogger(f"{ROOT}.{stage}")


def header(stage: str, message: str, *args) -> None:
    """Log a stage header at the left margin."""
    stage_logger(stage).info(message, *args, extra={"header": True})


def configure(log_file: Path | None, *, console: bool = True,
              level: int = logging.INFO) -> logging.Logger:
    """Send every udcomf message to ``log_file`` (appended) and the console.

    Calling it again replaces the handlers, so a wrapper may reconfigure
    without duplicating output.
    """
    root = logging.getLogger(ROOT)
    for handler in list(root.handlers):
        root.removeHandler(handler)
        handler.close()
    root.setLevel(level)
    root.propagate = False
    formatter = StageFormatter()
    if log_file is not None:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, mode="a", encoding="utf-8")
        file_handler.setFormatter(formatter)
        root.addHandler(file_handler)
    if console:
        stream = logging.StreamHandler(sys.stdout)
        stream.setFormatter(formatter)
        root.addHandler(stream)
    return root
