from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Type, Tuple
import time

SKIP = object()


class Section:
    def __init__(self, name: str, values: Dict[str, Any], sim: Any | None = None):
        self._name = name
        self.sim = sim
        for key, val in values.items():
            setattr(self, key, val)

    def run_steps(self, label: str, steps: List[Tuple[str, Callable[[], None]]]) -> None:
        for name, func in steps:
            start = time.perf_counter()
            print(f"[{label}] {name}...")
            func()
            elapsed = time.perf_counter() - start
            print(f"[{label}] {name} done in {elapsed:.3f} s")

    def write_changed_params(self) -> None:
        """
        Compare current values against defaults and write only changed parameters.

        This is a placeholder for logic that will diff section attributes against
        their default values and emit a minimal input file update.
        """

    def __repr__(self) -> str:
        lines = [f"{self._name}:"]
        for key, val in sorted(self.__dict__.items()):
            if key == "_name":
                continue
            lines.append(f"  {key}: {val}")
        return "\n".join(lines)


@dataclass(frozen=True)
class SectionSpec:
    name: str
    fields: List[str]
    defaults: Dict[str, Any | Callable[[Any], Any]]
    section_cls: Type[Section] = Section

    def get_default(self, key: str, ctx: Any) -> Any:
        default = self.defaults.get(key, SKIP)
        return default(ctx) if callable(default) else default
