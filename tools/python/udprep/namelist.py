from __future__ import annotations

from pathlib import Path
from typing import Any

from udbase import UDBase


def _format_value(value: Any) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, int):
        return f"{value:d}"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def update_namelist_value(sim: "UDBase", namelist: str, variable: str, value: Any) -> Path:
    """Update a namelist variable in namoptions.<id>.

    Parameters
    ----------
    sim : UDBase
        Populated UDBase instance; expnr and path are used to locate inputs.
    namelist : str
        Namelist block name (e.g., "TREES").
    variable : str
        Variable name to update (e.g., "ntrees").
    value : Any
        Value to write; bool/int/float/str are supported.
    """

    if sim is None:
        raise ValueError("sim (UDBase instance) must be provided")

    namelist_path = Path(sim.path) / f"namoptions.{sim.expnr}"
    if not namelist_path.is_file():
        raise FileNotFoundError(f"Missing {namelist_path}")

    lines = namelist_path.read_text(encoding="ascii").splitlines(keepends=True)
    nml_lower = namelist.strip().lower()
    var_lower = variable.strip().lower()
    value_str = _format_value(value)

    start_idx = None
    end_idx = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.lower().startswith("&") and stripped[1:].strip().lower() == nml_lower:
            start_idx = i
            break

    if start_idx is not None:
        for i in range(start_idx + 1, len(lines)):
            if lines[i].strip().startswith("/"):
                end_idx = i
                break
        if end_idx is None:
            raise ValueError(f"Namelist block '{namelist}' in {namelist_path} has no terminator '/'")

        updated = False
        for i in range(start_idx + 1, end_idx):
            line = lines[i]
            if "=" not in line:
                continue
            left, right = line.split("=", 1)
            if left.strip().lower() != var_lower:
                continue
            comment = ""
            if "!" in right:
                _, comment = right.split("!", 1)
                comment = "!" + comment.rstrip("\n")
            newline = f"{left.rstrip()} = {value_str}"
            if comment:
                newline = f"{newline} {comment}"
            lines[i] = newline + ("\n" if line.endswith("\n") else "")
            updated = True
            break

        if not updated:
            indent = "  "
            insert_line = f"{indent}{variable} = {value_str}\n"
            lines.insert(end_idx, insert_line)
    else:
        block_name = namelist.strip().upper()
        lines.append(f"&{block_name}\n")
        lines.append(f"  {variable} = {value_str}\n")
        lines.append("/\n")

    namelist_path.write_text("".join(lines), encoding="ascii", newline="\n")
    setattr(sim, variable.strip(), value)
    return namelist_path
