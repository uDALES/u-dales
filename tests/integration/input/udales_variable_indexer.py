#!/usr/bin/env python3
"""
uDALES input variable indexer.

This is adapted from the JSON branch tooling, but for this branch the
relevant contract is:
1) variable appears in a Fortran namelist
2) variable is broadcast across MPI ranks
3) variable is present in the editor/schema description
"""

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional


class UdalesVariableIndexer:
    def __init__(self, src_dir: str, schema_file: str):
        self.src_dir = Path(src_dir)
        self.schema_file = Path(schema_file)

        self.namelists = defaultdict(set)
        self.namelist_files = {}
        self.broadcasts = set()
        self.schema_vars = defaultdict(set)

        self.broadcast_counts = Counter()
        self.namelist_counts = Counter()
        self.current_aliases = {}

    def parse_fortran_files(self, file_path: Optional[str] = None):
        if file_path:
            files = [Path(file_path)]
        else:
            files = list(self.src_dir.glob("*.f90"))

        for f90_file in files:
            if f90_file.name.endswith(".backup") or ".backup." in f90_file.name:
                continue
            with open(f90_file, "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()
            self.current_aliases = {}
            self._extract_use_aliases(content)
            self._extract_namelists(content, f90_file.name)
            self._extract_broadcasts(content)

    def _extract_namelists(self, content: str, filename: str):
        keywords = {
            "integer", "real", "logical", "character", "complex", "double",
            "precision", "parameter", "intent", "optional", "pointer", "target",
            "save", "public", "private", "dimension", "kind", "len", "if", "then",
            "else", "elseif", "endif", "select", "case", "default", "endselect",
            "do", "enddo", "while", "exit", "cycle", "where", "elsewhere",
            "endwhere", "forall", "endforall", "goto", "continue", "stop",
            "return", "call", "function", "subroutine", "program", "module",
            "contains", "use", "only", "implicit", "none", "include", "data",
            "common", "equivalence", "external", "intrinsic", "sequence", "type",
            "endtype", "interface", "endinterface", "procedure", "true", "false",
            "read", "write", "print", "open", "close", "inquire", "format", "iostat",
            "err", "end", "file", "unit", "rec", "fmt", "nml", "advance", "size",
            "j", "k", "i", "n", "m", "l", "ii", "jj", "kk", "nn", "mm", "ll"
        }

        lines = content.split("\n")
        current_namelist = None
        collecting_vars = False

        for line in lines:
            raw_line = line.rstrip()
            proc_line = re.sub(r"!.*$", "", line).strip()
            if not proc_line:
                continue

            match = re.match(r"\s*namelist\s*/\s*(\w+)\s*/\s*&?\s*(.*)", line, re.IGNORECASE)
            if match:
                current_namelist = match.group(1).upper()
                self.namelist_files[current_namelist] = filename
                var_part = match.group(2)
                if var_part:
                    self._extract_vars_from_line(var_part, current_namelist, keywords)
                collecting_vars = raw_line.endswith("&")
                if not collecting_vars:
                    current_namelist = None
                continue

            if collecting_vars and current_namelist:
                var_line = re.sub(r"&\s*$", "", raw_line).strip()
                if var_line:
                    self._extract_vars_from_line(var_line, current_namelist, keywords)
                if not raw_line.endswith("&"):
                    collecting_vars = False
                    current_namelist = None

    def _extract_vars_from_line(self, line: str, namelist: str, keywords: set):
        for var in re.findall(r"\b([a-zA-Z_][a-zA-Z0-9_]*)\b", line):
            var_lower = self.current_aliases.get(var.lower(), var.lower())
            if var_lower not in keywords and len(var_lower) > 1 and not re.match(r"^[a-z]$", var_lower):
                self.namelists[namelist].add(var_lower)
                self.namelist_counts[var_lower] += 1

    def _extract_broadcasts(self, content: str):
        for var in re.findall(r"(?:call\s+)?MPI_BCAST\s*\(\s*([a-zA-Z_]\w*)", content, re.IGNORECASE):
            var_lower = self.current_aliases.get(var.lower(), var.lower())
            if var_lower not in {"comm3d", "mpierr", "mpi_integer", "mpi_real", "my_real", "mpi_logical", "mpi_character"}:
                self.broadcasts.add(var_lower)
                self.broadcast_counts[var_lower] += 1

    def _extract_use_aliases(self, content: str):
        content_no_comments = re.sub(r"!.*$", "", content, flags=re.MULTILINE)
        for line in content_no_comments.split("\n"):
            match = re.match(r"\s*use\s+\w+\s*,\s*(.*)", line.strip(), re.IGNORECASE)
            if not match:
                continue
            for part in [p.strip() for p in match.group(1).split(",") if p.strip()]:
                if part.lower().startswith("only:"):
                    continue
                alias_match = re.match(r"([a-zA-Z_]\w*)\s*=>\s*([a-zA-Z_]\w*)", part)
                if alias_match:
                    self.current_aliases[alias_match.group(1).lower()] = alias_match.group(2).lower()

    def parse_json_schema(self):
        with open(self.schema_file, "r") as f:
            schema = json.load(f)
        for namelist_name, namelist_def in schema.get("properties", {}).items():
            if isinstance(namelist_def, dict) and "properties" in namelist_def:
                for var_name in namelist_def["properties"].keys():
                    self.schema_vars[namelist_name.upper()].add(var_name.lower())

    def run_analysis(self):
        self.parse_fortran_files()
        self.parse_json_schema()


def main():
    parser = argparse.ArgumentParser(description="Analyze uDALES namelist/broadcast/schema consistency")
    parser.add_argument("--src-dir", default="../../../src")
    parser.add_argument("--schema-file", default="../../../docs/schemas/udales_input_schema.json")
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    src_dir = (script_dir / args.src_dir).resolve()
    schema_file = (script_dir / args.schema_file).resolve()

    if not src_dir.exists() or not schema_file.exists():
        print("required input paths not found")
        sys.exit(1)

    indexer = UdalesVariableIndexer(str(src_dir), str(schema_file))
    indexer.run_analysis()
    print("Indexed", len(indexer.namelists), "namelists")


if __name__ == "__main__":
    main()
