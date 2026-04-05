#!/usr/bin/env bash
set -euo pipefail

echo "hostname: $(hostname -f 2>/dev/null || hostname)"
echo "uname: $(uname -a)"

if command -v module >/dev/null 2>&1; then
  echo "module_version: $(module --version 2>&1 | head -n 1)"
else
  echo "module_version: <none>"
fi

echo "mpiexec: $(command -v mpiexec 2>/dev/null || echo '<none>')"
if command -v mpiexec >/dev/null 2>&1; then
  echo "mpiexec_version: $(mpiexec --version 2>&1 | head -n 1)"
fi

echo "python3: $(command -v python3 2>/dev/null || echo '<none>')"
if command -v python3 >/dev/null 2>&1; then
  echo "python3_version: $(python3 --version 2>&1)"
  python3 -c "import sysconfig, pathlib; inc = sysconfig.get_paths().get('include',''); print('python_include:', inc); print('python_h:', pathlib.Path(inc)/'Python.h')" 2>/dev/null || true
fi

