#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
MODULE_NAME="${MODULE_NAME:-directshortwave_f2py}"
SOURCE_FILE="${SOURCE_FILE:-tools/python/fortran/directShortwave_f2py.f90}"
TARGET_DIR="${TARGET_DIR:-${ROOT_DIR}/tools/python/udprep}"

cd "${ROOT_DIR}"

python - <<'PY'
import pathlib
import sysconfig

include_dir = pathlib.Path(sysconfig.get_paths().get("include", ""))
header = include_dir / "Python.h"
if not header.is_file():
    raise SystemExit(
        f"Python development headers not found at {header}. "
        "Use an interpreter with Python.h available."
    )
PY

python -m numpy.f2py -c -m "${MODULE_NAME}" "${SOURCE_FILE}" --f90flags='-ffree-line-length-none'

artifact="$(find . -maxdepth 1 -name "${MODULE_NAME}*.so" | head -n 1)"
if [ -z "${artifact}" ]; then
    echo "Build succeeded but no shared object was found for ${MODULE_NAME}" >&2
    exit 1
fi

mkdir -p "${TARGET_DIR}"
rm -f "${TARGET_DIR}/${MODULE_NAME}"*.so
mv -f "${artifact}" "${TARGET_DIR}/"
echo "Moved ${artifact} to ${TARGET_DIR}/"
