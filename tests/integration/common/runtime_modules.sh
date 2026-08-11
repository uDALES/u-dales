#!/usr/bin/env bash

# Load a caller-selected runtime module stack. Integration tests otherwise
# inherit the active environment and never infer a cluster from the mere
# presence of the `module` command.
load_udales_runtime_modules() {
    local module_spec="${UDALES_RUNTIME_MODULES:-}"
    local -a requested_modules

    if [[ -z "${module_spec//[[:space:]]/}" ]]; then
        return 0
    fi

    if ! command -v module >/dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
        # shellcheck disable=SC1091
        source /etc/profile.d/modules.sh || true
    fi

    if ! command -v module >/dev/null 2>&1; then
        echo "ERROR: UDALES_RUNTIME_MODULES was set, but the module command is unavailable." >&2
        return 1
    fi

    read -r -a requested_modules <<< "$module_spec"
    module load "${requested_modules[@]}"
}
