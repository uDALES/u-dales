#!/usr/bin/env bash
#
# Summarise compiler warnings from a build log into the GitHub step summary.
#
# Reporting only: this never fails the build, by design. CI does not pin its
# compilers (`apt install gfortran`, `brew install gcc` on rolling runner
# images), so each runner legitimately reports a different set of warnings, and
# a newer image can introduce new ones with no change to this repository.
# Gating on that would turn CI red on an unrelated PR. Surfacing the counts
# keeps the signal without the brittleness.
#
# Classes listed in BENIGN below were reviewed in #334 and are correct as they
# stand; everything else is flagged for review. If a class here starts hiding
# real problems, drop it from the list rather than muting the report.
#
# Usage: summarise_warnings.sh <build.log> <label>

set -uo pipefail

log="${1:?usage: summarise_warnings.sh <build.log> <label>}"
label="${2:-build}"
out="${GITHUB_STEP_SUMMARY:-/dev/stdout}"

if [ ! -f "$log" ]; then
    # Loud on stdout: a missing log means this reported nothing, and a silent
    # exit 0 would be indistinguishable from a clean build.
    echo "WARNING: no build log at '${log}' -- no warnings were summarised." >&2
    echo "### Compiler warnings — ${label}" >> "$out"
    echo >> "$out"
    echo "_No build log at \`${log}\` — nothing to report._" >> "$out"
    exit 0
fi

# gfortran/clang warning classes are reported as a trailing [-Wname] tag.
BENIGN_RE='^-W(compare-reals|unused-value)$'

counts="$(grep -ohE '\[-W[a-z-]+\]' "$log" | tr -d '[]' | sort | uniq -c | sort -rn || true)"

total=$(printf '%s' "$counts" | awk '{s+=$1} END {print s+0}')
actionable=$(printf '%s' "$counts" | awk -v re="$BENIGN_RE" '$2 !~ re {s+=$1} END {print s+0}')

# Echoed to the job log as well as the summary, so the log alone shows whether
# this actually ran against a real build.
echo "${label}: ${total} warning(s), ${actionable} actionable, from '${log}'"

{
    echo "### Compiler warnings — ${label}"
    echo
    if [ "$total" -eq 0 ]; then
        echo "None."
        exit 0
    fi
    if [ "$actionable" -eq 0 ]; then
        echo "**0 actionable** (${total} known-benign)."
    else
        echo "**${actionable} actionable** of ${total} — see below."
    fi
    echo
    echo "| count | class | |"
    echo "|---|---|---|"
    printf '%s\n' "$counts" | while read -r n c; do
        [ -n "${n:-}" ] || continue
        if printf '%s' "$c" | grep -qE "$BENIGN_RE"; then
            tag="known-benign (#334)"
        else
            tag="**review**"
        fi
        echo "| ${n} | \`${c}\` | ${tag} |"
    done
} >> "$out"

[ "$actionable" -eq 0 ] && exit 0

# Pair each warning with the file:line: header that precedes it, so the report
# points at source rather than just counting.
{
    echo
    echo "<details><summary>Actionable warning sites</summary>"
    echo
    echo '```'
    # Two shapes to handle:
    #   gfortran: a "file:line:col:" header line, then "Warning: msg [-Wflag]"
    #   clang:    "file:line:col: warning: msg [-Wflag]" all on one line
    awk -v re="$BENIGN_RE" '
        function emit(loc, msg, flag,   p, n) {
            if (flag ~ re) return
            sub(/:[0-9]+:?$/, "", loc)          # drop the column, keep file:line
            n = split(loc, p, "/")
            printf "%s: %s [%s]\n", p[n], msg, flag
        }
        # clang single-line form
        /:[0-9]+:[0-9]+: *warning:/ && /\[-W[a-z-]+\]/ {
            match($0, /\[-W[a-z-]+\]/)
            flag = substr($0, RSTART+1, RLENGTH-2)
            loc = $0; sub(/: *warning:.*/, "", loc)
            msg = $0; sub(/.*warning: /, "", msg); sub(/ \[-W.*/, "", msg)
            emit(loc, msg, flag)
            next
        }
        # gfortran location header
        /^[^ ].*:[0-9]+:[0-9]+:[[:space:]]*$/ { loc = $0; next }
        # gfortran message, refers back to the last header
        /^[[:space:]]*Warning:/ && /\[-W[a-z-]+\]/ {
            match($0, /\[-W[a-z-]+\]/)
            flag = substr($0, RSTART+1, RLENGTH-2)
            msg = $0; sub(/.*Warning: /, "", msg); sub(/ \[-W.*/, "", msg)
            if (loc != "") emit(loc, msg, flag)
        }
    ' "$log" | sort -u -t: -k1,1 -k2,2n | head -50
    echo '```'
    echo "</details>"
} >> "$out"

exit 0
