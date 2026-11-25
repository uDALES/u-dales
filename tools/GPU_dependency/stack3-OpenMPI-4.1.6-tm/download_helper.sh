#!/usr/bin/env bash
# Sourced by the *_s0_download.sh scripts. Provides
#
#   fetch_archive URL FILE
#
# which downloads FILE from URL unless FILE already exists. It tries curl with
# the system CA store first. On HX1 (RHEL 8) that store lacks the SSL.com root
# that downloads.unidata.ucar.edu's certificate chain ends in, so on a
# certificate failure it retries with newer CA bundles: $UDALES_CA_BUNDLE if
# set, then the certifi bundle of any python3 on PATH, then the certifi bundle
# of the uDALES Python environment (tools/python/.venv) if this file lives
# inside a uDALES checkout. Certificate checking is never switched off.

_ca_bundles() {
    [[ -n "${UDALES_CA_BUNDLE:-}" ]] && echo "$UDALES_CA_BUNDLE"
    python3 -m certifi 2>/dev/null || true
    local d
    d=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
    while [[ "$d" != / ]]; do
        if [[ -x "$d/tools/python/.venv/bin/python" ]]; then
            "$d/tools/python/.venv/bin/python" -m certifi 2>/dev/null || true
            break
        fi
        d=$(dirname -- "$d")
    done
}

fetch_archive() {
    local url=$1 out=$2 tmp=$2.part bundle
    [[ -f "$out" ]] && return 0
    if command -v curl >/dev/null 2>&1; then
        if curl -fL -o "$tmp" "$url"; then mv "$tmp" "$out"; return 0; fi
        for bundle in $(_ca_bundles); do
            [[ -r "$bundle" ]] || continue
            echo "Retrying with the CA bundle $bundle" >&2
            if curl -fL --cacert "$bundle" -o "$tmp" "$url"; then mv "$tmp" "$out"; return 0; fi
        done
    elif command -v wget >/dev/null 2>&1; then
        if wget -O "$tmp" "$url"; then mv "$tmp" "$out"; return 0; fi
        for bundle in $(_ca_bundles); do
            [[ -r "$bundle" ]] || continue
            echo "Retrying with the CA bundle $bundle" >&2
            if wget --ca-certificate="$bundle" -O "$tmp" "$url"; then mv "$tmp" "$out"; return 0; fi
        done
    else
        echo "Error: neither curl nor wget is available." >&2
        return 1
    fi
    rm -f "$tmp"
    echo "Error: could not download $url" >&2
    echo "If this is a certificate problem, point UDALES_CA_BUNDLE at a current CA bundle." >&2
    return 1
}
