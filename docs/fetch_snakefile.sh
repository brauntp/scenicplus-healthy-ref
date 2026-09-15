#!/usr/bin/env bash
# Fetch the pinned SCENIC+ Snakefile for READING, and verify it is the exact
# file this repo's documentation was written against.
#
# WHY THIS EXISTS
# ---------------
# The Snakefile is the computational core of this pipeline, and it is not in
# this repository. It ships inside the installed scenicplus package, at
#   <site-packages>/scenicplus/snakemake/Snakefile
# resolved at runtime by 03_pipeline/run_pipeline.sh via importlib.resources.
#
# That means an analyst cannot read the DAG without first building a ~900-package
# conda env -- so every statement this repo makes about the Snakefile (which
# rules run, what they consume, which config keys matter) is unverifiable at
# review time. This script closes that gap: it pulls the pinned file from the
# upstream tag and checks it byte-for-byte against a recorded SHA256.
#
# WHY IT IS FETCHED AND NOT VENDORED
# ----------------------------------
# scenicplus v1.0a2 declares "License :: Other/Proprietary License" in its
# pyproject.toml and ships no LICENCE.txt at that tag, so this repo does not
# redistribute the source. Fetching keeps upstream as the single source of
# truth and makes the pin tamper-evident instead of making a silent copy.
#
# The fetched file is for READING ONLY. It is written to a gitignored path and
# is never what runs -- run_pipeline.sh always resolves the INSTALLED package,
# so a stale or edited local copy cannot silently become the executed workflow.
#
#   Usage:  bash docs/fetch_snakefile.sh [--out DIR]
#   Exit:   0 = fetched and checksum matched
#           1 = download failed, or checksum MISMATCH (pin is stale or upstream moved)
set -euo pipefail

SCENICPLUS_TAG="v1.0a2"
SNAKEFILE_URL="https://raw.githubusercontent.com/aertslab/scenicplus/${SCENICPLUS_TAG}/src/scenicplus/snakemake/Snakefile"

# SHA256 of the Snakefile at tag v1.0a2, recorded when docs/PIPELINE_DAG.md was
# written. If this stops matching, the DAG documentation may no longer describe
# the file -- re-verify PIPELINE_DAG.md before trusting it.
SNAKEFILE_SHA256="b0f2e1553838068bb10d151c93e2d950deb2dc72681bc727644fd2fe6a21517a"

OUT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/reference"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --out) OUT_DIR="${2:?--out requires a directory}"; shift 2 ;;
        --out=*) OUT_DIR="${1#*=}"; shift ;;
        -h|--help) sed -n '2,31p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unrecognised option: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "$OUT_DIR"
OUT="${OUT_DIR}/Snakefile.${SCENICPLUS_TAG}"

echo "fetching scenicplus ${SCENICPLUS_TAG} Snakefile ..."
if ! curl -sS -L --fail --max-time 120 -o "$OUT" "$SNAKEFILE_URL"; then
    echo "ERROR: download failed: $SNAKEFILE_URL" >&2
    echo "  This needs outbound HTTPS. On a cluster login node behind a proxy," >&2
    echo "  fetch it elsewhere and copy it in, then re-run to check the sum." >&2
    exit 1
fi

got="$(sha256sum "$OUT" | cut -d' ' -f1)"
if [[ "$got" != "$SNAKEFILE_SHA256" ]]; then
    echo "ERROR: checksum MISMATCH -- this is not the file the docs describe." >&2
    echo "  expected: $SNAKEFILE_SHA256" >&2
    echo "  got:      $got" >&2
    echo "  Upstream may have re-tagged. Re-verify docs/PIPELINE_DAG.md against" >&2
    echo "  the new file before trusting either." >&2
    exit 1
fi

echo "ok: ${OUT}"
echo "    sha256 ${got} (matches the pin in docs/PIPELINE_DAG.md)"
echo ""
echo "This is a READ-ONLY reference copy. What actually runs is the INSTALLED"
echo "package, which run_pipeline.sh resolves independently. To confirm the"
echo "installed copy is the same file, with the scenicplus env active:"
echo ""
echo "  sha256sum \"\$(python -c 'from importlib.resources import files; print(files(\"scenicplus.snakemake\").joinpath(\"Snakefile\"))')\""
