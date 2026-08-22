#!/usr/bin/env bash
# =============================================================================
# Preflight for 03_pipeline/environment.yml -- run BEFORE `mamba env create`.
# =============================================================================
# `mamba env create` on this environment installs a ~900-line pinned stack and
# builds six packages from source. When it fails, it fails minutes in, with a
# solver trace that points at a symptom rather than a cause. Two real examples
# this script catches in seconds:
#
#   * `mallet` is packaged in conda-forge for linux-64 ONLY (single build,
#     2.0.8). On macOS the solve dies with "mallet =* * does not exist (perhaps
#     a typo or a missing channel)", which reads like a channel misconfiguration
#     and is not.
#   * `python=3.11.8` is a load-bearing PATCH pin (scenicplus declares
#     requires-python ">=3.8,<=3.11.8", inclusive). conda-forge has no
#     pandas 1.5.0 build for 3.11 -- which is fine, because pandas is installed
#     by pip where cp311 wheels exist -- but moving pandas into the conda
#     section would make the solve unsatisfiable.
#
# Exit 0 = safe to create. Exit 1 = fix what is printed first.
# =============================================================================
set -uo pipefail

YML="${1:-$(dirname "$0")/environment.yml}"
fail=0
note() { printf '  %s\n' "$*"; }
bad()  { printf '  FAIL: %s\n' "$*"; fail=1; }
ok()   { printf '  ok  : %s\n' "$*"; }

echo "=============================================================="
echo "preflight: $YML"
echo "=============================================================="

[[ -f "$YML" ]] || { echo "ERROR: no such file: $YML" >&2; exit 1; }

# ---------------------------------------------------------------- 1. platform
echo
echo "[1] platform"
UNAME_S="$(uname -s)"; UNAME_M="$(uname -m)"
case "$UNAME_S/$UNAME_M" in
    Linux/x86_64)  SUBDIR=linux-64  ;;
    Linux/aarch64) SUBDIR=linux-aarch64 ;;
    Darwin/arm64)  SUBDIR=osx-arm64 ;;
    Darwin/x86_64) SUBDIR=osx-64    ;;
    *)             SUBDIR=unknown   ;;
esac
note "$UNAME_S $UNAME_M  ->  conda subdir: $SUBDIR"
if grep -qE '^\s*-\s*mallet' "$YML"; then
    if [[ "$SUBDIR" == "linux-64" ]]; then
        ok "mallet is requested and this IS linux-64 (the only platform it is built for)"
    else
        bad "environment.yml requests 'mallet', which conda-forge builds ONLY for linux-64."
        note "     This platform is $SUBDIR, so the solve WILL fail with a misleading"
        note "     \"does not exist (perhaps a typo or a missing channel)\"."
        note "     Fix: comment out the 'mallet' line and install upstream Mallet:"
        note "       wget https://mimno.github.io/Mallet/dist/Mallet-202108-bin.tar.gz"
        note "       tar xzf Mallet-202108-bin.tar.gz"
        note "     then pass --mallet-path <dir>/bin/mallet to run_cistopic.py."
    fi
fi

# ---------------------------------------------------------------- 2. solver
echo
echo "[2] solver"
SOLVER=""
for c in mamba micromamba conda; do
    if command -v "$c" >/dev/null 2>&1; then SOLVER="$c"; break; fi
done
if [[ -z "$SOLVER" ]]; then
    bad "no mamba/micromamba/conda on PATH"
else
    ok "using $SOLVER ($(command -v "$SOLVER"))"
    [[ "$SOLVER" == "conda" ]] && note "     mamba is much faster on a lock this size; consider mamba/miniforge."
fi

# ------------------------------------------------------- 3. python patch pin
echo
echo "[3] python pin"
PYPIN="$(grep -oE '^\s*-\s*python=[0-9.]+' "$YML" | head -1 | grep -oE '[0-9.]+')"
if [[ -z "$PYPIN" ]]; then
    bad "no explicit python pin found"
elif [[ "$PYPIN" == *.*.* ]]; then
    ok "python=$PYPIN (patch-level pin -- required: scenicplus requires-python is <=3.11.8, inclusive)"
else
    bad "python=$PYPIN is a MINOR pin. scenicplus declares requires-python \">=3.8,<=3.11.8\";"
    note "     conda resolves python=3.11 to the newest 3.11.x, which is past 3.11.8,"
    note "     and pip then refuses with \"requires a different Python\". Pin the patch."
fi

# ------------------------------------------- 4. pandas must stay in the pip section
echo
echo "[4] pandas placement"
PIP_LINE="$(grep -n '^\s*-\s*pip:' "$YML" | head -1 | cut -d: -f1)"
if [[ -z "$PIP_LINE" ]]; then
    bad "no pip: section found"
else
    CONDA_SEC="$(head -n "$PIP_LINE" "$YML")"
    if grep -qE '^\s*-\s*pandas\s*=' <<<"$CONDA_SEC"; then
        bad "pandas is pinned in the CONDA section. conda-forge has no pandas 1.5.0 build"
        note "     for python 3.11, so the solve becomes unsatisfiable. Keep pandas in pip"
        note "     (cp311 wheels exist there)."
    else
        ok "pandas is not in the conda section (correct -- it is a pip pin)"
    fi
fi

# --------------------------------------------- 5. unresolved placeholders
echo
echo "[5] placeholders"
if grep -nE '<COMMIT_SHA>|<[A-Z_]{3,}>' "$YML" | grep -vE '^\s*[0-9]+:\s*#' | grep -q .; then
    bad "unresolved placeholders remain -- pip will fail on the git refs:"
    grep -nE '<COMMIT_SHA>|<[A-Z_]{3,}>' "$YML" | grep -vE ':\s*#' | sed 's/^/       /'
else
    ok "no unresolved placeholders outside comments"
fi

# --------------------------------------------- 6. git deps pinned to commits
echo
echo "[6] git dependency pins"
while IFS= read -r line; do
    [[ "$line" =~ ^[[:space:]]*# ]] && continue
    pkg="$(sed -E 's/^[[:space:]]*-[[:space:]]*([A-Za-z0-9_.-]+)[[:space:]]*@.*/\1/' <<<"$line")"
    ref="$(sed -E 's/.*@([A-Za-z0-9._-]+)[[:space:]]*$/\1/' <<<"$line")"
    if [[ "$ref" =~ ^[0-9a-f]{40}$ ]]; then
        ok "$pkg pinned to commit ${ref:0:12}"
    elif [[ "$ref" =~ ^v?[0-9] ]]; then
        ok "$pkg pinned to tag $ref"
    else
        bad "$pkg tracks moving ref '$ref' -- the env is not reproducible."
        note "     Install once, then: pip list --format=freeze | grep -i $pkg"
    fi
done < <(grep -E '^\s*-\s*[A-Za-z0-9_.-]+\s*@\s*git\+' "$YML")

# --------------------------------------------- 7. disk for the build
echo
echo "[7] disk"
avail_kb="$(df -Pk "${CONDA_PKGS_DIRS:-$HOME}" 2>/dev/null | awk 'NR==2{print $4}')"
if [[ -n "${avail_kb:-}" ]]; then
    avail_gb=$(( avail_kb / 1024 / 1024 ))
    if (( avail_gb < 15 )); then
        bad "only ${avail_gb} GB free on ${CONDA_PKGS_DIRS:-$HOME}; the package cache plus"
        note "     six source builds need ~15 GB. Set CONDA_PKGS_DIRS somewhere with space."
    else
        ok "${avail_gb} GB free on ${CONDA_PKGS_DIRS:-$HOME}"
    fi
fi

echo
echo "=============================================================="
if (( fail )); then
    echo "PREFLIGHT FAILED -- fix the items above before \`mamba env create\`."
    exit 1
fi
echo "PREFLIGHT PASSED -- safe to run:"
echo "    ${SOLVER:-mamba} env create -f $YML"
echo "=============================================================="
