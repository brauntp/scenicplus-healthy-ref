#!/usr/bin/env bash
# =============================================================================
# Preflight for the PAIRING step. Run before `mamba env create`, and again
# after activating, to confirm you are about to run in the right interpreter.
# =============================================================================
#     bash 03_pipeline/preflight_pairing.sh
#
# Answers one question: can the current python run 02_pair/aggregate_atac_sparse.py?
# The failure this exists to prevent is submitting a job that dies on
# `import numpy` after waiting in the queue.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
YML="$HERE/03_pipeline/pairing_env.yml"
FAIL=0

echo "=============================================================="
echo "pairing preflight"
echo "=============================================================="

# ---------------------------------------------------------------- [1] platform
echo "[1] platform"
UNAME="$(uname -s)-$(uname -m)"
echo "    $UNAME"
case "$UNAME" in
    Linux-x86_64) echo "    OK -- linux-64, the platform these pins were solved for" ;;
    *) echo "    NOTE: not linux-64. The pins in pairing_env.yml were verified for"
       echo "          linux-64; other platforms may solve differently." ;;
esac

# ---------------------------------------------------------------- [2] interpreter
echo
echo "[2] current python"
PY="$(command -v python || command -v python3 || true)"
if [[ -z "$PY" ]]; then
    echo "    FAIL: no python on PATH"
    FAIL=1
else
    echo "    $PY  ($("$PY" -c 'import sys;print(".".join(map(str,sys.version_info[:3])))' 2>/dev/null || echo '?'))"
    if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
        echo "    active conda env: $CONDA_DEFAULT_ENV"
        if [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
            echo "    WARNING: this is 'base'. The pairing step needs its own env;"
            echo "             installing into base is how base gets broken."
        fi
    else
        echo "    no conda env active"
    fi
fi

# ---------------------------------------------------------------- [3] imports
echo
echo "[3] required libraries in THIS python"
MISSING=()
if [[ -n "$PY" ]]; then
    for m in numpy scipy pandas h5py anndata mudata; do
        if v=$("$PY" -c "import $m,sys;print(getattr($m,'__version__','?'))" 2>/dev/null); then
            printf "    OK      %-8s %s\n" "$m" "$v"
        else
            printf "    MISSING %-8s\n" "$m"
            MISSING+=("$m")
        fi
    done
fi
if (( ${#MISSING[@]} )); then
    FAIL=1
    echo
    echo "    -> ${#MISSING[@]} library(ies) missing: ${MISSING[*]}"
    echo "    -> Create and activate the env FIRST:"
    echo "           mamba env create -f 03_pipeline/pairing_env.yml"
    echo "           conda activate scplus-pairing"
    echo "       (--dry-run needs these too -- it is cheap in time and memory,"
    echo "        not in dependencies.)"
fi

# ------------------------------------------------------- [4] does the CLI load?
echo
echo "[4] script loads and parses arguments"
if [[ -n "$PY" ]] && "$PY" "$HERE/02_pair/aggregate_atac_sparse.py" --help >/dev/null 2>&1; then
    echo "    OK -- aggregate_atac_sparse.py --help works"
else
    echo "    FAIL -- the script cannot even print its help in this python"
    FAIL=1
fi

# ------------------------------------------------------------- [5] solver present
echo
echo "[5] conda/mamba available to build the env"
SOLVER=""
for t in mamba micromamba conda; do
    command -v "$t" >/dev/null 2>&1 && { SOLVER="$t"; break; }
done
if [[ -n "$SOLVER" ]]; then
    echo "    OK -- $SOLVER ($(command -v $SOLVER))"
elif (( ${#MISSING[@]} == 0 )); then
    # Everything already imports, so there is nothing left to build. A solver
    # on PATH is only needed to CREATE the env, not to run in it.
    echo "    not on PATH -- but every library already imports, so the env"
    echo "    is already usable and nothing needs building. Not a problem."
else
    echo "    FAIL: libraries are missing AND no mamba/micromamba/conda is on"
    echo "          PATH to build the env with. Load the conda module first."
    FAIL=1
fi
[[ -f "$YML" ]] && echo "    env file: $YML" || { echo "    FAIL: $YML missing"; FAIL=1; }

# ----------------------------------------------------------------- [6] scratch
echo
echo "[6] disk for the output"
AVAIL=$(df -Pk . 2>/dev/null | awk 'NR==2{printf "%.0f", $4/1048576}')
echo "    $(pwd): ${AVAIL:-?} GB free"
if [[ -n "${AVAIL:-}" && "$AVAIL" -lt 60 ]]; then
    echo "    WARNING: the .h5mu is tens of GB at oversample 8. Write it"
    echo "             somewhere with room, not a quota'd \$HOME."
fi

echo
echo "=============================================================="
if (( FAIL )); then
    echo "PREFLIGHT FAILED -- fix the items above before submitting."
    echo "Most likely you have not created/activated the env yet:"
    echo "    mamba env create -f 03_pipeline/pairing_env.yml"
    echo "    conda activate scplus-pairing"
    echo "    bash 03_pipeline/preflight_pairing.sh     # re-run to confirm"
    exit 1
fi
echo "PREFLIGHT PASSED -- safe to run --dry-run, then sbatch."
echo "=============================================================="
