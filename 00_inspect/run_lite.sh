#!/usr/bin/env bash
# =============================================================================
# Run inspect_h5ad_lite.py with whatever python on this machine has h5py.
# =============================================================================
# The lite inspector needs h5py and nothing else. Rather than make you find an
# interpreter that has it, this scans every conda env's site-packages for h5py
# and re-executes the inspector with the first one that qualifies.
#
# All arguments are passed straight through:
#
#     bash 00_inspect/run_lite.sh "$REF/rna.h5ad" "$REF/atac.h5ad" \
#         --tsv "$REF/combined_glue_embeddings.tsv" --out report_lite
#
# If nothing has h5py it says so and prints the one-line fix, instead of
# failing with a ModuleNotFoundError from whichever python happened to be first
# on PATH.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SCRIPT="$HERE/inspect_h5ad_lite.py"
[[ -f "$SCRIPT" ]] || { echo "ERROR: $SCRIPT not found" >&2; exit 1; }

candidates() {
    # Currently-active and on-PATH pythons first -- cheapest if they work.
    command -v python3 2>/dev/null
    command -v python  2>/dev/null
    [[ -n "${CONDA_PREFIX:-}" ]] && printf '%s\n' "$CONDA_PREFIX/bin/python"
    # Then every env we can find on disk.
    local roots=()
    [[ -n "${MAMBA_ROOT_PREFIX:-}" ]] && roots+=("$MAMBA_ROOT_PREFIX/envs")
    [[ -n "${CONDA_PREFIX:-}"      ]] && roots+=("$(dirname "$CONDA_PREFIX")")
    [[ -n "${CONDA_ENVS_PATH:-}"   ]] && roots+=(${CONDA_ENVS_PATH//:/ })
    roots+=("$HOME/.conda/envs" "$HOME/micromamba/envs" "$HOME/miniconda3/envs"
            "$HOME/miniforge3/envs" "$HOME/anaconda3/envs" "$HOME/mambaforge/envs")
    for t in micromamba mamba conda; do
        if command -v "$t" >/dev/null 2>&1; then
            roots+=("$(dirname "$(dirname "$(command -v "$t")")")/envs")
        fi
    done
    for r in "${roots[@]}"; do
        [[ -d "$r" ]] || continue
        for e in "$r"/*/; do [[ -x "${e}bin/python" ]] && printf '%s\n' "${e}bin/python"; done
    done
}

# A python qualifies if h5py is present in ITS site-packages -- checked on the
# filesystem, so a broken interpreter cannot masquerade as a missing package.
PICK=""
while read -r py; do
    [[ -z "$py" || ! -x "$py" ]] && continue
    prefix="$(dirname "$(dirname "$py")")"
    if compgen -G "$prefix/lib/python*/site-packages/h5py" >/dev/null 2>&1; then
        # Confirm it actually imports before committing to it.
        if "$py" -c "import h5py" >/dev/null 2>&1; then PICK="$py"; break; fi
    fi
done < <(candidates | awk 'NF && !seen[$0]++')

if [[ -z "$PICK" ]]; then
    cat <<'EOF'
ERROR: no python with h5py found on this machine.

Fix (any one of these, seconds -- h5py ships as a binary wheel):

    pip install --user h5py
    # or, without touching base:
    micromamba create -y -n h5 -c conda-forge python=3.11 h5py
    micromamba run -n h5 python 00_inspect/inspect_h5ad_lite.py ...

EOF
    exit 1
fi

echo "[run_lite] using $PICK"
exec "$PICK" "$SCRIPT" "$@"
