#!/usr/bin/env bash
# =============================================================================
# Find (or build) environments able to run the stage-1 inspectors.
# =============================================================================
# Stage 1 needs two things that the base conda env does not have:
#
#   00_inspect/inspect_anndata.py  ->  python with anndata + mudata
#   00_inspect/inspect_archr.R     ->  Rscript with ArchR
#
# You almost certainly already have both somewhere: the GLUE integration was
# run in a python env with anndata, and the ArchRProject was built in an R with
# ArchR. This script finds them instead of building anything new.
#
# READ-ONLY by default. Pass --create to build a minimal inspection env only if
# nothing suitable is found.
#
#     bash 00_inspect/find_inspect_env.sh
#     bash 00_inspect/find_inspect_env.sh --create
# =============================================================================
set -uo pipefail

CREATE=0
[[ "${1:-}" == "--create" ]] && CREATE=1
HERE="$(cd "$(dirname "$0")" && pwd)"

echo "=============================================================="
echo "stage-1 prerequisite probe"
echo "=============================================================="

# ---------------------------------------------------------------- conda tool
CONDA_TOOL=""
for c in micromamba mamba conda; do
    command -v "$c" >/dev/null 2>&1 && { CONDA_TOOL="$c"; break; }
done
if [[ -z "$CONDA_TOOL" ]]; then
    echo "ERROR: no micromamba/mamba/conda on PATH." >&2
    exit 1
fi
echo "conda tool: $CONDA_TOOL ($(command -v "$CONDA_TOOL"))"

# List env prefixes. `<tool> env list` can fail for reasons unrelated to the
# envs themselves (an unreadable ~/.condarc is enough), and a silent empty list
# would make this probe report "nothing found" when everything is fine. So:
# take the tool's answer, ALWAYS also walk the standard envs directories, and
# say so when the tool contributed nothing.
env_prefixes() {
    local from_tool=""
    case "$CONDA_TOOL" in
        micromamba) from_tool="$(micromamba env list 2>/dev/null | awk '/\// {print $NF}')" ;;
        *)          from_tool="$("$CONDA_TOOL" env list 2>/dev/null | grep -v '^#' | awk 'NF{print $NF}')" ;;
    esac
    [[ -n "$from_tool" ]] && printf '%s\n' "$from_tool"

    # Filesystem fallback: every envs dir we can infer, plus the active env.
    local roots=()
    [[ -n "${MAMBA_ROOT_PREFIX:-}" ]] && roots+=("$MAMBA_ROOT_PREFIX/envs")
    [[ -n "${CONDA_PREFIX:-}"      ]] && roots+=("$(dirname "$CONDA_PREFIX")")
    [[ -n "${CONDA_ENVS_PATH:-}"   ]] && roots+=(${CONDA_ENVS_PATH//:/ })
    roots+=("$HOME/.conda/envs" "$HOME/micromamba/envs" "$HOME/miniconda3/envs"
            "$HOME/miniforge3/envs" "$HOME/anaconda3/envs" "$HOME/mambaforge/envs")
    local base
    base="$(dirname "$(dirname "$(command -v "$CONDA_TOOL")")")"
    roots+=("$base/envs")
    for r in "${roots[@]}"; do
        [[ -d "$r" ]] || continue
        for e in "$r"/*/; do [[ -d "$e" ]] && printf '%s\n' "${e%/}"; done
    done
    # The currently-active env may live outside all of the above.
    [[ -n "${CONDA_PREFIX:-}" ]] && printf '%s\n' "$CONDA_PREFIX"
}

# Deduplicate while preserving order.
env_prefixes_uniq() { env_prefixes | awk 'NF && !seen[$0]++'; }

N_ENVS="$(env_prefixes_uniq | wc -l | tr -d ' ')"
echo "environments discovered: $N_ENVS"
if [[ "$N_ENVS" == "0" ]]; then
    echo "  WARNING: found no environment prefixes at all. '$CONDA_TOOL env list'"
    echo "  may be failing (an unreadable ~/.condarc will do it). Try:"
    echo "      $CONDA_TOOL env list"
    echo "  and if that errors, pass the python you want directly:"
    echo "      /path/to/env/bin/python 00_inspect/inspect_anndata.py --help"
fi

# ------------------------------------------------- 1. python with anndata
echo
echo "[1] python envs with anndata + mudata"
PY_OK=()
while read -r prefix; do
    [[ -z "$prefix" || ! -x "$prefix/bin/python" ]] && continue
    out="$("$prefix/bin/python" - <<'PY' 2>/dev/null
try:
    import anndata, importlib.metadata as md
    a = md.version("anndata")
    try:    m = md.version("mudata")
    except Exception: m = None
    print(f"anndata={a} mudata={m or 'MISSING'}")
except Exception:
    pass
PY
)"
    if [[ -n "$out" ]]; then
        echo "  FOUND  $prefix"
        echo "         $out"
        [[ "$out" == *"mudata=MISSING"* ]] || PY_OK+=("$prefix")
    fi
done < <(env_prefixes_uniq)

if ((${#PY_OK[@]})); then
    echo
    echo "  -> run the python inspector with:"
    echo "     ${PY_OK[0]}/bin/python 00_inspect/inspect_anndata.py --rna ... --atac ... --out report_py"
else
    echo "  none found with BOTH anndata and mudata."
    if ((CREATE)); then
        echo
        echo "  creating a minimal inspection env (small; ~1-2 min) ..."
        "$CONDA_TOOL" create -y -n scplus-inspect -c conda-forge \
            python=3.11 anndata mudata pandas numpy scipy h5py pyarrow \
            && echo "  -> $CONDA_TOOL activate scplus-inspect" \
            || echo "  ERROR: env creation failed; see output above."
    else
        echo "  re-run with --create to build one, or activate the env you ran GLUE in"
        echo "  (mudata is the only piece GLUE envs sometimes lack:"
        echo "   pip install mudata  -- it is pure python, no build)."
    fi
fi

# ------------------------------------------------- 2. Rscript with ArchR
echo
echo "[2] Rscript with ArchR"
if ! command -v Rscript >/dev/null 2>&1; then
    echo "  Rscript is NOT on PATH."
    if command -v module >/dev/null 2>&1 || [[ -n "${MODULEPATH:-}" ]]; then
        echo "  candidate modules (load one, then re-run this probe):"
        (module avail R 2>&1 || true) | tr ' ' '\n' \
            | grep -iE '^R(/|-[0-9])' | sort -u | sed 's/^/    module load /' | head -12
        echo "    (nothing listed above means 'module avail R' found no R modules)"
    fi
    # R inside conda envs is a common alternative on clusters.
    echo "  conda envs containing Rscript:"
    found_r=0
    while read -r prefix; do
        [[ -n "$prefix" && -x "$prefix/bin/Rscript" ]] && { echo "    $prefix/bin/Rscript"; found_r=1; }
    done < <(env_prefixes_uniq)
    ((found_r)) || echo "    none"
else
    echo "  Rscript: $(command -v Rscript)"
    if Rscript -e 'quit(status = !requireNamespace("ArchR", quietly=TRUE))' 2>/dev/null; then
        echo "  ArchR: available"
        echo "  -> Rscript 00_inspect/inspect_archr.R --archr-project <proj> --out report_archr"
    else
        echo "  ArchR: NOT installed in this R."
        echo "     The ArchRProject was built somewhere -- use that R. Check for an"
        echo "     R module or conda env, or look at the project's own logs."
        echo "     ArchR is not on CRAN; installing it fresh is a detour you do not"
        echo "     need if the original R is still around."
    fi
fi

echo
echo "=============================================================="
echo "Note: the PYTHON inspector answers most of what stage 1 needs"
echo "(cell counts, latent location, label columns, id format, matrix"
echo "provenance). If the R side is a hassle, run the python one first"
echo "and send that report -- do not block on ArchR."
echo "=============================================================="
