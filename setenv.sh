#!/usr/bin/env bash
# =============================================================================
# Paths for this project. Source it in every new shell:
#
#     source setenv.sh
#
# Why this file exists: `export REF=...` typed at a prompt is gone the moment
# you log out, and the next login node gets nothing. An unset $REF makes
# "$REF/rna.h5ad" expand to "/rna.h5ad", which fails 20 stack frames deep in
# h5py. Keeping the paths in the repo makes them survive sessions, and the
# checks below fail loudly instead of silently producing "/".
#
# Safe to source repeatedly. Does not activate a conda env -- that is deliberate,
# since the h5py-only tools and the pairing step want different interpreters.
# =============================================================================

# --- the GLUE-integrated reference (EDIT IF IT MOVES) ------------------------
export REF=/home/groups/MaxsonBraunScratch/worme/projects/scATAC/251112_hematopoiesis_ref/integration/output/02

# --- pairing parameters, so the sbatch and manual runs cannot drift ----------
# Rationale for each value: docs/GROUPS.md
export GROUP_KEY=predicted_CellType_Broad
export LATENT_KEY=X_glue
export K=50
export OVERSAMPLE=8
export FLOOR=60
export MIN_CELLS=50

# --- outputs -----------------------------------------------------------------
export LABELS=atac_labels.tsv
export PAIRED=ACC_GEX.h5mu

# --- verify, don't assume ----------------------------------------------------
_ok=1
if [[ ! -d "$REF" ]]; then
    echo "setenv.sh: REF does not exist: $REF" >&2
    _ok=0
else
    for f in rna.h5ad atac.h5ad atac_metadata_with_transferred_labels.tsv; do
        [[ -f "$REF/$f" ]] || { echo "setenv.sh: missing $REF/$f" >&2; _ok=0; }
    done
fi
if (( _ok )); then
    echo "setenv.sh: REF=$REF"
    echo "setenv.sh: rna.h5ad, atac.h5ad and the label TSV are all present"
    echo "setenv.sh: k=$K oversample=$OVERSAMPLE floor=$FLOOR group=$GROUP_KEY"
else
    echo "setenv.sh: FIX THE PATHS ABOVE before running anything." >&2
fi
unset _ok
