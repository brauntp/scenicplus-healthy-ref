#!/usr/bin/env bash
# =============================================================================
# Size motif_enrichment_cistarget's memory -- as a batch job.
#
# WHY THIS EXISTS
# ---------------
# I documented size_cistarget_memory.py as "runs in seconds on a login node."
# Its first version counted motifs with
#
#     pa.ipc.open_file(src).get_batch(i).num_rows
#
# with `read_table(columns=[])` as a "safe" fallback. Measured in a fresh
# process on a 3.04 GB / 400k-column fixture:
#
#     get_batch(i).num_rows      1.43 s   peak RSS 3341 MB
#     read_table(columns=[])     0.82 s   peak RSS 3508 MB
#     read_table(columns=[one])  0.40 s   peak RSS  251 MB
#
# Both of those scale with FILE SIZE. On the real 33 GB rankings database that
# is tens of GB resident on a shared login node -- which is what made the run
# slow. The script now reads exactly one column (~131 KB of data) and only when
# --db is passed at all.
#
# So the login node is fine for the default invocation. This job exists for the
# case where you want --db against a database on a slow or cold filesystem, or
# simply do not want to find out on a shared node. The default path (counting
# .bed lines) measures 0.16 s for 21 sets totalling 1.69 M lines -- that has
# never needed a batch job.
#
#     sbatch slurm/size_cistarget.sbatch                    # sets only
#     sbatch slurm/size_cistarget.sbatch --with-db          # + confirm motifs
#
# Account and partition are not hardcoded: pass them on the command line, which
# overrides any directive here anyway.
#     sbatch --partition=... slurm/size_cistarget.sbatch
#SBATCH --job-name=size_ctx
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
# One bounded column read plus line counting. 8G is generous; it measured
# 283 MB peak on a 3.04 GB fixture after the fix.
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=jobs/%x-%j.out
#SBATCH --error=jobs/%x-%j.err

# No `set -e`: the report and the closing guidance should print even if the
# sizing script exits nonzero.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR" || exit 1

REGION_SETS="${REGION_SETS:-01_atac/region_sets}"
DB="${DB:-resources/cistarget_db/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather}"

WITH_DB=0
[[ "${1:-}" == "--with-db" ]] && WITH_DB=1

echo "==============================================================="
echo "job         : $SLURM_JOB_NAME ($SLURM_JOB_ID)"
echo "node        : $(hostname -s)"
echo "region sets : $REGION_SETS"
echo "read db     : $([[ $WITH_DB -eq 1 ]] && echo "yes -- $DB" || echo "no (motif count from --n-motifs)")"
echo "==============================================================="

if [[ ! -d "$REGION_SETS" ]]; then
    echo "ERROR: no region sets at $REGION_SETS" >&2
    echo "       sbatch slurm/region_sets.sbatch writes them." >&2
    exit 1
fi

ARGS=(--region-set-folder "$REGION_SETS")
if (( WITH_DB )); then
    if [[ ! -f "$DB" ]]; then
        echo "ERROR: --with-db given but no database at $DB" >&2
        exit 1
    fi
    ARGS+=(--db "$DB")
fi

/usr/bin/time -v python 03_pipeline/size_cistarget_memory.py "${ARGS[@]}" 2>&1
RC=$?

echo
echo "==============================================================="
echo "sizing exit: $RC"
echo
echo "Read the 'concurrent peak' table against the partition's per-node"
echo "memory limit. Two levers, and they are not equivalent:"
echo
echo "  --cpus-per-task N   halves the peak when halved, and costs wall"
echo "                      time in THIS rule only. Nothing scientific"
echo "                      changes."
echo "  --top-n N --write   makes every region set the same size, so the"
echo "                      peak becomes N x n_cpu rather than a sum over"
echo "                      the largest few. Cheaper AND predictable, but"
echo "                      membership becomes rank-based -- a region"
echo "                      admitted for one cell type can be rejected for"
echo "                      another. That is a scientific choice; see"
echo "                      docs/REGION_SETS.md before taking it."
echo "==============================================================="
exit "$RC"
