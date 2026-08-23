#!/usr/bin/env bash
# =============================================================================
# QC the paired object -- SLURM job.
#
#     conda activate scplus-pairing
#     mkdir -p jobs
#     sbatch slurm/qc_paired.sbatch
#
# NOT for a login node. Even with the column subsetting in qc_paired.py this
# reads the RNA block whole (~3 GB dense) and streams all 394k ATAC columns for
# the NaN scan. The full validator (validate_h5mu.py -> mudata.read) is worse: it
# materialises every matrix, ~40 GB, which is what hung an interactive session.
#
# Follows maxsonBraunLab/cutTag-pipeline: partition batch, logs in jobs/.
# =============================================================================
#SBATCH --job-name=qc_paired
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
# 96G, not 32G. The backed/subsetting path in qc_paired.py peaks near 6 GB, but
# validate_h5mu.py at the end calls mudata.read() unbacked and then scans every
# matrix. Derivation on this reference (25,323 metacells):
#     RNA   25,323 x  34,112 x 4 B =  3.2 GB
#     ATAC  25,323 x 393,832 x 4 B = 37.2 GB
#     resident after read          = 40.4 GB
#     + one transient copy in scan = 80.7 GB   <- the peak
#     x1.15 allocator headroom     = 92.9 GB   <- hence 96G
# Sizing this to the cheap step would OOM on the expensive one.
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=jobs/qc_paired_%j.out
#SBATCH --error=jobs/qc_paired_%j.err

# NO `set -e`. The checks below are deliberately allowed to fail independently:
# a FAIL verdict from qc_paired.py exits 2, and under `set -e` that would abort
# the script before the full validator ran -- skipping the expensive check
# exactly when a problem has been found and it is most wanted. Each step
# captures its own status instead, and the script exits on the aggregate.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR" || exit 1
[[ -f setenv.sh ]] || { echo "ERROR: setenv.sh not found" >&2; exit 1; }
# shellcheck disable=SC1091
source setenv.sh

OUT="${PAIRED:-ACC_GEX.h5mu}"
[[ -f "$OUT" ]] || { echo "ERROR: $OUT does not exist" >&2; exit 1; }

PAIRING_ENV="scplus-pairing"
if [[ "${CONDA_DEFAULT_ENV:-}" != "$PAIRING_ENV" ]]; then
    CONDA_BASE="${CONDA_BASE:-$HOME/miniconda3}"
    if [[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
        # shellcheck disable=SC1091
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        conda activate "$PAIRING_ENV"
    else
        echo "ERROR: '$PAIRING_ENV' not active and no conda.sh under $CONDA_BASE" >&2
        exit 1
    fi
fi
echo "env: ${CONDA_DEFAULT_ENV:-?}  |  object: $OUT ($(du -h "$OUT" | cut -f1))"

echo "=== [1/3] structure (metadata only, seconds) ==="
python 03_pipeline/peek_h5mu.py "$OUT" \
    --expect-metacells 25323 --expect-groups 24 --group-key "$GROUP_KEY"
PEEK_RC=$?
echo "peek exit: $PEEK_RC"

echo
echo "=== [2/3] biology: marker specificity and locus linkage ==="
/usr/bin/time -v python 02_pair/qc_paired.py \
    --h5mu             "$OUT" \
    --markers          00_inspect/marker_tss_hg38.tsv \
    --group-key        "$GROUP_KEY" \
    --expect-groups    24 \
    --expect-metacells 25323 \
    --out              qc_paired
QC_RC=$?
echo "qc exit: $QC_RC  (2 = FAIL verdict; the report still written)"

echo
echo "=== [3/3] full NaN / structural validation (materialises every matrix) ==="
# Runs unconditionally, INCLUDING after a FAIL above -- when something is wrong
# is exactly when the exhaustive scan is most informative. This is the step the
# 96G allocation exists for.
python 03_pipeline/validate_h5mu.py "$OUT"
VAL_RC=$?
echo "validate exit: $VAL_RC"

echo
echo "=== summary ==="
printf '  peek     : %s\n' "$( [[ $PEEK_RC -eq 0 ]] && echo PASS || echo "PROBLEM (rc=$PEEK_RC)" )"
printf '  qc       : %s\n' "$( [[ $QC_RC   -eq 0 ]] && echo PASS || echo "FAIL (rc=$QC_RC)" )"
printf '  validate : %s\n' "$( [[ $VAL_RC  -eq 0 ]] && echo PASS || echo "PROBLEM (rc=$VAL_RC)" )"
echo
echo "Read qc_paired.md -- the verdict line, then section 3 (locus linkage)."
echo "A FAIL on linkage means region-to-gene inference will not work on this"
echo "object and SCENIC+ should not be run until it is understood."
# Nonzero if ANY check failed, so sacct shows FAILED rather than a green job
# with a bad object.
if [[ $PEEK_RC -ne 0 || $QC_RC -ne 0 || $VAL_RC -ne 0 ]]; then
    exit 2
fi
exit 0
