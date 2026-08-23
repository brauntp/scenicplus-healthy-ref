#!/usr/bin/env bash
# =============================================================================
# Region sets for motif enrichment, from metacell DARs -- SLURM job.
#
# Fills input_data.region_set_folder, which is a REQUIRED config entry that
# nothing else in this pipeline produced. See
# 01_cistopic/region_sets_from_metacells.py for why this replaces pycisTopic
# topic modelling here, and what that costs you.
#
#   sbatch slurm/region_sets.sbatch
#
# Memory is flat: the ATAC matrix is streamed in blocks of --block regions from
# a backed object, so peak RSS is set by the block (1.9 GB at 20,000 regions x
# 25,323 metacells), not by the 37 GB matrix. 16G is generous.
#
# Time is the uncertainty: Mann-Whitney per region per group is the cost centre,
# 20 blocks x 24 groups. Tens of minutes to a couple of hours. 4h is padded --
# a walltime kill loses everything, since BEDs are written only at the end.
# =============================================================================
#SBATCH --job-name=region_sets
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=jobs/%x-%j.out
#SBATCH --error=jobs/%x-%j.err

# No `set -e`: the generator and the config check each capture their own status,
# so a nonzero generator still reaches the closing guidance instead of aborting
# silently. Input checks below exit explicitly.
set -uo pipefail

cd "${SLURM_SUBMIT_DIR:-$(pwd)}" || exit 1
mkdir -p jobs

if [[ -f setenv.sh ]]; then
    # shellcheck disable=SC1091
    source setenv.sh
else
    echo "ERROR: setenv.sh not found in $(pwd)" >&2
    exit 1
fi

PAIRED="${PAIRED:-ACC_GEX.h5mu}"
GROUP_KEY="${GROUP_KEY:-predicted_CellType_Broad}"
OUTDIR="${REGION_SET_DIR:-01_atac/region_sets}"

# Validate BEFORE importing anything: an empty variable expanding to a bare
# path is how an earlier job wrote to the filesystem root.
if [[ -z "${PAIRED}" || ! -f "${PAIRED}" ]]; then
    echo "ERROR: paired object not found: '${PAIRED}'" >&2
    echo "       set PAIRED in setenv.sh, or run slurm/pairing.sbatch first" >&2
    exit 1
fi

CONDA_BASE="${SCPLUS_CONDA_BASE:-/home/groups/MaxsonLab/braun/miniconda3}"
if [[ ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    echo "ERROR: no conda hook at ${CONDA_BASE}/etc/profile.d/conda.sh" >&2
    exit 1
fi
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${PAIRING_ENV:-scplus-pairing}" || {
    echo "ERROR: cannot activate ${PAIRING_ENV:-scplus-pairing}" >&2; exit 1; }

echo "host        : $(hostname)"
echo "job         : ${SLURM_JOB_ID:-none}"
echo "python      : $(which python)"
echo "paired      : ${PAIRED}"
echo "group key   : ${GROUP_KEY}"
echo "out dir     : ${OUTDIR}"
echo

echo "=== [1/2] region sets from metacell DARs ==="
python 01_cistopic/region_sets_from_metacells.py \
    --h5mu "${PAIRED}" \
    --group-key "${GROUP_KEY}" \
    --out-dir "${OUTDIR}" \
    --min-log2fc "${MIN_LOG2FC:-0.25}" \
    --fdr "${DAR_FDR:-0.05}" \
    --max-regions "${MAX_REGIONS:-20000}" \
    --min-metacells "${MIN_METACELLS:-5}"
RS_RC=$?
echo "generator exit: ${RS_RC}"

echo
echo "=== [2/2] does the config validate now? ==="
bash 03_pipeline/make_config.sh --check-only
CFG_RC=$?
echo "config check exit: ${CFG_RC}"

echo
echo "=== done ==="
if [[ "${RS_RC}" -eq 0 && "${CFG_RC}" -eq 0 ]]; then
    echo "Region sets written and the config validates. Next:"
    echo "  bash 03_pipeline/make_config.sh          # write config.yaml"
    echo "  sbatch slurm/scenicplus.sbatch"
    echo
    echo "Read ${OUTDIR}/DARs_cell_type.summary.json first. A group with zero"
    echo "regions contributes no eRegulon; a capped group hit --max-regions and"
    echo "its set is the top slice by effect size, not everything significant."
else
    echo "Something failed above -- read the two exit codes. A nonzero config"
    echo "check with a zero generator exit usually means a DIFFERENT config"
    echo "input is missing (the database, or the paired object)."
fi
exit $(( RS_RC != 0 ? RS_RC : CFG_RC ))
