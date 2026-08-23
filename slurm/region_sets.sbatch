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
# Memory is set by --block, but the peak is 5-6.5x the block's own size, not 1x.
# Job 10714844 was OOM-KILLED (exit 137) at 16G because I sized this from the
# block alone: a 20,000-region block is 1.89 GB, but ranks + a stable argsort +
# the gathered sorted copy took it to 9.59 GB MEASURED, and scipy's rankdata
# returned float64 (2x a float32 block) with an int64 argsort inside.
#
# Fixed on two fronts: ranks and tie terms now come from ONE argsort and stay in
# the block's dtype, and --block defaults to 4,000 (0.38 GB block, 2.45 GB peak
# measured, a 6.5x factor; the 20,000 block measured 5.1x -- the ratio falls
# with block size because the per-column python loop's temporaries do not scale).
# 32G gives ~13x headroom over 2.45 GB, deliberately generous: this is the second
# sizing attempt and the first was killed.
#
# Time is MEASURED: 2.14 ms per region for ranks+ties+all 24 groups at
# 25,323 metacells -> ~14 min for 393,832 regions IN ISOLATION. MEASURED on the
# real object: 42m43s (job 10714858), 3x my projection -- the per-block timing
# omitted reading each block from the backed object and the BH/selection pass.
# 2h still holds, at a 2.8x margin rather than the 8x I claimed
# (the combined routine is ~30% slower per block than rankdata alone, but it
# replaces rankdata AND a second np.sort). 2h is an ~8x margin, which matters
# because BEDs are written only at the end, so a walltime kill loses everything.
#
# The first version of this script called scipy.mannwhitneyu once per group,
# which re-ranks the same block 24 times: 1.84 ms per region per group, i.e.
# 4.8 h against the 4 h walltime I had set. It would have been killed with
# nothing written. Ranking once per block and deriving each group's U from its
# rank sum is 20x faster and matches scipy to 5e-08.
# =============================================================================
#SBATCH --job-name=region_sets
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
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
# Keep the effect-size/q-value dump by default. 63 MB uncompressed buys every
# future threshold question offline; without it, changing --min-log2fc costs
# another 43-minute pass over the matrix. See docs/REGION_SETS.md -- the first
# real run capped 20 of 21 groups, so the threshold WILL need revisiting.
DUMP_STATS="${DUMP_STATS:-01_atac/dar_stats.npz}"
# The pairing job's diagnostics carry independent_metacell_equiv, the honest
# per-group denominator. Without it the script falls back to
# metacells/oversample, which is APPROXIMATE and skips fewer groups.
DIAGNOSTICS="${DIAGNOSTICS:-}"
if [[ -z "${DIAGNOSTICS}" ]]; then
    for cand in pairing_diagnostics.csv pairing_diagnostics/pairing_diagnostics.csv; do
        [[ -f "$cand" ]] && { DIAGNOSTICS="$cand"; break; }
    done
fi
if [[ -n "${DIAGNOSTICS}" ]]; then
    echo "diagnostics : ${DIAGNOSTICS} (authoritative independence counts)"
else
    echo "diagnostics : NONE FOUND -- falling back to metacells/oversample."
    echo "              That is approximate and skips FEWER thin groups."
    echo "              Set DIAGNOSTICS=<path> if the CSV is elsewhere."
fi

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

echo "=== [1/3] self-test (seconds) ==="
python 01_cistopic/region_sets_from_metacells.py --self-test
ST_RC=$?
if [[ "${ST_RC}" -ne 0 ]]; then
    echo "ERROR: the script's own worked example disagrees with its arithmetic." >&2
    echo "       Refusing to run: a drifted example means the documented filter" >&2
    echo "       behaviour is not the implemented one." >&2
    exit 1
fi

echo
echo "=== [2/3] region sets from metacell DARs ==="
python 01_cistopic/region_sets_from_metacells.py \
    --h5mu "${PAIRED}" \
    --group-key "${GROUP_KEY}" \
    --out-dir "${OUTDIR}" \
    --min-log2fc "${MIN_LOG2FC:-0.25}" \
    --fdr "${DAR_FDR:-0.05}" \
    --max-regions "${MAX_REGIONS:-20000}" \
    --min-independent "${MIN_INDEPENDENT:-5}" \
    ${DUMP_STATS:+--dump-stats "${DUMP_STATS}"} \
    ${DIAGNOSTICS:+--diagnostics "${DIAGNOSTICS}"}
RS_RC=$?
echo "generator exit: ${RS_RC}"

echo
echo "=== [3/3] does the config validate now? ==="
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
