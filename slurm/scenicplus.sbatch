#!/usr/bin/env bash
# =============================================================================
# SCENIC+ v1.0a2 -- SLURM submit script (ARC)
# =============================================================================
# Submit:
#     sbatch --account=<ARC_ACCOUNT> --partition=<ARC_PARTITION> \
#            slurm/scenicplus.sbatch <path/to/config.yaml>
#
# Or fill the placeholders below and just:
#     sbatch slurm/scenicplus.sbatch 03_pipeline/config.yaml
#
# Anything passed on the sbatch command line overrides the #SBATCH directives,
# so the placeholders are a default, not a lock-in.
#
# This runs the ENTIRE pipeline as ONE fat job rather than having snakemake
# submit a job per rule. Rationale: the pipeline is a mostly-linear chain of
# 13 jobs (see the dry-run plan), the two expensive steps want all the cores on
# one node anyway, and the cisTarget rankings database is loaded from scratch by
# every rule that touches it -- so per-rule submission would pay that multi-GB
# read repeatedly and spend more time queueing than computing. If you do want
# per-rule submission later, that is a snakemake executor plugin
# (snakemake-executor-plugin-slurm) plus a profile, not an edit to this file.
# =============================================================================

#SBATCH --job-name=scenicplus
#SBATCH --partition=batch

# --- one node, all cores in one task -----------------------------------------
# Everything parallel in SCENIC+ is shared-memory within a single process tree:
# joblib (region_to_gene, tf_to_gene), the arboreto/dask local cluster, and
# pycistarget's multiprocessing pool. None of it is MPI, so asking for more than
# one node buys nothing and just makes the job harder to schedule.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<N_CPU>          # e.g. 32 or 48; feeds --cores AND config n_cpu

# --- memory ------------------------------------------------------------------
# Ask for a lot, and understand which term dominates. Two independent reasons:
#
# (1) THE cisTARGET DATABASE IS MEMORY-RESIDENT. input_data.ctx_db_fname is a
#     regions-vs-motifs RANKINGS feather that pycistarget reads whole into a
#     dense in-memory table -- there is no windowed or memory-mapped path. For a
#     human v10 SCREEN database that is order 1.1M regions x ~34k motifs of
#     int32, i.e. tens of GB resident for as long as the rule runs. It is read
#     by motif_enrichment_cistarget AND again by BOTH eGRN rules (the Snakefile
#     passes ctx_db_fname as --ranking_db_fname to eGRN_direct and
#     eGRN_extended), so this is not a one-off spike at the start. The DEM
#     scores database is loaded the same way in motif_enrichment_dem.
#     This term is set by your database, not by your data, and it is usually
#     the binding constraint.
#
# (2) REGION-TO-GENE GBM OVER N METACELLS. region_to_gene fits one gradient
#     boosting model per gene over that gene's candidate regions
#     (enhancer_to_gene.py:214-227). The work is spread across n_cpu joblib
#     workers, and each worker holds its own slice of the accessibility and
#     expression matrices plus the sklearn ensemble it is currently fitting.
#     Peak RSS therefore scales roughly with
#         n_cpu x (metacells x regions_in_search_space_for_that_gene)
#     -- so raising --cpus-per-task raises the memory floor too. Our GLUE
#     metacells keep N small (hundreds, not the ~10^4-10^5 of single cells),
#     which is what makes this affordable at all; a run on raw cells would not
#     fit. joblib spills its memmaps to TMPDIR (below), which converts some of
#     this into disk rather than RAM, but only some.
#
# Start at 240G. If the job dies with no python traceback -- just "Killed" or a
# bare nonzero exit -- that is the OOM killer, not a bug: check
# `seff $SLURM_JOB_ID` / `sacct -j $SLURM_JOB_ID -o MaxRSS` and raise.
#SBATCH --mem=<MEM>                      # e.g. 240G

# --- walltime ----------------------------------------------------------------
# region_to_gene and the two eGRN/GSEA rules dominate. With a few hundred
# metacells, a 150 kb search space and n_cpu in the 32-48 range, budget on the
# order of 12-24 h and set the ceiling generously -- a walltime kill mid-GBM
# wastes the whole run, and `run_pipeline.sh` can resume from completed outputs
# (that is part of why it pins --rerun-triggers mtime).
#SBATCH --time=<HH:MM:SS>                # e.g. 24:00:00

# --- logs --------------------------------------------------------------------
# %x = job name, %j = job id. Directory must EXIST before submit or slurm
# silently discards output.
#SBATCH --output=jobs/%x-%j.out
#SBATCH --error=jobs/%x-%j.err

# --- notifications (optional) ------------------------------------------------
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=<EMAIL>

set -Eeuo pipefail

CONFIG="${1:?usage: sbatch slurm/scenicplus.sbatch <config.yaml>}"

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$REPO_ROOT"

echo "==============================================================="
echo "job            : ${SLURM_JOB_NAME:-scenicplus} (${SLURM_JOB_ID:-nojob})"
echo "node           : $(hostname)"
echo "submit dir     : ${REPO_ROOT}"
echo "config         : ${CONFIG}"
echo "cpus-per-task  : ${SLURM_CPUS_PER_TASK:-unset}"
echo "mem            : ${SLURM_MEM_PER_NODE:-unset} MB"
echo "started        : $(date -Is)"
echo "==============================================================="

# --- conda env ---------------------------------------------------------------
# No rule in the Snakefile declares a `conda:` directive, so snakemake does not
# manage environments -- every rule shells out to whatever `scenicplus` is on
# PATH. The env therefore MUST be active in this shell before snakemake starts.
# `conda activate` needs the hook; `conda run` does not work here because we
# need the env live for the child shells snakemake spawns.
CONDA_BASE="<CONDA_BASE>"                # e.g. /home/<user>/miniforge3  (output of `conda info --base`)
CONDA_ENV="<CONDA_ENV_NAME>"             # e.g. scenicplus (see 03_pipeline/environment.yml)

# ARC-style module setup, if your site needs it. Delete if not applicable.
# module purge
# module load <MODULE_FOR_CONDA>

if [[ ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    echo "ERROR: no conda hook at ${CONDA_BASE}/etc/profile.d/conda.sh" >&2
    echo "       Set CONDA_BASE to the output of \`conda info --base\`." >&2
    exit 1
fi
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}" || { echo "ERROR: cannot activate env '${CONDA_ENV}'" >&2; exit 1; }
echo "env            : $(which python)"
command -v scenicplus >/dev/null || { echo "ERROR: 'scenicplus' not on PATH after activating ${CONDA_ENV}" >&2; exit 1; }

# --- TMPDIR on node-local scratch -------------------------------------------
# SCENIC+ writes LARGE temp files and the pipeline is explicitly built around
# that: params_general.temp_dir is passed as --temp_dir to cisTarget, DEM,
# tf_to_gene, region_to_gene and both eGRN rules, and inside region_to_gene it
# becomes joblib's `temp_folder` (enhancer_to_gene.py:221), where joblib
# memory-maps the per-worker accessibility/expression slices. On a shared
# parallel filesystem that traffic is both slow and antisocial; in $HOME it will
# hit a quota mid-run and kill a job hours in.
#
# Point everything at node-local scratch. ARC provides per-job scratch under one
# of the paths below -- CHECK WHICH ONE YOUR SITE USES and set SCRATCH_ROOT
# accordingly. If none exists, the fallback creates a job-scoped directory under
# the repo, which works but gives up the locality benefit.
SCRATCH_ROOT="<NODE_SCRATCH_ROOT>"       # e.g. /scratch/$SLURM_JOB_ID  or  /local/$USER  or  $TMPDIR
if [[ -z "${SCRATCH_ROOT}" || "${SCRATCH_ROOT}" == "<NODE_SCRATCH_ROOT>" || ! -d "${SCRATCH_ROOT}" ]]; then
    for cand in "/scratch/${SLURM_JOB_ID:-}" "/lscratch/${SLURM_JOB_ID:-}" "/local/${USER}" "${TMPDIR:-}"; do
        if [[ -n "$cand" && -d "$cand" && -w "$cand" ]]; then SCRATCH_ROOT="$cand"; break; fi
    done
fi
if [[ -z "${SCRATCH_ROOT}" || ! -d "${SCRATCH_ROOT}" ]]; then
    SCRATCH_ROOT="${REPO_ROOT}/.scratch"
    echo "WARNING: no node-local scratch found; falling back to ${SCRATCH_ROOT}" >&2
    echo "         This is on shared storage -- expect slower region_to_gene." >&2
fi

export TMPDIR="${SCRATCH_ROOT}/scenicplus.${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"
# Keep every temp-dir convention pointing at the same place: python's
# tempfile, joblib, numba's cache, and MPLCONFIGDIR (matplotlib writes a font
# cache on first import and warns loudly if HOME is read-only on the node).
export TMP="$TMPDIR" TEMP="$TMPDIR"
export JOBLIB_TEMP_FOLDER="$TMPDIR"
export NUMBA_CACHE_DIR="$TMPDIR/numba"
export MPLCONFIGDIR="$TMPDIR/mpl"
mkdir -p "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR"

echo "TMPDIR         : ${TMPDIR}"
df -h "$TMPDIR" | tail -1

# Clean up node-local scratch on exit -- SLURM does not always do it, and a
# leftover multi-GB temp tree will make the NEXT job on this node fail on space.
# NOTE: this deliberately removes only our own $TMPDIR subtree, never
# SCRATCH_ROOT itself.
cleanup() {
    rc=$?
    echo "--- cleanup (exit ${rc}) ---"
    if [[ -n "${TMPDIR:-}" && "$TMPDIR" == *"scenicplus."* && -d "$TMPDIR" ]]; then
        du -sh "$TMPDIR" 2>/dev/null || true
        rm -rf "$TMPDIR"
    fi
    echo "finished       : $(date -Is)"
    exit $rc
}
trap cleanup EXIT

# --- thread oversubscription guard ------------------------------------------
# joblib/arboreto fan out to n_cpu WORKER PROCESSES, and each worker's numpy
# will by default start its own OpenMP/BLAS thread pool sized to the whole node.
# That is n_cpu^2 threads fighting for n_cpu cores, which on a 48-core node
# turns a compute-bound run into a scheduler-bound one. Pin the inner pools to
# 1 thread and let joblib own the parallelism.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# Mallet is a JVM app and defaults to a small heap; pycisTopic's LDA step needs
# far more. Harmless if this run does no topic modelling.
export MALLET_MEMORY="${MALLET_MEMORY:-64g}"

# --- network note ------------------------------------------------------------
# rule download_genome_annotations calls pybiomart against
# params_data_preparation.biomart_host, so THIS COMPUTE NODE needs outbound
# HTTP. On clusters where compute nodes have no internet, run that one rule on
# a login node first and let --rerun-triggers mtime keep the result:
#     03_pipeline/run_pipeline.sh --config <cfg> -- --until download_genome_annotations
# then submit this job.

# --- run ---------------------------------------------------------------------
# Fall back to a number, not a placeholder: "<N_CPU>" would reach snakemake
# --cores as a literal string and fail with an unhelpful parse error.
CORES="${SLURM_CPUS_PER_TASK:-8}"

# run_pipeline.sh does the config/input validation, gates on validate_h5mu.py,
# writes the resolved config with n_cpu=CORES, and invokes snakemake with
# --rerun-triggers mtime. Do not call snakemake directly here -- the mtime
# trigger is what stops snakemake from deciding our hand-built paired MuData is
# stale and regenerating it via prepare_GEX_ACC (which would replace correctly
# paired GLUE metacells with randomly paired ones, silently).
srun --cpu-bind=cores \
    03_pipeline/run_pipeline.sh \
        --config "$CONFIG" \
        --cores "$CORES" \
        "${@:2}"

echo "==============================================================="
echo "pipeline completed"
echo "==============================================================="
