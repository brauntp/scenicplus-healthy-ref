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
# 16, not 32+: worker count scales with cores and each joblib worker adds to
# the 43.2 GB parent, so cores and --mem trade off. Raise with --mem together.
# Override at submit: sbatch --cpus-per-task=32 --mem=96G slurm/scenicplus.sbatch
#SBATCH --cpus-per-task=16               # feeds snakemake --cores AND config n_cpu

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
# `seff $SLURM_JOB_ID` / `sacct -j $SLURM_JOB_ID -o JobID,MaxRSS` and raise --
# reading the `.0` (run_pipeline) step, since `.batch` is the empty shell.
# MEASURED, not projected: 54G. The paired object is 40.4 GB dense and the
# probe measured a 1.07x overhead through both .to_df() calls (peak 43.2 GB,
# +25% headroom). The QC job independently hit 43.2 GB MaxRSS on the same
# object. An earlier projection said 96G on the assumption that .to_df()
# copies; it does not -- the frame shares memory with X.
#
# 54G was WRONG and OOM-killed the run at motif_enrichment_cistarget. That
# figure came from the region_to_gene probe (43.2 GB measured, +25%) -- a
# DIFFERENT rule. The header above says the database term is usually the binding
# constraint, and then I sized the job to the other one.
#
# The database term, computed from what the code actually loads:
#   pycistarget's cisTargetDatabase.load_db() with region_sets supplied calls
#   db.load(subset), NOT load_full(), and ctxcore pushes the column list into
#   pyarrow read_table(columns=...). So it reads only the database regions that
#   MATCH our peaks -- not all 1.84M. Good news, but at our geometry the subset
#   is still most of the file:
#
#     393,832 matched regions x 32,765 motifs x int32   = 48.1 GB
#     subset_to_pandas builds a pyarrow Table then a DataFrame, both briefly
#     live                                              = ~96 GB peak
#
#   (load_full() on the whole database would be ~224 GB. That path is not taken
#   here, and the peak_overlap_audit showed near-total query-side matching, so
#   assume most peaks match rather than hoping for a small subset.)
#
# CORRECTION, after 128G was ALSO OOM-killed at motif_enrichment_cistarget:
# the figure above is one database load. cisTarget does not do one.
#
#   run_motif_enrichment_cistarget() calls
#     joblib.Parallel(n_jobs=n_cpu)(delayed(_run_cistarget_single_region_set)(
#         ..., cistarget_db_fname=..., ...) for key in region_set_dict)
#   and the worker builds its OWN cisTargetDatabase from that filename. joblib's
#   default backend is loky -- separate PROCESSES -- so n_cpu region sets are in
#   flight at once, each holding its own slice.
#
# So the peak is the sum over the n_cpu LARGEST region sets, and
# --cpus-per-task is a memory lever for this rule, not just a speed knob.
# DEM survived the same allocation because its per-worker slice is a foreground
# set plus a CAPPED background, not a whole cell type's DAR set.
#
# Rules DO still serialise with each other (each declares threads = n_cpu;
# verified directly against snakemake), so --mem covers the largest single rule
# -- but for cisTarget that rule's own peak is n_cpu-fold.
#
# Size it from the real region sets before submitting:
#   python 03_pipeline/size_cistarget_memory.py \
#       --region-set-folder 01_atac/region_sets --db <rankings.feather>
#
# 128G is retained because it is right for every OTHER rule. For cisTarget,
# either lower --cpus-per-task or make the region sets equal-sized with
# choose_dar_threshold.py --top-n 5000 --write, which fits all 21 sets in 128G
# with every worker running. See docs/MEMORY.md.
# See docs/MEMORY.md and 03_pipeline/probe_region_to_gene_memory.py.
# 192G, raised from 128G after AUCell_direct and AUCell_extended both failed
# with OSError: [Errno 12] Cannot allocate memory at os.fork().
#
# That is NOT an OOM kill -- there is no oom_kill event and the process produced
# a full traceback. The allocation was REFUSED because the parent was already
# too large for the cgroup to admit another address space.
#
# Where the ~118 GB goes, in order (see docs/MEMORY.md):
#   mudata.read(ACC_GEX.h5mu)                          40.4 GB
#   rank_data() on BOTH modalities, both held at once  40.4 GB
#     (`ranking = np.zeros_like(mtx)` inherits float32 from the source, so
#      this is NOT the uint32 the ranking values are cast to)
#   aucell4r RawArray(c_uint32) for the region matrix  37.2 GB
#
# CORRECTION (after the 192G run failed IDENTICALLY): raising --mem does not fix
# it and --cpus-per-task IS the lever. The evidence:
#
#   1. It failed on gex_AUC -- the GENE matrix (3.2 GB), the small one. The
#      37 GB region RawArray had not been allocated yet.
#   2. 128G -> 192G produced a byte-identical traceback.
#   3. cnode-08-26 is on `batch`, which has 400 GB+ nodes.
#
# fork() duplicates page tables and, under strict overcommit accounting, the
# kernel charges the child's worst case rather than trusting copy-on-write. A
# parent holding ~81 GB (paired object + BOTH rank matrices) forking 16 children
# implies a commit charge no --mem on this cluster can satisfy.
#
# aucell4r takes a SERIAL branch at num_workers == 1 -- verified in the pySCENIC
# source: no RawArray, no Process, no fork. That is the fix, and it is why
# AUCell must run as a --target subset at --cpus-per-task 1 rather than by
# raising memory here. See RUNBOOK.md.
#SBATCH --mem=192G

# --- walltime ----------------------------------------------------------------
# region_to_gene and the two eGRN/GSEA rules dominate. With a few hundred
# metacells, a 150 kb search space and n_cpu in the 32-48 range, budget on the
# order of 12-24 h and set the ceiling generously -- a walltime kill mid-GBM
# wastes the whole run, and `run_pipeline.sh` can resume from completed outputs
# (that is part of why it pins --rerun-triggers mtime).
# 24h. The pairing job took 19 min and QC 13 min, but those are single steps;
# this DAG runs cisTarget and DEM over 1.8M database regions plus a GBM fit per
# gene. A walltime kill mid-GBM loses the whole rule, so err long.
#SBATCH --time=24:00:00

# --- logs --------------------------------------------------------------------
# %x = job name, %j = job id. Directory must EXIST before submit or slurm
# silently discards output.
#SBATCH --output=jobs/%x-%j.out
#SBATCH --error=jobs/%x-%j.err

# --- notifications (optional) ------------------------------------------------
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=<EMAIL>

# No `set -e` past the input checks: snakemake's exit status is captured so a
# failed run still reaches the epilogue, which is where the "what to read next"
# guidance lives. Same fix as pairing.sbatch, qc_paired.sbatch and
# region_sets.sbatch -- this script was the one that still had it.
set -Euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$REPO_ROOT" || exit 1

# Default to what make_config.sh writes. Requiring the path as an argument made
# job 10715268 die at submit-time-plus-queue-time on a usage error, for a file
# whose location is already known.
# $1 is the config path ONLY when it does not look like a flag. Otherwise every
# argument is a pass-through and the config falls back to the default.
#
# Without this guard, `sbatch slurm/scenicplus.sbatch --target ctx_results.hdf5`
# consumed "--target" as the config path and forwarded a bare "ctx_results.hdf5"
# -- so the config was wrong AND the flag was lost, with no error until the
# path check failed on a nonsense filename.
if [[ -n "${1:-}" && "${1}" != -* ]]; then
    CONFIG="$1"
    PASSTHRU=("${@:2}")
else
    CONFIG="03_pipeline/config.yaml"
    PASSTHRU=("$@")
fi
if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config not found: $CONFIG" >&2
    echo "       Write it first:  bash 03_pipeline/make_config.sh" >&2
    echo "       Or pass one:     sbatch slurm/scenicplus.sbatch <config.yaml>" >&2
    exit 1
fi
echo "config         : $CONFIG"

# A batch job does not inherit what you exported at a prompt, so read the
# project's settings from the file. Same pattern as pairing.sbatch.
if [[ -f setenv.sh ]]; then
    # shellcheck disable=SC1091
    source setenv.sh
else
    echo "ERROR: setenv.sh not found in $REPO_ROOT" >&2
    exit 1
fi

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
# House style in cutTag-pipeline is to `conda activate` BEFORE sbatch and let
# the job inherit it (sbatch exports the environment by default). pairing.sbatch
# and qc_paired.sbatch -- the two scripts that have actually run on this cluster
# -- honour that first and fall back to activating here. This does the same, so
# the three behave identically.
#
# NOTE: setenv.sh does NOT define these. An earlier comment here claimed they
# were "read from setenv.sh (sourced above)"; they are not, and the fallback
# below is a hardcoded default, not a project setting. Override by exporting
# CONDA_ENV or CONDA_BASE before sbatch.
CONDA_ENV="${CONDA_ENV:-scenicplus}"     # see 03_pipeline/environment.yml

if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    echo "using inherited conda env: $CONDA_DEFAULT_ENV"
    # The activate.d hook fires on `conda activate`, so an INHERITED env brings
    # LD_LIBRARY_PATH along only if it was set at submit time. Assert it rather
    # than assume: without the env's own lib dir ahead of /lib64, extensions
    # whose RPATH does not reach $CONDA_PREFIX/lib resolve the SYSTEM
    # libstdc++, which is too old for the CXXABI they need. Measured on this
    # cluster: 7 of 7 affected modules fail as-is, 0 of 7 with the lib dir
    # prepended, and every cisTarget/DEM/region-to-gene rule dies at import --
    # hours into a run, not here.
    case ":${LD_LIBRARY_PATH:-}:" in
        *":${CONDA_PREFIX}/lib:"*) ;;
        *) export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
           echo "prepended \$CONDA_PREFIX/lib to LD_LIBRARY_PATH (inherited env" \
                "did not carry it)" ;;
    esac
else
    CONDA_BASE="${CONDA_BASE:-/home/groups/MaxsonLab/braun/miniconda3}"
    if [[ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
        # shellcheck disable=SC1091
        source "${CONDA_BASE}/etc/profile.d/conda.sh"
        if ! conda activate "${CONDA_ENV}"; then
            echo "ERROR: cannot activate env '${CONDA_ENV}'." >&2
            echo >&2
            if ! conda env list 2>/dev/null | grep -qE "^${CONDA_ENV}\s"; then
                echo "       It does not exist. This is the ~900-package env" >&2
                echo "       with source builds -- allow an hour or more, and" >&2
                echo "       run it under screen/tmux, not a bare login shell:" >&2
                echo >&2
                echo "         screen -S spenv" >&2
                echo "         bash 03_pipeline/create_env.sh" >&2
                echo >&2
                echo "       Use create_env.sh, NOT mamba env create directly:" >&2
                echo "       it exports PIP_CONSTRAINT so pip's isolated build" >&2
                echo "       env gets setuptools<81. Without that the build" >&2
                echo "       dies on pybedtools AFTER the conda solve, because" >&2
                echo "       setuptools 84 dropped pkg_resources. The wrapper" >&2
                echo "       runs preflight_env.sh first." >&2
                echo >&2
                echo "       preflight_env.sh checks the platform, solver," >&2
                echo "       python pin, pandas placement, git pins and disk" >&2
                echo "       BEFORE the solve, so a bad spec fails in seconds" >&2
                echo "       rather than an hour in." >&2
            else
                echo "       It exists but will not activate -- check" >&2
                echo "       ${CONDA_BASE} is the right conda root." >&2
            fi
            exit 1
        fi
        echo "activated conda env: ${CONDA_ENV} (from ${CONDA_BASE})"
    else
        echo "ERROR: '${CONDA_ENV}' is not active and no conda.sh under" >&2
        echo "       CONDA_BASE=${CONDA_BASE}." >&2
        echo "       Either 'conda activate ${CONDA_ENV}' before sbatch (house" >&2
        echo "       style), or set CONDA_BASE when submitting." >&2
        exit 1
    fi
fi
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
        "${PASSTHRU[@]+"${PASSTHRU[@]}"}"
RC=$?

echo "==============================================================="
if [[ "$RC" -eq 0 ]]; then
    echo "pipeline completed"
    echo
    echo "Read first:"
    echo "  - the eRegulon table: which groups yielded regulons at all. Three"
    echo "    groups were dropped before region sets (Pro-Monocyte, cDC, Late"
    echo "    GMP) and docs/QC_RESULT.md flags both erythroid groups as the"
    echo "    worst-aligned in the pairing -- thin erythroid results are the"
    echo "    cross-modal gap surfacing downstream, not a pipeline failure."
    # Read the directive rather than naming a number: this line said "the 54G
    # here" long after --mem had moved to 128G and then 192G, so it advised
    # checking a figure the job never requested.
    # NOT escaped: these must EXPAND. An earlier version wrote \$SLURM_JOB_ID,
    # which printed the literal text for the user to retype -- the same defect
    # already fixed once in fix_cxxabi.sh's closing advice. My local test used
    # the unescaped form, so it passed while the shipped file did not expand.
    # MaxRSS lives on the `.0` (run_pipeline) srun step, NOT on `.batch` --
    # this script launches the pipeline under srun, so the batch shell itself
    # reads ~0.01G. Grep for the step that did the work, or the number looks
    # like ~10 MB.
    echo "  - sacct -j ${SLURM_JOB_ID:-<job-id>} --units=G \\"
    echo "        --format=JobID%-16,JobName%-14,State,Elapsed,ReqMem,MaxRSS,AllocCPUS \\"
    echo "      | grep -E 'JobID|run_pipeline|^-'"
    echo "    (the .0 step carries MaxRSS; .batch is the empty shell)"
    echo "    against the --mem this job actually requested"
    echo "    ($(( ${SLURM_MEM_PER_NODE:-0} / 1024 ))G). Every figure in"
    echo "    docs/MEMORY.md for the rules that ran here is derived from source"
    echo "    reading, not measured at this scale -- record the real MaxRSS."
else
    echo "pipeline FAILED (exit $RC)"
    echo
    echo "snakemake reports the failing RULE and its log path; read that before"
    echo "resubmitting. A rule that fails on a missing input usually means the"
    echo "config points somewhere the DAG does not expect --"
    echo "  bash 03_pipeline/make_config.sh --check-only"
    echo "revalidates every path in seconds."
fi
echo "==============================================================="
exit "$RC"
