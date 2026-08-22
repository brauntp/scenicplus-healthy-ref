#!/usr/bin/env bash
#SBATCH --job-name=ctdb_build
#SBATCH --output=logs/ctdb_build_%j.out
#SBATCH --error=logs/ctdb_build_%j.err
#
# ---- EDIT THESE FOUR (cluster-specific; there are no safe defaults) -------- #
#SBATCH --partition=CHANGEME_PARTITION
#SBATCH --account=CHANGEME_ACCOUNT
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
# --------------------------------------------------------------------------- #
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=CHANGEME@example.edu     # uncomment + set to get email
##SBATCH --qos=CHANGEME_QOS                   # uncomment only if your site needs it
##SBATCH --constraint=CHANGEME_FEATURE        # e.g. avx512, or a node generation
#
# ============================================================================
# slurm_build_db.sbatch -- submit build_cistarget_db.sh on a SLURM cluster (ARC)
# ============================================================================
#
# QUICK START
#   1. mkdir -p logs
#   2. Discover your partition/account (see "HOW TO FIND YOUR VALUES" below).
#   3. Either edit the four #SBATCH lines above, or override at submit time --
#      command-line flags beat #SBATCH directives, so no edit is strictly
#      required:
#
#      sbatch --partition=compute --account=my_alloc \
#             --time=24:00:00 --cpus-per-task=32 \
#             --export=ALL,BED=/path/peaks.bed,GENOME_FA=/path/hg38.fa,\
#MOTIFS_DIR=/path/v10nr_clust_public/singletons,\
#OUT_PREFIX=/scratch/$USER/db/hg38_myAML_v10nr_clust,\
#CTDB_DIR=/path/src/ctdb \
#             04_db/slurm_build_db.sbatch
#
#      Or simply export the variables in your shell and submit with --export=ALL
#      (the default on most sites), which is far less error-prone:
#
#        export BED=... GENOME_FA=... MOTIFS_DIR=... OUT_PREFIX=... CTDB_DIR=...
#        sbatch -p compute -A my_alloc 04_db/slurm_build_db.sbatch
#
#   4. ALWAYS dry-run first -- it validates every input and prints the resource
#      estimate in seconds, without burning an allocation:
#        DRY_RUN=1 bash 04_db/slurm_build_db.sbatch
#
# ============================================================================
# HOW TO FIND YOUR VALUES (partition/account are site-specific)
# ============================================================================
#   Partitions you may submit to, with their limits:
#     sinfo -o '%20P %5a %10l %6D %10m %20f'
#     scontrol show partition | egrep 'PartitionName|MaxTime|DefMemPerCPU|TRES'
#   Your account(s)/allocations:
#     sacctmgr show assoc user=$USER format=account,partition,qos,maxtres
#     sshare -U
#   Largest job you can actually get scheduled:
#     sinfo -o '%20P %10c %10m'        # cores and MB per node, per partition
#   QOS names (only if your site requires one):
#     sacctmgr show qos format=name,maxwall,maxtresperuser
#
#   If --partition or --account is wrong, sbatch rejects the job immediately
#   with "Invalid partition" / "Invalid account" -- cheap to discover, so just
#   try it.
#
# ============================================================================
# WHY THESE RESOURCES
# ============================================================================
# --cpus-per-task=32
#   Passed straight through as -t to create_cistarget_motif_databases.py, which
#   runs that many PARALLEL cbust processes (each cbust is invoked with its own
#   -t 1). Scaling is close to linear in cores until I/O bound, so more cores
#   is the main lever on wall time:
#       16 cpus -> ~20-35 h ;  32 -> ~10-18 h ;  64 -> ~6-10 h
#   for 300k x 501 bp regions against the ~5.9k-motif v10nr_clust collection.
#   Do not request more cores than one node has -- this job is single-node
#   (shared-memory), it cannot span nodes.
#
# --mem=128G
#   MEMORY IS THE BINDING CONSTRAINT, not CPU. The upstream README's formula:
#       scores matrix   = 4 bytes x n_regions x n_motifs
#       rankings matrix = 4 bytes x n_regions x n_motifs   (>32768 regions)
#   and it warns you need roughly 3x the database size in RAM while running
#   (scores db + its transpose + rankings db coexist during the transpose).
#       200k regions x 5,876 motifs x 4 B =  4.7 GB/matrix -> ~14 GB   -> 32G ok
#       300k regions x 5,876 motifs x 4 B =  7.1 GB/matrix -> ~21 GB   -> 64G
#       400k regions x 5,876 motifs x 4 B =  9.4 GB/matrix -> ~28 GB   -> 64G
#       400k regions x 10,249 motifs x 4 B = 16.4 GB/matrix -> ~49 GB  -> 128G
#   128G is generous headroom for the whole 200k-400k peak range even against
#   the full 10,249-file singleton collection. If your partition caps memory
#   lower, drop to --mem=64G for <=400k peaks against the clustered collection,
#   and watch for OOM (below).
#
# --time=24:00:00
#   Sized for 32 cpus with margin. If your partition caps wall time below the
#   time you need, DO NOT just raise --time -- either raise --cpus-per-task, or
#   split the motif list with -p/--partial (see "IF IT DOES NOT FIT" below).
#
# --nodes=1 --ntasks=1
#   Required. This is one multithreaded process, not an MPI job.
#
# ============================================================================
# ADJUSTING
# ============================================================================
#   Peak set        Fewer peaks scales BOTH memory and time linearly. 200k
#                   peaks needs ~2/3 the resources of 300k.
#   Motif set       Using only clustered metaclusters instead of all 10,249
#                   singletons cuts time and memory proportionally.
#   OOM             Symptom: job dies with "Out Of Memory" / exit 137, or
#                   `seff <jobid>` shows Memory Efficiency near 100%. Fix:
#                   double --mem, or use --partial (below).
#   TIMEOUT         Symptom: "DUE TO TIME LIMIT", exit 124-ish, partial output.
#                   Fix: more --cpus-per-task, longer --time, or --partial.
#   Check usage     seff <jobid>                     # after it finishes
#                   sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,ReqMem,AllocCPUS,ExitCode
#                   sstat -j <jobid> --format=JobID,MaxRSS,AveCPU   # while running
#   Watch live      squeue -u $USER ; tail -f logs/ctdb_build_<jobid>.out
#
# IF IT DOES NOT FIT (memory or wall-time capped)
#   Build the scores database in parts using the upstream --partial flag, then
#   combine. This script does NOT orchestrate that; do it manually per the
#   repo README:
#     1. For k in 1..N:  create_cistarget_motif_databases.py ... -p k N
#          -> <prefix>.part_000k_of_000N.regions_vs_motifs.scores.feather
#     2. combine_partial_regions_or_genes_vs_motifs_scores_cistarget_dbs.py
#     3. convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py
#   Each part needs ~1/N the memory. Note that -p changes the output filenames.
#
# ============================================================================
# SCRATCH / FILESYSTEM
# ============================================================================
#   Point OUT_PREFIX at fast parallel scratch, NOT $HOME. You need ~25-30 GB
#   free for a 300k-peak build (scores + rankings + the transpose intermediate
#   + the ~1 GB region FASTA), and $HOME quotas are usually far smaller. If
#   your site purges scratch, copy the two .feather files somewhere durable as
#   soon as the job finishes.
#
set -o errexit
set -o nounset
set -o pipefail

# --------------------------------------------------------------------------- #
# REQUIRED INPUTS -- set via environment (export before sbatch, or --export=)
# or edit the defaults here.
# --------------------------------------------------------------------------- #
BED="${BED:-/CHANGEME/path/to/consensus_peaks.bed}"
GENOME_FA="${GENOME_FA:-/CHANGEME/path/to/hg38.fa}"
MOTIFS_DIR="${MOTIFS_DIR:-/CHANGEME/path/to/v10nr_clust_public/singletons}"
OUT_PREFIX="${OUT_PREFIX:-/CHANGEME/scratch/db/hg38_custom_v10nr_clust}"
CTDB_DIR="${CTDB_DIR:-/CHANGEME/path/to/src/ctdb}"

# Optional
CONDA_ENV="${CONDA_ENV:-create_cistarget_databases}"
CONDA_SH="${CONDA_SH:-}"          # e.g. /opt/conda/etc/profile.d/conda.sh
MODULES="${MODULES:-}"            # e.g. "bedtools/2.31.0 python/3.10"
CBUST="${CBUST:-cbust}"
BG_PADDING="${BG_PADDING:-0}"
SEED="${SEED:-42}"
MASK="${MASK:-0}"
DRY_RUN="${DRY_RUN:-0}"

BUILD_SCRIPT="${BUILD_SCRIPT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/build_cistarget_db.sh}"

# Threads = the cores SLURM actually gave us (not what we asked for).
THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-8}}"

# --------------------------------------------------------------------------- #
# Banner
# --------------------------------------------------------------------------- #
echo "========================================================================"
echo "cisTarget custom database build"
echo "========================================================================"
echo "job id       : ${SLURM_JOB_ID:-<not under slurm>}"
echo "node         : $(hostname)"
echo "partition    : ${SLURM_JOB_PARTITION:-n/a}"
echo "account      : ${SLURM_JOB_ACCOUNT:-n/a}"
echo "cpus alloc   : ${SLURM_CPUS_PER_TASK:-n/a}  -> passing -t ${THREADS}"
echo "mem alloc    : ${SLURM_MEM_PER_NODE:-${SLURM_MEM_PER_CPU:-n/a}} (MB if numeric)"
echo "time limit   : $(command -v squeue >/dev/null 2>&1 && [ -n "${SLURM_JOB_ID:-}" ] && squeue -h -j "${SLURM_JOB_ID}" -o '%l' 2>/dev/null || echo n/a)"
echo "started      : $(date)"
echo "workdir      : $(pwd)"
echo "========================================================================"

# --------------------------------------------------------------------------- #
# Fail fast on unedited placeholders
# --------------------------------------------------------------------------- #
BAD=0
for v in BED GENOME_FA MOTIFS_DIR OUT_PREFIX CTDB_DIR; do
    case "${!v}" in
        */CHANGEME/*|CHANGEME*)
            echo "ERROR: ${v} is still a placeholder: ${!v}" >&2
            BAD=1 ;;
    esac
done
if [ "${BAD}" -ne 0 ]; then
    cat >&2 <<'EOM'

Set the required inputs, either by editing this file or by exporting them:

  export BED=/path/to/consensus_peaks.bed
  export GENOME_FA=/path/to/hg38.fa
  export MOTIFS_DIR=/path/to/v10nr_clust_public/singletons
  export OUT_PREFIX=/scratch/$USER/db/hg38_myAML_v10nr_clust
  export CTDB_DIR=/path/to/src/ctdb
  sbatch -p <partition> -A <account> 04_db/slurm_build_db.sbatch
EOM
    exit 1
fi

# Also catch unedited SLURM directives when running under a scheduler.
if [ -n "${SLURM_JOB_ID:-}" ] && [ "${SLURM_JOB_PARTITION:-}" = "CHANGEME_PARTITION" ]; then
    echo "ERROR: --partition is still CHANGEME_PARTITION." >&2
    echo "       Edit the #SBATCH line or pass --partition=<name> to sbatch." >&2
    echo "       List partitions with: sinfo -o '%20P %5a %10l %6D %10m'" >&2
    exit 1
fi

[ -f "${BUILD_SCRIPT}" ] || {
    echo "ERROR: build script not found: ${BUILD_SCRIPT}" >&2
    echo "       set BUILD_SCRIPT=/path/to/build_cistarget_db.sh" >&2
    exit 1
}

mkdir -p logs

# --------------------------------------------------------------------------- #
# Environment
# --------------------------------------------------------------------------- #
if [ -n "${MODULES}" ]; then
    echo "[env] module load ${MODULES}"
    # shellcheck disable=SC1091
    if command -v module >/dev/null 2>&1; then
        # shellcheck disable=SC2086
        module load ${MODULES}
    else
        echo "[env] WARNING: 'module' not available; skipping" >&2
    fi
fi

if [ -n "${CONDA_ENV}" ]; then
    if [ -n "${CONDA_SH}" ] && [ -f "${CONDA_SH}" ]; then
        # shellcheck disable=SC1090
        source "${CONDA_SH}"
    else
        for c in "${HOME}/miniconda3" "${HOME}/miniforge3" "${HOME}/anaconda3" \
                 /opt/conda /usr/local/conda; do
            if [ -f "${c}/etc/profile.d/conda.sh" ]; then
                # shellcheck disable=SC1091
                source "${c}/etc/profile.d/conda.sh"
                break
            fi
        done
    fi
    if command -v conda >/dev/null 2>&1; then
        echo "[env] conda activate ${CONDA_ENV}"
        conda activate "${CONDA_ENV}" || {
            echo "ERROR: could not activate conda env '${CONDA_ENV}'." >&2
            echo "       Create it (per the upstream README):" >&2
            echo "         conda create -n ${CONDA_ENV} 'python=3.10' 'numpy=1.21' \\" >&2
            echo "             'pandas>=1.4.1' 'pyarrow>=7.0.0' 'numba>=0.55.1' \\" >&2
            echo "             python-flatbuffers" >&2
            echo "       Or set CONDA_SH=/path/to/etc/profile.d/conda.sh," >&2
            echo "       or CONDA_ENV='' if your modules already provide Python." >&2
            exit 1
        }
    else
        echo "[env] WARNING: conda not found; relying on whatever Python is on \$PATH" >&2
    fi
fi

echo "[env] python : $(command -v python3 || echo MISSING)"
echo "[env] bedtools: $(command -v bedtools || echo MISSING)"
echo "[env] cbust  : $(command -v "${CBUST}" || echo "MISSING (${CBUST})")"

# Keep BLAS/OpenMP from oversubscribing: parallelism here is process-level
# (one cbust per thread slot), so nested threading only causes contention.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Node-local scratch for temporaries, if the site provides it.
if [ -n "${SLURM_JOB_ID:-}" ] && [ -d "/tmp" ]; then
    export TMPDIR="${TMPDIR:-/tmp}"
fi
echo "[env] TMPDIR : ${TMPDIR:-<unset>}"

# --------------------------------------------------------------------------- #
# Run
# --------------------------------------------------------------------------- #
ARGS=(
    --bed        "${BED}"
    --genome-fa  "${GENOME_FA}"
    --motifs-dir "${MOTIFS_DIR}"
    --out-prefix "${OUT_PREFIX}"
    --ctdb-dir   "${CTDB_DIR}"
    --cbust      "${CBUST}"
    --threads    "${THREADS}"
    --bgpadding  "${BG_PADDING}"
    --seed       "${SEED}"
)
[ "${MASK}" = "1" ] && ARGS+=(--mask)
[ "${DRY_RUN}" = "1" ] && ARGS+=(--dry-run)

echo
echo "[run] bash ${BUILD_SCRIPT} $(printf '%s ' "${ARGS[@]}")"
echo

START=$(date +%s)
set +o errexit
bash "${BUILD_SCRIPT}" "${ARGS[@]}"
RC=$?
set -o errexit
END=$(date +%s)

echo
echo "========================================================================"
echo "finished     : $(date)"
echo "elapsed      : $(( (END-START)/3600 ))h $(( ((END-START)%3600)/60 ))m"
echo "exit code    : ${RC}"
if [ -n "${SLURM_JOB_ID:-}" ] && command -v seff >/dev/null 2>&1; then
    echo "--- seff ${SLURM_JOB_ID} ---"
    seff "${SLURM_JOB_ID}" 2>/dev/null || true
elif [ -n "${SLURM_JOB_ID:-}" ] && command -v sacct >/dev/null 2>&1; then
    sacct -j "${SLURM_JOB_ID}" \
        --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem,AllocCPUS,ExitCode 2>/dev/null || true
fi
echo "========================================================================"

if [ "${RC}" -ne 0 ]; then
    cat >&2 <<EOM

BUILD FAILED (exit ${RC}). Triage:
  * exit 137 / "Out Of Memory"  -> raise --mem (see the memory table above)
                                   or split the motif list with -p/--partial
  * "DUE TO TIME LIMIT"         -> raise --time or --cpus-per-task
  * cbust not found             -> the build script prints the exact fix;
                                   wget https://resources.aertslab.org/cistarget/programs/cbust
                                   chmod a+x cbust
  * ModuleNotFoundError:
    'cistarget_db'              -> incomplete CTDB_DIR clone, or PYTHONSAFEPATH
                                   is set in this environment (the build script
                                   sets PYTHONPATH to work around it)
  * bedtools getfasta errors    -> contig-name mismatch (chr1 vs 1) or wrong
                                   genome build for these peaks
  Check MaxRSS vs ReqMem:  sacct -j ${SLURM_JOB_ID:-<jobid>} --format=JobID,MaxRSS,ReqMem,State,Elapsed
EOM
    exit "${RC}"
fi

cat <<EOM

NEXT STEP -- verify the database you just built actually covers your peaks:

  python 04_db/peak_overlap_audit.py \\
      --peaks '${BED}' \\
      --db    '${OUT_PREFIX}.regions_vs_motifs.rankings.feather' \\
      --out-prefix '${OUT_PREFIX}.audit'

A custom DB built from these same peaks should report ~0% dropped. If it does
not, the region IDs in the database do not match your ATAC var_names and
SCENIC+ will silently discard peaks.
EOM
