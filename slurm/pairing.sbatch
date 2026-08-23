#!/usr/bin/env bash
# =============================================================================
# GLUE-anchored metacell pairing -- SLURM submit script (ARC)
# =============================================================================
# Submit:
#     mkdir -p jobs                      # once; SLURM will not create it
#     sbatch slurm/pairing.sbatch
#
# Reads rna.h5ad + atac.h5ad + atac_labels.tsv, writes ACC_GEX.h5mu.
# NOT a login-node job: it holds ~10 GB of sparse input and allocates ~40 GB of
# dense metacell matrices.
#
# ORDER OF OPERATIONS -- the env comes FIRST. --dry-run imports numpy, scipy
# and anndata like the real run does, so it needs the env too; it is cheap in
# TIME and MEMORY, not in dependencies.
#
#   1.  mamba env create -f 03_pipeline/pairing_env.yml      # once, ~1 min
#   2.  conda activate scplus-pairing
#   3.  python 02_pair/aggregate_atac_sparse.py ... --dry-run    # login node OK
#   4.  sbatch --account=... --partition=... slurm/pairing.sbatch
#
# Step 3 prints the exact metacell count, the dense footprint and a --mem
# recommendation read from your actual data. If it disagrees with the #SBATCH
# --mem below, trust the dry run.
# =============================================================================

# Partition and log layout follow maxsonBraunLab/cutTag-pipeline house style:
# partition "batch", logs under jobs/<name>_%j.{out,err}. No --account: neither
# cutTag-pipeline's wrappers nor the lab's snakemake SLURM profile set one, and
# the default account (maxsonlab) applies.
#     sbatch slurm/pairing.sbatch
#SBATCH --job-name=glue_pairing
#SBATCH --partition=batch

# --- one node, one task ------------------------------------------------------
# The work is a sequence of sparse matrix products (scipy) plus kNN in a 50-dim
# latent space (numpy). Both are single-process; the parallelism that exists is
# inside BLAS, so cores help the kNN and nothing else. This is a MEMORY job, not
# a core-count job -- do not request 64 cores expecting it to go faster.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=04:00:00

# jobs/ rather than slurm/logs/, matching cutTag-pipeline. SLURM will NOT
# create the directory and silently discards output if it is missing, so
# `mkdir -p jobs` before the first submit (the lab README says the same).
#SBATCH --output=jobs/pairing_%j.out
#SBATCH --error=jobs/pairing_%j.err
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=you@ohsu.edu

# No `set -e`: the aggregation and the validator each capture their own status,
# so a nonzero aggregation still reaches the validator and the closing guidance
# instead of aborting the script silently. Input checks below exit explicitly.
set -uo pipefail

# --- paths and parameters ----------------------------------------------------
# A batch job does NOT inherit the shell variables you exported at a prompt, so
# source the project's own settings. Edit setenv.sh, not this file.
cd "$SLURM_SUBMIT_DIR" || exit 1
if [[ -f setenv.sh ]]; then
    # shellcheck disable=SC1091
    source setenv.sh
else
    echo "ERROR: setenv.sh not found in $SLURM_SUBMIT_DIR" >&2
    exit 1
fi

OUT="${PAIRED:-ACC_GEX.h5mu}"

# Refuse before queue time is spent rather than failing deep in h5py.
for f in "$REF/rna.h5ad" "$REF/atac.h5ad" "$LABELS"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing input: $f" >&2
        [[ "$f" == /*.h5ad && "$(dirname "$f")" == "/" ]] && \
            echo "       (parent is '/', so REF was empty -- check setenv.sh)" >&2
        exit 1
    fi
done
echo "inputs verified under REF=$REF"

# --- environment -------------------------------------------------------------
# This step needs anndata + mudata + scipy, so the h5py-only interpreter that
# ran the inspectors will NOT work. Build the env once:
#     mamba env create -f 03_pipeline/pairing_env.yml
# Set these three, then the script parses and runs. They are shell variables
# rather than bare <PLACEHOLDER> tokens so that `bash -n slurm/pairing.sbatch`
# succeeds on an unedited checkout -- a literal <TOKEN> is a redirection to a
# file named TOKEN and makes the file a syntax error, which hides any REAL
# syntax error until submit time.
PAIRING_ENV="scplus-pairing"          # matches 03_pipeline/pairing_env.yml

# House style in cutTag-pipeline is to `conda activate` BEFORE sbatch and let
# the job inherit it (sbatch exports the environment by default). Honour that
# when it is already active, and fall back to activating here otherwise, so the
# script works either way.
if [[ "${CONDA_DEFAULT_ENV:-}" == "$PAIRING_ENV" ]]; then
    echo "using inherited conda env: $CONDA_DEFAULT_ENV"
else
    CONDA_BASE="${CONDA_BASE:-$HOME/miniconda3}"
    if [[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
        # shellcheck disable=SC1091
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        conda activate "$PAIRING_ENV"
        echo "activated conda env: $PAIRING_ENV (from $CONDA_BASE)"
    else
        echo "ERROR: '$PAIRING_ENV' is not active and no conda.sh under" >&2
        echo "       CONDA_BASE=$CONDA_BASE." >&2
        echo "       Either 'conda activate $PAIRING_ENV' before sbatch (house" >&2
        echo "       style), or set CONDA_BASE at the top of this script." >&2
        exit 1
    fi
fi

python -c "import anndata, mudata, scipy, numpy; \
print('env ok:', 'anndata', anndata.__version__, 'mudata', mudata.__version__, \
      'scipy', scipy.__version__, 'numpy', numpy.__version__)"

# --- guard rail --------------------------------------------------------------
# Refuse rather than get OOM-killed an hour in. 60% of --mem leaves room for the
# sparse inputs, the h5mu write buffer, and interpreter overhead.
# Parse the unit -- do not assume MB. 'mem=80G' previously had its 'G' stripped
# and was read as 80 MB, giving --max-dense-gb 0, which silently DISABLED the
# guard (0 means "no limit") instead of setting it to 48.
MEM_RAW=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null \
          | grep -oE 'mem=[0-9]+[KMGT]?' | head -1 | sed 's/mem=//')
MEM_GB=0
if [[ -n "$MEM_RAW" ]]; then
    NUM="${MEM_RAW//[!0-9]/}"
    UNIT="${MEM_RAW//[0-9]/}"
    case "${UNIT:-M}" in
        K) MEM_GB=$(( NUM / 1048576 )) ;;
        M) MEM_GB=$(( NUM / 1024 )) ;;
        G) MEM_GB=$NUM ;;
        T) MEM_GB=$(( NUM * 1024 )) ;;
    esac
fi
if [[ "$MEM_GB" -gt 0 ]]; then
    MAX_DENSE=$(( MEM_GB * 6 / 10 ))
    echo "job memory ${MEM_RAW} = ${MEM_GB} GB -> --max-dense-gb ${MAX_DENSE}"
else
    MAX_DENSE=48
    echo "could not read job memory (scontrol said '${MEM_RAW:-nothing}');" \
         "falling back to --max-dense-gb ${MAX_DENSE}"
fi

mkdir -p jobs

echo "=== pairing started $(date) on $(hostname) ==="
/usr/bin/time -v python 02_pair/aggregate_atac_sparse.py \
    --atac                "$REF/atac.h5ad" \
    --rna                 "$REF/rna.h5ad" \
    --obs-tsv             "$LABELS" \
    --latent-key          X_glue \
    --group-key           "$GROUP_KEY" \
    `# NO --use-raw: rna.h5ad X is already raw int32 counts and has no .raw` \
    --cells-per-metacell  "$K" \
    --oversample          "$OVERSAMPLE" \
    --min-metacells-per-group "$FLOOR" \
    --min-cells-per-group "$MIN_CELLS" \
    --max-dense-gb        "$MAX_DENSE" \
    --out                 "$OUT" \
    --diagnostics         pairing_diagnostics
AGG_RC=$?
echo "=== pairing finished $(date) ==="

# --- verify the handoff before anything downstream depends on it -------------
python 03_pipeline/validate_h5mu.py "$OUT"
VAL_RC=$?
echo "aggregate exit: $AGG_RC   validate exit: $VAL_RC"

echo
echo "Next: place $OUT at output_data.combined_GEX_ACC_mudata in the SCENIC+"
echo "config so snakemake treats prepare_GEX_ACC as already satisfied, then"
echo "read pairing_diagnostics.csv -- the median_crossmodal_gap column tells"
echo "you which cell types paired badly, and independent_metacell_equiv is the"
echo "honest denominator for how well any eRegulon is supported."

if [[ ${AGG_RC:-1} -ne 0 || ${VAL_RC:-1} -ne 0 ]]; then
    echo "one or more steps failed -- see the exits above" >&2
    exit 2
fi
exit 0
