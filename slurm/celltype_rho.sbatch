#!/usr/bin/env bash
# =============================================================================
# Per-cell-type region-to-gene correlation -- SLURM job.
#
# ORDER OF OPERATIONS -- the env comes FIRST. --dry-run imports numpy, scipy,
# pandas and mudata like the real run does, so it needs the env too; it is cheap
# in TIME and MEMORY, not in dependencies.
#
#   1.  mamba env create -f 03_pipeline/pairing_env.yml     (once)
#   2.  conda activate scplus-pairing
#   3.  python 05_report/celltype_rho.py ... --dry-run      (login node, seconds)
#   4.  sbatch slurm/celltype_rho.sbatch
#
# Account and partition are intentionally NOT hardcoded: pass them on the
# command line, which overrides any directive here anyway.
#     sbatch --partition=... slurm/celltype_rho.sbatch
#SBATCH --job-name=celltype_rho
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
# 16G against a measured peak of ~2.5 GB for a 4,000-region block (the block,
# its argsort order array, and the rank matrix are live together; the factor is
# 5-6.5x the block's float32 size, measured in this project). The margin is for
# the per-group gene-expression rank matrices held across blocks.
#SBATCH --mem=16G
# Measured, not guessed: argsort is 166 us per (region, group), which over
# 360,980 regions x 23 groups is ~23 min; the correlation half is ~48 min
# (18.01 s per 50,000 pairs at n=5,058, scaled by each group's n). ~1.2 h
# single-core total, so 6 h is ~5x margin.
#SBATCH --time=06:00:00
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --error=slurm/logs/%x-%j.err

# No `set -e` past the input checks: the run and the summary each capture their
# own status, so a nonzero run still reaches the closing guidance instead of
# aborting silently.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p slurm/logs

# Paths live in setenv.sh so they survive logout and reach batch jobs.
if [[ -f setenv.sh ]]; then
    # shellcheck disable=SC1091
    source setenv.sh
fi

H5MU="${H5MU:-ACC_GEX.h5mu}"
ADJ="${ADJ:-03_pipeline/region_to_gene_adj.tsv}"
OUT="${OUT:-05_report/celltype_rho.tsv}"
GROUP_KEY="${GROUP_KEY:-predicted_CellType_Broad}"
BLOCK="${BLOCK:-4000}"

echo "=============================================================="
echo "per-cell-type region-to-gene correlation"
echo "=============================================================="
echo "python      : $(command -v python)"
echo "env         : ${CONDA_DEFAULT_ENV:-<none>}"
echo "paired      : $H5MU"
echo "adjacencies : $ADJ"
echo "out         : $OUT"
echo "host        : $(hostname)"
echo "job         : ${SLURM_JOB_ID:-<none>}"
echo

# Validate BEFORE any long read: an absent input should fail in a second, not
# after twenty minutes of streaming.
for f in "$H5MU" "$ADJ"; do
    if [[ ! -e "$f" ]]; then
        echo "ERROR: not found: $f" >&2
        echo "       set H5MU / ADJ in setenv.sh, or pass them in the environment." >&2
        exit 1
    fi
done
python -c "
import importlib.util as u, sys
miss = [m for m in ('numpy','scipy','pandas','mudata','anndata','h5py')
        if u.find_spec(m) is None]
if miss:
    sys.exit('ERROR: missing ' + ', '.join(miss) + '\n'
             '       interpreter: ' + sys.executable + '\n'
             '       conda activate scplus-pairing')
print('  imports ok')
" || exit 1

echo
echo "=== [1/2] plan ==="
python 05_report/celltype_rho.py --h5mu "$H5MU" --adj "$ADJ" --out "$OUT" \
    --group-key "$GROUP_KEY" --block "$BLOCK" --dry-run
PLAN_RC=$?
echo "plan exit: $PLAN_RC"

echo
echo "=== [2/2] per-group correlation ==="
python 05_report/celltype_rho.py --h5mu "$H5MU" --adj "$ADJ" --out "$OUT" \
    --group-key "$GROUP_KEY" --block "$BLOCK"
RUN_RC=$?
echo "run exit: $RUN_RC"

echo
echo "=== done ==="
if [[ "$RUN_RC" -eq 0 ]]; then
    echo "Wrote $OUT and its .detection_floors.csv"
    echo
    echo "Rank on specificity_z, NOT on raw rho differences: raw |rho| is not"
    echo "comparable across group sizes (95th percentile under pure noise was"
    echo "0.054 at n=1200 but 0.202 at n=90), so a raw contrast favours the"
    echo "smallest group. specificity_z shares one null scale across sizes."
    echo
    echo "Positive control: celltype_specific_links.csv holds 2,162 links"
    echo "attributed by regulon activity alone. Those should be enriched among"
    echo "high specificity_z links for the same group -- if they are not, the"
    echo "two routes disagree and that needs explaining before either is used."
    echo
    echo "Record the real cost while the record exists:"
    echo "  sacct -j ${SLURM_JOB_ID:-<job-id>} --format=JobID,MaxRSS,Elapsed,State --units=G | grep -v '\.batch'"
    echo "  MaxRSS is on the .0 step, not .batch (.batch is the empty shell)."
    echo "  Requested: ${SLURM_MEM_PER_NODE:-0} MB."
else
    echo "Run failed (exit $RUN_RC). A nonzero plan exit with a zero run exit"
    echo "usually means the adjacency file lacks a required column -- the plan"
    echo "step names it."
fi
echo "=============================================================="
exit "$RUN_RC"
