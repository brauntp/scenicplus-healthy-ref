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
# 8G. MEASURED: job 10732007 peaked at 4.44 GB (sacct .batch row) against my
# pre-run estimate of ~2.5 GB -- 1.8x UNDER, because the estimate counted only
# the streaming block and ignored the two persistent gene-rank structures. The
# accounting: per-group gene rank dict 1.70 GB + the full gene matrix 1.70 GB
# live together in phase 1 (3.40 GB), then ranks 1.70 + block work 1.51 GB in
# phase 2 (3.21 GB); the remainder is interpreter overhead. 8G is ~1.8x the
# measured peak. The first run requested 16G, which was 3.6x headroom.
#SBATCH --mem=8G
# MEASURED: job 10732007 ran in 48m24s wall for all 1,590,176 pairs x 23
# groups. My pre-run projection was ~72 min (23 min argsort + 48 min
# correlation), so 1.5x conservative. 2 h keeps ~2.5x margin; 6 h was excessive
# and only delays scheduling.
#SBATCH --time=02:00:00
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
    echo "  sacct -j ${SLURM_JOB_ID:-<job-id>} --format=JobID,JobName%-16,MaxRSS,Elapsed,State --units=G"
    echo "  Do NOT filter by step name. This script does not use srun, so there"
    echo "  is no .0 step and MaxRSS is on .batch -- an earlier version of this"
    echo "  banner printed 'grep -v .batch', which removed the only useful row."
    echo "  Requested: ${SLURM_MEM_PER_NODE:-0} MB."
else
    echo "Run failed (exit $RUN_RC). A nonzero plan exit with a zero run exit"
    echo "usually means the adjacency file lacks a required column -- the plan"
    echo "step names it."
fi
echo "=============================================================="
exit "$RUN_RC"
