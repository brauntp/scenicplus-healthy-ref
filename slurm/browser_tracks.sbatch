#!/usr/bin/env bash
# =============================================================================
# Genome-browser tracks -- SLURM job.
#
# Two deliverables, and only ONE of them needs the cluster:
#
#   * E-P LINK tracks come from celltype_rho_usable.csv.gz (170 MB) and can be
#     built on a laptop. This job builds them anyway so both land together.
#   * ACCESSIBILITY tracks need ACC_GEX.h5mu (40.5 GB), which stays on ARC.
#     That is the reason this is a job rather than a local script.
#
# ORDER OF OPERATIONS -- the env comes FIRST. --dry-run imports numpy/pandas/
# mudata like the real run does; it is cheap in TIME and MEMORY, not in
# dependencies.
#
#   1.  mamba env create -f 03_pipeline/pairing_env.yml     (once)
#   2.  conda activate scplus-pairing
#   3.  python 03_pipeline/build_genome_annotation.py --out-dir 05_report/anno
#   4.  sbatch --partition=batch slurm/browser_tracks.sbatch
#
#SBATCH --job-name=browser_tracks
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
# 16G. The streaming block dominates: 8,000 peaks x 25,323 metacells float32 is
# 0.75 GB resident, and the measured overhead factor on this object was 5-6.5x
# the block (block + any densification + the group-sum writes), so ~5 GB. The
# per-group sums are float64 vectors totalling 39 MB. The remaining margin is
# because the last job on this object (celltype_rho, 10732007) came in at 4.44 GB
# against my 2.5 GB estimate -- I under-predicted by 1.8x by counting only the
# block and ignoring what persists across blocks.
#SBATCH --mem=16G
# One pass over 393,832 peaks in 50 blocks, with a mean per group per block --
# far cheaper than the per-group ranking in celltype_rho (48 min). Budgeting 2 h
# for ~5x margin on an unmeasured job; record the real figure from sacct.
#SBATCH --time=02:00:00
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --error=slurm/logs/%x-%j.err

# No `set -e` past the input checks: each step captures its own status so a
# nonzero step still reaches the closing guidance.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p slurm/logs

if [[ -f setenv.sh ]]; then
    # shellcheck disable=SC1091
    source setenv.sh
fi

H5MU="${H5MU:-ACC_GEX.h5mu}"
# The 13-group table. celltype_rho.sbatch writes the 23-GROUP table
# (05_report/celltype_rho.tsv); make_usable_subset.py derives this one from it.
# That derivation step was missing from the repo -- the table had been produced
# interactively and published without a script -- so this default pointed at a
# path that never existed on the cluster and stage 3 failed with "not found".
# Stage 0 below now builds it if absent.
RHO_FULL="${RHO_FULL:-05_report/celltype_rho.tsv}"
FLOORS="${FLOORS:-05_report/celltype_rho.detection_floors.csv}"
LINKS="${LINKS:-05_report/celltype_rho_usable.csv.gz}"
MIN_MC="${MIN_MC:-400}"
ANNO="${ANNO:-05_report/anno}"
OUT="${OUT:-05_report/browser_tracks}"
GROUP_KEY="${GROUP_KEY:-predicted_CellType_Broad}"
BLOCK="${BLOCK:-8000}"

echo "=============================================================="
echo "genome-browser tracks"
echo "=============================================================="
echo "python      : $(command -v python)"
echo "env         : ${CONDA_DEFAULT_ENV:-<none>}"
echo "paired      : $H5MU"
echo "links       : $LINKS"
echo "annotation  : $ANNO"
echo "out         : $OUT"
echo "host        : $(hostname)"
echo "job         : ${SLURM_JOB_ID:-<none>}"
echo

# Validate BEFORE any long read: a missing input should fail in a second, not
# after twenty minutes of streaming.
# Stage 0: derive the 13-group table if it is absent. Cheap (one pass over a
# 1.59M-row table), and it removes the ordering trap where the job assumed a
# file that only ever existed on a laptop.
if [[ ! -e "$LINKS" ]]; then
    echo "=== [0/4] $LINKS absent -- deriving it from $RHO_FULL ==="
    if [[ ! -e "$RHO_FULL" || ! -e "$FLOORS" ]]; then
        echo "ERROR: cannot derive $LINKS: need both" >&2
        echo "         $RHO_FULL" >&2
        echo "         $FLOORS" >&2
        echo "       Run the correlation job first:" >&2
        echo "         sbatch --partition=batch slurm/celltype_rho.sbatch" >&2
        exit 1
    fi
    python 05_report/make_usable_subset.py --in "$RHO_FULL" --floors "$FLOORS" \
        --out "$LINKS" --min-metacells "$MIN_MC"
    SUB_RC=$?
    echo "subset exit: $SUB_RC"
    if (( SUB_RC != 0 )); then
        echo "ERROR: could not build $LINKS; not proceeding." >&2
        exit "$SUB_RC"
    fi
    echo
fi

MISSING=0
for f in "$H5MU" "$LINKS" "$ANNO/genome_annotation.tsv" "$ANNO/chromsizes.tsv"; do
    if [[ ! -e "$f" ]]; then
        echo "ERROR: not found: $f" >&2
        MISSING=1
    fi
done
if (( MISSING )); then
    echo "       Build the annotation first:" >&2
    echo "         python 03_pipeline/build_genome_annotation.py --out-dir $ANNO" >&2
    echo "       and set H5MU / LINKS in setenv.sh if they live elsewhere." >&2
    exit 1
fi
python -c "
import importlib.util as u, sys
miss = [m for m in ('numpy','pandas','mudata','anndata','h5py')
        if u.find_spec(m) is None]
if miss:
    sys.exit('ERROR: missing ' + ', '.join(miss) + '\n'
             '       interpreter: ' + sys.executable + '\n'
             '       conda activate scplus-pairing')
print('  imports ok')
bw = u.find_spec('pyBigWig') is not None
print('  pyBigWig:', 'yes' if bw else 'no -- bedGraph only (still loads in IGV/UCSC)')
" || exit 1

echo
echo "=== [1/3] plan the accessibility tracks ==="
python 05_report/export_accessibility_tracks.py --h5mu "$H5MU" \
    --chromsizes "$ANNO/chromsizes.tsv" --out-dir "$OUT" \
    --group-key "$GROUP_KEY" --block "$BLOCK" --dry-run
PLAN_RC=$?
echo "plan exit: $PLAN_RC"

echo
echo "=== [2/3] accessibility tracks (streams the 40.5 GB object) ==="
python 05_report/export_accessibility_tracks.py --h5mu "$H5MU" \
    --chromsizes "$ANNO/chromsizes.tsv" --out-dir "$OUT" \
    --group-key "$GROUP_KEY" --block "$BLOCK"
ACC_RC=$?
echo "accessibility exit: $ACC_RC"

echo
echo "=== [3/3] E-P link tracks ==="
python 05_report/export_browser_tracks.py --links "$LINKS" \
    --annotation "$ANNO/genome_annotation.tsv" --out-dir "$OUT"
LNK_RC=$?
echo "links exit: $LNK_RC"

echo
echo "=== done ==="
if (( ACC_RC == 0 && LNK_RC == 0 )); then
    echo "Wrote $OUT:"
    echo "  <group>.accessibility.bedGraph / .bw   peak-level signal"
    echo "  <group>.interact                       E-P arcs (UCSC)"
    echo "  <group>.bedpe                          E-P arcs (IGV)"
    echo "  <group>.peaks.bed                      peak ends only"
    echo
    echo "TWO THINGS TO CARRY INTO ANY INTERPRETATION:"
    echo "  1. The accessibility tracks are PEAK-LEVEL step functions, not"
    echo "     fragment coverage -- flat within each peak, absent between."
    echo "     This reference has no fragments. Comparing cell types at a peak"
    echo "     is valid; reading a summit or a boundary off the track is not."
    echo "  2. Arc gene-end anchors are RECONSTRUCTED from the stored Distance,"
    echo "     not looked up. They reproduce |Distance| exactly, and sit a median"
    echo "     ~1.9 kb from the nearest annotated TSS because the pipeline's"
    echo "     transcript choice is not recoverable from a fresh annotation."
    echo
    echo "Record the real cost while the record exists:"
    echo "  sacct -j ${SLURM_JOB_ID:-<job-id>} --format=JobID,JobName%-16,MaxRSS,Elapsed,State --units=G"
    echo "  Do NOT filter by step name: this script does not use srun, so there"
    echo "  is no .0 step and MaxRSS is on .batch."
    echo "  Requested: ${SLURM_MEM_PER_NODE:-0} MB."
else
    echo "Something failed above -- read the two exit codes. A nonzero plan exit"
    echo "with a zero accessibility exit usually means --group-key is wrong."
fi
echo "=============================================================="
exit $(( ACC_RC != 0 ? ACC_RC : LNK_RC ))
