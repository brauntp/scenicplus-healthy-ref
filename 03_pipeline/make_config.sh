#!/usr/bin/env bash
# =============================================================================
# Fill config.template.yaml's placeholders and validate the result.
#
# WHY THIS EXISTS
# -----------------------------------------------------------------------------
# The template carries five <PLACEHOLDER> tokens spread through 25 KB of
# annotated YAML. Every value is already known -- from setenv.sh, from the
# download destination, from the sbatch --cpus-per-task -- so hand-editing five
# sites in a long file is a transcription error waiting to happen, and a wrong
# path here fails at rule time, after cisTarget has loaded a 33 GiB database.
#
# This substitutes them, then checks that every path the config names actually
# exists (or is a prepare_GEX_ACC input, whose absence is the intended guard),
# so the failure is
# reported now rather than hours in.
#
# ORDER OF OPERATIONS
# -----------------------------------------------------------------------------
#   1. source setenv.sh                       (REF, PAIRED, GROUP_KEY, ...)
#   2. bash 04_db/download_precomputed_db.sh  (the database must be present)
#   3. bash 03_pipeline/make_config.sh        <- here
#   4. sbatch slurm/scenicplus.sbatch
#
# Usage
# -----
#   bash 03_pipeline/make_config.sh                     # writes 03_pipeline/config.yaml
#   bash 03_pipeline/make_config.sh --out /tmp/c.yaml   # elsewhere
#   bash 03_pipeline/make_config.sh --check-only        # validate, write nothing
#
# UNREACHED ENTRIES: the template names two inputs the DAG does not open in our
# situation -- the cisTopic object and the separate RNA h5ad, both inputs to
# prepare_GEX_ACC. This script reports them rather than erroring, because a
# missing file there is expected and in fact desirable.
#
# The REASON is not is_multiome. Both branches of the Snakefile conditional
# (prepare_GEX_ACC_multiome / prepare_GEX_ACC_non_multiome) declare the SAME two
# inputs; the flag only selects which variant is defined. What keeps the rule out
# of the DAG is that its OUTPUT -- the paired MuData at
# output_data.combined_GEX_ACC_mudata -- already exists, so snakemake has no
# reason to build it and never evaluates its inputs. Measured on a real dry run:
# with the mudata present, 13 jobs are planned and prepare_GEX_ACC is absent;
# remove the mudata and snakemake raises MissingInputException naming both files.
#
# So absent is the SAFE state, not a gap: if snakemake ever judged the paired
# object stale, a missing input aborts the run instead of quietly regenerating
# GLUE-paired metacells as randomly-paired ones. Do not create placeholders.
# =============================================================================
set -o errexit
set -o nounset
set -o pipefail

SELF="$(basename "$0")"
log()  { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
die()  { printf 'ERROR (%s): %s\n' "$SELF" "$*" >&2; exit 1; }

# Resolve the repo root from this script's own location, not from `pwd`. The
# shell builtin reports $PWD, which a parent can hand down as a relative path,
# and a relative path written into the config breaks the moment snakemake runs a
# rule from a different directory.
#
# abspath, deliberately NOT realpath: resolving symlinks can land back on a
# mount-relative name (it does in the dev sandbox), and there is no reason here
# to prefer the physical path over the one the user invoked.
REPO_ROOT="$(python3 -c 'import os,sys; print(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[1]))))' "${BASH_SOURCE[0]}")"
case "$REPO_ROOT" in
    /*) : ;;
    *)  printf 'ERROR (%s): could not resolve an absolute repo root (got %s)\n' \
            "$(basename "$0")" "$REPO_ROOT" >&2; exit 1 ;;
esac

TEMPLATE="${REPO_ROOT}/03_pipeline/config.template.yaml"
OUT="${REPO_ROOT}/03_pipeline/config.yaml"
CHECK_ONLY=0
FORCE=0

while [ $# -gt 0 ]; do
    case "$1" in
        --out)        OUT="${2:?--out needs a path}"; shift 2 ;;
        --out=*)      OUT="${1#*=}"; shift ;;
        --check-only) CHECK_ONLY=1; shift ;;
        --force)      FORCE=1; shift ;;
        -h|--help)    sed -n '2,34p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *)            die "unrecognised option: $1" ;;
    esac
done

[ -f "$TEMPLATE" ] || die "template not found: $TEMPLATE"

# --------------------------------------------------------------------------- #
# Settings come from setenv.sh, which is the single source of truth for this
# project. Sourcing rather than re-deriving keeps this in step with the sbatch
# scripts, which do the same.
# --------------------------------------------------------------------------- #
if [ -f "${REPO_ROOT}/setenv.sh" ]; then
    # shellcheck disable=SC1091
    . "${REPO_ROOT}/setenv.sh" >/dev/null 2>&1 || true
else
    die "setenv.sh not found in ${REPO_ROOT} -- it defines REF, PAIRED, GROUP_KEY"
fi

: "${GROUP_KEY:?setenv.sh did not define GROUP_KEY}"
: "${PAIRED:?setenv.sh did not define PAIRED}"
CISTARGET_DB="${CISTARGET_DB:-${REPO_ROOT}/resources/cistarget_db}"
N_CPU="${N_CPU:-16}"          # matches slurm/scenicplus.sbatch --cpus-per-task

# Derive DB_PREFIX from what is actually on disk rather than hardcoding it: the
# prefix differs between the SCREEN build and a custom build, and guessing it
# produces a missing-file error late in the run.
RANKINGS="$(ls "${CISTARGET_DB}"/*.regions_vs_motifs.rankings.feather 2>/dev/null | head -1 || true)"
if [ -n "$RANKINGS" ]; then
    DB_PREFIX="$(basename "$RANKINGS" .regions_vs_motifs.rankings.feather)"
else
    DB_PREFIX="hg38_screen_v10_clust"
    log "WARNING: no rankings feather found in ${CISTARGET_DB}"
    log "         assuming DB_PREFIX=${DB_PREFIX}; run the download first"
fi

# NOTE: there is deliberately no SUBFOLDER substitution. The template's only
# occurrence of it was inside a comment describing the layout
# (region_set_folder/<family>/*.bed), so filling it with a guess made that
# comment name a family that may not exist. SCENIC+ scans EVERY subdirectory of
# region_set_folder, so no single family needs naming in the config: whatever
# 01_cistopic/region_sets_from_metacells.py wrote (DARs_cell_type) is picked up,
# and adding topics later is a second subdirectory with no config change.

log "substitutions:"
log "  ABS_PATH          = ${REPO_ROOT}"
log "  DB_PREFIX         = ${DB_PREFIX}"
log "  N_CPU             = ${N_CPU}"
log "  CELLTYPE_OBS_KEY  = ${GROUP_KEY}"

# --------------------------------------------------------------------------- #
# Substitute. sed with a delimiter that cannot appear in a path.
# --------------------------------------------------------------------------- #
RENDERED="$(mktemp)"
trap 'rm -f "$RENDERED"' EXIT
sed -e "s|<ABS_PATH>|${REPO_ROOT}|g" \
    -e "s|<DB_PREFIX>|${DB_PREFIX}|g" \
    -e "s|<N_CPU>|${N_CPU}|g" \
    -e "s|<CELLTYPE_OBS_KEY>|${GROUP_KEY}|g" \
    "$TEMPLATE" > "$RENDERED"

LEFT="$(grep -oE '<[A-Z_]{3,}>' "$RENDERED" | sort -u || true)"
if [ -n "$LEFT" ]; then
    die "unsubstituted placeholders remain: $(echo "$LEFT" | tr '\n' ' ')
       Every token must be filled -- a literal <TOKEN> reaching snakemake fails
       with an unhelpful path error."
fi

# --------------------------------------------------------------------------- #
# Validate: does every referenced path exist? prepare_GEX_ACC's two inputs are
# exempt -- absent is correct there, and PRESENT is what gets flagged.
# --------------------------------------------------------------------------- #
python3 - "$RENDERED" "$REPO_ROOT" "$PAIRED" <<'PY'
import os, sys
try:
    import yaml
except ImportError:
    sys.exit("ERROR: this python has no 'pyyaml', which the config validation "
             "needs.\n"
             f"       interpreter: {sys.executable}\n"
             "\n"
             "       Do NOT just activate scplus-pairing -- older builds of "
             "that env\n"
             "       lack pyyaml too, which is how this failed in job 10714858. "
             "Either:\n"
             "         mamba install -n scplus-pairing pyyaml     # one pure "
             "wheel, seconds\n"
             "       or run this from the scenicplus env, which carries it via "
             "snakemake.")

rendered, root, paired = sys.argv[1], sys.argv[2], sys.argv[3]
cfg = yaml.safe_load(open(rendered))
inp = cfg.get("input_data", {})
out = cfg.get("output_data", {})

# Inputs to prepare_GEX_ACC, the one rule that reads them. It stays out of the
# DAG because its output (the paired MuData) already exists -- NOT because of
# is_multiome, which merely picks between two rule variants declaring identical
# inputs. Verified on a real snakemake dry run: mudata present -> 13 jobs, rule
# absent, these paths never evaluated; mudata absent -> MissingInputException
# naming both. Absent is therefore the intended state.
BRANCH_DEAD = {"cisTopic_obj_fname", "GEX_anndata_fname"}

# is_multiome lives under params_data_preparation (Snakefile:10 branches on it),
# NOT params_general -- reading the wrong key returned None. It is still read and
# reported here because it selects which prepare_GEX_ACC variant is defined and
# therefore what a rebuild WOULD do, but it does not decide whether these two
# inputs get opened; the presence of the paired mudata does.
multi = cfg.get("params_data_preparation", {}).get("is_multiome")
if multi is None:
    sys.exit("ERROR: params_data_preparation.is_multiome not found in the "
             "config -- cannot tell which inputs the DAG will open.")
# Whether these two inputs are required depends on the paired mudata, NOT on
# is_multiome. This condition previously read `if not multi`, which contradicted
# the reasoning above: on the non-multiome branch it forced both keys back to
# required even when the mudata was present -- and the measured dry run shows the
# mudata's presence alone keeps prepare_GEX_ACC out of the DAG, under EITHER
# branch (both variants declare the same inputs).
_mu = out.get("combined_GEX_ACC_mudata")
_mu_abs = (_mu if (_mu and os.path.isabs(_mu))
           else (os.path.join(root, _mu) if _mu else None))
MUDATA_PRESENT = bool(_mu_abs and os.path.exists(_mu_abs))
if not MUDATA_PRESENT:
    # No paired object, so prepare_GEX_ACC IS in the DAG and genuinely reads
    # both files. They stop being guards and become required inputs.
    BRANCH_DEAD = set()

fail = []
# is_multiome picks WHICH prepare_GEX_ACC variant is defined; it does not decide
# whether the rule runs. The presence of the paired mudata does that, and it is
# checked separately below.
print(f"  is_multiome            : {multi}"
      + ("  (selects prepare_GEX_ACC_multiome)" if multi
         else "  (selects prepare_GEX_ACC_non_multiome -- SCENIC+'s own"
              " label-random pairing)"))
for k, v in inp.items():
    if not isinstance(v, str):
        continue
    exists = os.path.exists(v)
    # "guard" not "dead": these are prepare_GEX_ACC's inputs, and their ABSENCE
    # is what makes an unwanted rebuild fail loudly rather than silently replace
    # GLUE-paired metacells. Present is the state worth flagging.
    if k in BRANCH_DEAD:
        tag = "guard" if not exists else "PRESENT"
    else:
        tag = "OK   " if exists else "MISS "
    print(f"  [{tag:<5}] {k:<26} {v}")
    if (k in ("cisTopic_obj_fname", "GEX_anndata_fname")
            and not MUDATA_PRESENT and not exists):
        # Not a guard any more: without the paired object prepare_GEX_ACC is in
        # the DAG and genuinely needs these. Say why, so the fix is obvious.
        print("           ^ required because the paired MuData is MISSING, so"
              " prepare_GEX_ACC would have to run. Build the paired object"
              " first (02_pair/) rather than supplying these.")
    if k in BRANCH_DEAD:
        if exists:
            print(f"           ^ prepare_GEX_ACC input is PRESENT. Absent is safer:"
                  f" a rebuild would then abort instead of silently"
                  f" re-pairing at random.")
        continue
    if not exists:
        fail.append((k, v))

# The paired object is named relative to the working directory in output_data.
mu = out.get("combined_GEX_ACC_mudata")
if mu:
    cand = mu if os.path.isabs(mu) else os.path.join(root, mu)
    ok = os.path.exists(cand)
    print(f"  [{'OK   ' if ok else 'MISS '}] combined_GEX_ACC_mudata    {cand}")
    if not ok:
        fail.append(("combined_GEX_ACC_mudata", cand))

if fail:
    print()
    print("MISSING INPUTS -- the DAG would fail on these:")
    for k, v in fail:
        print(f"  {k}: {v}")
    sys.exit(1)
print()
print("all required inputs present; prepare_GEX_ACC inputs absent as intended")
PY
VALIDATED=$?

if [ "$VALIDATED" -ne 0 ]; then
    die "validation failed -- config not written"
fi

if [ "$CHECK_ONLY" -eq 1 ]; then
    log "--check-only: validated, nothing written"
    exit 0
fi

if [ -e "$OUT" ] && [ "$FORCE" -eq 0 ]; then
    die "$OUT exists -- pass --force to overwrite"
fi
cp "$RENDERED" "$OUT"
log "wrote ${OUT}"
log ""
log "NEXT: sbatch slurm/scenicplus.sbatch"
log "      (or dry-run first: 03_pipeline/run_pipeline.sh --config ${OUT} -n)"
