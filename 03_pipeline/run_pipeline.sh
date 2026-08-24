#!/usr/bin/env bash
# =============================================================================
# SCENIC+ v1.0a2 pipeline driver
# =============================================================================
# Validates the config and every input it references, prints a resolved plan,
# then hands off to snakemake. Everything fails loud: no step is allowed to
# proceed on a missing or unreadable input, because the cheapest failure is the
# one that happens before the job is queued.
#
# Usage:
#   03_pipeline/run_pipeline.sh --config <config.yaml> [options]
#
# Options:
#   --config FILE      Config YAML (required). Copy config.template.yaml and
#                      fill in the <ANGLE_BRACKET> placeholders.
#   --cores N          Cores for snakemake AND the n_cpu passed to the CLI
#                      (default: $SLURM_CPUS_PER_TASK, else nproc, else 8).
#   --workdir DIR      Snakemake working directory; all relative output_data
#                      paths resolve here (default: dirname of --config).
#   --snakefile FILE   Snakefile to run. Default: resolved from the INSTALLED
#                      scenicplus package via importlib.resources, i.e.
#                      <site-packages>/scenicplus/snakemake/Snakefile -- the
#                      same lookup scenicplus's own init_snakemake uses.
#   -n | --dry-run     Passthrough to snakemake -n (also prints the plan).
#   --unlock           Run `snakemake --unlock` and exit. Use after a job is
#                      killed by SLURM (OOM / walltime) and leaves the .snakemake
#                      directory locked.
#   --validate-only    Run all checks, print the plan, do not invoke snakemake.
#   --                 Everything after this is forwarded verbatim to snakemake.
#
# Examples:
#   03_pipeline/run_pipeline.sh --config 03_pipeline/config.yaml -n
#   03_pipeline/run_pipeline.sh --config 03_pipeline/config.yaml --cores 32
#   03_pipeline/run_pipeline.sh --config 03_pipeline/config.yaml -- --until region_to_gene
# =============================================================================

set -Eeuo pipefail

# --- diagnostics -------------------------------------------------------------
RED=''; GRN=''; YLW=''; BLD=''; RST=''
if [[ -t 2 ]]; then RED=$'\033[31m'; GRN=$'\033[32m'; YLW=$'\033[33m'; BLD=$'\033[1m'; RST=$'\033[0m'; fi

die() {
    echo "" >&2
    echo "${RED}${BLD}ERROR${RST} $1" >&2
    shift || true
    for line in "$@"; do echo "      $line" >&2; done
    echo "" >&2
    exit 1
}
warn() { echo "${YLW}WARN${RST}  $*" >&2; }
ok()   { echo "${GRN}ok${RST}    $*"; }
info() { echo "      $*"; }

trap 'code=$?; [[ $code -ne 0 ]] && echo "${RED}run_pipeline.sh aborted (exit $code) at line $LINENO${RST}" >&2' ERR

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- defaults ----------------------------------------------------------------
CONFIG=""
CORES=""
WORKDIR=""
# Empty by default: resolved from the INSTALLED scenicplus package below, not
# from a path relative to this repo. An earlier default pointed at
#   ${REPO_ROOT}/src/scenicplus/src/scenicplus/snakemake/Snakefile
# -- a doubled src/scenicplus/ AND a location this repo never contains, so the
# job died in its first second on a clone instruction that was never needed.
# The Snakefile ships inside the package: <site-packages>/scenicplus/snakemake/.
SNAKEFILE=""
DRY_RUN=0
UNLOCK=0
VALIDATE_ONLY=0
EXTRA_ARGS=()

# --- argument parsing --------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)        CONFIG="${2:?--config requires a path}"; shift 2 ;;
        --config=*)      CONFIG="${1#*=}"; shift ;;
        --cores)         CORES="${2:?--cores requires an integer}"; shift 2 ;;
        --cores=*)       CORES="${1#*=}"; shift ;;
        --workdir)       WORKDIR="${2:?--workdir requires a path}"; shift 2 ;;
        --workdir=*)     WORKDIR="${1#*=}"; shift ;;
        --snakefile)     SNAKEFILE="${2:?--snakefile requires a path}"; shift 2 ;;
        --snakefile=*)   SNAKEFILE="${1#*=}"; shift ;;
        -n|--dry-run|--dryrun) DRY_RUN=1; shift ;;
        --unlock)        UNLOCK=1; shift ;;
        --validate-only) VALIDATE_ONLY=1; shift ;;
        -h|--help)      sed -n '2,34p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        --)             shift; EXTRA_ARGS+=("$@"); break ;;
        *)              die "unrecognised option: $1" \
                            "Run with --help. To forward flags to snakemake, put them after '--'." ;;
    esac
done

[[ -n "$CONFIG" ]] || die "--config is required." \
    "Copy 03_pipeline/config.template.yaml, fill in every <ANGLE_BRACKET> placeholder," \
    "and pass it here."
[[ -f "$CONFIG" ]] || die "config file not found: $CONFIG"
[[ -r "$CONFIG" ]] || die "config file not readable: $CONFIG"
CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"

# Resolve the Snakefile from the installed package unless --snakefile was
# given. importlib.resources is how scenicplus's own `init_snakemake` command
# locates it (cli/commands.py: files("scenicplus.snakemake").joinpath(...)), so
# this asks the same question the upstream code does.
if [[ -z "$SNAKEFILE" ]]; then
    SNAKEFILE="$(python - <<'PY' 2>/dev/null
import sys
try:
    from importlib.resources import files
    p = files("scenicplus.snakemake").joinpath("Snakefile")
    print(p if p.is_file() else "", end="")
except Exception:
    print("", end="")
PY
)"
fi

[[ -n "$SNAKEFILE" && -f "$SNAKEFILE" ]] || die \
    "cannot locate the scenicplus Snakefile." \
    "It ships inside the installed package, at" \
    "  <site-packages>/scenicplus/snakemake/Snakefile" \
    "and is found via importlib.resources -- the same lookup scenicplus's own" \
    "init_snakemake command uses. Getting nothing back means scenicplus is not" \
    "importable by the active python, or was installed without its package" \
    "data. Check:" \
    "  python -c 'from importlib.resources import files; print(files(\"scenicplus.snakemake\"))'" \
    "  bash 03_pipeline/smoke_test.sh" \
    "Or pass --snakefile /path/to/Snakefile explicitly."
SNAKEFILE="$(cd "$(dirname "$SNAKEFILE")" && pwd)/$(basename "$SNAKEFILE")"

if [[ -z "$WORKDIR" ]]; then WORKDIR="$(dirname "$CONFIG")"; fi
mkdir -p "$WORKDIR" || die "cannot create --workdir: $WORKDIR"
WORKDIR="$(cd "$WORKDIR" && pwd)"

# Core count: explicit > SLURM allocation > nproc > 8. Never guess high, since
# n_cpu is also the joblib/arboreto worker count and oversubscription on a
# shared node thrashes.
if [[ -z "$CORES" ]]; then
    if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
        CORES="$SLURM_CPUS_PER_TASK"
    elif command -v nproc >/dev/null 2>&1; then
        CORES="$(nproc)"
    else
        CORES=8
    fi
fi
[[ "$CORES" =~ ^[0-9]+$ && "$CORES" -ge 1 ]] || die "--cores must be a positive integer, got: $CORES"

# --- tool availability -------------------------------------------------------
command -v snakemake  >/dev/null 2>&1 || die "'snakemake' not on PATH." \
    "Activate the pipeline env:  conda activate scenicplus  (see 03_pipeline/environment.yml)"
command -v scenicplus >/dev/null 2>&1 || die "'scenicplus' CLI not on PATH." \
    "Every rule in the Snakefile shells out to 'scenicplus'. Activate the env, or" \
    "reinstall with:  pip install 'scenicplus @ git+https://github.com/aertslab/scenicplus@v1.0a2'"
command -v python     >/dev/null 2>&1 || die "'python' not on PATH."

# --- config reader -----------------------------------------------------------
# Read the config through Python/PyYAML rather than grep/sed. The config
# legitimately contains values that break naive shell parsing -- notably
# bc_transform_func, whose value is a nested-quoted lambda, and dem_promoter_space,
# written as 1_000 (a YAML int with an underscore separator).
cfg() {
    python - "$CONFIG" "$1" <<'PY'
import sys, yaml
path, key = sys.argv[1], sys.argv[2]
with open(path) as fh:
    cfg = yaml.safe_load(fh) or {}
node = cfg
for part in key.split("."):
    if not isinstance(node, dict) or part not in node:
        sys.exit(3)
    node = node[part]
if node is None:
    print("")
else:
    print(node)
PY
}

require_key() {
    local key="$1" val
    if ! val="$(cfg "$key")"; then
        die "config key missing: ${key}" \
            "The pinned Snakefile reads this key unconditionally; snakemake will die" \
            "with a KeyError while building the DAG. Restore it from" \
            "03_pipeline/config.template.yaml."
    fi
    printf '%s' "$val"
}

python -c 'import yaml' 2>/dev/null || die "PyYAML not importable in the active python." \
    "It is a snakemake dependency, so this usually means the wrong env is active."

python - "$CONFIG" <<'PY' || die "config is not valid YAML: $CONFIG" "See the parser message above."
import sys, yaml
try:
    with open(sys.argv[1]) as fh:
        yaml.safe_load(fh)
except yaml.YAMLError as e:
    print(e, file=sys.stderr); sys.exit(1)
PY

echo ""
echo "${BLD}=== [1/5] config keys ===${RST}"

# Every two-level key the pinned Snakefile reads. Enumerated from the source
# (151 config[...][...] lookups over 67 distinct keys) rather than trusted from
# the template, so a stale config is caught here instead of inside snakemake.
REQUIRED_KEYS=(
  input_data.cisTopic_obj_fname
  input_data.GEX_anndata_fname
  input_data.region_set_folder
  input_data.ctx_db_fname
  input_data.dem_db_fname
  input_data.path_to_motif_annotations
  output_data.combined_GEX_ACC_mudata
  output_data.dem_result_fname
  output_data.ctx_result_fname
  output_data.output_fname_dem_html
  output_data.output_fname_ctx_html
  output_data.cistromes_direct
  output_data.cistromes_extended
  output_data.tf_names
  output_data.genome_annotation
  output_data.chromsizes
  output_data.search_space
  output_data.tf_to_gene_adjacencies
  output_data.region_to_gene_adjacencies
  output_data.eRegulons_direct
  output_data.eRegulons_extended
  output_data.AUCell_direct
  output_data.AUCell_extended
  output_data.scplus_mdata
  params_general.temp_dir
  params_general.n_cpu
  params_general.seed
  params_data_preparation.bc_transform_func
  params_data_preparation.is_multiome
  params_data_preparation.key_to_group_by
  params_data_preparation.nr_cells_per_metacells
  params_data_preparation.direct_annotation
  params_data_preparation.extended_annotation
  params_data_preparation.species
  params_data_preparation.biomart_host
  params_data_preparation.search_space_upstream
  params_data_preparation.search_space_downstream
  params_data_preparation.search_space_extend_tss
  params_motif_enrichment.species
  params_motif_enrichment.annotation_version
  params_motif_enrichment.motif_similarity_fdr
  params_motif_enrichment.orthologous_identity_threshold
  params_motif_enrichment.annotations_to_use
  params_motif_enrichment.fraction_overlap_w_dem_database
  params_motif_enrichment.dem_max_bg_regions
  params_motif_enrichment.dem_balance_number_of_promoters
  params_motif_enrichment.dem_promoter_space
  params_motif_enrichment.dem_adj_pval_thr
  params_motif_enrichment.dem_log2fc_thr
  params_motif_enrichment.dem_mean_fg_thr
  params_motif_enrichment.dem_motif_hit_thr
  params_motif_enrichment.fraction_overlap_w_ctx_database
  params_motif_enrichment.ctx_auc_threshold
  params_motif_enrichment.ctx_nes_threshold
  params_motif_enrichment.ctx_rank_threshold
  params_inference.tf_to_gene_importance_method
  params_inference.region_to_gene_importance_method
  params_inference.region_to_gene_correlation_method
  params_inference.order_regions_to_genes_by
  params_inference.order_TFs_to_genes_by
  params_inference.gsea_n_perm
  params_inference.quantile_thresholds_region_to_gene
  params_inference.top_n_regionTogenes_per_gene
  params_inference.top_n_regionTogenes_per_region
  params_inference.min_regions_per_gene
  params_inference.rho_threshold
  params_inference.min_target_genes
)
for k in "${REQUIRED_KEYS[@]}"; do require_key "$k" >/dev/null; done
ok "all ${#REQUIRED_KEYS[@]} keys read by the Snakefile are present"

# Placeholder sweep -- an unfilled <ANGLE_BRACKET> is the single most common way
# this pipeline gets submitted broken. Scan the PARSED VALUES, not the raw text:
# the template's own comments legitimately mention tokens like <SUBFOLDER> when
# documenting directory layout, and a grep over the file flags those as errors.
# <N_CPU> is exempt -- this driver fills it from --cores when it writes the
# resolved config below.
placeholder_report="$(python - "$CONFIG" <<'PY'
import re, sys, yaml
with open(sys.argv[1]) as fh:
    cfg = yaml.safe_load(fh) or {}
pat = re.compile(r"<[A-Z0-9_]+>")
bad = []
def walk(node, path):
    if isinstance(node, dict):
        for k, v in node.items():
            walk(v, f"{path}.{k}" if path else str(k))
    elif isinstance(node, list):
        for i, v in enumerate(node):
            walk(v, f"{path}[{i}]")
    elif isinstance(node, str):
        for m in pat.findall(node):
            if m != "<N_CPU>":
                bad.append(f"{path}: {m}   (value: {node!r})")
for line in bad:
    print(line)
PY
)" || die "placeholder scan failed"

if [[ -n "$placeholder_report" ]]; then
    echo "" >&2
    echo "$placeholder_report" >&2
    die "config still contains unfilled <ANGLE_BRACKET> placeholders (listed above)." \
        "Replace every one before running. Only <N_CPU> may be left in place --" \
        "this driver substitutes it from --cores."
fi
ok "no unfilled placeholders in config values (<N_CPU> is filled by this driver)"

# --- resolved config ---------------------------------------------------------
# n_cpu has to be injected into the config rather than passed with --config.
# snakemake/cli.py:217 validates every --config key against `[a-zA-Z_]\w*$`, so
# only FLAT keys are accepted -- there is no command-line syntax for setting a
# nested key like params_general.n_cpu, and attempting one raises
# "Config entry must start with a valid identifier". So we materialise a resolved
# copy of the config with n_cpu set to --cores and hand THAT to --configfile.
# The copy also serves as the run's provenance record.
# params_general.temp_dir is resolved the same way. An empty value in the source
# config means "use the scratch this job was given", i.e. $TMPDIR -- which
# slurm/scenicplus.sbatch points at per-job node-local disk. That keeps a single
# source of truth for scratch instead of letting the config and the sbatch script
# disagree about where joblib spills.
RESOLVED_TEMP_DIR="$(cfg params_general.temp_dir)"
if [[ -z "$RESOLVED_TEMP_DIR" ]]; then
    if [[ -n "${TMPDIR:-}" ]]; then
        RESOLVED_TEMP_DIR="${TMPDIR%/}/scenicplus_tmp"
        info "params_general.temp_dir is empty -> using \$TMPDIR: ${RESOLVED_TEMP_DIR}"
    else
        die "params_general.temp_dir is empty and \$TMPDIR is not set." \
            "Either export TMPDIR (slurm/scenicplus.sbatch does this, pointing at" \
            "per-job node-local scratch) or set an explicit path in the config." \
            "This directory is joblib's temp_folder for region_to_gene -- it needs" \
            "tens of GB of fast local disk, so it must not silently default to /tmp."
    fi
fi

RESOLVED_CONFIG="${WORKDIR}/.resolved_config.yaml"
python - "$CONFIG" "$RESOLVED_CONFIG" "$CORES" "$RESOLVED_TEMP_DIR" <<'PY' \
    || die "failed to write the resolved config to ${WORKDIR}/.resolved_config.yaml"
import sys, yaml, datetime
src, dst, cores, temp_dir = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
with open(src) as fh:
    cfg = yaml.safe_load(fh)
pg = cfg.setdefault("params_general", {})
pg["n_cpu"] = cores
pg["temp_dir"] = temp_dir
with open(dst, "w") as fh:
    fh.write(f"# Resolved by run_pipeline.sh at {datetime.datetime.now().isoformat(timespec='seconds')}\n")
    fh.write(f"# Source: {src}\n")
    fh.write(f"# params_general.n_cpu    -> {cores}      (from --cores)\n")
    fh.write(f"# params_general.temp_dir -> {temp_dir}\n")
    fh.write("# Do not edit -- regenerated on every run.\n")
    yaml.safe_dump(cfg, fh, default_flow_style=False, sort_keys=False)
PY
ok "resolved config written  $RESOLVED_CONFIG  (n_cpu=${CORES}, temp_dir=${RESOLVED_TEMP_DIR})"

# --- resolve values ----------------------------------------------------------
CISTOPIC_OBJ="$(cfg input_data.cisTopic_obj_fname)"
GEX_ANNDATA="$(cfg input_data.GEX_anndata_fname)"
REGION_SET_FOLDER="$(cfg input_data.region_set_folder)"
CTX_DB="$(cfg input_data.ctx_db_fname)"
DEM_DB="$(cfg input_data.dem_db_fname)"
MOTIF_ANNOT="$(cfg input_data.path_to_motif_annotations)"
H5MU="$(cfg output_data.combined_GEX_ACC_mudata)"
FINAL_TARGET="$(cfg output_data.scplus_mdata)"
TEMP_DIR="$RESOLVED_TEMP_DIR"
CFG_N_CPU="$(cfg params_general.n_cpu)"
IS_MULTIOME="$(cfg params_data_preparation.is_multiome)"
DEM_BALANCE="$(cfg params_motif_enrichment.dem_balance_number_of_promoters)"
SS_UP="$(cfg params_data_preparation.search_space_upstream)"
SS_DOWN="$(cfg params_data_preparation.search_space_downstream)"
QUANTILES="$(cfg params_inference.quantile_thresholds_region_to_gene)"
TOPN_GENE="$(cfg params_inference.top_n_regionTogenes_per_gene)"
R2G_METHOD="$(cfg params_inference.region_to_gene_importance_method)"

# Relative output paths resolve against the snakemake working directory.
abspath_in_workdir() {
    case "$1" in /*) printf '%s' "$1" ;; *) printf '%s/%s' "$WORKDIR" "$1" ;; esac
}
H5MU_ABS="$(abspath_in_workdir "$H5MU")"
FINAL_TARGET_ABS="$(abspath_in_workdir "$FINAL_TARGET")"

# --- unlock shortcut ---------------------------------------------------------
if [[ "$UNLOCK" -eq 1 ]]; then
    echo ""
    echo "${BLD}=== unlocking ===${RST}"
    info "workdir:   $WORKDIR"
    info "snakefile: $SNAKEFILE"
    warn "Only do this if no snakemake process is still running on this workdir."
    # Must use the RESOLVED config: the raw one may still carry the <N_CPU>
    # placeholder, and the Snakefile's `threads: config["params_general"]["n_cpu"]`
    # rejects a non-numeric value with "Threads value has to be an integer,
    # float, or a callable" before it ever gets to unlocking.
    exec snakemake --snakefile "$SNAKEFILE" --directory "$WORKDIR" \
                   --configfile "$RESOLVED_CONFIG" --unlock
fi

# --- input existence checks --------------------------------------------------
echo ""
echo "${BLD}=== [2/5] referenced inputs exist ===${RST}"

check_file() {
    local path="$1" label="$2"; shift 2
    [[ -n "$path" ]] || die "$label is empty in the config."
    [[ -e "$path" ]] || die "$label does not exist: $path" "$@"
    [[ -f "$path" ]] || die "$label is not a regular file: $path"
    [[ -r "$path" ]] || die "$label exists but is not readable: $path" \
        "Check permissions and, on ARC, that the filesystem holding it is mounted on compute nodes."
    ok "$label  ($(du -h "$path" 2>/dev/null | cut -f1))  $path"
}

# --- prepare_GEX_ACC's inputs: absent is CORRECT, and load-bearing -----------
#
# I previously required these to exist, on the belief that snakemake evaluates
# every rule's inputs while building the DAG. That is wrong, and measured:
#
#   paired mudata EXISTS   -> prepare_GEX_ACC is not in the DAG at all
#                             (13 jobs planned, rule absent), and the missing
#                             cisTopic/GEX paths are never evaluated.
#   paired mudata ABSENT   -> the rule is needed and snakemake raises
#                             MissingInputException naming both files.
#
# Note also that `is_multiome` does NOT skip the rule -- both branches of the
# Snakefile conditional (prepare_GEX_ACC_multiome / _non_multiome) declare the
# SAME two inputs. The flag only selects which variant exists. What keeps the
# rule out of the DAG is that its output, the paired mudata, is already there.
#
# So these files being absent is exactly the state we want, and it is a SAFETY
# GUARD rather than a gap: if snakemake ever concludes the paired object is
# stale and tries to rebuild it, the run dies with MissingInputException instead
# of silently replacing GLUE-paired metacells with randomly paired ones. Do not
# create placeholder files to satisfy a checker -- that would disarm the guard.
#
# What must be verified instead is the OUTPUT they would have produced -- the
# paired MuData. That check already lives below, against H5MU_ABS (resolved
# against --workdir), so it is not repeated here.
for _p in "$CISTOPIC_OBJ" "$GEX_ANNDATA"; do
    if [[ -e "$_p" ]]; then
        warn "prepare_GEX_ACC input EXISTS: $_p"
        warn "  Absent is the safer state. With both inputs present, a snakemake"
        warn "  rerun that considers the paired object stale would regenerate it"
        warn "  by SCENIC+'s own label-random pairing, discarding the GLUE"
        warn "  metacells without error. run_pipeline.sh passes"
        warn "  --rerun-triggers mtime to make that unlikely; a missing input"
        warn "  makes it impossible."
    else
        ok "prepare_GEX_ACC input absent (intended guard): $_p"
    fi
done

check_file "$CTX_DB" "cisTarget rankings DB" \
    "Expected <prefix>.regions_vs_motifs.rankings.feather." \
    "NOTE: cisTarget needs the RANKINGS file and DEM needs the SCORES file." \
    "They are different files; swapping them fails late and confusingly."
check_file "$DEM_DB" "DEM scores DB" \
    "Expected <prefix>.regions_vs_motifs.scores.feather (scores, not rankings)."
check_file "$MOTIF_ANNOT" "motif annotation table" \
    "Expected e.g. motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl, matching" \
    "params_motif_enrichment.annotation_version."

# Cheap sanity check on the ctx/dem swap, which is otherwise only detected hours in.
case "$CTX_DB" in *rankings*) : ;; *) warn "ctx_db_fname does not contain 'rankings' -- verify it is the rankings DB." ;; esac
case "$DEM_DB" in *scores*)   : ;; *) warn "dem_db_fname does not contain 'scores' -- verify it is the scores DB." ;; esac
if [[ "$CTX_DB" == "$DEM_DB" ]]; then
    die "ctx_db_fname and dem_db_fname are the same file." \
        "cisTarget consumes the rankings database; DEM consumes the scores database."
fi

# region_set_folder: cli/commands.py:216-224 iterates the folder and only
# descends into SUBDIRECTORIES, reading *.bed inside them. Loose .bed files at
# the top level are ignored, which produces an empty enrichment rather than an
# error -- so check the structure, not just existence.
[[ -n "$REGION_SET_FOLDER" ]] || die "input_data.region_set_folder is empty."
[[ -d "$REGION_SET_FOLDER" ]] || die "region_set_folder is not a directory: $REGION_SET_FOLDER"
[[ -r "$REGION_SET_FOLDER" ]] || die "region_set_folder is not readable: $REGION_SET_FOLDER"
n_subdir=0; n_bed=0
while IFS= read -r d; do
    [[ -z "$d" ]] && continue
    c=$(find "$d" -maxdepth 1 -name '*.bed' -type f 2>/dev/null | wc -l | tr -d ' ')
    if [[ "$c" -gt 0 ]]; then n_subdir=$((n_subdir + 1)); n_bed=$((n_bed + c)); fi
done < <(find "$REGION_SET_FOLDER" -mindepth 1 -maxdepth 1 -type d 2>/dev/null)
if [[ "$n_subdir" -eq 0 ]]; then
    loose=$(find "$REGION_SET_FOLDER" -maxdepth 1 -name '*.bed' -type f 2>/dev/null | wc -l | tr -d ' ')
    die "region_set_folder contains no subdirectory holding .bed files: $REGION_SET_FOLDER" \
        "scenicplus/cli/commands.py:216-224 iterates the folder, keeps only entries that" \
        "are DIRECTORIES, and reads *.bed inside each. Required layout:" \
        "    region_set_folder/<set_family>/<region_set>.bed" \
        "Found ${loose} loose .bed file(s) at the top level -- those are silently ignored." \
        "Nest them one level deeper, e.g. region_set_folder/topics_otsu/Topic_1.bed"
fi
ok "region_set_folder  (${n_subdir} subfolder(s), ${n_bed} .bed file(s))  $REGION_SET_FOLDER"

# --- the hand-built paired MuData -------------------------------------------
echo ""
echo "${BLD}=== [3/5] paired MuData hand-off ===${RST}"
info "Our data are unpaired scRNA + scATAC co-embedded with scGLUE, so"
info "02_pair/glue_metacells.py builds the paired MuData itself. We do NOT use"
info "SCENIC+'s non-multiome path: process_non_multiome_data draws metacells"
info "uniformly at random within a label INDEPENDENTLY for RNA and ATAC"
info "(generate_pseudocells_for_numpy), which destroys the within-cell-type"
info "covariation region_to_gene measures. Placing the paired file at"
info "output_data.combined_GEX_ACC_mudata makes snakemake treat prepare_GEX_ACC"
info "as already satisfied, so neither prepare_GEX_ACC branch ever runs."
echo ""

if [[ ! -e "$H5MU_ABS" ]]; then
    # Distinguish "the file does not exist" from "the config named it relative
    # to the wrong base". The second is what actually happened once: the value
    # was 'ACC_GEX.h5mu', snakemake resolves relative paths against --directory
    # (03_pipeline/), and the object lives at the repo root -- so it looked
    # absent while sitting one level up.
    _hint=()
    if [[ ! "$H5MU" = /* ]]; then
        _alt="${REPO_ROOT}/${H5MU}"
        if [[ -e "$_alt" ]]; then
            _hint=(
              ""
              "FOUND IT one level up: $_alt"
              "The config value '${H5MU}' is RELATIVE, and snakemake resolves"
              "relative paths against --directory (${WORKDIR}) -- not the repo"
              "root. Make this key ABSOLUTE in 03_pipeline/config.template.yaml:"
              "    combined_GEX_ACC_mudata: \"<ABS_PATH>/${H5MU}\""
              "then re-run 03_pipeline/make_config.sh. It is the one output_data"
              "entry that names a PRE-EXISTING file, so it is the one that must"
              "not float with the working directory."
            )
        fi
    fi
    die "paired MuData not found at combined_GEX_ACC_mudata: $H5MU_ABS" \
        "(config value: '${H5MU}', resolved against workdir '${WORKDIR}')" \
        "${_hint[@]}" \
        "" \
        "Run 02_pair/glue_metacells.py first and write its output to exactly that path." \
        "Without it, snakemake will try to RUN prepare_GEX_ACC -- which for" \
        "is_multiome=${IS_MULTIOME} would either fail on our unpaired inputs or, worse," \
        "silently build randomly-paired metacells."
fi
[[ -f "$H5MU_ABS" ]] || die "combined_GEX_ACC_mudata is not a regular file: $H5MU_ABS"
[[ -r "$H5MU_ABS" ]] || die "combined_GEX_ACC_mudata is not readable: $H5MU_ABS"
ok "paired MuData present  ($(du -h "$H5MU_ABS" | cut -f1))  $H5MU_ABS"

# Structural validation. This is the gate before a multi-hour job: modality keys
# must be literally "scRNA"/"scATAC", obs_names must align positionally, region
# names must be coordinate-parseable.
VALIDATOR="${REPO_ROOT}/03_pipeline/validate_h5mu.py"
if [[ -f "$VALIDATOR" ]]; then
    echo ""
    if ! python "$VALIDATOR" "$H5MU_ABS"; then
        die "validate_h5mu.py rejected the paired MuData (see the report above)." \
            "Fix 02_pair/glue_metacells.py output before submitting; every failure it" \
            "reports would otherwise surface hours into the cluster job, or not at all."
    fi
else
    warn "validate_h5mu.py not found at $VALIDATOR -- skipping structural validation."
    warn "This is the gate before a multi-hour job; restore it."
fi

# --- temp dir ----------------------------------------------------------------
echo ""
echo "${BLD}=== [4/5] scratch / temp_dir ===${RST}"
if [[ ! -d "$TEMP_DIR" ]]; then
    mkdir -p "$TEMP_DIR" || die "cannot create params_general.temp_dir: $TEMP_DIR"
    ok "created temp_dir  $TEMP_DIR"
else
    ok "temp_dir exists  $TEMP_DIR"
fi
[[ -w "$TEMP_DIR" ]] || die "params_general.temp_dir is not writable: $TEMP_DIR"
if command -v df >/dev/null 2>&1; then
    avail_kb=$(df -Pk "$TEMP_DIR" | awk 'NR==2{print $4}')
    avail_gb=$(( avail_kb / 1024 / 1024 ))
    info "free space on temp_dir: ${avail_gb} GiB"
    if [[ "$avail_gb" -lt 50 ]]; then
        warn "under 50 GiB free on temp_dir. joblib memory-maps the region_to_gene"
        warn "arrays through here; running out mid-run kills the job late."
    fi
fi
case "$TEMP_DIR" in
    "$HOME"|"$HOME"/*)
        warn "temp_dir is under \$HOME. On ARC that is a quota'd shared filesystem;"
        warn "joblib memory-maps large arrays through temp_dir, so prefer node-local"
        warn "scratch (slurm/scenicplus.sbatch exports one)."
        ;;
esac

# --- resolved plan -----------------------------------------------------------
echo ""
echo "${BLD}=== [5/5] resolved plan ===${RST}"
cat <<PLAN
  snakefile              : ${SNAKEFILE}
  configfile (source)    : ${CONFIG}
  configfile (passed)    : ${RESOLVED_CONFIG}
  working directory      : ${WORKDIR}
  cores (--cores)        : ${CORES}
  n_cpu in source config : ${CFG_N_CPU}   -> resolved to ${CORES}
  temp_dir               : ${TEMP_DIR}

  paired MuData (input)  : ${H5MU_ABS}
  final target           : ${FINAL_TARGET_ABS}

  prepare_GEX_ACC        : BYPASSED (output pre-satisfied by 02_pair/glue_metacells.py)
  is_multiome            : ${IS_MULTIOME}  (selects which bypassed rule is defined)
  dem_balance_promoters  : ${DEM_BALANCE}  (True adds a DAG edge DEM <- genome_annotation)

  search space           : upstream ${SS_UP} / downstream ${SS_DOWN} bp (min max)
  region_to_gene method  : ${R2G_METHOD}
  r2g thresholds         : quantiles [${QUANTILES}], top-N per gene [${TOPN_GENE}]

  rules that will run    : motif_enrichment_cistarget, motif_enrichment_dem,
                           download_genome_annotations, prepare_menr,
                           get_search_space, tf_to_gene, region_to_gene,
                           eGRN_direct, eGRN_extended, AUCell_direct,
                           AUCell_extended, scplus_mudata
PLAN

if [[ "$VALIDATE_ONLY" -eq 1 ]]; then
    echo ""
    ok "--validate-only: all checks passed, not invoking snakemake."
    exit 0
fi


# --- invoke snakemake --------------------------------------------------------
SNAKE_ARGS=(
    --snakefile "$SNAKEFILE"
    --directory "$WORKDIR"

    # The upstream Snakefile hard-codes `configfile: "config/config.yaml"` on
    # line 3. Passing --configfile is still correct and is what we do:
    # snakemake/workflow.py:1454-1473 raises on the hard-coded path only when it
    # is missing AND no --configfile was given; with --configfile supplied it
    # falls through to update_config(config, overwrite_config). Caveat: if a file
    # DOES exist at ${WORKDIR}/config/config.yaml it is loaded first and then
    # overlaid by ours, so keep the run directory free of stray configs.
    --configfile "$RESOLVED_CONFIG"
    --cores "$CORES"

    # --- why --rerun-triggers mtime -----------------------------------------
    # Snakemake 8's default rerun triggers are {mtime, params, input, code,
    # software-env} (snakemake/settings.py:175, RerunTrigger.all()). Three of
    # those are actively harmful to this workflow:
    #
    #   code / params : every rule's `params:` block is a lambda closing over
    #     config, and we override params_general/n_cpu on the command line. So
    #     changing --cores, or touching any unrelated comment or default in the
    #     config, changes the recorded params fingerprint of nearly every rule
    #     and triggers a full rerun of a multi-hour pipeline for no scientific
    #     reason.
    #
    #   input : this is the one that matters most for us. We deliberately
    #     PRE-SATISFY output_data.combined_GEX_ACC_mudata with the h5mu built by
    #     02_pair/glue_metacells.py, so that prepare_GEX_ACC never runs. Under
    #     the default triggers, snakemake compares the recorded input set of that
    #     output against the rule's declared inputs (cisTopic_obj_fname,
    #     GEX_anndata_fname). Our file has no recorded provenance at all, so the
    #     input trigger fires and snakemake decides the h5mu is out of date and
    #     must be regenerated -- which would run prepare_GEX_ACC and overwrite
    #     our correctly-paired metacells with randomly-paired ones. That is
    #     exactly the failure mode this whole design exists to avoid, and it
    #     would be silent: the pipeline would complete successfully with wrong
    #     data.
    #
    #   software-env : no rule declares a conda/container directive, so this
    #     trigger contributes nothing but noise across env rebuilds.
    #
    # Restricting to mtime restores plain make semantics: rerun a rule only when
    # an input is newer than its output. That is the correct contract for a
    # pipeline whose first artifact is produced out-of-band, and it makes
    # resuming after a SLURM walltime kill actually work.
    #
    # Consequence to be aware of: if you edit a scientific parameter in the
    # config (search_space_upstream, rho_threshold, ...), snakemake will NOT
    # notice. Delete the affected outputs by hand, or run in a fresh --workdir.
    --rerun-triggers mtime

    --printshellcmds
    --keep-going
    --rerun-incomplete
)

if [[ "$DRY_RUN" -eq 1 ]]; then
    # Snakemake 8 removed --reason (it prints per-job reasons unconditionally in
    # a dry run), so passing it here is an argparse error. Verified against
    # snakemake 8.5.5 --help.
    SNAKE_ARGS+=(--dry-run)
    echo ""
    echo "${BLD}=== dry run ===${RST}"
else
    echo ""
    echo "${BLD}=== running ===${RST}"
fi

if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
    SNAKE_ARGS+=("${EXTRA_ARGS[@]}")
    info "forwarding to snakemake: ${EXTRA_ARGS[*]}"
fi

printf '      snakemake'; printf ' %q' "${SNAKE_ARGS[@]}"; printf '\n\n'

set +e
snakemake "${SNAKE_ARGS[@]}"
rc=$?
set -e

echo ""
if [[ $rc -ne 0 ]]; then
    echo "${RED}${BLD}snakemake exited ${rc}${RST}" >&2
    echo "  Triage:" >&2
    echo "    - 'Directory cannot be locked'  -> a previous run was killed; rerun with --unlock" >&2
    echo "    - killed with no python traceback -> almost always SLURM OOM. The cisTarget" >&2
    echo "      rankings DB is memory-resident; raise --mem in slurm/scenicplus.sbatch." >&2
    echo "    - 'No space left on device'     -> params_general.temp_dir filled up." >&2
    echo "    - logs are under ${WORKDIR}/.snakemake/log/" >&2
    exit $rc
fi

if [[ "$DRY_RUN" -eq 0 ]]; then
    if [[ -f "$FINAL_TARGET_ABS" ]]; then
        ok "final target written: $FINAL_TARGET_ABS ($(du -h "$FINAL_TARGET_ABS" | cut -f1))"
    else
        die "snakemake exited 0 but the final target is missing: $FINAL_TARGET_ABS" \
            "Check ${WORKDIR}/.snakemake/log/ -- --keep-going can mask a failed branch."
    fi
fi
ok "done"
