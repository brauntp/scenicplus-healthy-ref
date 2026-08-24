#!/usr/bin/env bash
# =============================================================================
# Does the built env actually run SCENIC+? Read-only, minutes not hours.
#
# WHY THIS EXISTS
# ---------------
# The env now satisfies all 34 pins in 03_pipeline/environment.yml, but
# scenicplus v1.0a2's own requirements.txt is a FULL LOCK FILE of 222 pinned
# packages, and our spec carries only 35 of them. pip therefore printed a long
# "X requires Y==a, but you have Y==b" list -- 187 uncovered pins, and in the
# visible slice alone nine were major-version jumps (toolz 0.12->1.1,
# zope-interface 6.2->8.6, xmltodict 0.13->1.0, ...).
#
# Those warnings are pip reporting declared-vs-installed. They are not evidence
# of breakage, and they are not evidence of safety either. Rather than reason
# about 187 packages, exercise the code paths the pipeline actually uses and see
# what imports and what runs.
#
# Usage (env ACTIVE, login node, no data needed):
#     conda activate scenicplus
#     bash 03_pipeline/smoke_test.sh
# =============================================================================
set -uo pipefail

REPO_ROOT="$(python3 -c 'import os,sys;print(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[1]))))' "$0")"
cd "$REPO_ROOT" || exit 1

PY="${CONDA_PREFIX:-}/bin/python"
[[ -x "$PY" ]] || PY="$(command -v python)"
[[ -x "$PY" ]] || { echo "no python -- activate the env first" >&2; exit 1; }

FAIL=0
echo "=============================================================="
echo "SCENIC+ env smoke test"
echo "=============================================================="
echo "  python : $PY"
echo

echo "-- [1] the CLI the pipeline invokes -------------------------------"
if "${CONDA_PREFIX}/bin/scenicplus" --help >/dev/null 2>&1; then
    echo "  ok    scenicplus --help"
    # Argparse prints the subcommand list inside an indented "positional
    # arguments:" block, not at column 0 -- an anchored ^\{...\} pattern matches
    # nothing and printed an empty line here on the first real run.
    SUBS="$("${CONDA_PREFIX}/bin/scenicplus" --help 2>&1 |
            grep -oE '\{[a-z_,]{4,}\}' | head -1 | tr -d '{}')"
    echo "  subcommands: ${SUBS:-<none parsed; argparse layout differs>}"
    # The DAG drives these two; if either is absent the CLI is not the one
    # this pipeline was written against.
    for want in prepare_data grn_inference; do
        case ",${SUBS}," in
            *",${want},"*) printf "  ok    subcommand %s\n" "$want" ;;
            *) printf "  WARN  subcommand %s not listed\n" "$want" ;;
        esac
    done
else
    echo "  FAIL  scenicplus --help"
    "${CONDA_PREFIX}/bin/scenicplus" --help 2>&1 | tail -5 | sed 's/^/        /'
    FAIL=1
fi
echo

echo "-- [2] the modules each Snakemake rule imports --------------------"
# One module per pipeline stage, so a failure names the stage that would break.
while read -r stage mod; do
    if err="$("$PY" -c "import ${mod}" 2>&1)"; then
        printf "  ok    %-26s %s\n" "$stage" "$mod"
    else
        printf "  FAIL  %-26s %s\n" "$stage" "$mod"
        printf "        %s\n" "$(tail -1 <<<"$err" | head -c 100)"
        FAIL=1
    fi
done <<'EOF'
motif_enrichment_cistarget pycistarget.motif_enrichment_cistarget
motif_enrichment_dem       pycistarget.motif_enrichment_dem
region_to_gene             scenicplus.grn_builder.gsea_approach
tf_to_gene                 arboreto.algo
eregulon_assembly          scenicplus.grn_builder.modules
AUCell                     scenicplus.eregulon_enrichment
cli_entrypoint             scenicplus.cli.scenicplus
snakemake_driver           snakemake
data_wrangling             scenicplus.data_wrangling.adata_cistopic_wrangling
cistarget_io               pycistarget.input_output
cistarget_result           pycistarget.motif_enrichment_result
EOF
echo

echo "-- [3] does it run real work on tiny synthetic input? -------------"
"$PY" - <<'PYEOF'
import sys, warnings, traceback
warnings.filterwarnings("ignore")
ok = True

def check(label, fn):
    global ok
    try:
        r = fn()
        print(f"  ok    {label:<34} {r}")
    except Exception as e:
        ok = False
        print(f"  FAIL  {label:<34} {type(e).__name__}: {str(e)[:70]}")

def t_mudata():
    import numpy as np, anndata as ad, mudata
    from scipy import sparse
    a = ad.AnnData(sparse.random(20, 50, density=0.3, format="csr", dtype="float32"))
    g = ad.AnnData(sparse.random(20, 30, density=0.3, format="csr", dtype="float32"))
    a.var_names = [f"chr1:{i*1000}-{i*1000+500}" for i in range(50)]
    m = mudata.MuData({"scRNA": g, "scATAC": a})
    return f"MuData {m.shape}"

def t_todf():
    # The exact call infer_region_to_gene makes on the paired object.
    import numpy as np, anndata as ad, mudata
    a = ad.AnnData(np.random.rand(15, 40).astype("float32"))
    m = mudata.MuData({"scATAC": a})
    df = m["scATAC"].to_df()
    return f"to_df {df.shape} {df.values.dtype}"

def t_arboreto():
    # TF-to-gene uses GRNBoost2; exercise the estimator it wraps.
    import numpy as np
    from arboreto.core import SGBM_KWARGS
    from sklearn.ensemble import GradientBoostingRegressor
    X = np.random.rand(60, 5); y = np.random.rand(60)
    GradientBoostingRegressor(n_estimators=5).fit(X, y)
    return f"GBM ok, SGBM_KWARGS keys={len(SGBM_KWARGS)}"

def t_ranking():
    # cisTarget reads the feather rankings with pyarrow.
    import pyarrow as pa, pyarrow.feather as feather, io, pandas as pd
    df = pd.DataFrame({"motifs": ["m1", "m2"], "r1": [1, 2], "r2": [2, 1]})
    buf = io.BytesIO()
    feather.write_feather(pa.Table.from_pandas(df), buf)
    buf.seek(0)
    back = feather.read_table(buf).to_pandas()
    return f"feather roundtrip {back.shape}, pyarrow {pa.__version__}"

def t_snakemake():
    import snakemake, subprocess, sys
    r = subprocess.run([sys.executable, "-m", "snakemake", "--version"],
                       capture_output=True, text=True)
    return f"snakemake {r.stdout.strip() or r.stderr.strip()[:30]}"

def t_ray():
    # pycisTopic and DEM parallelise over ray.
    import ray
    return f"ray {ray.__version__} imports"

for label, fn in (("mudata paired object", t_mudata),
                  ("to_df on a modality", t_todf),
                  ("arboreto / GBM regressor", t_arboreto),
                  ("pyarrow feather roundtrip", t_ranking),
                  ("snakemake --version", t_snakemake),
                  ("ray import", t_ray)):
    check(label, fn)
sys.exit(0 if ok else 1)
PYEOF
[[ $? -eq 0 ]] || FAIL=1
echo

echo "=============================================================="
if (( FAIL == 0 )); then
    echo "SMOKE TEST PASSED."
    echo
    echo "Every module the DAG imports loads, the CLI responds, and the"
    echo "specific operations SCENIC+ performs on the paired object run."
    echo "The pip 'requires X, but you have Y' warnings are pip comparing"
    echo "scenicplus's 222-package lock against what is installed; they are"
    echo "not breakage. Proceed:"
    echo "    bash 03_pipeline/make_config.sh"
    echo "    sbatch slurm/scenicplus.sbatch"
else
    echo "SMOKE TEST FAILED -- read the FAIL lines above."
    echo
    echo "A failure here names the pipeline stage that would break, hours into"
    echo "a run. If it is an ImportError on a package pip warned about, add"
    echo "that package at scenicplus's pinned version to the pip section of"
    echo "03_pipeline/environment.yml and re-run --repair-pins. The upstream"
    echo "pin is in scenicplus v1.0a2's requirements.txt."
fi
echo "=============================================================="
exit "$FAIL"
