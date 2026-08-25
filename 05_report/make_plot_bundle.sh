#!/usr/bin/env bash
# =============================================================================
# Build the local-analysis bundle ON THE CLUSTER, then print the rsync line.
#
# WHY: the final object is 41 GB and ~99% of it is the input matrices carried
# through. Everything a local plotting session needs is a few hundred MB, and the
# one piece that lives inside the 41 GB file (TF expression) is extracted here as
# ~22 MB. Do not transfer scplusmdata.h5mu.
#
#   conda activate scplus-pairing
#   bash 05_report/make_plot_bundle.sh
#
# Login-node safe: both extractors read in blocks/chunks and never materialise a
# whole matrix. Stated because this project has already shipped one tool wrongly
# labelled cheap.
# =============================================================================

# No `set -e`: each step captures its own status so a failure in one extractor
# still reaches the other and the closing summary.
set -uo pipefail

cd "$(dirname "$0")/.." || exit 1
OUT="05_report/plot_bundle"

echo "=============================================================="
echo "plot bundle"
echo "=============================================================="
echo "  repo    : $PWD"
echo "  python  : $(command -v python)"
echo "  out dir : $OUT"
echo

# -- interpreter check, before any long read ----------------------------------
# `python` may not exist at all (some systems ship only python3), in which case
# the heredoc below would die with a bare "command not found" instead of saying
# what to do. Check for the interpreter first.
if ! command -v python >/dev/null 2>&1; then
    echo "ERROR: no 'python' on PATH." >&2
    echo "       conda activate scplus-pairing" >&2
    echo "       (that env provides python, h5py, pandas and pyarrow)" >&2
    exit 1
fi
python - <<'PY' || exit 1
import sys
need = ["h5py", "numpy", "pandas", "pyarrow"]
missing = [m for m in need if __import__("importlib.util", fromlist=["x"])
           .find_spec(m) is None]
if missing:
    sys.exit(f"ERROR: missing {', '.join(missing)}\n"
             f"       interpreter: {sys.executable}\n"
             "       conda activate scplus-pairing\n"
             "       (pyarrow is needed for the .parquet outputs)")
print(f"  imports ok ({sys.executable})")
PY
echo

echo "=== [1/2] peak-to-gene links (the reusable resource) ==="
python 05_report/export_peak_gene_links.py --out-dir "$OUT"
LINK_RC=$?
echo "links exit: $LINK_RC"
echo

echo "=== [2/2] AUC objects + TF expression slice ==="
python 05_report/extract_for_plots.py --out-dir "$OUT"
EXTR_RC=$?
echo "extract exit: $EXTR_RC"
echo

echo "=== bundle contents ==="
if [[ -d "$OUT" ]]; then
    ls -lh "$OUT" | awk 'NR>1 {printf "  %-34s %s\n", $9, $5}'
    echo "  ------------------------------------------------"
    echo "  total: $(du -sh "$OUT" | cut -f1)"
else
    echo "  NOTHING WRITTEN -- read the two exit codes above"
fi
echo

if (( LINK_RC == 0 && EXTR_RC == 0 )); then
    echo "=== pull it down (run this on the LAPTOP, not here) ==="
    echo "    mkdir -p ~/comp_ws/scenicplus_bundle"
    echo "    rsync -avhP --stats \\"
    echo "      arc:${PWD}/${OUT}/ \\"
    echo "      ~/comp_ws/scenicplus_bundle/"
    echo
    echo "  'arc' is the ssh config Host alias, NOT this node's hostname."
    echo "  rsync shells out to ssh, which reads ~/.ssh/config, so the alias"
    echo "  resolves and any ProxyJump/User in it applies -- the two-hop login"
    echo "  (acc.ohsu.edu then arc.ohsu.edu) collapses to one step exactly as"
    echo "  'ssh arc' does. Do not substitute \$(hostname) here: the login"
    echo "  node's internal name is not reachable from off-cluster."
    echo "  Trailing slash on the source is load-bearing: with it you get the"
    echo "  CONTENTS; without it a nested plot_bundle/ directory."
    echo
    echo "  Do NOT rsync scplusmdata.h5mu (41 GB, ~99% input matrices)."
    echo "  peak_gene_links.README.md travels with the data and carries the"
    echo "  caveats -- read it before using the links in another project."
else
    echo "One of the extractors failed; the bundle is incomplete."
    echo "A nonzero links exit with a zero extract exit usually means"
    echo "region_to_gene_adj.tsv is not where the config says it is."
fi
echo "=============================================================="
exit $(( LINK_RC != 0 ? LINK_RC : EXTR_RC ))
