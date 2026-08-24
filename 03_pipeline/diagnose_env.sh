#!/usr/bin/env bash
# =============================================================================
# What is actually installed in the active env, and does it match the spec?
#
# WHY THIS EXISTS
# ---------------
# `create_env.sh --repair-pins` reported "STILL DRIFTED: 33 of 34 ... installed
# MISSING" while all five git packages imported and the scenicplus console
# script existed. Those two facts cannot both be true: scenicplus imports
# pandas. So one of the two reports was wrong, and the useful move is to ask
# the environment directly rather than add more inference.
#
# This is read-only. It installs, removes and modifies nothing.
#
# Usage (with the env ACTIVE):
#     conda activate scenicplus
#     bash 03_pipeline/diagnose_env.sh
# =============================================================================
set -uo pipefail

REPO_ROOT="$(python3 -c 'import os,sys;print(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[1]))))' "$0")"
cd "$REPO_ROOT" || exit 1
YML="03_pipeline/environment.yml"

PY="${CONDA_PREFIX:-}/bin/python"
[[ -x "$PY" ]] || PY="$(command -v python)"
[[ -x "$PY" ]] || { echo "no python found -- activate the env first" >&2; exit 1; }

echo "=============================================================="
echo "environment diagnosis"
echo "=============================================================="
echo "  CONDA_DEFAULT_ENV : ${CONDA_DEFAULT_ENV:-<unset>}"
echo "  CONDA_PREFIX      : ${CONDA_PREFIX:-<unset>}"
echo "  python            : $PY"
echo "  version           : $("$PY" -c 'import sys;print(".".join(map(str,sys.version_info[:3])))' 2>&1)"
echo "  pip               : $("$PY" -m pip --version 2>&1 | head -1)"
echo

echo "-- [1] how many distributions does each mechanism see? ---------"
echo "  importlib.metadata : $("$PY" -c '
from importlib.metadata import distributions
print(sum(1 for _ in distributions()))' 2>&1 | tail -1)"
echo "  pip list           : $("$PY" -m pip list --format=freeze \
    --disable-pip-version-check 2>/dev/null | grep -c '==')"
echo "  site-packages dirs : $("$PY" -c '
import site,sys
for p in site.getsitepackages(): print(p)' 2>/dev/null | head -2 | tr '\n' ' ')"
echo
echo "  If these disagree wildly, the metadata query is the problem, not the"
echo "  environment. A --force-reinstall interrupted mid-flight can leave"
echo "  partially-written .dist-info directories that one reader tolerates and"
echo "  another does not."
echo

echo "-- [2] can the packages the pipeline needs be imported? --------"
for mod in numpy pandas scipy h5py anndata mudata sklearn matplotlib \
           pycisTopic pycistarget scenicplus loomxpy pyscenic; do
    v="$("$PY" -c "
import ${mod} as m
print(getattr(m, '__version__', 'no __version__'))" 2>&1 | tail -1)"
    if "$PY" -c "import ${mod}" >/dev/null 2>&1; then
        printf "  ok      %-14s %s\n" "$mod" "$v"
    else
        printf "  FAIL    %-14s %s\n" "$mod" "$(head -c 60 <<<"$v")"
    fi
done
echo

echo "-- [3] the pinned versions, as the metadata reports them -------"
python3 - "$YML" /tmp/.diag_pins.$$ <<'PYEOF'
import sys, yaml
d = yaml.safe_load(open(sys.argv[1]))
pip = [x for y in d["dependencies"] if isinstance(y, dict) for x in y.get("pip", [])]
open(sys.argv[2], "w").write("\n".join(pip) + "\n")
PYEOF
"$PY" - "/tmp/.diag_pins.$$" <<'PYEOF'
import sys
from importlib.metadata import distributions, version, PackageNotFoundError
want = {}
for line in open(sys.argv[1]):
    line = line.strip()
    if "==" in line and "git+" not in line:
        n, v = line.split("==", 1)
        want[n.strip()] = v.strip()
norm = lambda s: s.strip().lower().replace("_", "-")
meta = {}
for d in distributions():
    n = d.metadata["Name"]
    if n:
        meta[norm(n)] = d.version
ok = bad = gone = 0
for name, w in sorted(want.items()):
    key = norm(name)
    got = meta.get(key)
    if got is None:
        # Second opinion: version() resolves differently from iterating
        # distributions() in some broken-metadata cases.
        try:
            got = version(name)
        except PackageNotFoundError:
            got = None
    if got is None:
        print(f"  ABSENT   {name:<20} spec {w}")
        gone += 1
    elif got != w:
        print(f"  DIFFERS  {name:<20} spec {w:<15} installed {got}")
        bad += 1
    else:
        ok += 1
print(f"\n  {ok} match, {bad} differ, {gone} absent  (of {len(want)} pins)")
PYEOF
rm -f "/tmp/.diag_pins.$$"
echo
echo "=============================================================="
echo "Send the whole output. What matters: whether section 1's counts"
echo "agree, and whether section 3's ABSENT list overlaps section 2's"
echo "FAIL list. A package that imports but reports ABSENT means"
echo "damaged metadata -- the code is there, the .dist-info is not."
echo "=============================================================="
