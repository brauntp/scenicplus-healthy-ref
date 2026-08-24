#!/usr/bin/env bash
# =============================================================================
# Build the scenicplus env. Use this INSTEAD of `mamba env create` directly.
#
# WHY THIS WRAPPER EXISTS
# -----------------------
# `mamba env create -f 03_pipeline/environment.yml` fails with:
#
#     pybedtools uses setuptools (...) for installation but setuptools was
#     not found
#     CondaEnvException: Pip failed
#
# and it fails AFTER the 234-package conda solve has already succeeded.
#
# THE MECHANISM (measured, not inferred)
# --------------------------------------
# setuptools 84.0.0 (2026-08-08) removed the bundled `pkg_resources` module.
# pybedtools 0.9.1's setup.py opens with `import pkg_resources` inside a
# try/except ImportError whose handler prints exactly that message.
#
# Which setuptools the build sees depends on WHICH BUILD PATH pip takes, and
# that depends on whether `wheel` is installed in the target environment:
#
#   wheel ABSENT  -> pip uses the PEP 517 path with an ISOLATED build env,
#                    which installs `setuptools>=40.8.0` unbounded, i.e. 84.
#                    ("Installing build dependencies ... Getting requirements
#                    to build wheel: error")   -> THIS IS THE FAILURE
#   wheel PRESENT -> pip uses the legacy setup.py path, which runs in the
#                    TARGET env against the setuptools installed there.
#                    ("Preparing metadata (setup.py) ... done")
#
# Full matrix, measured on pybedtools 0.9.1 with setuptools 80.9.0 in the env:
#
#   no wheel, isolation ON   -> FAILED   <- what happened on the cluster
#   no wheel, isolation OFF  -> BUILT
#   wheel present, ON        -> BUILT
#   wheel present, OFF       -> BUILT
#
# THE FIX is therefore in environment.yml, not here: `setuptools<81` and
# `wheel` are now conda-section dependencies, so the legacy path is taken and
# it uses a setuptools that still has pkg_resources. conda-forge resolves the
# bound to 80.10.2 (verified for linux-64/py3.11; the solve is still 234
# packages).
#
# PIP_CONSTRAINT is still exported below, as a second line of defence for pip
# versions that DO apply constraints to build dependencies. It is not the
# primary mechanism: it did not prevent the failure on the cluster, so do not
# rely on it.
#
# Also worth knowing: three of the five source-built sdists ship no
# pyproject.toml (pybedtools 0.9.1, pyranges 0.0.111, tspex 0.6.3), so all
# three take whichever path pip chooses. pybedtools is simply the one pip
# reached first -- fixing it alone would have moved the failure, not removed
# it. All five build cleanly under the fix: those three plus pyrle 0.0.39 and
# MACS2 2.2.9.1.
#
# Usage:
#     bash 03_pipeline/create_env.sh                # build it
#     bash 03_pipeline/create_env.sh --dry-run      # show what would run
#
# Run under screen/tmux: the pip section compiles from source and takes an
# hour or more.
# =============================================================================
set -uo pipefail

REPO_ROOT="$(python3 -c 'import os,sys;print(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[1]))))' "$0")"
cd "$REPO_ROOT" || exit 1

YML="03_pipeline/environment.yml"
CONS="03_pipeline/pip-build-constraints.txt"
DRY=0
[[ "${1:-}" == "--dry-run" ]] && DRY=1

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

[[ -f "$YML" ]]  || die "$YML not found (run from the repo, or let this script cd)"
[[ -f "$CONS" ]] || die "$CONS not found -- it carries the setuptools bound"

SOLVER=""
for c in mamba micromamba conda; do
    command -v "$c" >/dev/null 2>&1 && { SOLVER="$c"; break; }
done
[[ -n "$SOLVER" ]] || die "no mamba/micromamba/conda on PATH"

ENV_NAME="$(python3 -c "
import sys
for line in open(sys.argv[1]):
    if line.startswith('name:'):
        print(line.split(':', 1)[1].strip()); break
" "$YML")"

echo "=============================================================="
echo "building the ${ENV_NAME} env"
echo "=============================================================="
echo "  solver      : $SOLVER ($(command -v $SOLVER))"
echo "  spec        : $YML"
echo "  constraints : $CONS"
grep -vE '^[[:space:]]*(#|$)' "$CONS" | sed 's/^/                /'
echo
echo "  The real fix is in the spec: setuptools<81 and wheel are conda"
echo "  dependencies, so pip takes the legacy setup.py build path in the"
echo "  TARGET env instead of an isolated env that fetches setuptools 84"
echo "  (which dropped pkg_resources, breaking pybedtools 0.9.1's setup.py)."
echo "  PIP_CONSTRAINT below is a second line of defence, not the mechanism."
echo

if "$SOLVER" env list 2>/dev/null | grep -qE "^${ENV_NAME}[[:space:]]"; then
    die "env '${ENV_NAME}' already exists. Remove it first:
       ${SOLVER} env remove -n ${ENV_NAME}
     or build under a different name:
       PIP_CONSTRAINT=${CONS} ${SOLVER} env create -f ${YML} -n <other>"
fi

for req in "setuptools<81" "wheel"; do
    if ! grep -qE "^\s+- ${req}\s*$" "$YML"; then
        die "$YML is missing the conda dependency '${req}'.
     That bound is what makes pip take the legacy setup.py build path.
     Without it the build fails on pybedtools after the conda solve --
     see the header of this script."
    fi
done
log "spec carries the build-path fix (setuptools<81, wheel)"

log "validating the spec first (seconds)"
if ! bash 03_pipeline/preflight_env.sh; then
    die "preflight failed -- fix the items above before spending an hour here"
fi
echo

CMD=("$SOLVER" env create -f "$YML")
if (( DRY )); then
    echo "--dry-run: would run"
    echo "    PIP_CONSTRAINT=$CONS ${CMD[*]}"
    exit 0
fi

log "starting the build -- expect an hour or more"
export PIP_CONSTRAINT="$REPO_ROOT/$CONS"
"${CMD[@]}"
RC=$?

echo
echo "=============================================================="
if [[ "$RC" -eq 0 ]]; then
    echo "env created. Smoke test it before submitting:"
    echo "    conda activate ${ENV_NAME}"
    echo "    scenicplus --help"
    echo "    python -c \"import pycisTopic, pycistarget, scenicplus; print('ok')\""
    echo
    echo "Then:  sbatch slurm/scenicplus.sbatch"
else
    echo "build FAILED (exit $RC)"
    echo
    echo "If the failure is a pip BUILD error naming a package, check whether"
    echo "its sdist has a pyproject.toml. Without one, pip synthesizes an"
    echo "unbounded setuptools requirement and picks up whatever is newest --"
    echo "which is how the pybedtools break happened. Add a bound to"
    echo "$CONS rather than editing $YML: the constraint is about the"
    echo "BUILD environment, which the env file cannot reach."
fi
echo "=============================================================="
exit "$RC"
