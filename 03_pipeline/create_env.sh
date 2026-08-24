#!/usr/bin/env bash
# =============================================================================
# Build the scenicplus env. Use this INSTEAD of `mamba env create` directly.
#
# WHY THIS WRAPPER EXISTS
# -----------------------
# `mamba env create -f 03_pipeline/environment.yml` fails on 2026-08-24 with:
#
#     pybedtools uses setuptools (...) for installation but setuptools was
#     not found
#     CondaEnvException: Pip failed
#
# and it fails AFTER the 234-package conda solve has already succeeded.
#
# The cause is not pybedtools and not this environment file. setuptools 84.0.0
# (2026-08-08) REMOVED the bundled `pkg_resources` module. pybedtools 0.9.1's
# setup.py opens with `import pkg_resources` inside a try/except ImportError
# whose handler prints exactly that message. And because the 0.9.1 sdist ships
# no pyproject.toml, pip synthesizes a build requirement of
# `setuptools>=40.8.0` -- unbounded -- so the isolated build environment gets
# the newest setuptools, which no longer has pkg_resources.
#
# MEASURED, not assumed:
#   * setuptools 84.0.0 -> ModuleNotFoundError: No module named 'pkg_resources'
#   * setuptools 80.9.0 -> imports (with a deprecation warning)
#   * pinning setuptools in the CONDA section does NOT help: pip's build
#     isolation creates a fresh environment and ignores what is installed in
#     the target env. Verified by installing setuptools 80.9.0 into a venv and
#     watching the pybedtools build fail anyway.
#   * PIP_CONSTRAINT is honoured when pip populates the isolated build env, so
#     it is the one mechanism that reaches the build. With
#     `setuptools<81`, all four source-built packages build cleanly:
#     pybedtools 0.9.1, pyranges 0.0.111, tspex 0.6.3, pyrle 0.0.39,
#     MACS2 2.2.9.1.
#
# Three of the source-built sdists carry NO pyproject.toml (pybedtools,
# pyranges, tspex), so all three are exposed to the same break -- pybedtools is
# simply the one pip reached first.
#
# A requirements/environment file cannot express this: `--no-build-isolation`
# is rejected inside a requirements file, and a conda env file has nowhere to
# put a pip constraint. Hence a wrapper.
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
echo "  PIP_CONSTRAINT is exported so it reaches pip's ISOLATED build"
echo "  environment. Without it the build dies on pybedtools after the"
echo "  conda solve has already succeeded (setuptools 84 dropped"
echo "  pkg_resources; see the header of this script)."
echo

if "$SOLVER" env list 2>/dev/null | grep -qE "^${ENV_NAME}[[:space:]]"; then
    die "env '${ENV_NAME}' already exists. Remove it first:
       ${SOLVER} env remove -n ${ENV_NAME}
     or build under a different name:
       PIP_CONSTRAINT=${CONS} ${SOLVER} env create -f ${YML} -n <other>"
fi

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
