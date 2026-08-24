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
# THE MECHANISM (measured across pip versions)
# -------------------------------------------
# setuptools 84.0.0 (2026-08-08) removed the bundled `pkg_resources` module.
# pybedtools 0.9.1's setup.py opens with `import pkg_resources` inside a
# try/except ImportError whose handler prints exactly that message.
#
# Whether the build sees a setuptools that still HAS pkg_resources depends on
# the build path, and the available paths depend on the PIP VERSION:
#
#   pip 24.0   wheel in env -> legacy setup.py path, TARGET env  -> BUILT
#   pip 25.2   wheel in env -> legacy setup.py path, TARGET env  -> BUILT
#   pip 26.2.1 wheel in env -> legacy path GONE; PEP 517 isolated
#                              build env installs setuptools>=40.8.0
#                              unbounded, i.e. 84               -> FAILED
#
# So two earlier fixes were both dead ends, and both looked right in a venv
# running an older pip:
#   * PIP_CONSTRAINT=setuptools<81      -- not applied to build deps here
#   * setuptools<81 + wheel in conda    -- only helps where the legacy path
#                                          still exists, i.e. pip < 26
# PIP_NO_BUILD_ISOLATION as an ENV VAR is inert on pip 26 at any value
# (tested 1, true, yes). Only the COMMAND-LINE flag works.
#
# THE FIX: two pip passes, which is why this wrapper drives pip itself instead
# of letting `mamba env create` run the yml's pip section.
#
#   PASS A  the five sdists with no usable wheel, with --no-build-isolation
#           --no-deps. Their build deps (setuptools<81, wheel, cython, numpy)
#           come from the conda section, already installed.
#             pybedtools 0.9.1, pyranges 0.0.111, tspex 0.6.3,
#             pyrle 0.0.39, MACS2 2.2.9.1
#   PASS B  the full pinned pip list, WITH isolation. Three of the five git
#           dependencies declare their own PEP 518 backends (hatchling for
#           pycistarget, poetry for LoomXpy, setuptools-scm for pycisTopic and
#           scenicplus), so isolation must stay ON for them --
#           --no-build-isolation cannot be applied globally.
#
# Pass A's installs satisfy pass B's pins exactly, so pip reports "Requirement
# already satisfied" and does not rebuild them (verified on pip 26.2.1).
#
# Three of the five sdists ship no pyproject.toml (pybedtools, pyranges,
# tspex); polars 0.20.13 is NOT in the list because its cp38-abi3 manylinux
# wheel works on 3.11.
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
echo "  Five stages. pip 26 removed the legacy setup.py build path, so the"
echo "  five wheel-less sdists are built --no-build-isolation against the"
echo "  conda-installed setuptools<81 (84 dropped pkg_resources, which"
echo "  pybedtools 0.9.1 imports). The pinned PyPI list then installs WITH"
echo "  isolation, the git packages --no-deps (scenicplus pins LoomXpy@main"
echo "  against our commit pin), and pip check drives a measured repair."
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
     Pass A runs pip with --no-build-isolation, which installs NO build
     dependencies -- so setuptools and wheel must already be in the env,
     and setuptools must be < 81 because 84 dropped pkg_resources.
     See the header of this script."
    fi
done
log "spec provides pass A's build deps (setuptools<81, wheel)"

log "validating the spec first (seconds)"
if ! bash 03_pipeline/preflight_env.sh; then
    die "preflight failed -- fix the items above before spending an hour here"
fi
echo


# The five sdists with no usable wheel, built WITHOUT isolation. Kept here
# rather than parsed out of the yml so the list is explicit and reviewable;
# create_env.sh cross-checks each against the yml's pip section below, so the
# two cannot drift silently.
NOISO=(
    "pybedtools==0.9.1"
    "pyranges==0.0.111"
    "pyrle==0.0.39"
    "tspex==0.6.3"
    "MACS2==2.2.9.1"
)

# Drift check: every NOISO entry must appear verbatim in the yml's pip section,
# otherwise pass A would pin a version pass B then rebuilds.
for spec in "${NOISO[@]}"; do
    grep -qE "^[[:space:]]+- ${spec//./\\.}[[:space:]]*$" "$YML" || \
        die "'${spec}' is in this script's NOISO list but not in ${YML}'s pip
     section. Pass A would install a version pass B does not pin, and pip
     would rebuild it under isolation -- the failure this wrapper exists to
     avoid."
done
log "NOISO list matches the spec (${#NOISO[@]} packages)"

# Write a pip requirements file holding ONLY the yml's pip section, so pass B
# installs exactly what the spec pins.
PIPREQ="$(mktemp -t scplus-pip.XXXXXX)"
trap 'rm -f "$PIPREQ"' EXIT
python3 - "$YML" "$PIPREQ" <<'PYEOF'
import sys, yaml
spec, out = sys.argv[1], sys.argv[2]
d = yaml.safe_load(open(spec))
pip = [x for y in d["dependencies"] if isinstance(y, dict) for x in y.get("pip", [])]
if not pip:
    sys.exit("no pip section found in " + spec)
open(out, "w").write("\n".join(pip) + "\n")
print(f"pip section: {len(pip)} requirements")
PYEOF
[[ -s "$PIPREQ" ]] || die "failed to extract the pip section from $YML"

# Build the conda side WITHOUT the pip section: passing the yml directly would
# let mamba run pip once, with isolation, which is the failure.
CONDAYML="$(mktemp -t scplus-conda.XXXXXX).yml"
trap 'rm -f "$PIPREQ" "${PIPREQ}.pypi" "${PIPREQ}.git" "$CONDAYML"' EXIT
python3 - "$YML" "$CONDAYML" <<'PYEOF'
import sys, yaml
spec, out = sys.argv[1], sys.argv[2]
d = yaml.safe_load(open(spec))
d["dependencies"] = [x for x in d["dependencies"] if isinstance(x, str)]
yaml.safe_dump(d, open(out, "w"), default_flow_style=False, sort_keys=False)
PYEOF

if (( DRY )); then
    echo "--dry-run: would run, in this order"
    echo
    echo "  [1/5]  $SOLVER env create -f <conda-only copy of $YML>"
    echo "         (the pip section is stripped: letting mamba run it would"
    echo "          install everything under isolation, which is the failure)"
    echo "  [2/5]  pip install --no-deps --no-build-isolation \\"
    echo "             ${NOISO[*]}"
    echo "  [3/5]  pip install -r <PyPI lines of the pip section>"
    echo "         ($(grep -vc 'git+' "$PIPREQ") pinned requirements, isolation ON)"
    echo "  [4/5]  pip install --no-deps -r <git+ lines of the pip section>"
    echo "         ($(grep -c 'git+' "$PIPREQ") packages; --no-deps because"
    echo "          scenicplus pins LoomXpy@main against our commit pin, which"
    echo "          pip calls ResolutionImpossible, and because resolving them"
    echo "          freely overrides our own pins)"
    echo "  [5/5]  pip check -> install whatever it reports missing, then"
    echo "         verify the pinned versions survived"
    echo
    echo "  PIP_CONSTRAINT=$CONS is exported for every pip pass."
    exit 0
fi

log "[1/5] conda section (${ENV_NAME})"
"$SOLVER" env create -f "$CONDAYML" || die "conda solve failed"

PY="$("$SOLVER" run -n "$ENV_NAME" python -c 'import sys; print(sys.executable)' 2>/dev/null)"
[[ -n "$PY" ]] || die "cannot locate python in the new env"
PIPV="$("$PY" -m pip --version | awk '{print $2}')"
log "env python: $PY  (pip $PIPV)"
# The pip major version decides which build paths exist. Recorded rather than
# gated on: pass A works on both sides of the boundary, so this is diagnostic.
case "$PIPV" in
    2[0-5].*) log "  pip < 26: the legacy setup.py path still exists here." ;;
    *)        log "  pip >= 26: legacy setup.py path removed; pass A is what" \
                  "makes the five sdists build." ;;
esac

log "[2/5] pass A -- ${#NOISO[@]} sdists, --no-build-isolation --no-deps"
export PIP_CONSTRAINT="$REPO_ROOT/$CONS"
"$PY" -m pip install --no-cache-dir --no-deps --no-build-isolation \
    --disable-pip-version-check "${NOISO[@]}"
A_RC=$?

# Pass B is split, because the five git dependencies cannot be resolved
# together with our pins. scenicplus v1.0a2's requirements.txt declares
#     loomxpy @ git+https://github.com/aertslab/LoomXpy@main
# and pycisTopic declares the same, while we pin a specific COMMIT. pip treats
# two different URLs for one name as an irreconcilable conflict:
#     ResolutionImpossible ... The user requested loomxpy 0.4.2 (from ...@<sha>)
#     scenicplus 1.0a2 depends on loomxpy 0.4.2 (from ...@main)
# There is no version range that resolves this -- the URLs differ, so pip
# refuses regardless of the versions matching.
#
# Blanket --no-deps is NOT the answer: the five git packages declare 196
# requirements that our own pinned list does not cover, and skipping them
# leaves an env that imports but breaks at use.
#
# There is a second, independent reason not to let the resolver near the git
# group. Resolving pycisTopic WITH its declared dependencies silently overrides
# our pins -- measured locally on pip 26.2.1:
#     scanpy       1.8.2         -> 1.11.0
#     anndata      0.10.5.post1  -> 0.11.4
#     scikit-learn 1.3.2         -> 1.5.2
#     scipy        1.12.0        -> 1.17.1
#     numba        0.59.0        -> 0.67.0
#     polars       0.20.13       -> 1.44.0
#     pyranges     0.0.111       -> 0.0.127
#     mudata       0.2.3         -> dropped entirely
# Those pins are transcribed from scenicplus's own requirements.txt lock, so
# losing them is losing the thing this env file exists to reproduce.
#
# So: pinned PyPI requirements first (isolation ON, resolver free within our
# pins), then the git packages with --no-deps -- their pinned refs ARE the
# specification, and the moving @main refs are exactly what we refuse to honour
# -- then ask pip what is missing and install only that. `pip check` names
# every unsatisfied requirement, so the repair is measured rather than guessed.
#
# Note which conflict is real: pycisTopic + our pinned loomxpy resolves fine
# (pycisTopic declares loomxpy with no URL). Only scenicplus pins
# LoomXpy@main, and pip treats two URLs for one name as irreconcilable
# regardless of the versions matching.
log "[3/5] pass B1 -- pinned PyPI requirements (isolation ON)"
grep -v "git+" "$PIPREQ" > "${PIPREQ}.pypi"
"$PY" -m pip install --no-cache-dir --disable-pip-version-check -r "${PIPREQ}.pypi"
B1_RC=$?

log "[4/5] pass B2 -- the git packages, --no-deps (pinned refs are the spec)"
grep "git+" "$PIPREQ" > "${PIPREQ}.git"
"$PY" -m pip install --no-cache-dir --no-deps --disable-pip-version-check \
    -r "${PIPREQ}.git"
B2_RC=$?

log "[5/5] repair -- install what pip check reports missing"
MISSING="$("$PY" -m pip check 2>&1 |
    sed -nE 's/^.* requires ([A-Za-z0-9._-]+), which is not installed\.$/\1/p' |
    sort -u | tr '\n' ' ')"
if [[ -n "${MISSING// /}" ]]; then
    log "  missing: ${MISSING}"
    # No pins here: these are transitive requirements the upstream packages
    # declare without pinning, so pip resolves them. The pinned list from B1 is
    # already installed and constrains what pip may pick.
    "$PY" -m pip install --no-cache-dir --disable-pip-version-check ${MISSING}
    R_RC=$?
    log "  re-checking"
    "$PY" -m pip check || log "  pip check still reports issues -- read them above"
else
    log "  nothing missing"
    R_RC=0
fi

# Did anything override the pins? The repair step resolves freely, and the git
# packages' own metadata wants newer versions of several of these. Report drift
# rather than assume it did not happen.
log "verifying the pinned versions survived"
"$PY" - "$PIPREQ" <<'PYEOF'
import subprocess, sys, re
want = {}
for line in open(sys.argv[1]):
    line = line.strip()
    if "==" in line and "git+" not in line:
        n, v = line.split("==", 1)
        want[n.strip().lower().replace("_", "-")] = v.strip()
have = {}
out = subprocess.run([sys.executable, "-m", "pip", "list", "--format=freeze",
                      "--disable-pip-version-check"],
                     capture_output=True, text=True).stdout
for line in out.splitlines():
    if "==" in line:
        n, v = line.split("==", 1)
        have[n.strip().lower().replace("_", "-")] = v.strip()
drift = [(n, w, have.get(n, "MISSING")) for n, w in sorted(want.items())
         if have.get(n) != w]
if drift:
    print(f"  {len(drift)} of {len(want)} pinned versions differ from the spec:")
    for n, w, h in drift:
        print(f"    {n:<22} spec {w:<16} installed {h}")
    print("  The env may still work, but it is no longer the pinned "
          "configuration.")
else:
    print(f"  all {len(want)} pinned versions match the spec")
PYEOF

B_RC=$(( B1_RC != 0 ? B1_RC : (B2_RC != 0 ? B2_RC : R_RC) ))

RC=$(( A_RC != 0 ? A_RC : B_RC ))

echo
echo "=============================================================="
if [[ "$RC" -eq 0 ]]; then
    echo "env created (pass A exit $A_RC, pass B exit $B_RC)."
    echo "Smoke test it before submitting:"
    echo "    conda activate ${ENV_NAME}"
    echo "    scenicplus --help"
    echo "    python -c \"import pycisTopic, pycistarget, scenicplus; print('ok')\""
    echo "    pip check"
    echo
    echo "Then:  sbatch slurm/scenicplus.sbatch"
else
    echo "build FAILED (pass A exit $A_RC, pass B exit $B_RC)"
    echo
    if [[ "$A_RC" -ne 0 ]]; then
        echo "Pass A failed, so a package built WITHOUT isolation is missing a"
        echo "build dependency from the conda section. --no-build-isolation"
        echo "installs nothing: whatever its setup.py imports must already be"
        echo "in the env. Add it to the conda section of $YML."
    else
        echo "Pass B failed (B1 exit $B1_RC, B2 exit $B2_RC, repair exit $R_RC)."
        echo
        echo "B1 -- 'setuptools was not found' for a package NOT in the NOISO"
        echo "array means that package also needs building without isolation:"
        echo "add its exact pin to NOISO. The tell in the log is 'Installing"
        echo "build dependencies' (isolated) versus 'Requirement already"
        echo "satisfied' (pass A handled it)."
        echo
        echo "B2 -- ResolutionImpossible here would be surprising, since it runs"
        echo "--no-deps. A git ref that no longer resolves is the likelier"
        echo "cause; 03_pipeline/preflight_env.sh checks reachability."
        echo
        echo "repair -- a package pip check named could not be installed. Read"
        echo "the name: if it is a conda-forge package, adding it to the conda"
        echo "section of $YML is usually better than letting pip build it."
    fi
fi
echo "=============================================================="
exit "$RC"
