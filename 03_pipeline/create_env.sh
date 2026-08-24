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
REPAIR=0
case "${1:-}" in
    --dry-run)     DRY=1 ;;
    --repair-pins) REPAIR=1 ;;
    "")            ;;
    *)             echo "usage: $0 [--dry-run|--repair-pins]" >&2; exit 2 ;;
esac

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

[[ -f "$YML" ]]  || die "$YML not found (run from the repo, or let this script cd)"
[[ -f "$CONS" ]] || die "$CONS not found -- it carries the setuptools bound"

SOLVER=""
for c in mamba micromamba conda; do
    command -v "$c" >/dev/null 2>&1 && { SOLVER="$c"; break; }
done
[[ -n "$SOLVER" ]] || die "no mamba/micromamba/conda on PATH"


# Locate an env's python WITHOUT parsing `<solver> env list`.
#
# That listing can fail for reasons unrelated to the envs themselves -- an
# unreadable ~/.condarc is enough -- and it fails by printing nothing useful to
# stdout while still exiting 0. Grepping it then reports "env does not exist"
# for an env that is not only present but ACTIVE, which is exactly what
# happened: `bash create_env.sh --repair-pins` refused inside an activated
# scenicplus prompt.
#
# The filesystem cannot lie about this. Search, in order: the active env, the
# solver's own envs dir, and the standard roots.
env_python() {
    local name="$1" c
    # Already activated? Then $CONDA_PREFIX is authoritative.
    if [[ "${CONDA_DEFAULT_ENV:-}" == "$name" && -x "${CONDA_PREFIX:-}/bin/python" ]]; then
        printf '%s\n' "${CONDA_PREFIX}/bin/python"; return 0
    fi
    local roots=()
    [[ -n "${MAMBA_ROOT_PREFIX:-}" ]] && roots+=("${MAMBA_ROOT_PREFIX}/envs")
    [[ -n "${CONDA_ENVS_PATH:-}" ]]   && roots+=("${CONDA_ENVS_PATH}")
    # `conda info --base` is a different code path from `env list` and usually
    # survives a bad rc file; fall back to the usual install locations.
    local base
    base="$("$SOLVER" info --base 2>/dev/null | tail -1)"
    [[ -n "$base" && -d "$base" ]] && roots+=("${base}/envs")
    roots+=("$HOME/miniconda3/envs" "$HOME/miniforge3/envs" "$HOME/mambaforge/envs"
            "$HOME/anaconda3/envs" "$HOME/micromamba/envs" "$HOME/.conda/envs")
    for c in "${roots[@]}"; do
        [[ -x "${c}/${name}/bin/python" ]] && { printf '%s\n' "${c}/${name}/bin/python"; return 0; }
    done
    return 1
}

ENV_NAME="$(python3 -c "
import sys
for line in open(sys.argv[1]):
    if line.startswith('name:'):
        print(line.split(':', 1)[1].strip()); break
" "$YML")"

if (( REPAIR )); then
    # Repair an EXISTING env in place, rather than an hour-long rebuild. Fixes
    # the two things the first real build got wrong: pinned versions that the
    # unconstrained repair step upgraded, and git packages that never landed.
    PY="$(env_python "$ENV_NAME")" || die "cannot find a python for env '${ENV_NAME}'.
     Searched \$CONDA_PREFIX (if active), \$MAMBA_ROOT_PREFIX/envs,
     \`${SOLVER} info --base\`/envs and the usual ~/ install roots.
     If the env lives elsewhere, activate it and re-run:
       conda activate ${ENV_NAME}
       bash 03_pipeline/create_env.sh --repair-pins"
    echo "=============================================================="
    echo "repairing ${ENV_NAME} in place"
    echo "=============================================================="
    echo "  python: $PY"
    echo
    # Confirm this really is the SCENIC+ env before force-reinstalling 34
    # pinned packages into it. env_python() resolves by NAME, and a name can
    # point somewhere unexpected -- during development this branch started
    # rewriting an unrelated env that happened to match. Two cheap signals:
    # the interpreter version the spec pins, and at least one package only this
    # env would carry.
    PYV="$("$PY" -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])' 2>/dev/null)"
    WANTV="$(grep -oE '^[[:space:]]+- python=[0-9.]+' "$YML" | head -1 |
             grep -oE '[0-9]+\.[0-9]+\.[0-9]+')"
    echo "  python version: ${PYV}  (spec pins ${WANTV})"
    if [[ -n "$WANTV" && "$PYV" != "$WANTV" ]]; then
        die "that interpreter is ${PYV}, but ${YML} pins python=${WANTV}.
     Refusing to force-reinstall 34 pinned packages into an env that is
     probably not the one this spec built. If it IS the right env, the
     python pin has drifted and a rebuild is the honest fix."
    fi

    PIPREQ="$(mktemp -t scplus-pip.XXXXXX)"
    trap 'rm -f "$PIPREQ" "${PIPREQ}-pypi" "${PIPREQ}-gitreqs" "${PIPREQ}-cons"' EXIT
    python3 - "$YML" "$PIPREQ" <<'PYEOF'
import sys, yaml
d = yaml.safe_load(open(sys.argv[1]))
pip = [x for y in d["dependencies"] if isinstance(y, dict) for x in y.get("pip", [])]
open(sys.argv[2], "w").write("\n".join(pip) + "\n")
PYEOF
    grep -v "git+" "$PIPREQ" > "${PIPREQ}-pypi"
    grep    "git+" "$PIPREQ" > "${PIPREQ}-gitreqs"
    cat "${PIPREQ}-pypi" "$REPO_ROOT/$CONS" > "${PIPREQ}-cons"

    # WHY THIS IS NO LONGER --force-reinstall ON ALL 34
    # -------------------------------------------------
    # The first version ran `pip install --force-reinstall --no-deps -r <all 34
    # pins>`. That was wrong twice over:
    #
    #   * It omitted --no-build-isolation, and five of those 34 are wheel-less
    #     sdists (pybedtools, pyranges, pyrle, tspex, MACS2). On pip 26 they
    #     cannot build with isolation on -- the same failure this whole wrapper
    #     exists to avoid -- so the step could not succeed as written.
    #   * --force-reinstall uninstalls each package before reinstalling it. A
    #     failure partway through a 34-package list therefore leaves the
    #     already-processed ones REMOVED. After that run, mudata, sklearn and
    #     matplotlib no longer imported and 33 of 34 pins reported absent.
    #
    # A repair should touch only what is actually wrong. So: read the current
    # state, act on the difference, and give the sdists the flag they need.
    log "[1/3] finding which pins actually drifted"
    DRIFTED="$("$PY" - "$PIPREQ" <<'PYEOF'
import sys
from importlib.metadata import distributions, version, PackageNotFoundError
want = {}
for line in open(sys.argv[1]):
    line = line.strip()
    if "==" in line and "git+" not in line:
        n, v = line.split("==", 1)
        want[n.strip()] = v.strip()
norm = lambda t: t.strip().lower().replace("_", "-")
meta = {}
for d in distributions():
    nm = d.metadata["Name"]
    if nm:
        meta[norm(nm)] = d.version
out = []
for name, w in want.items():
    got = meta.get(norm(name))
    if got is None:
        try:
            got = version(name)
        except PackageNotFoundError:
            got = None
    if got != w:
        out.append(f"{name}=={w}")
print(" ".join(out))
PYEOF
)"
    if [[ -z "${DRIFTED// /}" ]]; then
        log "  nothing drifted -- the pins already match the spec"
    else
        n_drift=$(wc -w <<<"$DRIFTED")
        log "  ${n_drift} to reinstall"
        # Split them: the wheel-less sdists need --no-build-isolation, the rest
        # do not care. No --force-reinstall -- a plain pinned install downgrades
        # an installed package and leaves everything else alone (verified).
        SD=(); WH=()
        for spec in $DRIFTED; do
            base="$(cut -d= -f1 <<<"$spec" | tr '[:upper:]' '[:lower:]')"
            case "$base" in
                pybedtools|pyranges|pyrle|tspex|macs2) SD+=("$spec") ;;
                *)                                    WH+=("$spec") ;;
            esac
        done
        export PIP_CONSTRAINT="${PIPREQ}-cons"
        P1=0
        if (( ${#WH[@]} )); then
            log "  wheels (${#WH[@]}): ${WH[*]}"
            "$PY" -m pip install --no-cache-dir --disable-pip-version-check \
                "${WH[@]}" || P1=$?
        fi
        if (( ${#SD[@]} )); then
            log "  sdists (${#SD[@]}), --no-build-isolation: ${SD[*]}"
            "$PY" -m pip install --no-cache-dir --no-build-isolation \
                --disable-pip-version-check "${SD[@]}" || P1=$?
        fi
        if (( P1 != 0 )); then
            die "pin reinstall failed (exit ${P1}). Nothing was force-removed,
     so the env is no worse than before this command. Read the pip error:
     if a pin genuinely cannot be satisfied, that pin has to move, and a
     clean rebuild is the honest fix."
        fi
    fi
    P1=0

    log "[2/3] git packages -- reinstalling only those that do not import"
    MISSING_GIT=()
    for mod in pycisTopic pycistarget scenicplus loomxpy pyscenic; do
        "$PY" -c "import $mod" >/dev/null 2>&1 || MISSING_GIT+=("$mod")
    done
    if (( ${#MISSING_GIT[@]} == 0 )); then
        log "  all five import -- leaving them alone"
        P2=0
    else
        log "  not importing: ${MISSING_GIT[*]} -- reinstalling all five"
        "$PY" -m pip install --no-cache-dir --no-deps \
            --disable-pip-version-check -r "${PIPREQ}-gitreqs"
        P2=$?
    fi

    log "[3/3] repairing what pip check now reports, WITHIN the pins"
    M="$("$PY" -m pip check 2>&1 |
        sed -nE 's/^.* requires ([A-Za-z0-9._-]+), which is not installed\.$/\1/p' |
        sort -u | tr '\n' ' ')"
    if [[ -n "${M// /}" ]]; then
        log "  missing: ${M}"
        "$PY" -m pip install --no-cache-dir --disable-pip-version-check ${M}
    else
        log "  nothing missing"
    fi
    P3=$?

    echo
    log "re-verifying"
    "$PY" 03_pipeline/_check_pins.py "$PIPREQ"
    D=$?
    for mod in pycisTopic pycistarget scenicplus loomxpy pyscenic; do
        "$PY" -c "import $mod" 2>/dev/null && echo "  ok      import ${mod}" ||
            { echo "  MISSING import ${mod}"; D=1; }
    done
    if [[ -x "${PY%/python}/scenicplus" ]]; then
        echo "  ok      ${PY%/python}/scenicplus exists"
    else
        echo "  MISSING ${PY%/python}/scenicplus"; D=1
    fi
    RC=$(( P1 != 0 ? P1 : (P2 != 0 ? P2 : (P3 != 0 ? P3 : D)) ))
    echo
    echo "=============================================================="
    (( RC == 0 )) && echo "repair complete -- smoke test, then sbatch." ||
        echo "repair INCOMPLETE (exit $RC) -- read the drift/import lines above."
    echo "=============================================================="
    exit "$RC"
fi

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

if env_python "$ENV_NAME" >/dev/null 2>&1; then
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
trap 'rm -f "$PIPREQ" "${PIPREQ}-pypi" "${PIPREQ}-gitreqs" "${PIPREQ}-cons" "$CONDAYML"' EXIT
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

PY="$(env_python "$ENV_NAME")" ||
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
grep -v "git+" "$PIPREQ" > "${PIPREQ}-pypi"
"$PY" -m pip install --no-cache-dir --disable-pip-version-check -r "${PIPREQ}-pypi"
B1_RC=$?

log "[4/5] pass B2 -- the git packages, --no-deps (pinned refs are the spec)"
grep "git+" "$PIPREQ" > "${PIPREQ}-gitreqs"
"$PY" -m pip install --no-cache-dir --no-deps --disable-pip-version-check \
    -r "${PIPREQ}-gitreqs"
B2_RC=$?

log "[5/5] repair -- install what pip check reports missing"
MISSING="$("$PY" -m pip check 2>&1 |
    sed -nE 's/^.* requires ([A-Za-z0-9._-]+), which is not installed\.$/\1/p' |
    sort -u | tr '\n' ' ')"
if [[ -n "${MISSING// /}" ]]; then
    log "  missing: ${MISSING}"
    # CONSTRAIN THE REPAIR TO OUR PINS. Installed-already is NOT a constraint:
    # pip will happily upgrade an installed package to satisfy a new one. On the
    # first real run this step took pandas 1.5.0 -> 3.0.5, matplotlib
    # 3.6.3 -> 3.11.1, polars 0.20.13 -> 1.44.0, statsmodels and tqdm too --
    # a pandas MAJOR version jump against a spec that pins 1.5.0.
    #
    # PIP_CONSTRAINT holds them: verified that an explicit `pip install
    # --upgrade pandas` under a `pandas==1.5.3` constraint leaves the installed
    # version alone. If a missing package genuinely cannot coexist with the
    # pins, pip now fails here and says so, instead of silently rewriting the
    # environment the spec describes.
    cat "${PIPREQ}-pypi" "$REPO_ROOT/$CONS" > "${PIPREQ}-cons"
    PIP_CONSTRAINT="${PIPREQ}-cons" \
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
"$PY" 03_pipeline/_check_pins.py "$PIPREQ"
DRIFT_RC=$?

# The git packages install --no-deps, so nothing verifies they arrived. And the
# console script is the thing the pipeline actually invokes: scenicplus's
# pyproject.toml declares [project.scripts] scenicplus =
# "scenicplus.cli.scenicplus:main", so `scenicplus --help` failing means the
# package did not install even if pass B exited 0.
log "verifying the git packages and the CLI"
CLI_RC=0
for mod in pycisTopic pycistarget scenicplus loomxpy pyscenic; do
    if "$PY" -c "import $mod" 2>/dev/null; then
        v="$("$PY" -c "import $mod,sys; print(getattr($mod,'__version__','?'))" 2>/dev/null)"
        echo "  ok      import ${mod}  (${v})"
    else
        echo "  MISSING import ${mod}"
        CLI_RC=1
    fi
done
if "$PY" -m scenicplus --help >/dev/null 2>&1 ||
   "${PY%/python}/scenicplus" --help >/dev/null 2>&1; then
    echo "  ok      scenicplus CLI responds"
else
    echo "  MISSING scenicplus CLI -- 'scenicplus --help' will not work."
    echo "          The entry point IS declared upstream, so this means the"
    echo "          package itself is absent or its console script was not"
    echo "          written. Check: ${PY%/python}/scenicplus"
    CLI_RC=1
fi

B_RC=$(( B1_RC != 0 ? B1_RC : (B2_RC != 0 ? B2_RC : R_RC) ))
# Drift and a missing CLI are build FAILURES, not notes: an env with pandas 3
# where the spec pins 1.5.0 will fail somewhere inside a multi-hour pipeline
# run instead of here.
(( B_RC == 0 )) && B_RC=$(( DRIFT_RC != 0 ? DRIFT_RC : CLI_RC ))

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
