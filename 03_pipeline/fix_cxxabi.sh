#!/usr/bin/env bash
# =============================================================================
# Install the runtime C++ library into an EXISTING env, so compiled extensions
# stop falling back to the system /lib64/libstdc++.so.6.
#
# THE PROBLEM
# -----------
# smoke_test.sh reported, for seven of eleven pipeline stages:
#     ImportError: /lib64/libstdc++.so.6: version `CXXABI_1.3.15' not found
#
# PyPI wheels for the compiled packages in this stack (sorted-nearest, ncls,
# pyrle and their kin, reached via pyranges) are built against newer GCC than a
# RHEL-family login node ships. The env's conda section listed c-compiler and
# cxx-compiler -- the toolchain for BUILDING -- but not libstdcxx-ng, the
# runtime those built extensions LOAD. With no libstdc++ under the env prefix,
# the dynamic linker uses the system one, and that one is too old.
#
# 03_pipeline/environment.yml now lists libstdcxx-ng and libgcc-ng, so a fresh
# build gets this right. This script fixes an env that already exists, in
# seconds, without the hour-long rebuild.
#
# Usage:
#     conda activate scenicplus
#     bash 03_pipeline/fix_cxxabi.sh
# =============================================================================
set -uo pipefail

PREFIX="${CONDA_PREFIX:-}"
[[ -n "$PREFIX" && -x "${PREFIX}/bin/python" ]] || {
    echo "ERROR: activate the env first (conda activate scenicplus)" >&2; exit 1; }
PY="${PREFIX}/bin/python"

SOLVER=""
for c in mamba micromamba conda; do command -v "$c" >/dev/null && { SOLVER="$c"; break; }; done
[[ -n "$SOLVER" ]] || { echo "ERROR: no mamba/micromamba/conda on PATH" >&2; exit 1; }

highest_cxxabi() {  # $1 = path to a libstdc++.so
    strings "$1" 2>/dev/null | grep -oE '^CXXABI_1\.3\.[0-9]+$' | sort -V | tail -1
}

echo "=============================================================="
echo "libstdc++ runtime fix"
echo "=============================================================="
echo "  env    : $PREFIX"
echo "  solver : $SOLVER"
echo

SYS="$(readlink -f /lib64/libstdc++.so.6 2>/dev/null ||
       readlink -f /usr/lib64/libstdc++.so.6 2>/dev/null)"
if [[ -n "$SYS" ]]; then
    echo "  system : $SYS"
    echo "           highest CXXABI $(highest_cxxabi "$SYS")"
fi

ENVLIB="$(ls -1 "${PREFIX}"/lib/libstdc++.so.6* 2>/dev/null | head -1)"
if [[ -n "$ENVLIB" ]]; then
    echo "  env    : $ENVLIB"
    echo "           highest CXXABI $(highest_cxxabi "$ENVLIB")"
    echo
    echo "  The env already ships its own libstdc++. If imports still fail, the"
    echo "  linker is not preferring it -- that is a different problem from a"
    echo "  missing library, and \$LD_LIBRARY_PATH is the lever. Check:"
    echo "      echo \"\$LD_LIBRARY_PATH\""
    echo "  Anything ahead of ${PREFIX}/lib there will win."
else
    echo "  env    : NO libstdc++ under the env prefix -- this is the cause."
fi
echo

MODS=(pycistarget.motif_enrichment_cistarget
      pycistarget.motif_enrichment_dem
      pycistarget.input_output
      pycistarget.motif_enrichment_result
      scenicplus.grn_builder.gsea_approach
      scenicplus.grn_builder.modules
      scenicplus.eregulon_enrichment)

# Import the seven affected modules under a given LD_LIBRARY_PATH, print a
# per-module verdict, and return the failure count.
try_imports() {   # $1 = LD_LIBRARY_PATH value ("" for unset), $2 = label
    local llp="$1" label="$2" n=0 mod err
    echo "  [${label}]"
    for mod in "${MODS[@]}"; do
        if [[ -z "$llp" ]]; then
            err="$(env -u LD_LIBRARY_PATH "$PY" -c "import ${mod}" 2>&1)"
        else
            err="$(env LD_LIBRARY_PATH="$llp" "$PY" -c "import ${mod}" 2>&1)"
        fi
        if [[ $? -eq 0 ]]; then
            printf "    ok    %s\n" "$mod"
        else
            printf "    FAIL  %s\n" "$mod"
            printf "          %s\n" "$(tail -1 <<<"$err" | head -c 160)"
            n=$((n + 1))
        fi
    done
    return "$n"
}

if [[ -z "$ENVLIB" ]]; then
    echo "-- installing libstdcxx-ng + libgcc-ng ---------------------------"
    "$SOLVER" install -y -p "$PREFIX" -c conda-forge libstdcxx-ng libgcc-ng
    RC=$?
    if (( RC != 0 )); then
        echo
        echo "install FAILED (exit $RC). Nothing was removed. If the solver"
        echo "refuses on a conflict, a fresh build with the corrected"
        echo "03_pipeline/environment.yml is the fallback."
        exit "$RC"
    fi
    ENVLIB="$(ls -1 "${PREFIX}"/lib/libstdc++.so.6* 2>/dev/null | head -1)"
    if [[ -n "$ENVLIB" ]]; then
        echo "  now present: $ENVLIB  (highest CXXABI $(highest_cxxabi "$ENVLIB"))"
    else
        echo "  WARNING: still no libstdc++ under ${PREFIX}/lib after install."
    fi
    echo
else
    echo "-- libstdcxx-ng already present, skipping the install -------------"
    echo
fi

# THE SECOND HALF OF THE PROBLEM
# ------------------------------
# Having the library is necessary, not sufficient. A conda-BUILT extension
# carries RPATH $ORIGIN/../../.. and finds $CONDA_PREFIX/lib by itself. A
# manylinux WHEEL from PyPI carries no such RPATH, so the dynamic linker falls
# through to ld.so.cache and /lib64 -- the old one. The failing modules here all
# come from PyPI wheels (sorted-nearest, ncls, pyrle, reached via pyranges),
# which is why installing the library alone can leave the error unchanged.
#
# So test both, and let the result decide rather than assuming.
echo "-- does it import? testing with and without the env lib dir -------"
try_imports "" "as-is"
N_PLAIN=$?
echo
try_imports "${PREFIX}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "LD_LIBRARY_PATH=\$CONDA_PREFIX/lib"
N_LLP=$?
echo
echo "  as-is: ${N_PLAIN}/${#MODS[@]} failed    with env lib dir: ${N_LLP}/${#MODS[@]} failed"
echo

FAIL=0
if (( N_PLAIN == 0 )); then
    echo "=============================================================="
    echo "FIXED, and no LD_LIBRARY_PATH needed. Proceed:"
    echo "    bash 03_pipeline/smoke_test.sh"
    echo "    bash 03_pipeline/make_config.sh && sbatch slurm/scenicplus.sbatch"
    echo "=============================================================="
    exit 0
elif (( N_LLP == 0 )); then
    # Make it permanent. A hook under the env's activate.d fires for every
    # `conda activate` INCLUDING inside a batch job, which is the case that
    # matters -- typing the export in a login shell would not reach sbatch.
    HOOK="${PREFIX}/etc/conda/activate.d/zzz_libstdcxx.sh"
    mkdir -p "$(dirname "$HOOK")"
    cat > "$HOOK" <<'HOOKEOF'
# Prepend this env's lib dir so PyPI manylinux wheels find the conda
# libstdc++ instead of the system /lib64 one, which is too old for the
# CXXABI these wheels were built against. Conda-built extensions do not
# need this (they carry an RPATH); pip-installed ones do.
#
# Guarded: LD_LIBRARY_PATH is inherited by child processes and by batch
# jobs, so an unguarded prepend accumulates duplicates every time the env
# is activated inside an already-activated shell.
case ":${LD_LIBRARY_PATH:-}:" in
    *":${CONDA_PREFIX}/lib:"*) ;;
    *) export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}" ;;
esac
HOOKEOF
    echo "  wrote activation hook: $HOOK"
    echo "  (fires on every `conda activate`, including inside batch jobs)"
    echo
    echo "=============================================================="
    echo "FIXED via LD_LIBRARY_PATH, and made permanent."
    echo
    echo "The library was present but the linker preferred /lib64: PyPI"
    echo "manylinux wheels carry no RPATH into the env. The hook above"
    echo "prepends \$CONDA_PREFIX/lib on activation, so it applies in batch"
    echo "jobs too -- an export typed here would not have."
    echo
    echo "Re-activate to pick it up, then confirm:"
    echo "    conda deactivate && conda activate \$(basename $PREFIX)"
    echo "    bash 03_pipeline/smoke_test.sh"
    echo "=============================================================="
    exit 0
else
    FAIL=1
    echo "=============================================================="
    echo "STILL FAILING both ways (${N_PLAIN} and ${N_LLP} of ${#MODS[@]})."
    echo
    echo "That rules out the search-path explanation, so the wheel wants an"
    echo "ABI newer than this env's libstdc++ provides. Get the specifics:"
    echo "    bash 03_pipeline/diagnose_cxxabi.sh"
    echo "Its section 4 names the offending .so and its owning package, and"
    echo "section 5 shows that package against scenicplus's lock. Pinning it"
    echo "to the lock version is then both the fix and a step toward the"
    echo "upstream lock -- sorted-nearest is one of the 187 uncovered pins."
    echo "=============================================================="
fi
exit "$FAIL"
