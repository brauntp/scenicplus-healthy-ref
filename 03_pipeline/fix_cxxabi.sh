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

echo "-- installing libstdcxx-ng + libgcc-ng ---------------------------"
"$SOLVER" install -y -p "$PREFIX" -c conda-forge libstdcxx-ng libgcc-ng
RC=$?
echo

if (( RC != 0 )); then
    echo "=============================================================="
    echo "install FAILED (exit $RC). Nothing was removed."
    echo "If the solver refuses because of a conflict, a fresh build with the"
    echo "corrected 03_pipeline/environment.yml is the fallback."
    echo "=============================================================="
    exit "$RC"
fi

ENVLIB="$(ls -1 "${PREFIX}"/lib/libstdc++.so.6* 2>/dev/null | head -1)"
if [[ -n "$ENVLIB" ]]; then
    echo "  now present: $ENVLIB"
    echo "  highest CXXABI: $(highest_cxxabi "$ENVLIB")"
else
    echo "  WARNING: still no libstdc++ under ${PREFIX}/lib after install."
fi
echo

echo "-- re-testing the seven stages that failed ----------------------"
FAIL=0
for mod in pycistarget.motif_enrichment_cistarget \
           pycistarget.motif_enrichment_dem \
           pycistarget.input_output \
           pycistarget.motif_enrichment_result \
           scenicplus.grn_builder.gsea_approach \
           scenicplus.grn_builder.modules \
           scenicplus.eregulon_enrichment; do
    if err="$("$PY" -c "import ${mod}" 2>&1)"; then
        printf "  ok    %s\n" "$mod"
    else
        printf "  FAIL  %s\n" "$mod"
        printf "        %s\n" "$(tail -1 <<<"$err" | head -c 140)"
        FAIL=1
    fi
done
echo

echo "=============================================================="
if (( FAIL == 0 )); then
    echo "FIXED. Confirm the whole thing, then proceed:"
    echo "    bash 03_pipeline/smoke_test.sh"
    echo "    bash 03_pipeline/make_config.sh"
    echo "    sbatch slurm/scenicplus.sbatch"
else
    echo "STILL FAILING. Two things left to distinguish, and"
    echo "03_pipeline/diagnose_cxxabi.sh reports both:"
    echo "  * the env now HAS a new enough libstdc++ but the linker prefers"
    echo "    the system one -- \$LD_LIBRARY_PATH, fixable per-job"
    echo "  * the offending wheel wants an ABI newer than even conda-forge's"
    echo "    libstdcxx provides -- then pin that package back to the version"
    echo "    scenicplus's requirements.txt specifies, which is also a step"
    echo "    toward the upstream lock"
fi
echo "=============================================================="
exit "$FAIL"
