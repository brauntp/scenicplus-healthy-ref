#!/usr/bin/env bash
# =============================================================================
# Which compiled extension needs a newer libstdc++ than the system provides?
#
# smoke_test.sh reported, for six of eleven pipeline stages:
#     ImportError: /lib64/libstdc++.so.6: version `CXXABI_1.3.15' not found
#                  (required by /home/groups/Maxs...   <- path truncated
#
# This is NOT a version-pin mismatch in the Python sense. A wheel was built
# against a newer GCC than this machine's system libstdc++ supports:
#     CXXABI_1.3.13 -> GCC 11    CXXABI_1.3.14 -> GCC 13
#     CXXABI_1.3.15 -> GCC 14    (what the failing wheel wants)
#
# Two candidate fixes, and which one is right depends on WHICH library it is:
#   A. the env's own libstdc++ (conda-forge ships a newer one) is present but
#      not being found first -- a pip-installed wheel links against whatever
#      the dynamic linker finds, and /lib64 wins unless told otherwise.
#   B. the wheel is a NEWER release than scenicplus's lock pins, and the
#      lock-pinned version predates the new ABI -- in which case pinning it is
#      both the fix and a step toward the lock.
#
# Read-only. Run with the env ACTIVE.
# =============================================================================
set -uo pipefail

PY="${CONDA_PREFIX:-}/bin/python"
[[ -x "$PY" ]] || { echo "activate the env first" >&2; exit 1; }
PREFIX="${CONDA_PREFIX}"
SITE="$("$PY" -c 'import site;print(site.getsitepackages()[0])')"

echo "=============================================================="
echo "CXXABI diagnosis"
echo "=============================================================="
echo "  env    : $PREFIX"
echo "  site   : $SITE"
echo

echo "-- [1] the full error, untruncated -------------------------------"
"$PY" -c "import pycistarget.motif_enrichment_cistarget" 2>&1 | tail -3 | fold -w 200
echo

echo "-- [2] what does the SYSTEM libstdc++ provide? -------------------"
for lib in /lib64/libstdc++.so.6 /usr/lib64/libstdc++.so.6; do
    [[ -e "$lib" ]] || continue
    echo "  $lib -> $(readlink -f "$lib")"
    echo "    highest CXXABI: $(strings "$(readlink -f "$lib")" 2>/dev/null |
        grep -oE '^CXXABI_1\.3\.[0-9]+$' | sort -V | tail -1)"
    echo "    highest GLIBCXX: $(strings "$(readlink -f "$lib")" 2>/dev/null |
        grep -oE '^GLIBCXX_3\.4\.[0-9]+$' | sort -V | tail -1)"
done
echo

echo "-- [3] does the ENV ship its own, and is it newer? ---------------"
found=0
for lib in "$PREFIX"/lib/libstdc++.so.6*; do
    [[ -e "$lib" ]] || continue
    found=1
    echo "  $lib"
    echo "    highest CXXABI: $(strings "$lib" 2>/dev/null |
        grep -oE '^CXXABI_1\.3\.[0-9]+$' | sort -V | tail -1)"
done
if (( found == 0 )); then
    echo "  NONE -- the env has no libstdc++ of its own."
    echo "  Fix A is unavailable until conda-forge's libstdcxx-ng is installed:"
    echo "      conda install -p ${PREFIX} -c conda-forge libstdcxx-ng"
fi
echo

echo "-- [4] which .so files demand CXXABI_1.3.15, and who owns them? --"
# Scan compiled extensions for the specific symbol version. Cheap: only .so
# files under site-packages, and `strings` not a full ldd walk.
hits=0
while IFS= read -r so; do
    if strings "$so" 2>/dev/null | grep -qx "CXXABI_1.3.15"; then
        hits=$((hits + 1))
        rel="${so#"$SITE"/}"
        dist="${rel%%/*}"
        ver="$("$PY" -c "
from importlib.metadata import distributions
t='${dist}'.lower().replace('_','-')
for d in distributions():
    n=(d.metadata['Name'] or '').lower().replace('_','-')
    if n==t or n.replace('-','')==t.replace('-',''):
        print(f'{d.metadata[\"Name\"]}=={d.version}'); break
else: print('(owner not in metadata)')" 2>/dev/null)"
    printf "  %-52s %s\n" "$rel" "$ver"
    fi
done < <(find "$SITE" -name "*.so" -type f 2>/dev/null)
if (( hits == 0 )); then
    echo "  none found by symbol scan -- the requirement may come from a"
    echo "  transitively loaded library outside site-packages. Section 1's"
    echo "  untruncated path names it."
fi
echo

echo "-- [5] those packages vs scenicplus's lock ------------------------"
"$PY" - <<'PYEOF'
import urllib.request
from importlib.metadata import distributions
lock = {}
try:
    t = urllib.request.urlopen(
        "https://raw.githubusercontent.com/aertslab/scenicplus/v1.0a2/"
        "requirements.txt", timeout=20).read().decode()
    for line in t.splitlines():
        line = line.strip()
        if "==" in line and "git+" not in line and not line.startswith("#"):
            n, v = line.split("==", 1)
            lock[n.strip().lower().replace("_", "-")] = v.strip()
except Exception as e:
    print(f"  could not fetch the upstream lock ({type(e).__name__}) -- "
          f"compare by hand against scenicplus v1.0a2 requirements.txt")
    raise SystemExit(0)
# The compiled packages in this stack, whatever the scan found.
watch = ["sorted-nearest", "ncls", "pyrle", "pyranges", "pybedtools",
         "pyarrow", "numba", "llvmlite", "ray", "tables", "scikit-image"]
have = {}
for d in distributions():
    n = (d.metadata["Name"] or "").lower().replace("_", "-")
    if n:
        have[n] = d.version
print(f"  {'package':<18}{'lock':<14}{'installed':<14}status")
for p in watch:
    l, h = lock.get(p, "-"), have.get(p, "ABSENT")
    st = "match" if l == h else ("NEWER THAN LOCK" if l != "-" and h != "ABSENT"
                                 else "")
    print(f"  {p:<18}{l:<14}{h:<14}{st}")
PYEOF
echo
echo "=============================================================="
echo "Send sections 1, 3 and 4. Section 1 names the library, 3 says"
echo "whether the env can supply a newer libstdc++ at all, and 4 says"
echo "whether the offender is a package we can simply pin back."
echo "=============================================================="
