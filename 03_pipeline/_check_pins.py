import subprocess, sys

def installed():
    """Installed distributions, by importlib.metadata -- the same mechanism
    `import` uses. Preferred over parsing `pip list`: on the first real repair
    run pip list reported 33 of 34 packages absent while every dependent
    package still imported, which is impossible. Whatever caused that, asking
    the metadata directly removes the middleman."""
    out = {}
    try:
        from importlib.metadata import distributions
        for d in distributions():
            n = (d.metadata["Name"] or "").strip().lower().replace("_", "-")
            if n:
                out[n] = (d.version or "?").strip()
    except Exception as e:
        print(f"  importlib.metadata failed: {e}", file=sys.stderr)
    if len(out) < 20:
        # A real env here has hundreds. A near-empty answer is a broken query,
        # not a broken environment -- say so instead of listing every pin as
        # missing.
        r = subprocess.run([sys.executable, "-m", "pip", "list",
                            "--format=freeze", "--disable-pip-version-check"],
                           capture_output=True, text=True)
        for line in r.stdout.splitlines():
            if "==" in line:
                n, v = line.split("==", 1)
                out.setdefault(n.strip().lower().replace("_", "-"), v.strip())
    return out

want = {}
for line in open(sys.argv[1]):
    line = line.strip()
    if "==" in line and "git+" not in line:
        n, v = line.split("==", 1)
        want[n.strip().lower().replace("_", "-")] = v.strip()

have = installed()
print(f"  {len(have)} distributions visible to {sys.executable}")
if len(have) < 20:
    print("  THAT COUNT IS IMPLAUSIBLE for this env -- treating it as a failed")
    print("  query, not 34 missing packages. Check by hand:")
    print(f"    {sys.executable} -m pip list | wc -l")
    print(f"    {sys.executable} -c 'import pandas; print(pandas.__version__)'")
    sys.exit(4)

drift = [(n, w, have.get(n, "MISSING")) for n, w in sorted(want.items())
         if have.get(n) != w]
missing = [d for d in drift if d[2] == "MISSING"]
if len(missing) > len(want) // 2:
    print(f"  {len(missing)} of {len(want)} pins report MISSING while "
          f"{len(have)} distributions are installed.")
    print("  That combination is suspicious rather than informative -- a real")
    print("  env with the git packages importing cannot be missing pandas.")
    print("  Sample of what IS installed:")
    for n in sorted(have)[:8]:
        print(f"    {n}=={have[n]}")
if drift:
    print(f"  {len(drift)} of {len(want)} pinned versions differ from the spec:")
    for n, w, h in drift:
        print(f"    {n:<22} spec {w:<16} installed {h}")
    sys.exit(3)
print(f"  all {len(want)} pinned versions match the spec")
