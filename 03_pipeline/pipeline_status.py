#!/usr/bin/env python3
"""
Which pipeline outputs exist, which rule runs next, and is there room to write.

WHY THIS EXISTS
---------------
A run failed after 3.5 hours, immediately after logging

    R2G  INFO  Done!
    SCENIC+  INFO  Saving region to gene adjacencies to region_to_gene_adj.tsv

and the failure banner said only "pipeline FAILED (exit 1)". The first question
after a long partial run is not "why" but "what survived" -- snakemake resumes
from completed outputs, so a rule that finished is 71 minutes you do not pay
again. Nothing in this repo answered that.

This reads the resolved config, reports every output_data target with size and
age, names the next rule the DAG would run, and checks free space on the
filesystem the outputs land on (a truncated TSV from a full disk looks exactly
like a crash mid-write).

Read-only. Safe on a login node -- it stats files, it does not open them.

Usage:
    python 03_pipeline/pipeline_status.py --config 03_pipeline/config.yaml
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
import time
from pathlib import Path

# Rule -> the output_data keys that rule produces, in DAG order. Taken from the
# pinned v1.0a2 Snakefile; the order is the topological one snakemake follows,
# so the first rule with a missing output is where a resume starts.
STAGES: list[tuple[str, list[str]]] = [
    ("prepare_GEX_ACC (bypassed: pre-built)", ["combined_GEX_ACC_mudata"]),
    ("motif_enrichment_cistarget", ["ctx_result_fname", "output_fname_ctx_html"]),
    ("motif_enrichment_dem", ["dem_result_fname", "output_fname_dem_html"]),
    ("prepare_menr", ["cistromes_direct", "cistromes_extended", "tf_names"]),
    ("download_genome_annotations", ["genome_annotation", "chromsizes"]),
    ("get_search_space", ["search_space"]),
    ("tf_to_gene", ["tf_to_gene_adjacencies"]),
    ("region_to_gene", ["region_to_gene_adjacencies"]),
    ("eGRN_direct", ["eRegulons_direct"]),
    ("eGRN_extended", ["eRegulons_extended"]),
    ("AUCell_direct", ["AUCell_direct"]),
    ("AUCell_extended", ["AUCell_extended"]),
    ("scplus_mudata", ["scplus_mdata"]),
]


def human(n: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if n < 1024 or unit == "TiB":
            return f"{n:,.1f} {unit}" if unit != "B" else f"{int(n)} B"
        n /= 1024
    return f"{n:.1f} TiB"


def age(mtime: float) -> str:
    d = time.time() - mtime
    if d < 3600:
        return f"{d/60:.0f}m ago"
    if d < 86400:
        return f"{d/3600:.1f}h ago"
    return f"{d/86400:.1f}d ago"



def last_failure(workdir: Path) -> dict | None:
    """Name the rule that failed, from snakemake's own log.

    The job's .err file shows the traceback but not always which RULE owned it,
    and the .out file shows neither (snakemake writes to stderr). Snakemake's own
    log records `Error in rule <name>:` plus the shell exit code, which is the
    only place the two are stated together.

    Exit code 137 = SIGKILL = the OOM killer in practice. Distinguishing that
    from an ordinary non-zero exit matters, because the remedy is different:
    more memory versus a real bug.
    """
    logs = sorted((workdir / ".snakemake" / "log").glob("*.snakemake.log"),
                  key=lambda p: p.stat().st_mtime, reverse=True)
    if not logs:
        return None
    txt = logs[0].read_text(errors="replace")
    # run_pipeline.sh passes --keep-going, so ONE run can fail several rules on
    # independent branches -- and a later success (region_to_gene finishing
    # after another rule died) makes the tail of the log look healthy. Collect
    # every failure, paired with the exit code reported just before it.
    lines = txt.splitlines()
    fails = []
    pending_code = None
    for line in lines:
        m = re.search(r"returned non-zero exit status (\d+)", line)
        if m:
            pending_code = int(m.group(1))
        if line.startswith("Error in rule "):
            fails.append({"rule": line[len("Error in rule "):].rstrip(":").strip(),
                          "exit": pending_code,
                          "oom": pending_code in (137, 9)})
            pending_code = None
    # Snakemake repeats each failure in an end-of-run summary, where no exit
    # code precedes it -- so the same rule appears twice, once with a code and
    # once with None. Keep one entry per rule, preferring the one that carries
    # the code.
    dedup: dict[str, dict] = {}
    for f in fails:
        prev = dedup.get(f["rule"])
        if prev is None or (prev["exit"] is None and f["exit"] is not None):
            dedup[f["rule"]] = f
    fails = list(dedup.values())

    # Snakemake's own tally, which counts only rules that COMPLETED.
    m = re.findall(r"(\d+) of (\d+) steps \((\d+)%\) done", txt)
    tally = m[-1] if m else None
    return {"log": logs[0], "fails": fails, "tally": tally}


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--workdir", type=Path, default=None,
                    help="Snakemake's --directory. Defaults to the config's "
                         "own directory, which is what run_pipeline.sh uses. "
                         "Relative output_data paths resolve against THIS, not "
                         "the repo root.")
    args = ap.parse_args()

    try:
        import yaml
    except ImportError:
        sys.exit("ERROR: no PyYAML in this python.\n"
                 "       conda activate scenicplus   (or scplus-pairing)")

    if not args.config.is_file():
        sys.exit(f"ERROR: config not found: {args.config}\n"
                 f"       Run 03_pipeline/make_config.sh first.")
    cfg = yaml.safe_load(open(args.config))
    out = cfg.get("output_data", {})
    workdir = (args.workdir or args.config.parent).resolve()

    def resolve(key: str) -> Path | None:
        v = out.get(key)
        if not isinstance(v, str):
            return None
        p = Path(v)
        return p if p.is_absolute() else workdir / p

    print("=" * 74)
    print("pipeline status")
    print("=" * 74)
    print(f"  config  : {args.config}")
    print(f"  workdir : {workdir}")
    print()

    print("-- outputs, in DAG order -------------------------------------------")
    next_rule = None
    done = partial = 0
    for rule, keys in STAGES:
        marks = []
        for k in keys:
            p = resolve(k)
            if p is None:
                marks.append((k, "NOT IN CONFIG", None))
            elif p.exists():
                st = p.stat()
                marks.append((k, human(st.st_size), age(st.st_mtime)))
            else:
                marks.append((k, None, None))
        have = [m for m in marks if m[1] and m[1] != "NOT IN CONFIG"]
        if len(have) == len(marks):
            tag, done = "done", done + 1
        elif have:
            tag, partial = "PARTIAL", partial + 1
            if next_rule is None:
                next_rule = rule
        else:
            tag = "pending"
            if next_rule is None:
                next_rule = rule
        print(f"  [{tag:<7}] {rule}")
        for k, size, when in marks:
            p = resolve(k)
            nm = p.name if p else "?"
            if size and size != "NOT IN CONFIG":
                print(f"              {nm:<34} {size:>12}  {when}")
            elif size == "NOT IN CONFIG":
                print(f"              {k:<34} {'--':>12}  not in config")
            else:
                print(f"              {nm:<34} {'absent':>12}")
    print()

    print("-- last failure, from snakemake's own log --------------------------")
    fail = last_failure(workdir)
    if fail is None:
        print("  no snakemake log under .snakemake/log/ -- no run has started here.")
    else:
        print(f"  log : {fail['log']}")
        if fail["tally"]:
            d, t, pct = fail["tally"]
            print(f"  snakemake's tally: {d} of {t} steps ({pct}%) COMPLETED in "
                  f"that run")
            print(f"  (this counts steps that ran, not outputs that exist -- a")
            print(f"   rule satisfied before the run started is not counted)")
        if not fail["fails"]:
            print("  no 'Error in rule' lines: the last run did not fail at a")
            print("  rule. It either completed or was killed outright.")
        else:
            print()
            for f in fail["fails"]:
                oom = "  <- SIGKILL / OOM" if f["oom"] else ""
                print(f"  FAILED RULE: {f['rule']}   exit {f['exit']}{oom}")
            if any(f["oom"] for f in fail["fails"]):
                print()
                print("  Exit 137 is SIGKILL -- in a SLURM job that is the OOM")
                print("  killer, and the rules above are what to size --mem for.")
                print("  Heavy rules are SERIALISED (each declares threads =")
                print("  n_cpu; verified against snakemake), so --mem covers the")
                print("  largest single rule, not a sum.")
            if len(fail["fails"]) > 1:
                print()
                print("  MORE THAN ONE rule failed. run_pipeline.sh passes")
                print("  --keep-going, so independent branches continue after a")
                print("  failure -- which means a later rule finishing does NOT")
                print("  mean the run was healthy, and the log's tail can look")
                print("  clean while an earlier branch is dead.")
    print()

    print("-- resume point ----------------------------------------------------")
    if next_rule is None:
        print("  Every output exists. The DAG is complete; snakemake will")
        print("  report 'Nothing to be done'.")
    else:
        print(f"  next rule: {next_rule}")
        print(f"  completed rules: {done}   partial: {partial}")
        print()
        print("  Snakemake resumes from completed outputs, so re-submitting")
        print("  does NOT redo the rules marked done above. run_pipeline.sh")
        print("  pins --rerun-triggers mtime precisely so a touched file does")
        print("  not invalidate finished work.")
        if partial:
            print()
            print("  A PARTIAL rule is the dangerous case: some outputs written,")
            print("  some not. Snakemake will rerun it, but a half-written file")
            print("  from a kill mid-write can also be a TRUNCATED file that")
            print("  looks complete. If the rule above is partial, delete its")
            print("  outputs before resuming.")
    print()

    print("-- free space where outputs land -----------------------------------")
    try:
        du = shutil.disk_usage(workdir)
        print(f"  {workdir}")
        print(f"    total {human(du.total)}   used {human(du.used)}   "
              f"free {human(du.free)}")
        if du.free < 20 * 1024**3:
            print("    WARNING: under 20 GiB free. The region-to-gene and")
            print("    search-space TSVs are large (tens of millions of rows),")
            print("    and a write that runs out of space fails in a way that")
            print("    reads like a crash rather than a disk error.")
    except OSError as exc:
        print(f"  could not stat {workdir}: {exc}")
    tmp = cfg.get("params_general", {}).get("temp_dir") or os.environ.get("TMPDIR")
    if tmp and Path(tmp).parent.exists():
        try:
            du = shutil.disk_usage(Path(tmp).parent)
            print(f"  temp_dir parent: {Path(tmp).parent}")
            print(f"    free {human(du.free)}")
        except OSError:
            pass
    print()
    print("=" * 74)


if __name__ == "__main__":
    main()
