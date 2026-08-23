#!/usr/bin/env bash
# =============================================================================
# Did the pairing job succeed? -- read-only postmortem
# =============================================================================
#     bash slurm/check_pairing.sh            # most recent glue_pairing job
#     bash slurm/check_pairing.sh 12345678   # a specific job id
#
# Answers, in order:
#   1. how SLURM says it ended, and what it actually used (MaxRSS vs the 80G
#      request -- the number that tells us whether the estimate was right)
#   2. the tail of the job log
#   3. whether ACC_GEX.h5mu exists and is complete
#   4. the pairing diagnostics: which cell types paired badly
#
# "FAILED" with an exit code of 1 in seconds is usually a path or env problem;
# "OUT_OF_MEMORY", or CANCELLED with MaxRSS near the request, means the memory
# estimate was low.
# =============================================================================
set -uo pipefail
JOB="${1:-}"

echo "=============================================================="
echo "pairing job postmortem"
echo "=============================================================="

# --- 1. how it ended --------------------------------------------------------
if command -v sacct >/dev/null 2>&1; then
    if [[ -z "$JOB" ]]; then
        JOB=$(sacct -u "$USER" --name=glue_pairing -n -X \
                    --starttime now-7days --format=JobID%-20 2>/dev/null \
              | awk 'NF{j=$1} END{print j}')
    fi
    if [[ -z "$JOB" ]]; then
        echo "[1] no job named 'glue_pairing' in the last 7 days."
        echo "    Either it was never submitted, or it ran under another name."
        echo "    Recent jobs:"
        sacct -u "$USER" --starttime now-2days \
              --format=JobID%-16,JobName%-18,State%-14,Elapsed,MaxRSS 2>/dev/null | head -12
    else
        echo "[1] job $JOB"
        sacct -j "$JOB" --format=JobID%-18,JobName%-16,State%-14,ExitCode%-8,Elapsed%-10,Timelimit%-10,ReqMem%-9,MaxRSS%-10,NodeList%-14 2>/dev/null
        echo
        STATE=$(sacct -j "$JOB" -n -X --format=State 2>/dev/null | head -1 | tr -d ' ')
        echo "    verdict: ${STATE:-unknown}"
        case "$STATE" in
            COMPLETED) echo "    -> the job exited 0. Check [3] that the output is real." ;;
            OUT_OF_MEMORY|OOM*) echo "    -> memory. Compare MaxRSS above with ReqMem; raise" ;
                                echo "       --mem in slurm/pairing.sbatch, or lower --oversample." ;;
            TIMEOUT)   echo "    -> hit the walltime. Raise --time (batch caps at 1-12:00:00)." ;;
            FAILED)    echo "    -> nonzero exit. The log in [2] has the reason; a failure in" ;
                       echo "       seconds is almost always a path or environment problem." ;;
            CANCELLED*) echo "    -> cancelled. If MaxRSS is near ReqMem this was likely the" ;
                        echo "       OOM killer rather than a manual scancel." ;;
        esac
        command -v seff >/dev/null 2>&1 && { echo; echo "    seff:"; seff "$JOB" 2>/dev/null | sed 's/^/      /'; }
    fi
else
    echo "[1] sacct unavailable"
fi

# --- 2. the log -------------------------------------------------------------
echo
echo "[2] job log"
# Read BOTH streams. The .sbatch sends stdout and stderr to separate files, so
# a python traceback lands in .err while .out just stops -- looking only at
# .out reports "no errors found" on a job that died with a traceback.
LOGS=()
for cand in "jobs/pairing_${JOB}.out" "jobs/pairing_${JOB}.err"; do
    [[ -f "$cand" ]] && LOGS+=("$cand")
done
if (( ${#LOGS[@]} == 0 )); then
    while IFS= read -r f; do [[ -n "$f" ]] && LOGS+=("$f"); done < <(
        ls -t jobs/pairing_*.out jobs/pairing_*.err 2>/dev/null | head -2)
fi
if (( ${#LOGS[@]} )); then
    for LOG in "${LOGS[@]}"; do
        echo "    $LOG  ($(wc -l < "$LOG") lines)"
    done
    STDERR_LOG=""
    for LOG in "${LOGS[@]}"; do [[ "$LOG" == *.err ]] && STDERR_LOG="$LOG"; done
    # Show the stderr tail FIRST when it has content -- that is where the cause is.
    if [[ -n "$STDERR_LOG" && -s "$STDERR_LOG" ]]; then
        echo
        echo "    --- STDERR tail (the failure cause lives here) ---"
        tail -30 "$STDERR_LOG" | sed 's/^/      /'
    elif [[ -n "$STDERR_LOG" ]]; then
        echo "    (stderr file is empty)"
    else
        echo "    NOTE: no .err file found. If the job failed with no error in"
        echo "          .out, the traceback went to a stderr file that is"
        echo "          missing -- check the --error path in pairing.sbatch."
    fi
    echo
    echo "    --- stdout tail ---"
    tail -20 "${LOGS[0]}" | sed 's/^/      /'
    echo
    echo "    --- error-shaped lines across BOTH streams ---"
    grep -niE "error|traceback|killed|oom|no such file|refus|cannot|exceeded|attributeerror|keyerror|valueerror" \
         "${LOGS[@]}" 2>/dev/null | head -15 | sed 's/^/      /' || echo "      none"
    grep -hiE "maximum resident set size" "${LOGS[@]}" 2>/dev/null \
        | sed 's/^/      /' || true
else
    echo "    no log found under jobs/. If the directory did not exist at submit"
    echo "    time, SLURM discarded the output silently -- mkdir -p jobs and"
    echo "    resubmit."
fi

# --- 3. the output ----------------------------------------------------------
echo
echo "[3] output file"
OUT="${PAIRED:-ACC_GEX.h5mu}"
if [[ -f "$OUT" ]]; then
    echo "    $OUT  $(du -h "$OUT" | cut -f1)  modified $(date -r "$OUT" 2>/dev/null)"
    echo "    -> verifying it is a complete, readable MuData:"
    python 03_pipeline/validate_h5mu.py "$OUT" 2>&1 | tail -12 | sed 's/^/      /'
else
    echo "    $OUT does NOT exist -- the job did not get to the write step."
fi

# --- 4. diagnostics ---------------------------------------------------------
echo
echo "[4] pairing diagnostics"
if [[ -f pairing_diagnostics.csv ]]; then
    python - <<'PY'
import csv
rows = list(csv.DictReader(open("pairing_diagnostics.csv")))
if not rows:
    print("      diagnostics file is empty")
else:
    def f(r, k):
        try: return float(r[k])
        except (KeyError, ValueError, TypeError): return float("nan")
    rows.sort(key=lambda r: -f(r, "median_crossmodal_gap"))
    print(f"      {len(rows)} groups. Worst cross-modal gaps -- RNA and ATAC do")
    print("      not occupy the same latent region here, so these links are the")
    print("      least trustworthy:")
    print(f"      {'group':<26}{'gap':>8}{'metacells':>11}{'indep':>7}")
    for r in rows[:8]:
        print(f"      {r.get('group','?'):<26}{f(r,'median_crossmodal_gap'):>8.4f}"
              f"{r.get('n_metacells','?'):>11}{r.get('independent_metacell_equiv','?'):>7}")
    tot = sum(int(r.get("n_metacells", 0) or 0) for r in rows)
    print(f"\n      total metacells: {tot:,}   (expected 25,323)")
PY
else
    echo "    pairing_diagnostics.csv not found -- the run did not reach the end."
fi

echo
echo "=============================================================="
