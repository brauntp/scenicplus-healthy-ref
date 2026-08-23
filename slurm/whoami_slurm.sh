#!/usr/bin/env bash
# =============================================================================
# What account and partition should I submit to?
# =============================================================================
#     bash slurm/whoami_slurm.sh
#
# Read-only. Asks SLURM four questions:
#   1. what have I successfully submitted before?   (most reliable answer)
#   2. what accounts am I associated with?
#   3. what partitions exist, and which is the default?
#   4. do I even need to specify them?
#
# If a default account and partition exist, `sbatch slurm/pairing.sbatch` works
# with no flags at all.
# =============================================================================
set -uo pipefail
command -v sinfo >/dev/null 2>&1 || {
    echo "ERROR: no SLURM commands on PATH. Are you on a login node?" >&2
    exit 1
}

echo "=============================================================="
echo "SLURM submission targets for ${USER:-$(whoami)}"
echo "=============================================================="

# --- 1. history: what already worked ----------------------------------------
echo
echo "[1] your recent jobs -- the account/partition here are known-good"
if command -v sacct >/dev/null 2>&1; then
    sacct -u "$USER" --starttime now-90days \
          --format=JobID%-14,Account%-18,Partition%-16,State%-12,Elapsed \
          2>/dev/null | head -25
    echo
    echo "    most-used account/partition pairs in the last 90 days:"
    sacct -u "$USER" --starttime now-90days -n -X \
          --format=Account%-30,Partition%-30 2>/dev/null \
        | awk 'NF==2{c[$1" "$2]++} END{for(k in c) printf "      %-45s %d jobs\n", k, c[k]}' \
        | sort -k2 -rn | head -8
else
    echo "    sacct not available (no accounting database configured)"
fi

# --- 2. entitlements --------------------------------------------------------
echo
echo "[2] accounts you are associated with"
if command -v sacctmgr >/dev/null 2>&1; then
    sacctmgr -n show associations user="$USER" \
             format=Account%-24,Partition%-20,QOS%-30 2>/dev/null | head -20 \
        || echo "    (sacctmgr returned nothing)"
    echo
    echo "    default account:"
    sacctmgr -n show user "$USER" format=DefaultAccount%-30 2>/dev/null \
        | sed 's/^/      /' || true
else
    echo "    sacctmgr not available"
fi
if command -v sshare >/dev/null 2>&1; then
    echo
    echo "    fairshare view (also lists your accounts):"
    sshare -U -u "$USER" 2>/dev/null | head -10 | sed 's/^/      /'
fi

# --- 3. partitions ----------------------------------------------------------
echo
echo "[3] partitions (* marks the default)"
sinfo -o "%20P %10a %12l %8D %10c %10m %N" 2>/dev/null | head -25

echo
echo "    partitions with enough memory for the pairing job (>=80 GB/node):"
sinfo -h -o "%P|%m|%c|%l" 2>/dev/null \
    | awk -F'|' '{gsub(/\+/,"",$2); if ($2+0 >= 80000) printf "      %-20s %8.0f GB  %3s cores  maxtime %s\n", $1, $2/1024, $3, $4}' \
    | sort -u | head -15

# --- 4. is anything even required? ------------------------------------------
echo
echo "[4] do you need to pass --account / --partition?"
DEF_PART=$(sinfo -h -o "%P" 2>/dev/null | grep '\*' | tr -d '*' | head -1)
DEF_ACCT=$(sacctmgr -n show user "$USER" format=DefaultAccount 2>/dev/null | awk 'NF{print $1; exit}')
[[ -n "$DEF_PART" ]] && echo "    default partition: $DEF_PART" \
                     || echo "    no default partition -- you must pass --partition"
[[ -n "$DEF_ACCT" && "$DEF_ACCT" != "(null)" ]] \
    && echo "    default account  : $DEF_ACCT" \
    || echo "    no default account -- you may need --account"

# The most-used partition from history, which beats the site default: a default
# is often an INTERACTIVE partition, and a multi-hour batch job does not belong
# in a pool shared with people waiting at a prompt.
HIST_PART=""
HIST_ACCT=""
if command -v sacct >/dev/null 2>&1; then
    read -r HIST_ACCT HIST_PART < <(
        sacct -u "$USER" --starttime now-90days -n -X \
              --format=Account%-30,Partition%-30 2>/dev/null \
        | awk 'NF==2{c[$1" "$2]++} END{m=0; for(k in c) if(c[k]>m){m=c[k]; b=k} print b}')
fi

echo
echo "=============================================================="
if [[ -n "$HIST_PART" && -n "$HIST_ACCT" ]]; then
    echo "RECOMMENDED -- what you actually use:"
    echo "    sbatch --account=$HIST_ACCT --partition=$HIST_PART slurm/pairing.sbatch"
    if [[ -n "$DEF_PART" && "$DEF_PART" != "$HIST_PART" ]]; then
        echo
        echo "Note: the site DEFAULT partition is '$DEF_PART', which is not the"
        echo "one you normally use. Passing --partition explicitly avoids"
        echo "landing somewhere unintended -- especially if the default is an"
        echo "interactive pool, which is the wrong home for a multi-hour job."
    fi
elif [[ -n "$DEF_PART" && -n "$DEF_ACCT" && "$DEF_ACCT" != "(null)" ]]; then
    echo "No job history to learn from. Defaults are:"
    echo "    sbatch --account=$DEF_ACCT --partition=$DEF_PART slurm/pairing.sbatch"
    echo "Check in [3] that '$DEF_PART' is a batch partition, not interactive."
else
    echo "Submit with whatever [1] shows worked before:"
    echo "    sbatch --account=<acct> --partition=<part> slurm/pairing.sbatch"
fi
echo
echo "The pairing job needs --mem=80G and about 4 hours; check the maxtime"
echo "column in [3] so you do not pick a partition that caps below that."
echo "=============================================================="
