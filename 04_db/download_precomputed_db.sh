#!/usr/bin/env bash
#
# download_precomputed_db.sh -- fetch the precomputed hg38 SCREEN region-based
#                               cisTarget database + motif-to-TF annotation.
#
# ============================================================================
# WHAT YOU GET
# ============================================================================
# Base URL (verified reachable; server sends Content-Length, so sizes below are
# the real byte counts, not estimates):
#
#   https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/
#
#   FILE                                                  BYTES         SIZE
#   ----------------------------------------------------  ------------  --------
#   hg38_screen_v10_clust.regions_vs_motifs.rankings.feather  35,192,958,114  32.8 GiB (35.2 GB)
#   hg38_screen_v10_clust.regions_vs_motifs.scores.feather    13,882,267,682  12.9 GiB (13.9 GB)
#   ------------------------------------------------------------------------
#   both feathers                                             49,075,225,796  45.7 GiB (49.1 GB)
#
#   Motif-to-TF annotation (note: it is under /cistarget/motif2tf/, NOT
#   /motif2tf/ at the site root -- the root path 404s):
#   https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
#                                                                 98,718,421   94 MiB
#
#   TOTAL DOWNLOAD                                            49,173,944,217  45.8 GiB
#
# BOTH feathers are required by the default SCENIC+ Snakemake DAG. cisTarget
# reads the RANKINGS file; DEM reads the SCORES file -- and DEM is not optional:
# the Snakefile lists dem_results.hdf5 and dem_results.html as unconditional
# `output:` entries, so snakemake will run the DEM rule and fail on a missing
# scores feather. --skip-scores exists for the case where you have edited the
# workflow to drop DEM; on a stock checkout it produces a pipeline that dies
# partway through, after cisTarget has already spent its hours.
#
# DATABASE CONTENT (decoded from the Arrow schema/footer of the scores file):
#   1,837,304 regions x 5,876 motif rows. Regions are SCREEN cCREs on the 24
#   primary contigs (chr1-22, X, Y) -- no scaffolds, no alts. Region widths are
#   VARIABLE: min 150 bp, median 272 bp, max 350 bp (mean 267).
#
#   That width distribution matters: every SCREEN region is NARROWER than an
#   ArchR 501 bp fixed-width peak. A DB region fully inside your peak gives
#   overlap/peak_len as low as 150/501 = 0.30, which FAILS a 0.4 threshold on
#   the peak side -- but passes on the DB side (overlap/db_len = 1.0), and
#   pycistarget accepts either. Run peak_overlap_audit.py to get the real
#   number for your peak set before trusting this database.
#
# ============================================================================
# CHECKSUM VERIFICATION -- READ THIS
# ============================================================================
# There is NO aggregate manifest at
#     https://resources.aertslab.org/cistarget/databases/sha256sum.txt
# That URL returns HTTP 404 (verified). Do not build a workflow around it.
#
# What the server DOES provide is a per-file SHA1 sidecar next to each feather:
#     <file>.sha1sum.txt      containing "<sha1>  <filename>"
# Confirmed present, with these values:
#   1688a925f22d312769798258d990f13866bb4924  hg38_screen_v10_clust.regions_vs_motifs.rankings.feather
#   07b5e527d2ed082e081e439e68dffa77b5f6129c  hg38_screen_v10_clust.regions_vs_motifs.scores.feather
#
# This script therefore verifies with SHA1 from the sidecars (fetched live, so
# it stays correct if upstream republishes), and falls back to a strict
# byte-size check when a sidecar is unavailable. The expected SHA1s above are
# pinned as a cross-check: if the sidecar disagrees with the pin, you are
# warned loudly, because that means the upstream file changed.
#
# The motif2tf .tbl has no published checksum; it is size-checked only.
#
# NOTE on hashing cost: sha1sum over 46 GB takes 3-8 minutes on a decent
# filesystem and is I/O bound. Use --no-verify to skip (not recommended).
#
# ============================================================================
# USAGE
# ============================================================================
#   ./download_precomputed_db.sh --dest /path/to/resources/cistarget
#   ./download_precomputed_db.sh --dest DIR --skip-scores      # rankings only
#   ./download_precomputed_db.sh --dest DIR --dry-run
#
# Resume is on by default (wget -c / curl -C -). Re-running after an interrupted
# transfer continues where it stopped; a completed file is skipped unless
# --force. Interrupted-then-resumed files are exactly why the checksum step
# exists -- verify before you trust.
#
set -o errexit
set -o nounset
set -o pipefail

SCRIPT_NAME="$(basename "${0}")"

DB_BASE_URL="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based"
MOTIF2TF_URL="https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

RANKINGS_FILE="hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
SCORES_FILE="hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
MOTIF2TF_FILE="motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

# Exact Content-Length values observed from the server.
RANKINGS_BYTES=35192958114
SCORES_BYTES=13882267682
MOTIF2TF_BYTES=98718421

# Pinned SHA1s from the published .sha1sum.txt sidecars (cross-check only).
RANKINGS_SHA1_PINNED="1688a925f22d312769798258d990f13866bb4924"
SCORES_SHA1_PINNED="07b5e527d2ed082e081e439e68dffa77b5f6129c"

DEST=""
SKIP_SCORES=0
SKIP_MOTIF2TF=0
VERIFY=1
FORCE=0
DRY_RUN=0
DISK_MARGIN_GB=5

log()  { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${*}" >&2; }
warn() { printf '[%s] WARNING: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${*}" >&2; }

die() {
    printf '\n' >&2
    printf '========================================================================\n' >&2
    printf 'ERROR (%s): %s\n' "${SCRIPT_NAME}" "${1}" >&2
    shift
    while [ "${#}" -gt 0 ]; do printf '  -> %s\n' "${1}" >&2; shift; done
    printf '========================================================================\n' >&2
    exit 1
}

usage() { sed -n '2,/^set -o errexit/p' "${0}" | sed 's/^# \{0,1\}//; $d'; exit "${1:-0}"; }

while [ "${#}" -gt 0 ]; do
    case "${1}" in
        --dest)          DEST="${2:?--dest needs a value}"; shift 2 ;;
        --skip-scores)   SKIP_SCORES=1; shift ;;
        --skip-motif2tf) SKIP_MOTIF2TF=1; shift ;;
        --no-verify)     VERIFY=0; shift ;;
        --force)         FORCE=1; shift ;;
        --dry-run)       DRY_RUN=1; shift ;;
        --margin-gb)     DISK_MARGIN_GB="${2:?--margin-gb needs a value}"; shift 2 ;;
        -h|--help)       usage 0 ;;
        *) die "unknown option: ${1}" "run '${SCRIPT_NAME} --help' for usage" ;;
    esac
done

[ -n "${DEST}" ] || die "--dest is required" \
    "target directory for the database files, e.g. --dest resources/cistarget" \
    "you need ~46 GB free there (~33 GB with --skip-scores)"

# --------------------------------------------------------------------------- #
# Downloader
# --------------------------------------------------------------------------- #
DOWNLOADER=""
if command -v wget >/dev/null 2>&1; then
    DOWNLOADER="wget"
elif command -v curl >/dev/null 2>&1; then
    DOWNLOADER="curl"
else
    die "neither wget nor curl found on \$PATH" \
        "conda install -c conda-forge wget"
fi

SHA1_CMD=""
if command -v sha1sum >/dev/null 2>&1; then
    SHA1_CMD="sha1sum"
elif command -v shasum >/dev/null 2>&1; then
    SHA1_CMD="shasum -a 1"
fi
if [ "${VERIFY}" -eq 1 ] && [ -z "${SHA1_CMD}" ]; then
    warn "no sha1sum/shasum available -- falling back to size-only verification"
fi

# --------------------------------------------------------------------------- #
# Build the work list
# --------------------------------------------------------------------------- #
FILES=()
URLS=()
SIZES=()
SHA1S=()

FILES+=("${RANKINGS_FILE}"); URLS+=("${DB_BASE_URL}/${RANKINGS_FILE}")
SIZES+=("${RANKINGS_BYTES}"); SHA1S+=("${RANKINGS_SHA1_PINNED}")

if [ "${SKIP_SCORES}" -eq 0 ]; then
    FILES+=("${SCORES_FILE}"); URLS+=("${DB_BASE_URL}/${SCORES_FILE}")
    SIZES+=("${SCORES_BYTES}"); SHA1S+=("${SCORES_SHA1_PINNED}")
fi

if [ "${SKIP_MOTIF2TF}" -eq 0 ]; then
    FILES+=("${MOTIF2TF_FILE}"); URLS+=("${MOTIF2TF_URL}")
    SIZES+=("${MOTIF2TF_BYTES}"); SHA1S+=("")
fi

TOTAL_BYTES=0
for s in "${SIZES[@]}"; do TOTAL_BYTES=$(( TOTAL_BYTES + s )); done
TOTAL_GIB="$(awk -v b="${TOTAL_BYTES}" 'BEGIN{printf "%.1f", b/1024/1024/1024}')"

# --------------------------------------------------------------------------- #
# Plan
# --------------------------------------------------------------------------- #
cat >&2 <<PLAN
========================================================================
download_precomputed_db.sh -- PLAN
========================================================================
destination : ${DEST}
downloader  : ${DOWNLOADER} (resume enabled)
verify      : $([ "${VERIFY}" -eq 1 ] && echo "SHA1 via published .sha1sum.txt sidecars + size" || echo "DISABLED (--no-verify)")

files:
PLAN
for i in "${!FILES[@]}"; do
    printf '  %-58s %6.1f GiB\n' "${FILES[$i]}" \
        "$(awk -v b="${SIZES[$i]}" 'BEGIN{print b/1024/1024/1024}')" >&2
done
cat >&2 <<PLAN

total       : ${TOTAL_GIB} GiB (${TOTAL_BYTES} bytes)
disk needed : ${TOTAL_GIB} GiB + ${DISK_MARGIN_GB} GiB margin

note: BOTH feathers are needed by the stock Snakemake DAG. cisTarget reads
      the rankings file, DEM reads the scores file, and DEM's outputs are
      unconditional targets in the Snakefile -- so --skip-scores yields a
      pipeline that fails at the DEM rule unless you have edited the workflow.
========================================================================
PLAN

# --------------------------------------------------------------------------- #
# Disk check
# --------------------------------------------------------------------------- #
if [ "${DRY_RUN}" -eq 0 ]; then
    mkdir -p "${DEST}" || die "cannot create destination: ${DEST}"
    [ -w "${DEST}" ]   || die "destination not writable: ${DEST}"
fi
[ -d "${DEST}" ] || { log "--dry-run: ${DEST} does not exist yet"; }

if [ -d "${DEST}" ]; then
    AVAIL_KB="$(df -Pk "${DEST}" | awk 'NR==2{print $4}')"
    AVAIL_GIB="$(awk -v k="${AVAIL_KB}" 'BEGIN{printf "%.1f", k/1024/1024}')"

    # Credit bytes already on disk so a resume is not blocked by its own size.
    ALREADY=0
    for i in "${!FILES[@]}"; do
        f="${DEST}/${FILES[$i]}"
        [ -f "${f}" ] && ALREADY=$(( ALREADY + $(wc -c < "${f}" | tr -d ' ') ))
    done
    NEED_BYTES=$(( TOTAL_BYTES - ALREADY ))
    [ "${NEED_BYTES}" -lt 0 ] && NEED_BYTES=0
    NEED_GIB="$(awk -v b="${NEED_BYTES}" -v m="${DISK_MARGIN_GB}" \
        'BEGIN{printf "%.1f", b/1024/1024/1024 + m}')"

    log "free space in ${DEST}: ${AVAIL_GIB} GiB"
    log "still to download    : $(awk -v b="${NEED_BYTES}" 'BEGIN{printf "%.1f", b/1024/1024/1024}') GiB (+${DISK_MARGIN_GB} GiB margin = ${NEED_GIB} GiB)"

    if awk -v a="${AVAIL_GIB}" -v n="${NEED_GIB}" 'BEGIN{exit !(a < n)}'; then
        die "not enough free disk in ${DEST}: ${AVAIL_GIB} GiB free, ${NEED_GIB} GiB needed" \
            "options:" \
            "  * --skip-scores        omit the 12.9 GiB scores feather (SCENIC+ needs only rankings)" \
            "  * --dest <bigger_fs>   point at scratch/project storage" \
            "  * --margin-gb 1        reduce the safety margin (default ${DISK_MARGIN_GB} GiB)" \
            "NOTE: on many clusters \$HOME has a quota far below 46 GiB -- use /scratch or a project dir."
    fi
fi

if [ "${DRY_RUN}" -eq 1 ]; then
    log "--dry-run: nothing downloaded."
    exit 0
fi

# --------------------------------------------------------------------------- #
# Fetch helpers
# --------------------------------------------------------------------------- #
fetch() {  # fetch <url> <outfile>
    local url="${1}" out="${2}"
    if [ "${DOWNLOADER}" = "wget" ]; then
        # -c resume, --tries, --progress=dot:giga keeps SLURM logs small
        wget -c --tries=5 --timeout=60 --waitretry=15 \
             --progress=dot:giga -O "${out}" "${url}"
    else
        # curl -C - resumes; --retry handles transient failures.
        # -# gives a single-line progress bar instead of the multi-line meter,
        # which keeps SLURM .out files from filling with carriage returns.
        curl -fL -C - --retry 5 --retry-delay 15 --connect-timeout 60 \
             -# -o "${out}" "${url}"
    fi
}

remote_size() {  # remote_size <url> -> bytes on stdout ("" if unknown)
    local url="${1}"
    if command -v curl >/dev/null 2>&1; then
        curl -sSIL "${url}" 2>/dev/null \
            | awk 'BEGIN{IGNORECASE=1} /^content-length:/{v=$2} END{gsub(/\r/,"",v); print v}'
    fi
}

fetch_sha1_sidecar() {  # fetch_sha1_sidecar <file_url> -> sha1 on stdout
    local url="${1}.sha1sum.txt"
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL "${url}" 2>/dev/null | awk 'NF{print $1; exit}'
    elif command -v wget >/dev/null 2>&1; then
        wget -qO- "${url}" 2>/dev/null | awk 'NF{print $1; exit}'
    fi
}

# --------------------------------------------------------------------------- #
# Download loop
# --------------------------------------------------------------------------- #
FAILED=0
SUMMARY=()

for i in "${!FILES[@]}"; do
    name="${FILES[$i]}"
    url="${URLS[$i]}"
    want_bytes="${SIZES[$i]}"
    pinned_sha1="${SHA1S[$i]}"
    out="${DEST}/${name}"

    printf '\n' >&2
    log "=== ${name} ==="

    # Skip if already complete and correctly sized.
    if [ -f "${out}" ] && [ "${FORCE}" -eq 0 ]; then
        have="$(wc -c < "${out}" | tr -d ' ')"
        if [ "${have}" -eq "${want_bytes}" ]; then
            log "already complete (${have} bytes) -- skipping download"
        else
            log "partial file present: ${have} / ${want_bytes} bytes ($(awk -v h="${have}" -v w="${want_bytes}" 'BEGIN{printf "%.1f", 100*h/w}')%) -- resuming"
        fi
    fi

    have=0
    [ -f "${out}" ] && have="$(wc -c < "${out}" | tr -d ' ')"

    if [ "${FORCE}" -eq 1 ]; then
        log "--force: removing existing file"
        rm -f "${out}"
        have=0
    fi

    if [ "${have}" -ne "${want_bytes}" ]; then
        # Sanity-check the published size against what we have pinned, so a
        # silently-republished upstream file is visible rather than a mystery
        # checksum failure.
        rsize="$(remote_size "${url}" || true)"
        if [ -n "${rsize}" ] && [ "${rsize}" != "${want_bytes}" ]; then
            warn "server reports ${rsize} bytes but this script expects ${want_bytes}."
            warn "upstream may have republished ${name}. Using the SERVER size."
            want_bytes="${rsize}"
        fi

        log "downloading -> ${out}"
        if ! fetch "${url}" "${out}"; then
            warn "download FAILED for ${name}"
            warn "re-run the script to resume from $( [ -f "${out}" ] && wc -c < "${out}" | tr -d ' ' || echo 0) bytes"
            SUMMARY+=("FAIL(download)  ${name}")
            FAILED=1
            continue
        fi
    fi

    # -- size check ------------------------------------------------------- #
    have="$(wc -c < "${out}" | tr -d ' ')"
    if [ "${have}" -ne "${want_bytes}" ]; then
        warn "SIZE MISMATCH for ${name}: got ${have}, expected ${want_bytes}"
        SUMMARY+=("FAIL(size)      ${name}")
        FAILED=1
        continue
    fi
    log "size OK: ${have} bytes"

    # -- checksum --------------------------------------------------------- #
    if [ "${VERIFY}" -eq 0 ] || [ -z "${SHA1_CMD}" ]; then
        SUMMARY+=("OK(size only)   ${name}")
        continue
    fi

    # There is no aggregate sha256sum.txt (that URL 404s); the real artifact is
    # a per-file .sha1sum.txt sidecar. Fetch it live so we track upstream.
    expected="$(fetch_sha1_sidecar "${url}" || true)"

    if [ -z "${expected}" ]; then
        if [ -n "${pinned_sha1}" ]; then
            warn "could not fetch ${name}.sha1sum.txt -- using the SHA1 pinned in this script"
            expected="${pinned_sha1}"
        else
            log "no published checksum for ${name} -- size verification only"
            SUMMARY+=("OK(size only)   ${name}")
            continue
        fi
    elif [ -n "${pinned_sha1}" ] && [ "${expected}" != "${pinned_sha1}" ]; then
        warn "published SHA1 for ${name} differs from the value pinned in this script!"
        warn "  published: ${expected}"
        warn "  pinned   : ${pinned_sha1}"
        warn "upstream has republished this file. Trusting the PUBLISHED value."
    fi

    log "computing SHA1 (46 GiB total takes several minutes; I/O bound)..."
    actual="$(${SHA1_CMD} "${out}" | awk '{print $1}')"

    if [ "${actual}" = "${expected}" ]; then
        log "SHA1 OK: ${actual}"
        SUMMARY+=("OK(sha1)        ${name}")
    else
        warn "SHA1 MISMATCH for ${name}"
        warn "  expected: ${expected}"
        warn "  actual  : ${actual}"
        warn "the file is corrupt (interrupted resume, bad disk, or truncated proxy)"
        warn "delete it and re-download: rm '${out}' && ${SCRIPT_NAME} --dest '${DEST}'"
        SUMMARY+=("FAIL(sha1)      ${name}")
        FAILED=1
    fi
done

# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
printf '\n' >&2
printf '========================================================================\n' >&2
printf 'DOWNLOAD SUMMARY\n' >&2
printf '========================================================================\n' >&2
for s in "${SUMMARY[@]}"; do printf '  %s\n' "${s}" >&2; done
printf '========================================================================\n' >&2

if [ "${FAILED}" -ne 0 ]; then
    die "one or more files failed" \
        "re-run the same command -- completed files are skipped and partial" \
        "files resume from where they stopped"
fi

RANKINGS_PATH="${DEST}/${RANKINGS_FILE}"
log "all files verified."
log ""
log "NEXT STEP -- do this BEFORE running SCENIC+:"
log "  python 04_db/peak_overlap_audit.py \\"
log "      --peaks <your_consensus_peaks.bed> \\"
log "      --db '${RANKINGS_PATH}' \\"
log "      --out-prefix 04_db/audit_screen"
log ""
log "SCREEN regions are cCREs of 150-350 bp, narrower than 501 bp ArchR peaks,"
log "and cover only chr1-22,X,Y. Peaks that fail the 0.4 overlap rule are"
log "dropped from motif enrichment with NO warning. The audit tells you how"
log "many, and whether a custom build is worth the CPU-hours."
