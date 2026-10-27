#!/usr/bin/env bash
# Re-run the jobs from the repository root and compare their logs with git diff.
# Usage: ./verify_all.sh [-j JOBS] [LOG_PATH_FILTER]
# A full run costs roughly 43 CPU hours (the archive-only logs are not rerun); start with a subset.
# Magma jobs run first (JOBS at a time), then the other jobs one at a time, in the order of
# jobs.txt, since they read the logs and data files written by earlier jobs.
# Jobs whose command starts with "archived:" are not run (see jobs.txt and the README).
# A log whose command is "written by the job of LOG" is checked after the job of LOG has run;
# selecting such a log selects that job.
# Exit 0 means all selected jobs passed their completion and result checks.
# Failed jobs give a nonzero status (their count, capped at 125).
set -uo pipefail

usage_error() {
    echo "ERROR: $*" >&2
    echo "Usage: ./verify_all.sh [-j JOBS] [LOG_PATH_FILTER]" >&2
    exit 2
}

JOBS=1
while getopts "j:" opt; do
    case $opt in j) JOBS=$OPTARG;; *) usage_error "invalid option";; esac
done
shift $((OPTIND-1))
[[ $JOBS =~ ^[1-9][0-9]*$ ]] || usage_error "-j requires a positive integer"
[ "$#" -le 1 ] || usage_error "expected at most one log-path filter"
FILTER="${1:-}"
[ -r jobs.txt ] || usage_error "jobs.txt is missing or unreadable; run from the repository root"

# Validate the whole manifest before overwriting any logs. The fourth column is optional.
JOB_RE=$'^([^\t]+)\t([^\t]+)\t([^\t]+)(\t([^\t]+))?$'
DEP_PREFIX='written by the job of '
declare -A SEEN=() KIND=() PRIMARY=() WANT=()
ALL=(); LINE=0
while IFS= read -r job || [ -n "$job" ]; do
    LINE=$((LINE+1))
    [[ $job =~ ^[[:space:]]*(#|$) ]] && continue
    [[ $job =~ $JOB_RE ]] || usage_error "jobs.txt:$LINE: expected three or four tab-separated fields"
    log=${BASH_REMATCH[1]}; cmd=${BASH_REMATCH[2]}
    expected=${BASH_REMATCH[5]:-}
    [ -z "${SEEN[$log]+present}" ] || usage_error "jobs.txt:$LINE: duplicate log path $log"
    SEEN[$log]=$LINE
    if [ -n "$expected" ]; then
        rc=0; grep -E -- "$expected" < /dev/null > /dev/null || rc=$?
        [ "$rc" -le 1 ] || usage_error "jobs.txt:$LINE: invalid result expression"
    fi
    ALL+=("$job")
    if [[ $cmd == archived:* ]]; then
        KIND[$log]=archived
    elif [[ $cmd == "$DEP_PREFIX"* ]]; then
        KIND[$log]=dependent; PRIMARY[$log]=${cmd#"$DEP_PREFIX"}
    elif [[ $cmd == magma\ * ]]; then
        KIND[$log]=magma
    else
        KIND[$log]=other
    fi
done < jobs.txt
[ "${#ALL[@]}" -gt 0 ] || usage_error "jobs.txt contains no jobs"
for log in "${!PRIMARY[@]}"; do
    p=${PRIMARY[$log]}
    [ -n "${KIND[$p]+present}" ] || usage_error "jobs.txt:${SEEN[$log]}: $log names an unknown job $p"
    [[ ${KIND[$p]} == magma || ${KIND[$p]} == other ]] || usage_error "jobs.txt:${SEEN[$log]}: $log names $p, which is not a runnable job"
done

# Selection: the logs matching the filter, the jobs that write them, and the logs those jobs write.
for job in "${ALL[@]}"; do
    log=${job%%$'\t'*}
    if [[ -z $FILTER || $log == *"$FILTER"* ]]; then
        WANT[$log]=1
        [ "${KIND[$log]}" = dependent ] && WANT[${PRIMARY[$log]}]=1
    fi
done
for log in "${!PRIMARY[@]}"; do
    [ -n "${WANT[${PRIMARY[$log]}]+present}" ] && WANT[$log]=1
done
MAGMA_JOBS=(); OTHER_JOBS=(); DEPENDENT=(); ARCHIVED=(); NSEL=0; NEED_MAGMA=0
for job in "${ALL[@]}"; do
    log=${job%%$'\t'*}
    [ -n "${WANT[$log]+present}" ] || continue
    case ${KIND[$log]} in
        archived) ARCHIVED+=("$log");;
        dependent) DEPENDENT+=("$job"); NSEL=$((NSEL+1));;
        magma) MAGMA_JOBS+=("$job"); NEED_MAGMA=1; NSEL=$((NSEL+1));;
        *) OTHER_JOBS+=("$job"); NSEL=$((NSEL+1));;
    esac
done
[ $((NSEL+${#ARCHIVED[@]})) -gt 0 ] || usage_error "no log path matches '$FILTER'"
[ "$NSEL" -gt 0 ] || usage_error "only archive-only logs match '$FILTER'; nothing to run"

PATCH=mdmagma-magma2.30.patch
if [ "$NEED_MAGMA" -eq 1 ]; then
    command -v magma > /dev/null || usage_error "magma is not on PATH"
    for d in mdmagma/v2/mdmagma.spec mdmagma/Magma/magma.spec gl2/magma.spec common/fast_hecke.m "$PATCH"; do
        [ -r "$d" ] && continue
        echo "Missing $d. Initialise the pinned packages from the repository root:" >&2
        echo "    git submodule update --init" >&2
        echo "    git -C mdmagma submodule update --init Magma" >&2
        exit 1
    done
fi

STATUSDIR=$(mktemp -d) || exit 1
trap 'rm -rf -- "$STATUSDIR"' EXIT
export STATUSDIR

# The closing line of a log: the last non-blank line, ignoring carriage returns, a final
# "Script done on ..." line written by script(1), and final lines holding only a CPU time.
closing_line() {
    tr -d '\r' < "$1" | grep -av '^[[:space:]]*$' | awk '
        { line[NR] = $0 }
        END {
            n = NR
            while (n > 0 && (line[n] ~ /^Script done on / || line[n] ~ /^[[:space:]]*[0-9]+\.[0-9]+[[:space:]]*$/)) n--
            if (n > 0) print line[n]
        }'
}
fail_job() {
    local record
    echo "!!! $2"
    record=$(mktemp "$STATUSDIR/fail.XXXXXX") || return 1
    printf '%s\n' "$1" > "$record"
    return 1
}
pass_job() {
    local record
    record=$(mktemp "$STATUSDIR/ok.XXXXXX") || return 1
    printf '%s\n' "$1" > "$record"
}
# The checks on a finished log: no Magma error, the closing marker, the expected result.
check_log() {
    local log=$1 marker=$2 expected=$3 last
    [ -f "$log" ] || { fail_job "$log" "MISSING: $log was not written"; return 1; }
    # Magma reports most errors on standard output, carries on and exits with status 0.
    if grep -aEq '(Runtime|User|System|Internal) error|Assertion failed' "$log"; then
        fail_job "$log" "MAGMA ERROR: $log reports an error"
        return 1
    fi
    last=$(closing_line "$log")
    if [[ $last != "$marker" && $last != "$marker "* ]]; then
        fail_job "$log" "INCOMPLETE: $log should end in '$marker', last line is '${last:-<empty>}'"
        return 1
    fi
    if [ -n "$expected" ] && ! grep -aEq -- "$expected" "$log"; then
        fail_job "$log" "UNEXPECTED RESULT: $log does not match '$expected'"
        return 1
    fi
    pass_job "$log"
}
run_one() {
    local log cmd marker expected rc
    IFS=$'\t' read -r log cmd marker expected <<< "$1"
    mkdir -p "$(dirname "$log")" || { fail_job "$log" "cannot create log directory: $log"; return 1; }
    echo "==> $log : $cmd"
    # Several scripts end without quit; Magma then stops at the end of its standard input.
    rc=0; (eval "$cmd") < /dev/null > "$log" 2>&1 || rc=$?
    [ "$rc" -eq 0 ] || { fail_job "$log" "FAILED (exit $rc): $log"; return 1; }
    check_log "$log" "$marker" "$expected"
}
check_dependent() {
    local log cmd marker expected primary
    IFS=$'\t' read -r log cmd marker expected <<< "$1"
    primary=${cmd#"$DEP_PREFIX"}
    if ! grep -qxF -- "$primary" /dev/null "$STATUSDIR"/ok.* 2> /dev/null; then
        fail_job "$log" "NOT CHECKED: $log is written by the job of $primary, which failed"
        return 1
    fi
    echo "==> $log : $cmd"
    check_log "$log" "$marker" "$expected"
}
export DEP_PREFIX
export -f closing_line fail_job pass_job check_log run_one check_dependent

# Magma V2.30 rejects two statements of the pinned mdmagma; apply the patch if needed.
patch_mdmagma() {
    if git -C mdmagma apply --reverse --check "../$PATCH" > /dev/null 2>&1; then
        echo "==> $PATCH is applied to mdmagma"
    elif git -C mdmagma apply "../$PATCH"; then
        echo "==> applied $PATCH to mdmagma"
    else
        echo "!!! could not apply $PATCH to mdmagma; see README, Getting the code" >&2
        return 1
    fi
}

# Compile every package once, in separate processes, before several jobs attach them at the
# same time: concurrent first-time attachment can cause Magma signature-cache errors.
# Packages that define the same intrinsics are attached in separate processes.
init_packages() {
    local set rc last i=0
    echo "==> initialising the Magma packages"
    for set in \
        "mdmagma/Magma/magma.spec" \
        "mdmagma/v2/mdmagma.spec C54/hasse_bsgs_count.m" \
        "mdmagma/v2/mdmagma.spec common/class_rational_reduction.m C54/C54_coordinate_place.m" \
        "mdmagma/v2/mdmagma.spec common/class_rational_reduction.m common/fast_diamond_hash_m2.m" \
        "mdmagma/v2/mdmagma.spec common/class_rational_reduction.m common/fast_full_aut_hash_m3_15.m" \
        "mdmagma/v2/mdmagma.spec common/class_rational_reduction.m common/fast_hecke.m common/fast_hecke_11.m common/fast_hecke_cyclic_commonzero.m C54/C54_coordinate_place.m" \
        "mdmagma/v2/mdmagma.spec common/fast_diamond_hash_m2.m C3x18/fast_full_aut_hash_m3_18_root.m" \
        "mdmagma/v2/mdmagma.spec common/fast_diamond_hash_m2.m common/fast_diamond_hash_cyclic_bc.m" \
        "mdmagma/v2/mdmagma.spec common/fast_diamond_hash_m2.m common/fast_full_aut_hash_m4_16.m" \
        "mdmagma/v2/mdmagma.spec common/fast_full_aut_hash_m3_15.m C3x15/fast_full_aut_hash_m3_15_native.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m C3x18/fast_full_aut_hash_m3_18_root.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m C6x12/model_6_12_root.m common/fast_diamond_hash_m2.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/class_rational_reduction.m common/fast_full_aut_hash_m3_15.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_diamond_hash_m2.m C2x28/fast_full_aut_hash_m2_28.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_diamond_hash_m2.m C2x32/fast_full_aut_hash_m2_32.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_diamond_hash_m2.m C4x24/fast_full_aut_hash_m4_24_root.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_diamond_hash_m2.m common/fast_diamond_hash_cyclic_bc.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_diamond_hash_m2.m common/fast_full_aut_hash_m4_16.m common/fast_full_aut_hash_m4_16_relative.m C4x20/fast_full_aut_hash_m4_20.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_diamond_hash_m2.m common/fast_full_aut_hash_m4_16.m common/fast_full_aut_hash_m4_16_relative.m common/class_rational_reduction.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_full_aut_hash_m3_15.m C3x15/fast_full_aut_hash_m3_15_native.m C3x15/fast_full_aut_hash_m3_15_keytable.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_full_aut_hash_m3_15.m C3x21/fast_full_aut_m3_21.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_full_aut_hash_m3_15.m C3x24/fast_full_aut_hash_m3_24.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_hecke_11.m C6x12/model_6_12_root.m common/class_rational_reduction.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_hecke_11.m common/fast_hecke_cyclic_commonzero.m common/class_rational_reduction.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_hecke_11.m common/fast_hecke_cyclic_commonzero.m common/fast_diamond_hash_m2.m common/fast_diamond_hash_cyclic_bc.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke.m common/fast_hecke_cyclic_commonzero.m common/class_rational_reduction.m" \
        "mdmagma/v2/mdmagma.spec common/fast_hecke_11.m common/fast_hecke.m" \
        "gl2/magma.spec"; do
        i=$((i+1))
        {
            printf 'SetColumns(0);\n'
            for f in $set; do
                case $f in *.spec) printf 'AttachSpec("%s");\n' "$f";; *) printf 'Attach("%s");\n' "$f";; esac
            done
            printf 'print "PACKAGE_READY";\nquit;\n'
        } > "$STATUSDIR/init.m" || return 1
        echo "    $set"
        rc=0; magma -n -b "$STATUSDIR/init.m" < /dev/null > "$STATUSDIR/init.log" 2>&1 || rc=$?
        last=$(closing_line "$STATUSDIR/init.log")
        if [ "$rc" -ne 0 ] || [ "$last" != PACKAGE_READY ]; then
            echo "!!! could not attach $set (exit $rc):" >&2
            tail -20 "$STATUSDIR/init.log" >&2
            return 1
        fi
    done
}

echo "$NSEL of $((${#ALL[@]}-${#ARCHIVED[@]})) jobs selected (${#MAGMA_JOBS[@]} Magma, ${#OTHER_JOBS[@]} other, ${#DEPENDENT[@]} logs written by these); $JOBS Magma jobs at a time"
for log in "${ARCHIVED[@]}"; do echo "--- skipping archive-only $log"; done
if [ "$NEED_MAGMA" -eq 1 ]; then
    patch_mdmagma || exit 1
    init_packages || exit 1
fi

RUN_STATUS=0
if [ "${#MAGMA_JOBS[@]}" -gt 0 ]; then
    if [ "$JOBS" = 1 ]; then
        for job in "${MAGMA_JOBS[@]}"; do run_one "$job" || RUN_STATUS=$?; done
    else
        printf '%s\n' "${MAGMA_JOBS[@]}" | xargs -d '\n' -P "$JOBS" -I{} bash -uo pipefail -c 'run_one "$@"' _ {} || RUN_STATUS=$?
    fi
fi
for job in "${OTHER_JOBS[@]}"; do run_one "$job" || true; done
for job in "${DEPENDENT[@]}"; do check_dependent "$job" || true; done

shopt -s nullglob
FAILED=("$STATUSDIR"/fail.*); PASSED=("$STATUSDIR"/ok.*)
NFAIL=${#FAILED[@]}
MISSING=$((NSEL-${#FAILED[@]}-${#PASSED[@]}))
if [ "$MISSING" -ne 0 ]; then
    echo "!!! execution did not record a result for every selected job" >&2
    [ "$MISSING" -gt 0 ] && NFAIL=$((NFAIL+MISSING)) || NFAIL=$((NFAIL+1))
elif [ "$RUN_STATUS" -ne 0 ] && [ "$NFAIL" -eq 0 ]; then
    echo "!!! job execution failed (exit $RUN_STATUS)" >&2
    NFAIL=1
fi
if [ "$NFAIL" -gt 0 ]; then
    echo "$NFAIL of $NSEL jobs failed:"
    [ "${#FAILED[@]}" -eq 0 ] || sort "${FAILED[@]}" | sed 's/^/    /'
    [ "$NFAIL" -le 125 ] || NFAIL=125
    exit "$NFAIL"
fi
echo "all $NSEL jobs completed (${#ARCHIVED[@]} archive-only logs not rerun); now run  git diff  to compare with the committed logs"
