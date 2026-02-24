#!/usr/bin/env bash
# run_tests.sh - Run all RVV test binaries and report results
#
# Usage:
#   ./scripts/run_tests.sh local              # Run locally (on riscv64)
#   ./scripts/run_tests.sh ssh HOST PORT      # Run via SSH
#
# Exit code: 0 if all tests pass, 1 if any fail

set -euo pipefail

MODE="${1:-local}"
SSH_HOST="${2:-}"
SSH_PORT="${3:-22}"

# Colors (if terminal supports it)
if [ -t 1 ]; then
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[0;33m'
    CYAN='\033[0;36m'
    NC='\033[0m'
else
    RED='' GREEN='' YELLOW='' CYAN='' NC=''
fi

# Find all test binaries via MANIFEST
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_DIR"

MANIFEST="tests/MANIFEST"
if [ ! -f "$MANIFEST" ]; then
    echo -e "${YELLOW}No MANIFEST found. Run 'python3 generate_tests.py' first.${NC}"
    exit 1
fi

# Parse MANIFEST: each line is "<path.S> <check_count>"
# Build parallel arrays of test binaries and their check counts
TESTS=()
CHECKS=()
while IFS=' ' read -r spath ncheck; do
    [ -z "$spath" ] && continue
    bin="${spath%.S}.bin"
    [ -f "$bin" ] || continue
    TESTS+=("$bin")
    CHECKS+=("${ncheck:-0}")
done < "$MANIFEST"

TOTAL=${#TESTS[@]}
if [ "$TOTAL" -eq 0 ]; then
    echo -e "${YELLOW}No test binaries found. Run 'make build' first.${NC}"
    exit 1
fi

# Sum total checks
TOTAL_CHECKS=0
for n in "${CHECKS[@]}"; do
    TOTAL_CHECKS=$((TOTAL_CHECKS + n))
done

echo -e "${CYAN}Running $TOTAL tests ($TOTAL_CHECKS checks)...${NC}"
echo "=========================================="

PASSED=0
FAILED=0
PASSED_CHECKS=0
FAILED_CHECKS=0
ERRORS=""

for i in "${!TESTS[@]}"; do
    test_bin="${TESTS[$i]}"
    ncheck="${CHECKS[$i]}"
    test_name="${test_bin%.bin}"
    test_name="${test_name#tests/}"

    # Run the test
    set +e
    if [ "$MODE" = "local" ]; then
        "./$test_bin" >/dev/null 2>&1
        rc=$?
    elif [ "$MODE" = "ssh" ]; then
        ssh -p "$SSH_PORT" "$SSH_HOST" \
            "cd ~/rvv-tests && ./$test_bin" >/dev/null 2>&1
        rc=$?
    else
        echo "Unknown mode: $MODE"
        exit 1
    fi
    set -e

    if [ "$rc" -eq 0 ]; then
        PASSED=$((PASSED + 1))
        PASSED_CHECKS=$((PASSED_CHECKS + ncheck))
        printf "${GREEN}  PASS${NC}  %-45s  ${CYAN}(%d checks)${NC}\n" "$test_name" "$ncheck"
    else
        FAILED=$((FAILED + 1))
        # On failure, rc is the check number that failed (1-based),
        # so checks before it passed and the rest are unknown.
        # We count failed_check=1, passed_from_this=rc-1 (checks before the failing one).
        if [ "$rc" -le "$ncheck" ] && [ "$rc" -gt 0 ]; then
            passed_before=$((rc - 1))
            PASSED_CHECKS=$((PASSED_CHECKS + passed_before))
            FAILED_CHECKS=$((FAILED_CHECKS + 1))
            # Remaining checks after the failed one are not executed
            skipped=$((ncheck - rc))
            printf "${RED}  FAIL${NC}  %-45s  ${RED}(check %d/%d failed)${NC}\n" \
                "$test_name" "$rc" "$ncheck"
        else
            # Unexpected exit code (signal, etc.) — all checks unknown
            FAILED_CHECKS=$((FAILED_CHECKS + ncheck))
            printf "${RED}  FAIL${NC}  %-45s  ${RED}(exit code: %d, %d checks)${NC}\n" \
                "$test_name" "$rc" "$ncheck"
        fi
        ERRORS="${ERRORS}\n  ${test_name} (exit code: ${rc})"
    fi
done

echo "=========================================="
printf "Results: ${GREEN}%d passed${NC}, ${RED}%d failed${NC} / %d tests\n" \
    "$PASSED" "$FAILED" "$TOTAL"
printf "Checks:  ${GREEN}%d passed${NC}, ${RED}%d failed${NC} / %d total\n" \
    "$PASSED_CHECKS" "$FAILED_CHECKS" "$TOTAL_CHECKS"

if [ "$FAILED" -gt 0 ]; then
    echo -e "\n${RED}Failed tests:${NC}${ERRORS}"
    exit 1
else
    echo -e "\n${GREEN}All tests passed!${NC}"
    exit 0
fi
