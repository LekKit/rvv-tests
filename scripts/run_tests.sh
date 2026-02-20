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

# Find all test binaries
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_DIR"

TESTS=$(find tests -name '*.bin' -type f | sort)
TOTAL=$(echo "$TESTS" | wc -l | tr -d ' ')

if [ "$TOTAL" -eq 0 ]; then
    echo -e "${YELLOW}No test binaries found. Run 'make build' first.${NC}"
    exit 1
fi

echo -e "${CYAN}Running $TOTAL RVV tests...${NC}"
echo "=========================================="

PASSED=0
FAILED=0
ERRORS=""

run_test() {
    local test_bin="$1"
    local test_name="${test_bin%.bin}"
    test_name="${test_name#tests/}"

    if [ "$MODE" = "local" ]; then
        "./$test_bin" 2>/dev/null
        return $?
    elif [ "$MODE" = "ssh" ]; then
        ssh -p "$SSH_PORT" "$SSH_HOST" "cd $PROJECT_DIR && ./$test_bin" 2>/dev/null
        return $?
    fi
}

for test_bin in $TESTS; do
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
        printf "${GREEN}  PASS${NC}  %s\n" "$test_name"
    else
        FAILED=$((FAILED + 1))
        printf "${RED}  FAIL${NC}  %s  (exit code: %d)\n" "$test_name" "$rc"
        ERRORS="${ERRORS}\n  ${test_name} (exit code: ${rc})"
    fi
done

echo "=========================================="
echo -e "Results: ${GREEN}${PASSED} passed${NC}, ${RED}${FAILED} failed${NC} / ${TOTAL} total"

if [ "$FAILED" -gt 0 ]; then
    echo -e "\n${RED}Failed tests:${NC}${ERRORS}"
    exit 1
else
    echo -e "\n${GREEN}All tests passed!${NC}"
    exit 0
fi
