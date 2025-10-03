#!/bin/bash

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color
BOLD='\033[1m'

# Function to print section headers
print_header() {
    echo -e "\n${BOLD}${CYAN}════════════════════════════════════════════════════════════════${NC}"
    echo -e "${BOLD}${CYAN}  $1${NC}"
    echo -e "${BOLD}${CYAN}════════════════════════════════════════════════════════════════${NC}\n"
}

# Function to print subsection headers
print_subsection() {
    echo -e "\n${BOLD}${MAGENTA}── $1 ──${NC}\n"
}

# Function to print success
print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

# Function to print error
print_error() {
    echo -e "${RED}✗ $1${NC}"
}

# Function to print info
print_info() {
    echo -e "${BLUE}ℹ $1${NC}"
}

# Function to format time
format_time() {
    local seconds=$1
    if (( $(echo "$seconds < 1" | bc -l) )); then
        printf "%.3fs" $seconds
    else
        printf "%.2fs" $seconds
    fi
}

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Initialize counters
TOTAL_START=$(date +%s.%N)
CODON_PASS=0
CODON_FAIL=0
PYTHON_PASS=0
PYTHON_FAIL=0

print_header "TRViz Test Suite"
print_info "Running all tests for Codon and Python implementations"
print_info "Working directory: $SCRIPT_DIR"

# ============================================================================
# CODON TESTS
# ============================================================================
print_header "CODON TESTS"

# Check if codon is installed
if ! command -v codon &> /dev/null; then
    print_error "Codon is not installed or not in PATH"
    print_info "Skipping Codon tests"
    CODON_SKIP=1
else
    print_success "Codon found: $(which codon)"
    
    CODON_TEST_DIR="$SCRIPT_DIR/codon/tests"
    
    if [ ! -d "$CODON_TEST_DIR" ]; then
        print_error "Codon test directory not found: $CODON_TEST_DIR"
        CODON_SKIP=1
    else
        cd "$CODON_TEST_DIR"
        
        # Find all test files
        CODON_TEST_FILES=(test_*.codon)
        
        print_info "Found ${#CODON_TEST_FILES[@]} test file(s)"
        
        for test_file in "${CODON_TEST_FILES[@]}"; do
            if [ -f "$test_file" ]; then
                print_subsection "Running $test_file"
                
                TEST_START=$(date +%s.%N)
                
                if codon run "$test_file"; then
                    TEST_END=$(date +%s.%N)
                    TEST_TIME=$(echo "$TEST_END - $TEST_START" | bc)
                    print_success "$test_file passed"
                    ((CODON_PASS++))
                else
                    TEST_END=$(date +%s.%N)
                    TEST_TIME=$(echo "$TEST_END - $TEST_START" | bc)
                    print_error "$test_file failed"
                    ((CODON_FAIL++))
                fi
            fi
        done
        
        cd "$SCRIPT_DIR"
    fi
fi

# ============================================================================
# PYTHON TESTS
# ============================================================================
print_header "PYTHON TESTS"

# Check if pytest is installed
if ! command -v pytest &> /dev/null; then
    print_error "pytest is not installed or not in PATH"
    print_info "Skipping Python tests"
    print_info "Install with: pip install pytest"
    PYTHON_SKIP=1
else
    print_success "pytest found: $(which pytest)"
    
    PYTHON_TEST_DIR="$SCRIPT_DIR/python/tests"
    
    if [ ! -d "$PYTHON_TEST_DIR" ]; then
        print_error "Python test directory not found: $PYTHON_TEST_DIR"
        PYTHON_SKIP=1
    else
        cd "$SCRIPT_DIR/python"
        
        # Find all test files
        PYTHON_TEST_FILES=(tests/test_*.py)
        
        print_info "Found ${#PYTHON_TEST_FILES[@]} test file(s)"
        
        for test_file in "${PYTHON_TEST_FILES[@]}"; do
            if [ -f "$test_file" ]; then
                print_subsection "Running $test_file"
                
                TEST_START=$(date +%s.%N)
                
                if pytest "$test_file" -v; then
                    TEST_END=$(date +%s.%N)
                    TEST_TIME=$(echo "$TEST_END - $TEST_START" | bc)
                    print_success "$test_file passed"
                    ((PYTHON_PASS++))
                else
                    TEST_END=$(date +%s.%N)
                    TEST_TIME=$(echo "$TEST_END - $TEST_START" | bc)
                    print_error "$test_file failed"
                    ((PYTHON_FAIL++))
                fi
            fi
        done
        
        cd "$SCRIPT_DIR"
    fi
fi

# ============================================================================
# SUMMARY
# ============================================================================
TOTAL_END=$(date +%s.%N)
TOTAL_TIME=$(echo "$TOTAL_END - $TOTAL_START" | bc)

print_header "TEST SUMMARY"

if [ -z "$CODON_SKIP" ]; then
    echo -e "${BOLD}Codon Tests:${NC}"
    echo -e "  ${GREEN}Passed: $CODON_PASS${NC}"
    if [ $CODON_FAIL -gt 0 ]; then
        echo -e "  ${RED}Failed: $CODON_FAIL${NC}"
    else
        echo -e "  Failed: $CODON_FAIL"
    fi
else
    echo -e "${BOLD}Codon Tests:${NC} ${YELLOW}SKIPPED${NC}"
fi

echo ""

if [ -z "$PYTHON_SKIP" ]; then
    echo -e "${BOLD}Python Tests:${NC}"
    echo -e "  ${GREEN}Passed: $PYTHON_PASS${NC}"
    if [ $PYTHON_FAIL -gt 0 ]; then
        echo -e "  ${RED}Failed: $PYTHON_FAIL${NC}"
    else
        echo -e "  Failed: $PYTHON_FAIL"
    fi
else
    echo -e "${BOLD}Python Tests:${NC} ${YELLOW}SKIPPED${NC}"
fi

echo ""
echo -e "${BOLD}Total Time:${NC} $(format_time $TOTAL_TIME)"

# Exit with appropriate code
TOTAL_FAIL=$((CODON_FAIL + PYTHON_FAIL))
if [ $TOTAL_FAIL -gt 0 ]; then
    echo ""
    print_error "Some tests failed!"
    exit 1
else
    echo ""
    print_success "All tests passed!"
    exit 0
fi

