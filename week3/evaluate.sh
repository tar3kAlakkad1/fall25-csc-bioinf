#!/bin/bash

# Navigate to the week3 directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Print the header
echo "Language    Runtime"
echo "-------------------"

# Run Python tests
cd python/test
PYTHON_OUTPUT=$(python3 run_timed_tests.py 2>&1)
PYTHON_EXIT_CODE=$?
cd "$SCRIPT_DIR"

if [ $PYTHON_EXIT_CODE -ne 0 ]; then
    echo "Python tests failed!" >&2
    echo "$PYTHON_OUTPUT" >&2
    exit 1
fi

echo "$PYTHON_OUTPUT"

# Run Codon tests
cd codon
CODON_OUTPUT=$(codon run test_phylo.py 2>&1)
CODON_EXIT_CODE=$?
cd "$SCRIPT_DIR"

if [ $CODON_EXIT_CODE -ne 0 ]; then
    echo "Codon tests failed!" >&2
    echo "$CODON_OUTPUT" >&2
    exit 1
fi

echo "$CODON_OUTPUT"

exit 0

