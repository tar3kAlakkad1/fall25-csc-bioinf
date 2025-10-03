# Testing Guide for TRViz

This document explains how to run tests for both Python and Codon implementations.

## Python Tests

To run the unified Python test file:

```bash
cd /path/to/fall25-csc-bioinf
python week2/test.py
```

The test file will automatically detect that it's running in Python and import from `week2/python/trviz/`.

## Codon Tests

Due to Codon's import resolution (which works relative to the file location, not the current working directory), the unified test file works differently:

### Option 1: Run Individual Codon Test Files

The recommended way to test the Codon implementation is to run the individual test files:

```bash
cd /path/to/fall25-csc-bioinf/week2/codon
codon run tests/test_decomposer.codon
codon run tests/test_motif_encoder.codon
codon run tests/test_motif_aligner.codon
codon run tests/test_utils.codon
codon run tests/test_refinement.codon
```

### Option 2: Run All Tests at Once

```bash
cd /path/to/fall25-csc-bioinf/week2/codon
for test in tests/*.codon; do echo "Running $test..."; codon run "$test"; done
```

## Unified Test File (`week2/test.py`)

The `week2/test.py` file contains a unified test suite that works for both Python and Codon:

- **Python**: Automatically imports from `week2/python/trviz/`
- **Codon**: Due to import path resolution,  use the individual test files in `week2/codon/tests/` instead

The file uses the `__codon__` built-in constant to detect which runtime it's in:

```python
try:
    if __codon__:
        IS_CODON = True
except NameError:
    IS_CODON = False
```

## Summary

- **For Python**: Run `python week2/test.py` from the project root
- **For Codon**: Run individual test files from `week2/codon/` directory

