# AGENTS.md

This file provides guidance to Agents when working with code in this repository.

## Repository Overview

This is a CSC 427 Bioinformatics course repository focusing on porting Python/Cython bioinformatics libraries to Codon (a Python compiler). Each weekly deliverable involves implementing and optimizing bioinformatics algorithms in both Python and Codon.

## Project Structure

The repository is organized by weekly deliverables:

- **week1/**: De Bruijn graph implementation for genome assembly
- **week2/**: TRViz - Tandem Repeat Visualizer library (Python → Codon port)
- **week3/**: Biotite phylogenetic tree module (Cython → Codon port)
- **week4/**: Sequence alignment algorithms (Python → Codon port)

Each week contains:
- `python/` - Reference Python implementation
- `codon/` - Codon port (when applicable)
- `Report.md` - Detailed implementation and porting notes
- `ai.md` - AI-assisted development notes (when applicable)
- `evaluate.sh` or test scripts for benchmarking

## Common Development Commands

### Week 2: TRViz (Tandem Repeat Visualization)

**Python setup and tests:**
```bash
cd week2/python
pip install -e .  # Install trviz package in development mode
pytest tests/     # Run all tests
```

**Codon tests:**
```bash
cd week2/codon
# Run individual test files
codon run tests/test_decomposer.codon
codon run tests/test_motif_encoder.codon
codon run tests/test_motif_aligner.codon

# Or run all tests
for test in tests/*.codon; do codon run "$test"; done
```

**Environment for Python interop (required for HMM mode):**
```bash
export CODON_PYTHON=/path/to/libpython3.11.dylib
# Find path with: python3-config --prefix
```

### Week 3: Phylogenetic Trees (UPGMA & Neighbor-Joining)

**Python tests:**
```bash
cd week3/python/test
python3 run_timed_tests.py  # Runs test_distances, test_upgma, test_neighbor_joining
```

**Codon tests:**
```bash
cd week3/codon
codon run test_phylo.py
```

**Evaluation (both languages):**
```bash
bash week3/evaluate.sh
```

### Week 4: Sequence Alignment Algorithms

**Build Codon executable:**
```bash
cd week4
codon build -release -o codon/align codon/main.py
```

**Run alignments manually:**
```bash
# Python
python3 week4/python/main.py run global \
  week4/python/data/MT-human.fa MT_human \
  week4/python/data/MT-orang.fa MT_orang

# Codon (after building)
week4/codon/align run global \
  week4/python/data/MT-human.fa MT_human \
  week4/python/data/MT-orang.fa MT_orang
```

**Full evaluation with timing:**
```bash
bash week4/evaluate.sh
```

This runs all four algorithms (global, local, semi-global, affine) on all test datasets and prints a formatted runtime comparison table.

### CI/CD

The repository uses GitHub Actions (`.github/workflows/actions.yml`):

**Current CI setup:**
- Installs Codon compiler
- Installs Python dependencies (NumPy, etc.)
- Runs `week4/evaluate.sh` for testing and benchmarking
- Handles macOS-specific libomp setup for Codon/NumPy interop

To run CI checks locally, execute the same evaluation scripts used in CI.

## Key Porting Patterns (Python/Cython → Codon)

### Type System Differences

**Import and type syntax:**
```python
# Python
from typing import Dict, List
def foo(x: List[str]) -> Dict[str, int]:
    result: Dict = {}

# Codon
# No imports needed - use native syntax
def foo(x: list[str]) -> dict[str, int]:
    result = dict[str, int]()  # Explicit initialization required
```

**No None initialization:**
```python
# Python
self.data = None  # Will assign later

# Codon
self.data = dict[str, str]()  # Must initialize with proper type
```

### Cython to Codon Conversions

**Remove Cython-specific constructs:**
```python
# Cython
cdef class TreeNode:
    cdef int value
    cdef bint is_leaf

def __cinit__(self):
    pass

# Codon
class TreeNode:
    value: int
    is_leaf: bool

def __init__(self):
    pass
```

**Memory views to NumPy arrays:**
```python
# Cython
def process(float[:,:] matrix):
    pass

# Codon
def process(matrix):  # Use standard NumPy array
    pass
```

### Testing Differences

**pytest → Plain functions:**
```python
# Python
import pytest

@pytest.fixture
def data():
    return [1, 2, 3]

@pytest.mark.parametrize("x,y", [(1,2), (3,4)])
def test_foo(data, x, y):
    assert x < y

# Codon
def test_foo():
    data = [1, 2, 3]
    # Inline test cases
    assert 1 < 2, "Test case 1 failed"
    assert 3 < 4, "Test case 2 failed"
```

### Python Interop in Codon

When external Python libraries (like pomegranate) are needed:

```python
@python
def use_external_library():
    from python import pomegranate as pg
    # Use pg as normal
```

Remember to set `CODON_PYTHON` environment variable for Python interop to work.

### String Formatting

**Only f-strings supported:**
```python
# Python - multiple formats work
"Value: {}".format(x)
"Value: %s" % x
f"Value: {x}"

# Codon - only f-strings
f"Value: {x}"
```

### Function Parameters

**No **kwargs:**
```python
# Python
def decompose(seq, motifs, **kwargs):
    score = kwargs.get('match_score', 1)

# Codon - explicit parameters
def decompose(seq: str, motifs: list[str],
              match_score: int = 1,
              mismatch_score: int = -1):
    pass
```

## Architecture Notes

### Week 2: TRViz Module Structure

The TRViz library has four main components:
1. **Decomposer** - Decomposes tandem repeat sequences into motif units using DP or HMM
2. **Motif Encoder** - Maps motifs to ASCII symbols for visualization
3. **Motif Aligner** - Aligns encoded motif sequences (uses restructured score matrix in Codon)
4. **Utils** - FASTA parsing and helper functions

**Key Codon adaptations:**
- HMM functionality separated into `hmm_helper.codon` to avoid loading Python runtime unless needed
- Backtrack arrays changed from `dtype=object` to explicit int64 arrays with extra dimension
- Score matrices restructured as nested dicts (uniform types required)

### Week 3: Phylogenetic Tree Implementation

**Core classes:**
- `TreeNode` - Individual nodes with parent/child relationships
- `Tree` - Wrapper managing root node with traversal methods

**Algorithms:**
- `upgma()` - UPGMA clustering (assumes molecular clock)
- `neighbor_joining()` - NJ clustering (no molecular clock assumption)

Both use distance matrices and build trees bottom-up by iteratively merging closest clusters.

### Week 4: Sequence Alignment Algorithms

Four DP-based alignment algorithms implemented:
1. **Global (Needleman-Wunsch)** - Aligns entire sequences
2. **Local (Smith-Waterman)** - Finds best local alignment
3. **Semi-global (fitting)** - No penalty for leading/trailing gaps
4. **Affine gap** - Uses three matrices (M, Ix, Iy) for gap-open vs gap-extend penalties

**Scoring schemes:**
- Match: +3, Mismatch: -3, Gap: -2
- Affine: gap-open = -5, gap-extend = -1

**Implementation differences:**
- Python: NumPy-backed matrices
- Codon: Explicit loops with typed arrays (loop-based DP)
- Codon affine uses float matrices with `-inf` for clean initialization

## macOS Development Notes

**OpenMP library for Codon:**
```bash
brew install libomp
export DYLD_LIBRARY_PATH="/opt/homebrew/opt/libomp/lib:$DYLD_LIBRARY_PATH"
```

The `evaluate.sh` scripts attempt to handle this automatically when Homebrew is available.

## Common Gotchas

1. **Mixed-type dictionaries not supported in Codon** - Restructure as nested dicts with uniform types
2. **Codon imports are relative to file location** - Not to current working directory
3. **NumPy dtype=object not supported** - Use explicit typed arrays
4. **Floating-point comparisons** - Use tolerance-based equality checks (≈1e-3) for tree distances
5. **Implicit string concatenation doesn't work** - Use explicit `+` operator
6. **TypeError bugs in original Python code** - Fixed during porting (missing `raise` keywords found)
7. **Test framework differences** - No pytest fixtures/parametrize in Codon

## Performance Expectations

Codon typically provides 2-10x speedup over pure Python for compute-intensive tasks. Performance gains come from:
- Native compilation (no interpreter overhead)
- Static type checking enabling optimizations
- Efficient loop execution
- Compiled NumPy operations

For I/O-bound or NumPy-heavy operations, speedups may be smaller since NumPy already uses C backends.
