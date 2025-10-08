# How this document was created
The prompt used to create this document is found below. The `ai.md` is being used in Cursor for this project and I usually break steps down in a Markdown file to start using AI in helping complete a task. 

## Prompt (to create my `ai.md`)
```
**Goal:** First, carefully read the task description. Then, read the Notes.md that I have made after reading the task description. Finally, create an ai.md that outlines the steps to complete the task below. Required source code files are attached to this message. Once done, answer the following question: do I need to also port the file tree.pyx to Codon? 

**Task:**
Week 3 Deliverable
Updated 2 days ago by Ibrahim Numanagić
Added to reading list
Preliminaries
Use the same setup as in Week 1.

Week 3 deliverable
The goals of this week’s deliverable are:

code analysis: get familiar with the unknown codebases;
evolutionary trees: port popular biotite’s phylo module to Codon;
Your job is to port biotite’s phylo package to Codon. Note that the source is written in Cython, not pure Python.

The source code can be found on Github. You only need to port the parts needed to run the two tests in test_phylo.py:

test_distances,
test_upgma,
test_neighbor_joining.
While you can use from python import for missing pieces, note that this is not needed here.

You must run Python and Codon tests and output the time of all tests in milliseconds. Use the time module in Python and Codon (do not time this from your evaluate.sh).

This assignment should be done alone.

What (and when) to submit?
Submit a single link (the repository URL pointing to the commit to be graded) via Brightspace before Sun, Oct 12, at 23:59. The final commit must be successfully checked by the GitHub CI.

Example:

https://github.com/Jong-hun-Park/trviz/commit/894be7b424d8bf8ad4ec474ab8caa97dccd3181e 
Your commit must have a working evaluate.sh (or any other extension) that is called from CI and that produces the following output:

Language    Runtime
-------------------
python      2000ms
codon       1000ms
The runtime is the time needed to run all tests (in milliseconds).

Important
Do not submit whole biotite. Assume that it is installed (and install it in CI). Just extract a few files that are needed for the assignment.
Grading
You will get one point for each ported test. Please use the following CI setup: https://github.com/inumanag/fall25-csc-bioinf/blob/main/.github/workflows/actions.yml. Don’t forget to install biotite in CI.

You are expected to understand the code. You might be called to defend your submission at random during office hours.

Please describe all steps and gotchas in report.md in week3/ directory. Also include the time needed to complete this assignment (an estimate is OK). And don’t forget ai.md as well.



**Notes:** Here's a Notes.md that I made when going through the task description. 

# Goal
Port `biotite`’s phylo package to Codon. The source code is written in Cython, not pure Python.

# Source code
The source can be found on [GitHub](https://github.com/biotite-dev/biotite/tree/v1.4.0/src/biotite/sequence/phylo).

# Steps
## Getting necessary files
### Source Code
There are four files found under `biotite/src/biotite/sequence/phylo`:

1) `__init__.py`: some necessary stuff for managing the Python package.

2) `nj.pyx`: a Cython file that only contains the function `neighbor_joining(np.ndarray distances)`. The function performs "hierarchical clustering using the *neighbor joining* algorithm." In contrast to UPGMA this algorithm does not assume a constant evolution rate. The resulting tree is considered to be unrooted.

3) `tree.pyx`: a Cython file that defines an object-oriented implementation of a `class Tree` with it's relevant functions. As described in the docstrings, a `class Tree` "represents a rooted tree (e.g. alignment guide tree or phylogenetic tree)."

4) `upgma.pyx`: a Cython file that only contains the function `upgma(np.ndarray distances)`. The function performs "hierarchical clustering using the *unweighted pair group method with arithmetic mean* (UPGMA)." This algorithm produces leaf nodes with the same distance to the root node. In the context of evolution this means a constant evolution rate (molecular clock).

### Tests
According to the deliverable instructions posted by the instructor on [piazza](https://piazza.com/class/mevafycmxgp28j/post/53), we only need the tests `test_distances`, `test_upgma`, and `test_neighbour_joining` from the test file `test_phylo.py`. 

The [test file](https://github.com/biotite-dev/biotite/blob/v1.4.0/tests/sequence/test_phylo.py) can be found on Github as well.

## Setting up Cython/Python
With the code files needed identified, main objective now is to get the Cython/Python tests running.

### Source code path
The Cython/Python src code will be in `fall25-csc-bioinf/week3/python/src` and the test file will be in `fall25-csc-bioinf/week3/python/test`.

### Getting tests running
No clue how to do this at the time of writing this tbh but we will figure it out lol.

## Setting up Codon
Once Cython/Python side of things is working, we will start porting the identified files to Codon (including the test file). 

### Source code path
The Codon src code will be in `fall25-csc-bioinf/week3/codon/src` and the test file will be in `fall25-csc-bioinf/week3/python/test`. 

## Setting up a CI/CD test
### Timing Python and Codon results
Must run Python and Codon tests and output the time of all tests in milliseconds. Use the time module in Python and Codon (do not time this from your evaluate.sh).

### Creating the `evaluate.sh`
Create a `evaluate.sh` bash script to use in `.github/workflows/actions.yml`

### Update `actions.yml`
Use `evaulate.sh` script in `.github/workflows/actions.yml` and change the code found in `.github/workflows/actions.yml` to trigger Week 3's evaluate script.
```

# ai.md - Week 3 Phylo Module Porting Plan

## Overview
This document outlines the steps to port biotite's phylo module from Cython to Codon, focusing on the three required tests: `test_distances`, `test_upgma`, and `test_neighbor_joining`.

## Files Analysis

### Source Files to Port
1. **tree.pyx** (MUST PORT - ~800 lines)
   - Contains `Tree` class
   - Contains `TreeNode` class (cdef class)
   - Contains helper functions: `_get_leaves`, `_get_leaf_count`, `_create_path_to_root`, `as_binary`, `_as_binary`
   - Contains `TreeError` exception
   - Dependencies: numpy, networkx (only for as_graph method - not needed), copy module

2. **upgma.pyx** (~150 lines)
   - Contains `upgma()` function
   - Dependencies: Tree, TreeNode from tree.pyx, numpy

3. **nj.pyx** (~170 lines)
   - Contains `neighbor_joining()` function  
   - Dependencies: Tree, TreeNode from tree.pyx, numpy

4. **__init__.py** (simple)
   - Package initialization
   - Just imports from the three modules

### Test Files
- **test_phylo.py**: Extract only the three required tests:
  - `test_distances` (lines 128-136)
  - `test_upgma` (lines 41-56)
  - `test_neighbor_joining` (lines 59-90)
  - Also need fixtures: `distances`, `upgma_newick`, `tree`

## Dependency Analysis

### What Each Test Needs

**test_upgma:**
- `phylo.upgma()` → needs upgma.pyx
- `phylo.Tree.from_newick()` → needs Tree class from tree.pyx
- `tree.get_distance()` → needs Tree.get_distance() method

**test_neighbor_joining:**
- `phylo.neighbor_joining()` → needs nj.pyx
- `phylo.TreeNode()` → needs TreeNode class from tree.pyx
- `phylo.Tree()` → needs Tree class from tree.pyx
- TreeNode equality comparison

**test_distances:**
- `tree.root` → needs Tree.root property
- `tree.leaves` → needs Tree.leaves property
- `tree.root.distance_to()` → needs TreeNode.distance_to() method
- `tree.get_distance()` → needs Tree.get_distance() method

### Critical Dependencies
All tests depend on:
- `Tree` class (from tree.pyx)
- `TreeNode` class (from tree.pyx)

## Porting Strategy

### Phase 1: Setup Environment (Week 3 Directory Structure)
```
fall25-csc-bioinf/week3/
├── python/
│   ├── src/
│   │   ├── __init__.py
│   │   ├── tree.pyx
│   │   ├── upgma.pyx
│   │   └── nj.pyx
│   └── test/
│       └── test_phylo.py (extracted tests only)
├── codon/
│   ├── src/
│   │   ├── tree.codon
│   │   ├── upgma.codon
│   │   └── nj.codon
│   └── test/
│       └── test_phylo.codon
├── evaluate.sh
├── report.md
├── ai.md
└── README.md
```

### Phase 2: Get Python/Cython Tests Working First

**Step 2.1: Copy source files**
- Copy tree.pyx, upgma.pyx, nj.pyx, __init__.py to `week3/python/src/`
- Copy test_phylo.py to `week3/python/test/`

**Step 2.2: Extract only required tests**
- From test_phylo.py, keep only:
  - Required fixtures: `distances()`, `upgma_newick()`, `tree()`
  - Required tests: `test_distances()`, `test_upgma()`, `test_neighbor_joining()`
- Remove all other tests and fixtures

**Step 2.3: Setup Cython build**
- Create setup.py for building .pyx files
- Build the Cython extensions
- Verify tests run successfully with pytest

**Step 2.4: Add timing to Python tests**
```python
import time
start = time.time()
# run tests
end = time.time()
print(f"python      {int((end-start)*1000)}ms")
```

### Phase 3: Port to Codon (File by File)

**Step 3.1: Port tree.pyx → tree.codon**

Key challenges:
- **Cython cdef classes** → Codon classes (regular Python-like classes)
- **Type annotations**: Convert Cython type hints to Codon type hints
  - `cdef int` → `int`
  - `cdef float` → `float`
  - `cdef TreeNode` → `TreeNode`
  - Memory views like `float32[:]` → use List[float] or similar
- **Properties**: Keep @property decorators, should work in Codon
- **Static methods**: Keep @staticmethod, should work in Codon
- **Recursion**: Codon handles recursion, but watch for stack depth
- **None handling**: Codon is stricter about None/Optional types
- **Exception handling**: Keep TreeError exception class

Order of porting within tree.codon:
1. `TreeError` exception class (trivial)
2. `TreeNode` class:
   - `__init__` method
   - Properties (parent, children, index, distance)
   - `_set_parent` method
   - `is_leaf()`, `is_root()`, `as_root()` methods
   - `copy()` method
   - `get_leaves()` method (needs helper)
   - `distance_to()` and `lowest_common_ancestor()` methods
   - `from_newick()` static method (complex parsing logic)
   - `__eq__` and `__hash__` methods
3. Helper functions:
   - `_get_leaves()`
   - `_create_path_to_root()`
4. `Tree` class:
   - `__init__` method
   - Properties (root, leaves)
   - `get_distance()` method
   - `from_newick()` static method
   - `__len__`, `__eq__`, `__hash__` methods

**Codon-specific considerations:**
- Use `Optional[T]` for nullable types
- Use `List[T]` instead of Python lists in type hints
- Use `Tuple[T, ...]` for tuples
- May need to use `@python` decorator for any complex Python interop
- Codon numpy support is good, but check compatibility

**Step 3.2: Port upgma.pyx → upgma.codon**

Key challenges:
- **Numpy arrays**: Codon has numpy support, but check:
  - Memory views (`float32[:,:]`) → may need regular numpy arrays
  - `np.loadtxt`, `np.allclose` - check if available
  - Array slicing and indexing
- **Type conversions**: `.astype(np.float32, copy=True)`
- **Cython decorators**: Remove `@cython.boundscheck(False)` and `@cython.wraparound(False)`

**Step 3.3: Port nj.pyx → nj.codon**

Similar challenges to upgma.pyx:
- Memory views and numpy arrays
- Algorithmic logic should port directly
- Watch for integer division (Python 3 behavior)

### Phase 4: Port Tests to Codon

**Step 4.1: Create test_phylo.codon**

Key considerations:
- Codon uses pytest-compatible testing
- Import statements: `from python import pytest` if needed
- Fixtures: Should work similarly
- Assertions: Standard Python assertions should work
- File I/O: `np.loadtxt` for loading distance matrix
- String operations: For loading newick strings

**Step 4.2: Add timing**
```python
import time
start = time.time()
# run tests  
end = time.time()
print(f"codon       {int((end-start)*1000)}ms")
```

### Phase 5: Testing and Debugging

**Step 5.1: Test each component individually**
- Test TreeNode creation and basic operations
- Test Tree creation from TreeNode
- Test from_newick parsing
- Test upgma algorithm
- Test neighbor_joining algorithm

**Step 5.2: Run full test suite**
- Ensure all three tests pass
- Compare results with Python version
- Verify timing output format

**Step 5.3: Debug common issues**
- Type mismatches (Codon is stricter)
- None/Optional handling
- Integer vs float division
- Array indexing and bounds
- Recursion depth (if any issues)

### Phase 6: CI/CD Setup

**Step 6.1: Create evaluate.sh**
```bash
#!/bin/bash

echo "Language    Runtime"
echo "-------------------"

# Run Python tests
cd week3/python
python_time=$(python test_phylo_timed.py)
echo "$python_time"

# Run Codon tests  
cd ../codon
codon_time=$(codon run test_phylo_timed.codon)
echo "$codon_time"
```

**Step 6.2: Update .github/workflows/actions.yml**
- Add biotite installation: `pip install biotite`
- Add pytest installation for Python tests
- Ensure Codon is installed
- Call `./week3/evaluate.sh`
- Verify output format matches requirements

**Step 6.3: Test CI locally**
- Use act or similar to test GitHub Actions locally
- Verify all dependencies install correctly
- Verify tests run and output is correct

## Potential Gotchas and Challenges

### Cython → Codon Translation Issues

1. **Memory Views**
   - Cython: `cdef float32[:] array_view`
   - Codon: May need to use regular numpy arrays or List[float]

2. **Typed Definitions**
   - Cython: `cdef int i = 0`
   - Codon: `i: int = 0` or just `i = 0` with type inference

3. **Extension Types (cdef class)**
   - Cython: `cdef class TreeNode:`
   - Codon: Regular `class TreeNode:` should work

4. **String Parsing (in from_newick)**
   - Complex nested logic with string manipulation
   - May need careful testing of edge cases
   - Watch for off-by-one errors

5. **Tuple vs List**
   - TreeNode uses tuples for immutable children
   - Ensure proper typing: `Tuple[TreeNode, ...]`

6. **None/Optional Types**
   - Codon is stricter about Optional types
   - Need explicit `Optional[TreeNode]` for nullable fields

7. **Numpy Compatibility**
   - Check which numpy functions are available in Codon
   - May need workarounds for unsupported functions
   - Test array operations carefully

8. **Recursion**
   - from_newick uses recursion
   - get_leaves uses recursion
   - Ensure no stack overflow issues

9. **Float Precision**
   - Watch for floating point comparison issues
   - May need pytest.approx equivalent in Codon

10. **File I/O**
    - Loading test data (distances.txt, newick files)
    - Ensure Codon can read these files

### Testing Gotchas

1. **Fixture Dependencies**
   - tree fixture depends on distances fixture
   - Ensure proper order in Codon tests

2. **Pytest Compatibility**
   - Check if Codon's pytest support covers all features used
   - @pytest.fixture decorator support
   - pytest.approx for floating point comparison

3. **Assertion Messages**
   - May need to adjust assertion format for Codon

## Minimal Viable Product (MVP) Scope

To get tests passing, focus on:

**From tree.pyx (REQUIRED):**
- TreeNode class: __init__, properties, is_leaf, is_root, as_root, _set_parent, copy, distance_to, lowest_common_ancestor, get_leaves, from_newick, __eq__
- Tree class: __init__, root, leaves, get_distance, from_newick, __len__, __eq__
- Helper functions: _get_leaves, _create_path_to_root
- TreeError exception

**NOT needed (can skip):**
- Tree.as_graph() method (uses networkx)
- Tree.to_newick() method (only reading, not writing)
- TreeNode.to_newick() method
- as_binary() function
- get_indices() method
- get_leaf_count() method
- __hash__ methods (unless needed for equality)

**From upgma.pyx (REQUIRED):**
- upgma() function

**From nj.pyx (REQUIRED):**
- neighbor_joining() function

## Verification Steps

1. ✅ Python tests pass (all 3)
2. ✅ Python timing output correct format
3. ✅ Codon code compiles without errors
4. ✅ Codon tests pass (all 3)
5. ✅ Codon timing output correct format
6. ✅ evaluate.sh produces correct output
7. ✅ CI workflow runs successfully
8. ✅ Timing comparison makes sense (Codon should be faster)

## Time Estimate

- Phase 1 (Setup): 30-45 minutes
- Phase 2 (Python working): 1-2 hours
- Phase 3 (Port to Codon): 6-10 hours
  - tree.pyx: 4-6 hours (most complex)
  - upgma.pyx: 1-2 hours
  - nj.pyx: 1-2 hours
- Phase 4 (Port tests): 1-2 hours
- Phase 5 (Testing/Debug): 2-4 hours
- Phase 6 (CI/CD): 1 hour

**Total Estimate: 12-20 hours**

## Success Criteria

1. All three tests pass in both Python and Codon
2. Timing output matches required format
3. CI runs successfully
4. Code is clean and well-documented
5. report.md explains challenges and solutions