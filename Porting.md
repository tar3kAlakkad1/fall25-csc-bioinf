# Porting Python/Cython Phylogenetic Tree Code to Codon

## Executive Summary

This document outlines the comprehensive steps needed to port the phylogenetic tree implementation from Python/Cython to Codon. The codebase consists of three main Cython source files (`tree.pyx`, `nj.pyx`, `upgma.pyx`) and one Python test file (`test_phylo.py`).

**Key Challenge**: The code heavily uses Cython-specific features (typed memoryviews, cimport, cdef) that don't exist in Codon. While Codon is similar to Python, it requires significant refactoring of Cython code.

---

## 1. Major Incompatibilities and Required Changes

### 1.1 Cython-Specific Syntax Removal

**Issue**: Codon does not support Cython syntax (`cimport`, `cdef`, `ctypedef`, memoryviews).

**Files Affected**: All `.pyx` files

**Changes Required**:

#### Remove Cython Imports and Type Definitions
- **Remove**: `cimport cython`, `cimport numpy as np`
- **Remove**: All `ctypedef` declarations (e.g., `ctypedef np.float32_t float32`)
- **Replace with**: Standard Python imports and Codon's native typing

**Before** (nj.pyx, lines 9-17):
```python
cimport cython
cimport numpy as np

from .tree import Tree, TreeNode
import numpy as np

ctypedef np.float32_t float32
ctypedef np.uint8_t uint8
ctypedef np.uint32_t uint32
```

**After**:
```python
import numpy as np
from tree import Tree, TreeNode

# Use Codon's native types directly
# float32 -> float (Codon uses f32 internally)
# uint8 -> u8
# uint32 -> u32
```

#### Replace Cython Decorators
- **Remove**: `@cython.boundscheck(False)`, `@cython.wraparound(False)`
- Codon performs these optimizations automatically

#### Replace `cdef` Functions and Variables

**Pattern**: `cdef TYPE variable_name` → Standard Python type annotations

**Before** (nj.pyx, lines 76-79):
```python
cdef int i=0, j=0, k=0, u=0
cdef int i_min=0, j_min=0
cdef float32 dist=0, dist_min, dist_sum=0
cdef float32 node_dist_i=0, node_dist_j=0, node_dist_k=0
```

**After**:
```python
i: int = 0
j: int = 0
k: int = 0
u: int = 0
i_min: int = 0
j_min: int = 0
dist: float = 0.0
dist_min: float
dist_sum: float = 0.0
node_dist_i: float = 0.0
node_dist_j: float = 0.0
node_dist_k: float = 0.0
```

#### Replace Typed Memoryviews

**Pattern**: `TYPE[:] variable` or `TYPE[:,:] variable` → Use NumPy arrays with type annotations

**Before** (nj.pyx, lines 101-115):
```python
cdef uint8[:] is_clustered_v = np.full(
    distances.shape[0], False, dtype=np.uint8
)
cdef float32[:] divergence_v = np.zeros(
    distances.shape[0], dtype=np.float32
)
cdef float32[:,:] corr_distances_v = np.zeros(
    (distances.shape[0],) * 2, dtype=np.float32
)
cdef float32[:,:] distances_v = distances.astype(np.float32, copy=True)
```

**After**:
```python
is_clustered_v: np.ndarray[u8] = np.full(
    distances.shape[0], False, dtype=np.uint8
)
divergence_v: np.ndarray[float] = np.zeros(
    distances.shape[0], dtype=np.float32
)
corr_distances_v: np.ndarray[float] = np.zeros(
    (distances.shape[0],) * 2, dtype=np.float32
)
distances_v: np.ndarray[float] = distances.astype(np.float32, copy=True)
```

### 1.2 Class Definition Changes (`tree.pyx`)

**Issue**: Cython's `cdef class` syntax is not supported in Codon.

**Changes Required**:

#### Convert `cdef class` to Standard Python Class

**Before** (tree.pyx, line 347):
```python
cdef class TreeNode:
    cdef int _index
    cdef float _distance
    cdef bint _is_root
    cdef TreeNode _parent
    cdef tuple _children
    
    def __cinit__(self, children=None, distances=None, index=None):
```

**After**:
```python
class TreeNode:
    _index: int
    _distance: float
    _is_root: bool
    _parent: Optional[TreeNode]
    _children: Optional[Tuple[TreeNode, ...]]
    
    def __init__(self, children=None, distances=None, index=None):
```

**Key Changes**:
- Replace `cdef class` with `class`
- Replace `cdef` member variables with type-annotated class attributes
- Replace `__cinit__` with `__init__`
- Replace `bint` with `bool`
- Add `Optional` type hints where appropriate

#### Convert `cdef` Methods to Regular Methods

**Pattern**: `cdef METHOD_NAME(...)` → `def METHOD_NAME(...) -> return_type:`

**Before** (tree.pyx, lines 484-488):
```python
def _set_parent(self, TreeNode parent not None, float distance):
    if self._parent is not None or self._is_root:
        raise TreeError("Node already has a parent")
    self._parent = parent
    self._distance = distance
```

**After**:
```python
def _set_parent(self, parent: TreeNode, distance: float) -> None:
    if self._parent is not None or self._is_root:
        raise TreeError("Node already has a parent")
    self._parent = parent
    self._distance = distance
```

### 1.3 Type Annotation Conversions

**Cython → Codon Type Mapping**:

| Cython Type | Codon Type |
|-------------|------------|
| `bint` | `bool` |
| `int` | `int` |
| `float` | `float` |
| `float32` | `float` (or `f32` if needed) |
| `uint8` | `u8` |
| `uint32` | `u32` |
| `str` | `str` |
| `list` | `List[T]` |
| `tuple` | `Tuple[T, ...]` |
| `dict` | `Dict[K, V]` |
| `np.ndarray` | `np.ndarray[dtype]` |

**Example Conversions**:

```python
# Function parameter annotations
def neighbor_joining(distances: np.ndarray[float]) -> Tree:
    ...

# Optional types (was: TreeNode or None)
_parent: Optional[TreeNode] = None

# Tuple types (was: tuple)
_children: Optional[Tuple[TreeNode, ...]] = None

# List types (was: list)
leaves: List[TreeNode] = []
```

### 1.4 Import Path Changes

**Issue**: Relative imports with `.module` syntax may need adjustment.

**Before**:
```python
from .tree import Tree, TreeNode
from ...file import InvalidFileError
from ...copyable import Copyable
```

**After** (Option 1 - Same directory):
```python
from tree import Tree, TreeNode
# Note: InvalidFileError and Copyable will need to be defined locally
# or imported from Python if they exist elsewhere
```

**After** (Option 2 - Use Python interop for missing modules):
```python
from python import biotite
InvalidFileError = biotite.file.InvalidFileError
Copyable = biotite.copyable.Copyable
```

---

## 2. File-by-File Conversion Plan

### 2.1 `tree.pyx` → `tree.codon` (or `tree.py`)

**Complexity**: High - Contains complex class definitions, recursion, graph operations

**Required Changes**:

1. **Remove Cython imports** (lines 9-10)
   - Remove `cimport cython`, `cimport numpy as np`

2. **Remove NetworkX dependency or use Python interop** (line 14)
   - **Option A**: Port graph operations to pure Codon (complex)
   - **Option B**: Use `from python import networkx as nx` for `as_graph()` method
   - **Recommendation**: Use Option B initially

3. **Remove external dependencies** (lines 15-16)
   - `InvalidFileError` and `Copyable` need to be:
     - Defined locally, OR
     - Imported from Python interop

4. **Convert Tree class** (lines 19-345)
   - Already uses standard Python class syntax
   - Update type hints in `__init__` parameter: `TreeNode root not None` → `root: TreeNode`
   - Update internal variable types

5. **Convert TreeNode cdef class** (lines 347-991)
   - **Critical**: Replace `cdef class TreeNode:` with `class TreeNode:`
   - Add type annotations for all class attributes (lines 424-428)
   - Replace `__cinit__` with `__init__` (line 430)
   - Remove all `cdef` from method definitions
   - Add proper type annotations to all methods

6. **Convert helper functions** (lines 993-1027)
   - Remove `cdef` prefix: `cdef _get_leaves(...)` → `def _get_leaves(...) -> None:`
   - Remove `cdef` prefix: `cdef int _get_leaf_count(...)` → `def _get_leaf_count(...) -> int:`
   - Remove `cdef` prefix: `cdef list _create_path_to_root(...)` → `def _create_path_to_root(...) -> List[TreeNode]:`
   - Add proper type annotations

7. **Update method signatures with type hints**:
   - `get_distance(self, index1, index2, bint topological=False)` → `get_distance(self, index1: int, index2: int, topological: bool = False) -> float`
   - `to_newick(self, labels=None, bint include_distance=True, round_distance=None)` → `to_newick(self, labels: Optional[List[str]] = None, include_distance: bool = True, round_distance: Optional[int] = None) -> str`
   - And so on for all methods

8. **String operations** (lines 860-861)
   - Codon uses ASCII strings, not Unicode
   - Verify string operations work correctly
   - Test with `.split()`, `.join()`, slicing, etc.

### 2.2 `nj.pyx` → `nj.codon` (or `nj.py`)

**Complexity**: Medium - Algorithmic code with heavy NumPy usage

**Required Changes**:

1. **Remove Cython preamble** (lines 9-17)
   - Remove all `cimport` statements
   - Remove all `ctypedef` declarations

2. **Update function signature** (line 25)
   - Remove `@cython` decorators (lines 23-24)
   - Add type annotations:
   
   **Before**:
   ```python
   @cython.boundscheck(False)
   @cython.wraparound(False)
   def neighbor_joining(np.ndarray distances):
   ```
   
   **After**:
   ```python
   def neighbor_joining(distances: np.ndarray[float]) -> Tree:
   ```

3. **Convert local variables** (lines 76-79)
   - Replace all `cdef TYPE var` with `var: TYPE`
   - Use appropriate Codon types

4. **Convert typed arrays** (lines 96-115)
   - Replace memoryviews with NumPy arrays
   - Add type annotations to array variables

5. **Update array indexing**
   - Memoryview syntax `distances_v[i,j]` should work in NumPy
   - Verify all array operations are compatible

6. **Handle MAX_FLOAT constant** (line 20)
   - Verify `np.finfo(np.float32).max` works in Codon-NumPy
   - May need to use Python interop or define constant directly

### 2.3 `upgma.pyx` → `upgma.codon` (or `upgma.py`)

**Complexity**: Medium - Similar structure to nj.pyx

**Required Changes**:

1. **Remove Cython preamble** (lines 9-17)
   - Same as nj.pyx

2. **Update function signature** (line 25)
   - Remove `@cython` decorators
   - Add type annotations:
   
   ```python
   def upgma(distances: np.ndarray[float]) -> Tree:
   ```

3. **Convert local variables** (lines 68-72)
   - Replace all `cdef TYPE var` with `var: TYPE`

4. **Convert typed arrays** (lines 86-103)
   - Replace memoryviews with NumPy arrays

5. **Update return statement** (line 164)
   - Should work as-is

### 2.4 `__init__.py` - Not Needed in Codon

**Complexity**: N/A - Remove this file

**Important Note**: Codon does **not** use Python's `__init__.py` package initialization system. 

**Changes Required**:

1. **Do not create `__init__.py` or `__init__.codon`**
   - Codon uses a simpler module system
   - Import directly from individual `.py` or `.codon` files
   - No package initialization file needed

2. **Import pattern changes**:
   
   **Python package style** (won't work in Codon):
   ```python
   import phylo  # imports from phylo/__init__.py
   phylo.Tree, phylo.upgma, phylo.neighbor_joining
   ```
   
   **Codon style** (direct imports):
   ```python
   from tree import Tree, TreeNode
   from upgma import upgma
   from nj import neighbor_joining
   ```

3. **Implications for project structure**:
   - Don't create nested package hierarchies with `__init__.py`
   - Keep modules flat or use explicit imports from each file
   - The `__all__` exports from `__init__.py` become irrelevant

### 2.5 `test_phylo.py` → `test_phylo.codon` or `test_phylo.py`

**Complexity**: High - pytest is not natively supported in Codon

**Required Changes**:

1. **Remove pytest dependency**
   - **Option A**: Rewrite tests without pytest (convert to regular functions)
   - **Option B**: Run tests using Python interop
   - **Option C**: Create a minimal test runner in Codon
   - **Recommendation**: Option A for native Codon, Option B for quick validation

2. **Convert fixtures to regular functions** (lines 13-33)
   
   **Before**:
   ```python
   @pytest.fixture
   def upgma_newick():
       with open(join(data_dir("sequence"), "newick_upgma.txt"), "r") as file:
           newick = file.read().strip()
       return newick
   ```
   
   **After**:
   ```python
   def get_upgma_newick() -> str:
       with open(join(data_dir("sequence"), "newick_upgma.txt"), "r") as file:
           newick = file.read().strip()
       return newick
   ```

3. **Convert test functions to regular functions**
   - Remove `@pytest.fixture` decorators
   - Explicitly call setup functions instead of dependency injection
   - Replace `pytest.approx()` with manual tolerance checks

4. **Replace assertions**
   
   **Before**:
   ```python
   assert tree.get_distance(i, j) == pytest.approx(
       ref_tree.get_distance(i, j), abs=1e-3
   )
   ```
   
   **After**:
   ```python
   distance_diff = abs(tree.get_distance(i, j) - ref_tree.get_distance(i, j))
   if distance_diff > 1e-3:
       raise AssertionError(f"Distance mismatch: {distance_diff}")
   ```

5. **Remove timing fixtures** (lines 36-56)
   - Can reimplement with Codon's timing if needed
   - Or use Python's `time` module (should be supported)

6. **Update import statement** (line 8)
   
   **Before**:
   ```python
   import biotite.sequence.phylo as phylo
   ```
   
   **After**:
   ```python
   import phylo  # or however the module is structured in Codon
   # Or use direct imports
   from tree import Tree, TreeNode
   from nj import neighbor_joining
   from upgma import upgma
   ```

### 2.6 `util.py` → Keep for testing or port

**Complexity**: Low-Medium

**Options**:
1. **Keep as Python** - Use Python interop to call it
2. **Port to Codon** - Will need to handle:
   - `os.path` operations (supported in Codon)
   - `urllib` (may need Python interop)
   - `importlib` (may need Python interop)
   - `shutil` (may need Python interop)

**Recommendation**: Keep as Python and use interop, or simplify to just the `data_dir()` function which is straightforward to port.

---

## 3. New Dependencies and Definitions Needed

### 3.1 Define Missing Exception Classes

Since `InvalidFileError` is from an external package, define it locally:

```python
class InvalidFileError(Exception):
    """Raised when a file format is invalid"""
    pass
```

### 3.2 Define or Remove Copyable Base Class

**Option A** - Define locally:
```python
class Copyable:
    """Base class for objects that support copying"""
    def copy(self):
        raise NotImplementedError("Subclasses must implement copy()")
```

**Option B** - Remove inheritance and keep just the copy method implementation.

**Recommendation**: Option B (simpler)

### 3.3 Handle NetworkX Dependency

For `Tree.as_graph()` method:

**Option A** - Use Python interop:
```python
from python import networkx as nx
```

**Option B** - Implement minimal graph structure in Codon:
```python
class DiGraph:
    def __init__(self):
        self.edges: List[Tuple[Any, Any, Dict[str, float]]] = []
    
    def add_edge(self, node1, node2, **attrs):
        self.edges.append((node1, node2, attrs))
```

**Recommendation**: Option A for initial port, Option B if graph operations are critical

---

## 4. Step-by-Step Porting Process

### Phase 1: Preparation
1. ✅ Analyze codebase structure
2. ✅ Identify all Cython-specific features
3. ✅ Identify external dependencies
4. ✅ Create porting strategy document

### Phase 2: Core Infrastructure (Start Here)
1. **Create new directory structure**
   ```
   week3/codon/
   ├── src/
   │   ├── tree.py       # Converted from tree.pyx
   │   ├── nj.py         # Converted from nj.pyx
   │   ├── upgma.py      # Converted from upgma.pyx
   │   └── common.py     # Helper file with exceptions
   └── test/
       ├── test_phylo.py # Converted test file
       ├── util.py       # Simplified or use Python interop
       └── sequence/
           └── data/     # Copy data files as-is
   ```
   
   **Note**: No `__init__.py` needed - Codon doesn't use package initialization files!

2. **Create helper/utility file** (`common.py`)
   ```python
   # Define exceptions
   class TreeError(Exception):
       """An exception that occurs in context of tree topology."""
       pass
   
   class InvalidFileError(Exception):
       """Raised when a file format is invalid"""
       pass
   ```

### Phase 3: Port Core Classes (Critical Path)

1. **Port `tree.py` first** (most complex, other modules depend on it)
   - [ ] Remove Cython syntax
   - [ ] Convert TreeNode class
   - [ ] Convert Tree class
   - [ ] Convert helper functions
   - [ ] Add type annotations
   - [ ] Test basic functionality

2. **Port `upgma.py`** (simpler algorithm)
   - [ ] Remove Cython syntax
   - [ ] Convert function signature
   - [ ] Convert local variables
   - [ ] Replace memoryviews with arrays
   - [ ] Test with sample data

3. **Port `nj.py`** (similar to upgma)
   - [ ] Remove Cython syntax
   - [ ] Convert function signature
   - [ ] Convert local variables
   - [ ] Replace memoryviews with arrays
   - [ ] Test with sample data

### Phase 4: Port Tests

1. **Create test runner** (if not using pytest via Python)
   ```python
   def run_test(test_name: str, test_func):
       try:
           test_func()
           print(f"✓ {test_name} passed")
           return True
       except Exception as e:
           print(f"✗ {test_name} failed: {e}")
           return False
   ```

2. **Convert test functions**
   - [ ] Convert fixtures to regular functions
   - [ ] Update imports
   - [ ] Replace pytest-specific assertions
   - [ ] Add manual test runner

3. **Run and debug tests**
   - [ ] Test with small datasets first
   - [ ] Verify numerical accuracy
   - [ ] Check edge cases

### Phase 5: Optimization and Polish

1. **Performance testing**
   - [ ] Compare execution time with Python version
   - [ ] Identify bottlenecks
   - [ ] Add parallel processing if beneficial

2. **Code cleanup**
   - [ ] Remove unused imports
   - [ ] Add docstrings where needed
   - [ ] Ensure consistent style

3. **Documentation**
   - [ ] Update README with Codon-specific instructions
   - [ ] Document any behavioral differences
   - [ ] Add build/run instructions

---

## 5. Potential Issues and Solutions

### 5.1 NumPy Compatibility

**Issue**: Not all NumPy functions may be implemented in Codon-NumPy.

**Solution**:
- Check Codon-NumPy documentation for supported functions
- For unsupported functions, use Python interop:
  ```python
  from python import numpy as pynp
  result = pynp.some_unsupported_function(data)
  ```

### 5.2 String Handling

**Issue**: Codon uses ASCII strings, not Unicode.

**Solution**:
- Test all string operations with ASCII data
- If Unicode is needed, use Python interop for those specific operations

### 5.3 Type Inference Issues

**Issue**: Codon may struggle with complex type inference in some cases.

**Solution**:
- Add explicit type annotations where compiler complains
- Use type comments if needed: `x = some_func()  # type: SomeType`

### 5.4 Optional/None Handling

**Issue**: Codon's Optional handling may differ from Python.

**Solution**:
- Use `Optional[T]` explicitly in type hints
- Always check `is None` before accessing optional values
- Initialize optional class attributes properly

### 5.5 Recursion Depth

**Issue**: Tree operations use recursion which may have limits.

**Solution**:
- Test with deep trees
- Consider iterative alternatives if recursion depth is a problem
- Codon should handle typical tree depths well

### 5.6 NetworkX Dependency

**Issue**: `as_graph()` method requires NetworkX.

**Solution**:
- Use Python interop: `from python import networkx as nx`
- Or implement minimal graph structure in Codon
- Or mark method as requiring Python interop in docs

### 5.7 File I/O

**Issue**: `from_newick()` reads Newick format strings.

**Solution**:
- Codon supports file operations natively
- Test string parsing operations carefully
- ASCII string handling should be sufficient for Newick format

### 5.8 Module Organization and Imports

**Issue**: Python's package system with `__init__.py` doesn't exist in Codon.

**Solution**:
- Do NOT create `__init__.py` or `__init__.codon` files
- Import directly from individual module files:
  ```python
  from tree import Tree, TreeNode
  from upgma import upgma
  from nj import neighbor_joining
  ```
- Keep module structure flat or ensure all files are in the same directory
- If you need to import from subdirectories, Codon may support relative paths to files, but avoid complex package hierarchies

---

## 6. Testing Strategy

### 6.1 Unit Testing Approach

Since pytest is not natively supported:

**Option 1: Simple Test Runner** (Recommended for initial port)
```python
def test_suite():
    passed = 0
    failed = 0
    
    tests = [
        ("distances", test_distances),
        ("upgma", test_upgma),
        ("neighbor_joining", test_neighbor_joining),
    ]
    
    for name, func in tests:
        try:
            func()
            print(f"✓ {name}")
            passed += 1
        except Exception as e:
            print(f"✗ {name}: {e}")
            failed += 1
    
    print(f"\n{passed} passed, {failed} failed")
    return failed == 0

if __name__ == "__main__":
    success = test_suite()
    exit(0 if success else 1)
```

**Option 2: Use Python pytest with Codon extension module**
- Compile Codon code as Python extension
- Call from pytest as before
- More complex but preserves test infrastructure

### 6.2 Test Data

- Copy test data files as-is:
  - `distances.txt`
  - `newick_upgma.txt`
- No changes needed to data files

### 6.3 Validation Approach

1. **Compare outputs with Python version**
   - Run same inputs through both implementations
   - Verify tree structures are identical
   - Check numerical values within tolerance

2. **Test edge cases**
   - Empty inputs
   - Single node trees
   - Large trees
   - Symmetric matrices
   - Various distance values

3. **Performance benchmarking**
   - Measure execution time
   - Compare with Python/Cython version
   - Should be faster than pure Python, competitive with Cython

---

## 7. Build and Deployment

### 7.1 Build Process

**Compile individual files:**
```bash
codon build --release tree.py -o tree.o
codon build --release nj.py -o nj.o
codon build --release upgma.py -o upgma.o
```

**Run tests:**
```bash
codon run test/test_phylo.py
```

**Create Python extension (optional):**
```bash
codon build --relocation-model=pic --plugin=python tree.py -o phylo.so
```

### 7.2 Integration Options

**Option A: Standalone Codon Application**
- Compile to native binary
- Run independently of Python
- Fastest execution

**Option B: Python Extension Module**
- Compile as `.so` that can be imported from Python
- Allows gradual migration
- Can use with existing Python tests

**Option C: Hybrid Approach**
- Core algorithms in Codon
- Tests and utilities in Python
- Best of both worlds

---

## 8. Expected Challenges by Difficulty

### Easy (Mechanical conversion)
- Removing `cimport` statements
- Removing `ctypedef` declarations
- Removing `@cython` decorators
- Changing file extensions
- Simple type annotation additions

### Medium (Requires careful attention)
- Converting memoryviews to NumPy arrays
- Converting `cdef` functions to regular functions
- Converting `cdef class` to regular class
- Updating all type annotations consistently
- Converting `bint` to `bool`
- Removing dependency on `InvalidFileError` and `Copyable`

### Hard (May require restructuring)
- Converting pytest fixtures to regular functions
- Handling NetworkX dependency in `as_graph()`
- Ensuring NumPy operations are all supported
- Converting recursive operations if depth is an issue
- Debugging type inference issues
- Handling `Optional` types correctly throughout

### Very Hard (May require significant work)
- Performance optimization if Codon version is slower
- Handling any NumPy functions not supported by Codon-NumPy
- Complex type relationships in Tree/TreeNode classes
- Ensuring exact numerical compatibility with original

---

## 9. Recommended Workflow

### Week 1: Foundation
- [ ] Set up Codon development environment
- [ ] Create project structure
- [ ] Port utility functions and exceptions
- [ ] Create basic test infrastructure

### Week 2: Core Implementation  
- [ ] Port TreeNode class completely
- [ ] Port Tree class completely
- [ ] Test basic tree operations
- [ ] Debug any type issues

### Week 3: Algorithms
- [ ] Port UPGMA algorithm
- [ ] Port Neighbor-Joining algorithm
- [ ] Test both algorithms with simple cases
- [ ] Verify numerical accuracy

### Week 4: Testing and Validation
- [ ] Convert all tests
- [ ] Run comprehensive test suite
- [ ] Compare outputs with Python version
- [ ] Performance benchmarking

### Week 5: Polish and Documentation
- [ ] Code cleanup
- [ ] Performance optimization
- [ ] Complete documentation
- [ ] Final validation

---

## 10. Success Criteria

The port is considered successful when:

1. ✅ All Cython syntax has been removed
2. ✅ Code compiles without errors in Codon
3. ✅ All tests pass with correct outputs
4. ✅ Numerical results match Python version (within tolerance)
5. ✅ Performance is equal to or better than Python version
6. ✅ Code is maintainable and well-documented
7. ✅ File I/O operations work correctly
8. ✅ Tree structure operations are correct
9. ✅ Both UPGMA and Neighbor-Joining produce correct trees
10. ✅ Edge cases are handled properly

---

## 11. Quick Reference: Common Pattern Conversions

### Cython → Codon Cheat Sheet

```python
# Variable declarations
cdef int x = 5                    →  x: int = 5
cdef float32 y                    →  y: float
cdef bint flag = True             →  flag: bool = True

# Function definitions
cdef int func(int x):             →  def func(x: int) -> int:
cdef float32[:] get_array():      →  def get_array() -> np.ndarray[float]:

# Class definitions
cdef class MyClass:               →  class MyClass:
    cdef int attr                 →      attr: int
    def __cinit__(self, ...):     →      def __init__(self, ...):

# Array declarations
cdef float32[:] arr               →  arr: np.ndarray[float]
cdef float32[:,:] matrix          →  matrix: np.ndarray[float]

# Typed parameters
def func(np.ndarray[float32] x):  →  def func(x: np.ndarray[float]) -> ReturnType:
def func(TreeNode node not None): →  def func(node: TreeNode) -> ReturnType:

# Imports
cimport numpy as np               →  import numpy as np
cimport cython                    →  # Remove
from .module import X             →  from module import X

# Decorators
@cython.boundscheck(False)        →  # Remove (automatic in Codon)
@cython.wraparound(False)         →  # Remove (automatic in Codon)
```

---

## 12. Additional Resources

### Codon Documentation
- **Main docs**: Check Codon documentation for supported features
- **NumPy compatibility**: Review Codon-NumPy docs for available functions
- **Python interop**: Reference guide for when you need Python features
- **Type system**: Understanding Codon's type inference and annotations

### Testing Approach
- Start with smallest components (TreeNode)
- Build up to larger structures (Tree)
- Test algorithms individually (UPGMA, Neighbor-Joining)
- Integration testing last

### Performance Expectations
- Should be **significantly faster** than pure Python
- May be comparable to or faster than Cython
- NumPy operations should be well-optimized
- Tree traversal operations will benefit from compilation

---

## Conclusion

Porting this phylogenetic tree codebase from Python/Cython to Codon is **feasible but non-trivial**. The main challenges are:

1. **Removing Cython-specific syntax** (mechanical but extensive)
2. **Converting class definitions** (requires careful type management)
3. **Handling external dependencies** (NetworkX, pytest, biotite modules)
4. **Converting tests** (pytest is not native to Codon)

**Estimated effort**: 40-60 hours for complete port with testing

**Risk level**: Medium - Most patterns have clear conversions, but type system and dependencies may require creative solutions

**Recommended approach**: Incremental port starting with core classes, extensive testing at each stage, use Python interop for complex dependencies initially, then optimize.

The resulting Codon code should be **significantly faster** than pure Python while maintaining correctness and readability.

