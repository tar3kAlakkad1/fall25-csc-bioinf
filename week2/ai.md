# TRviz to Codon Porting Guide

## Assignment Overview
Port the Python package **trviz** (Tandem Repeat Visualizer) to Codon, excluding the visualization components. TRviz is a bioinformatics library for decomposing, encoding, aligning, and visualizing tandem repeat DNA sequences.

**Repository:** https://github.com/Jong-hun-Park/trviz  
**Paper:** https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad058/7143378

---

## Project Structure Understanding

### TRviz Modules (4 core modules)
1. **Decomposer** - Decomposes TR sequences into constituent motifs using dynamic programming
2. **Encoder** - Encodes motifs into symbols (max 90 distinct symbols)
3. **Aligner** - Aligns encoded sequences using MAFFT with custom scoring matrices
4. **Visualizer** - Generates plots (**SKIP THIS - leave in Python**)

### Additional Components
- **utils.py** - Helper functions (e.g., reading FASTA files)
- **main.py** - TandemRepeatVizWorker class that orchestrates the pipeline

---

## Step-by-Step Implementation Plan

### Phase 1: Setup and Repository Analysis (Est. 2-3 hours)

#### Step 1.1: Clone and Explore the Repository
```bash
cd week2/
git clone https://github.com/Jong-hun-Park/trviz.git
cd trviz/
```

**Action Items:**
- Examine the directory structure
- Identify all Python files in the `trviz/` package directory
- Look at `tests/` directory to understand existing test structure
- Review `setup.py` or `pyproject.toml` for dependencies
- Read through the example notebook if available

**Use Cursor AI:** Ask Cursor to "Summarize the structure of the trviz package and list all classes and functions"

#### Step 1.2: Understand Core Algorithms
Read through each module to understand:
- **Decomposer**: Dynamic programming algorithm (similar to Needleman-Wunsch)
- **Encoder**: Motif-to-symbol mapping with frequency thresholds
- **Aligner**: MAFFT wrapper with custom scoring matrix generation
- **Utils**: File I/O operations

**Use Cursor AI:** For complex algorithms, ask "Explain how the decompose method works in decomposer.py"

#### Step 1.3: Document Dependencies
Create a list of what can/cannot be ported:
- ✅ Pure Python logic (decomposition algorithm, encoding)
- ✅ Basic string operations
- ✅ File I/O
- ⚠️ BioPython usage (use `from python import Bio` in Codon)
- ⚠️ subprocess calls to MAFFT (keep in Python)
- ❌ Matplotlib/plotting (skip completely)
- ❌ Complex parsing libraries (use `from python import`)

---

### Phase 2: Create Codon Port Structure (Est. 1-2 hours)

#### Step 2.1: Set Up Directory Structure
```bash
# In week2/
mkdir -p trviz_codon/
cd trviz_codon/
```

Create these files:
```
week2/
├── trviz_codon/
│   ├── __init__.py           # Codon package init
│   ├── decomposer.codon      # Core decomposition logic
│   ├── encoder.codon         # Encoding logic
│   ├── aligner.codon         # Alignment logic (may use python import)
│   ├── utils.codon           # Utility functions
│   └── main.codon            # Main worker class
├── test.py                   # Unified test file (Codon + Python)
├── report.md                 # Your assignment report
└── ai.md                     # This file
```

#### Step 2.2: Create test.py Template
```python
# test.py - Works with both Python and Codon
from python import unittest  # Only needed for imports, not for tests

# Import your modules
try:
    # Try Codon imports first
    from trviz_codon.decomposer import Decomposer
    from trviz_codon.encoder import Encoder
    # ... other imports
except:
    # Fall back to Python imports
    from trviz.decomposer import Decomposer
    from trviz.encoder import Encoder
    # ... other imports

# Test decorator for Codon
def test(func):
    """Decorator to mark test functions"""
    func._is_test = True
    return func

# ============= DECOMPOSER TESTS =============
@test
def test_decompose_basic():
    """Test basic sequence decomposition"""
    decomposer = Decomposer()
    sequence = "ACCTTGACCTTGACCTTGACCTTG"
    motifs = ["ACCTTG"]
    result = decomposer.decompose(sequence, motifs)
    expected = ["ACCTTG", "ACCTTG", "ACCTTG", "ACCTTG"]
    assert result == expected, f"Expected {expected}, got {result}"
    print("✓ test_decompose_basic passed")

@test
def test_decompose_multiple_motifs():
    """Test decomposition with multiple motifs"""
    # Add your test here
    pass

# ============= ENCODER TESTS =============
@test
def test_encode_motifs():
    """Test motif encoding"""
    # Add your test here
    pass

# ============= ALIGNER TESTS =============
# @test  # Comment out if not implementing alignment in Codon
# def test_align_sequences():
#     """Test sequence alignment"""
#     pass  # Skip if using python import for MAFFT

# ============= UTILS TESTS =============
@test
def test_read_fasta():
    """Test FASTA file reading"""
    # Add your test here
    pass

# Run all tests
def run_tests():
    """Execute all test functions"""
    import sys
    test_funcs = [obj for name, obj in globals().items() 
                  if callable(obj) and hasattr(obj, '_is_test')]
    
    passed = 0
    failed = 0
    
    for test_func in test_funcs:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"✗ {test_func.__name__} failed: {e}")
            failed += 1
    
    print(f"\n{'='*50}")
    print(f"Tests: {passed} passed, {failed} failed")
    print(f"{'='*50}")
    
    sys.exit(0 if failed == 0 else 1)

if __name__ == "__main__":
    run_tests()
```

---

### Phase 3: Port Core Modules (Est. 8-12 hours)

#### Step 3.1: Port Decomposer Module (Highest Priority)

**Use Codon MCP Server:** When you encounter Codon-specific syntax issues, use your MCP server to ask:
- "How do I declare types in Codon?"
- "What's the Codon equivalent of Python's defaultdict?"
- "How do I handle optional parameters in Codon?"

**Porting Strategy:**
1. Start with the core `decompose()` method
2. Implement the dynamic programming algorithm
3. Handle string operations (Codon has good string support)
4. Test incrementally

**Common Gotchas:**
- Codon requires type annotations
- List comprehensions may need explicit types
- Use `Optional[T]` for optional parameters
- Codon doesn't have all Python built-ins

**Example Conversion:**
```python
# Python
def decompose(self, sequence, motifs, gap_penalty=-1):
    result = []
    # ... algorithm
    return result

# Codon
def decompose(sequence: str, motifs: List[str], gap_penalty: int = -1) -> List[str]:
    result: List[str] = []
    # ... algorithm
    return result
```

**Use Cursor:** 
- Select the Python code and ask: "Convert this to Codon with type annotations"
- For complex sections: "Explain what this code does and suggest Codon-compatible alternatives"

#### Step 3.2: Port Encoder Module

**Key Components:**
- Motif frequency calculation
- Symbol assignment (90 max symbols)
- Motif-to-symbol mapping

**Considerations:**
- May need to use dictionaries (supported in Codon)
- String operations should work directly
- Counter functionality might need manual implementation or `from python import collections`

**Use MCP Server:** Ask about:
- "How to implement Counter in Codon?"
- "Best practices for dictionary operations in Codon"

#### Step 3.3: Port Aligner Module (Partial - MAFFT in Python)

**Strategy:**
- Port the scoring matrix generation to Codon
- Keep MAFFT system calls in Python using `from python import subprocess`
- Port sequence preparation and post-processing if possible

**Example:**
```python
# In aligner.codon
from python import subprocess

def run_mafft(sequences: List[str]) -> List[str]:
    """Wrapper for MAFFT - delegates to Python"""
    # Use Python's subprocess
    pass
```

**Use MCP Server:** Ask "How to call Python functions from Codon?"

#### Step 3.4: Port Utils Module

**Components to Port:**
- FASTA file parsing (may need BioPython)
- Basic file I/O
- Data structure utilities

**Strategy:**
```python
# For BioPython functionality
from python import Bio
from python.Bio import SeqIO

def get_sample_and_sequence_from_fasta(filepath: str):
    # Use Python's Bio.SeqIO
    pass
```

---

### Phase 4: Testing (Est. 3-4 hours)

#### Step 4.1: Write Comprehensive Tests

For each module, test:
1. **Basic functionality** - Simple inputs, expected outputs
2. **Edge cases** - Empty inputs, single element, very long sequences
3. **Multiple motifs** - Complex decomposition scenarios
4. **Error handling** - Invalid inputs

**Test Coverage Goals:**
- Decomposer: 70%+ of methods
- Encoder: 60%+ of methods
- Aligner: 40%+ (if partially ported)
- Utils: 50%+ of methods

#### Step 4.2: Run Tests in Both Environments

```bash
# Test with Python
python test.py

# Test with Codon
codon test.py
```

**Ensure same results!** Both should pass the same tests.

#### Step 4.3: Comment Out Unimplemented Tests

```python
# @test  # Commented because visualization not ported
# def test_generate_plot():
#     pass
```

---

### Phase 5: CI/CD Setup (Est. 1-2 hours)

#### Step 5.1: Update GitHub Actions Workflow

Modify `.github/workflows/actions.yml`:

```yaml
name: Github CI
on: [push]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - run: echo "🎉 The job was automatically triggered by a ${{ github.event_name }} event."
      
      - name: Install Codon
        run: |
          mkdir -p ${HOME}/.codon
          curl -L https://github.com/exaloop/codon/releases/download/v0.19.3/codon-linux-x86_64.tar.gz | tar zxvf - --strip-components=1 -C ${HOME}/.codon
          curl -L https://github.com/exaloop/seq/releases/download/v0.11.5/seq-linux-x86_64.tar.gz | tar zxvf - -C ${HOME}/.codon/lib/codon/plugins
          export PATH=${PATH}:${HOME}/.codon/bin
          echo "PATH=${PATH}" >> $GITHUB_ENV
      
      - name: Check out repository code
        uses: actions/checkout@v5
      
      - uses: actions/setup-python@v6
        with:
          python-version: '3.13'
      
      - name: Set up Codon Python bridge
        run: |
          pip install find_libpython
          export CODON_PYTHON=$(find_libpython)
          echo "Found Python at: ${CODON_PYTHON}"
          echo "CODON_PYTHON=${CODON_PYTHON}" >> $GITHUB_ENV
      
      - name: Install BioPython and dependencies for Week 2
        run: |
          pip install biopython
          # Install MAFFT if tests require it (optional)
          # sudo apt-get install -y mafft
      
      - name: Week 1
        run: bash week1/evaluate.sh
      
      - name: Week 2 - Python tests
        run: |
          cd week2
          python test.py
      
      - name: Week 2 - Codon tests
        run: |
          cd week2
          codon test.py
```

**Key Changes:**
- Added BioPython installation
- Added Week 2 test stages
- Consider adding MAFFT if alignment tests need it

---

### Phase 6: Documentation (Est. 1-2 hours)

#### Step 6.1: Create report.md

```markdown
# Week 2: TRviz to Codon Port - Report

## Summary
Successfully ported X% of trviz functionality to Codon, excluding visualization.

## Modules Ported

### Decomposer (100% ported)
- Core decomposition algorithm implemented
- Dynamic programming approach maintained
- All tests passing

### Encoder (95% ported)
- Motif encoding logic ported
- Used `from python import` for Counter
- Minor performance optimization applied

### Aligner (40% ported)
- Scoring matrix generation ported
- MAFFT calls kept in Python (external dependency)
- Sequence pre/post-processing ported

### Utils (80% ported)
- FASTA parsing uses BioPython via Python import
- File I/O operations ported

## Gotchas and Challenges

1. **Type Annotations Required**
   - Solution: Added explicit types to all functions
   - Time lost: ~1 hour debugging type errors

2. **No Built-in Counter**
   - Solution: Used `from python import collections`
   - Alternative: Could implement custom Counter

3. **MAFFT External Dependency**
   - Solution: Kept subprocess calls in Python
   - Rationale: External tool, not worth porting wrapper

4. **String Operations**
   - Most worked directly in Codon
   - Some advanced slicing needed adjustment

5. **List Comprehensions**
   - Required explicit type hints in some cases
   - Example: `[str(x) for x in items]` → needed type annotation

## Test Coverage
- Total methods in original: XX
- Methods ported: YY (ZZ%)
- Tests implemented: AA
- All tests passing: Yes/No

## Time Breakdown
- Repository analysis: X hours
- Decomposer porting: X hours
- Encoder porting: X hours
- Aligner porting: X hours
- Utils porting: X hours
- Testing: X hours
- CI setup: X hours
- Documentation: X hours
**Total: XX hours**

## Future Improvements
1. Implement custom Counter for better performance
2. Optimize dynamic programming algorithm
3. Add more edge case tests
4. Consider porting visualization to web-based solution
```

#### Step 6.2: Update ai.md (this file)
Add actual time spent and any deviations from the plan.

---

## Workflow Tips for Cursor + Codon MCP

### Using Cursor Effectively

1. **Code Explanation:**
   - Select Python code → Ask: "Explain this function's algorithm"
   - Helps understand before porting

2. **Codon Conversion:**
   - Select Python function → Ask: "Convert to Codon with type annotations"
   - Review and adjust the output

3. **Debugging:**
   - Select error-prone code → Ask: "What Codon syntax errors might be here?"

4. **Optimization:**
   - Select Codon code → Ask: "Can this be optimized for Codon?"

### Using Codon MCP Server

**When to use MCP:**
- Syntax questions: "What's the Codon syntax for X?"
- Best practices: "Best way to handle optional parameters in Codon?"
- Type system: "How do generics work in Codon?"
- Imports: "How to use Python libraries from Codon?"

**Common MCP Queries:**
```
1. "Show me how to define a generic function in Codon"
2. "What's the equivalent of Python's defaultdict in Codon?"
3. "How to handle None values in Codon?"
4. "Best practices for error handling in Codon"
5. "How to benchmark Codon code performance?"
```

### Debugging Strategy

1. **Incremental Testing:**
   ```bash
   # Test each function as you port it
   codon run -plugin seq -o test_decomp test_decomposer.py
   ```

2. **Type Errors:**
   - Most common issue in Codon
   - Use MCP to clarify type system
   - Add print statements with types: `print(type(variable))`

3. **Performance Comparison:**
   ```python
   import time
   start = time.time()
   # Your code
   print(f"Time: {time.time() - start}")
   ```

---

## Grading Rubric Alignment

### Code Coverage (Main Grade Factor)
- **Target:** Port 70%+ of non-visualization methods
- **Strategy:** Focus on Decomposer (most important) → Encoder → Utils → Aligner
- **Track:** Keep count of methods ported vs. total

### Tests Required
- ✅ Single test.py file
- ✅ Uses @test decorator (not unittest)
- ✅ Same results for `python test.py` and `codon test.py`
- ✅ Unimplemented tests commented out

### CI Expected
- ✅ Uses provided CI template
- ✅ Installs BioPython
- ✅ Runs both Python and Codon tests
- ✅ All tests pass

### Documentation
- ✅ report.md in week2/ directory
- ✅ Describes all steps and gotchas
- ✅ Includes time estimate
- ✅ ai.md uploaded

### Understanding Required
- Be prepared to explain:
  - How decomposition algorithm works
  - Why certain parts weren't ported
  - Trade-offs between Python imports vs. native Codon

---

## Quick Reference

### Essential Commands
```bash
# Clone repository
git clone https://github.com/Jong-hun-Park/trviz.git

# Run Python tests
python test.py

# Run Codon tests
codon test.py

# Run with Codon Python bridge
export CODON_PYTHON=$(find_libpython)
codon run -plugin seq test.py

# Check Codon compilation
codon build -plugin seq your_file.codon
```

### File Organization Checklist
- [ ] week2/trviz_codon/ directory created
- [ ] All .codon files created
- [ ] test.py with @test decorators
- [ ] report.md completed
- [ ] ai.md uploaded
- [ ] CI workflow updated
- [ ] Git commits with clear messages

### Pre-Submission Checklist
- [ ] Both Python and Codon tests pass
- [ ] CI passes on GitHub
- [ ] report.md includes time estimate
- [ ] All unimplemented tests commented out
- [ ] Code coverage documented
- [ ] No syntax errors in Codon files
- [ ] BioPython installed in CI
- [ ] No breaking changes to CI for week1

---

## Estimated Time Investment

| Phase | Estimated Time |
|-------|----------------|
| Setup & Analysis | 2-3 hours |
| Structure Creation | 1-2 hours |
| Module Porting | 8-12 hours |
| Testing | 3-4 hours |
| CI Setup | 1-2 hours |
| Documentation | 1-2 hours |
| **Total** | **16-25 hours** |

*Actual time may vary based on Python/Codon experience and complexity of the original code.*

---

## Additional Resources

- **TRviz Documentation:** https://trviz.readthedocs.io/
- **Codon Documentation:** https://docs.exaloop.io/codon/
- **Codon Examples:** https://github.com/exaloop/codon/tree/develop/test
- **BioPython Docs:** https://biopython.org/wiki/Documentation
- **MAFFT Manual:** https://mafft.cbrc.jp/alignment/software/manual/manual.html

---

*This guide is designed to be used with Cursor AI and a Codon MCP server for optimal efficiency. Update this file as you progress through the assignment.*