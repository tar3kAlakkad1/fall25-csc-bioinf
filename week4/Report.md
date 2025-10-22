### Week 4 — Edit Distance and DP Alignments

This deliverable implements four alignment algorithms in both Python and Codon, together with a reproducible evaluation script that prints the required runtime table.

### Implemented algorithms and files

- **Global alignment (Needleman–Wunsch)**
- **Local alignment (Smith–Waterman)**
- **Semi-global alignment (fitting)**
- **Affine gap penalty global alignment**

Python implementation:

- `week4/python/alignment.py` — algorithm implementations (NumPy for DP matrices)
- `week4/python/utils.py` — minimal FASTA parser
- `week4/python/main.py` — command-line wrapper for running a single alignment

Codon implementation:

- `week4/codon/alignment.py` — algorithm implementations (loop-based DP; Codon-friendly types)
- `week4/codon/main.py` — CLI for running/time-measuring alignments; compiled to `week4/codon/align`

Evaluation and helper:

- `week4/evaluate.sh` — builds Codon binary if needed and prints the runtime table
- `week4/test.py` — cross-language timing harness used by `evaluate.sh`

Data used for tests (bundled in repo):

- `week4/python/data/MT-human.fa`
- `week4/python/data/MT-orang.fa`
- `week4/python/data/q1.fa` (contains `q1..q5` records)
- `week4/python/data/t1.fa` (contains `t1..t5` records)

### Scoring schemes

- **Match**: +3
- **Mismatch**: -3
- **Gap (linear)**: -2
- **Affine gaps**: gap-open = -5, gap-extend = -1

### Algorithm notes

- **Global (NW)**: DP initialized with cumulative gap penalties on first row/col; returns `dp[m,n]`.
- **Local (SW)**: DP floored at 0; returns the maximum over all cells.
- **Semi-global (fitting)**: No penalties for leading gaps (first row/col initialized to 0). Objective is the best score that reaches the end of either sequence; we return the max over the last row or last column.
- **Affine global**: Uses three matrices: `M` (match/mismatch), `Ix` (gap in seq2; vertical), `Iy` (gap in seq1; horizontal). Transitions implement open vs extend with the specified penalties. In Python we use large negative sentinels for integer matrices; in Codon we use `-inf` with float matrices for clean initialization.

### How to run locally

Prerequisites:

- Python 3.11+ with NumPy
- Codon compiler on PATH (`codon --version`)
- macOS note: `libomp` may be required at runtime; the script tries to add it automatically when Homebrew is available.

Run the full evaluation (prints the required table):

```bash
bash week4/evaluate.sh
```

Example: run one alignment manually (Python):

```bash
python3 week4/python/main.py run global \
  week4/python/data/MT-human.fa MT_human \
  week4/python/data/MT-orang.fa MT_orang
```

Example: run/time one alignment (Codon binary):

```bash
# Build is handled by evaluate.sh, but you can build manually:
codon build -release -o week4/codon/align week4/codon/main.py

# Just run (no timing printed)
week4/codon/align run global \
  week4/python/data/MT-human.fa MT_human \
  week4/python/data/MT-orang.fa MT_orang

# Print runtime in ms for the same call
week4/codon/align time global \
  week4/python/data/MT-human.fa MT_human \
  week4/python/data/MT-orang.fa MT_orang
```

### What the evaluation prints

The output follows the requested format and includes both languages for each test:

```
Method          Language    Runtime
--------------------------------------
global-mt_human python      <ms>
global-mt_human codon       <ms>
... (also q1..q5 vs t1..t5 for global, local, semi-global, affine)
```

### CI setup

- We use the provided workflow (`.github/workflows/actions.yml`). CI calls `week4/evaluate.sh`, which:
  - Ensures NumPy is available for Python.
  - Builds the Codon binary if sources changed.
  - On macOS runners, attempts to install `libomp` via Homebrew and sets `DYLD_LIBRARY_PATH` if needed.

### Gotchas and notes

- **Local vs semi-global**: Local alignment floors at 0 (can start/end anywhere). Semi-global (fitting) removes penalties for starting gaps and scores the best alignment that reaches the end of one sequence; make sure to take the max over the last row/column.
- **Affine initialization**: For correctness, `M[0,0]=0` and the other two matrices start at negative infinity (or a large negative sentinel for integer DP) except along the first row/column where leading gaps are allowed according to affine costs.
- **Performance considerations**: Python uses NumPy-backed matrices for speed; Codon uses explicit loops and typed arrays. For large inputs, memory is `O(mn)`. No traceback is computed, only scores, to keep runtime/memory reasonable.
- **macOS OpenMP**: If Codon/Numpy fails to load `libomp.dylib`, run `brew install libomp` and re-run the evaluation; the script attempts this automatically when Homebrew exists.

### Time spent

- Approximately 6–8 hours (design, implementation, validation, and CI troubleshooting).
