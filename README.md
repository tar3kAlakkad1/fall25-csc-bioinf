# CSC 427 — Bioinformatics

[![CI](https://github.com/tar3kAlakkad1/fall25-csc-bioinf/actions/workflows/actions.yml/badge.svg)](https://github.com/tar3kAlakkad1/fall25-csc-bioinf/actions)

Course project for CSC 427 Bioinformatics (Fall 2025). Each weekly deliverable implements bioinformatics algorithms in Python and ports them to [Codon](https://github.com/exaloop/codon), a high-performance Python compiler, to study the trade-offs between productivity and runtime performance.

## Repository Structure

```
.
├── week1/             De Bruijn graph genome assembly
├── week2/             TRViz — Tandem Repeat Visualizer (Python → Codon)
├── week3/             Phylogenetic trees: UPGMA & Neighbor-Joining (Cython → Codon)
├── week4/             Sequence alignment: NW, SW, semi-global, affine (Python → Codon)
├── week5_and_6/       Variant calling & haplotype phasing pipeline
├── week6/             Single-cell RNA-seq analysis (alevin-fry + scanpy)
└── .github/workflows/ CI pipeline
```

Each deliverable directory contains:

| File / Dir | Purpose |
|---|---|
| `python/` | Reference Python (or Cython) implementation |
| `codon/` | Codon port |
| `Report.md` | Implementation write-up and porting notes |
| `ai.md` | AI-assisted development log |
| `evaluate.sh` | Benchmark script comparing Python vs Codon runtimes |

## Weekly Deliverables

### Week 1 — De Bruijn Graph Assembly

Builds de Bruijn graphs from k-mers and finds Eulerian paths to reconstruct genome sequences. Python-only reference implementation translated from the [Fudan macrogeomics project](https://github.com/fudan-generics/macrogeomics).

### Week 2 — TRViz (Tandem Repeat Visualizer)

Full port of the [TRViz](https://github.com/Jong-hun-Park/TRviz) library to Codon. Four modules: motif **Decomposer** (DP / HMM), **Motif Encoder**, **Motif Aligner**, and **Utils**. HMM functionality uses Codon's `@python` interop to call pomegranate.

### Week 3 — Phylogenetic Trees

Ports biotite's [phylo module](https://github.com/biotite-dev/biotite/tree/v1.4.0/src/biotite/sequence/phylo) from Cython to Codon. Implements **UPGMA** (assumes molecular clock) and **Neighbor-Joining** (no clock assumption) hierarchical clustering on distance matrices, with Newick format I/O.

### Week 4 — Sequence Alignment

Four dynamic-programming alignment algorithms, each implemented in both NumPy-backed Python and loop-based Codon:

| Algorithm | Description |
|---|---|
| **Global** (Needleman-Wunsch) | End-to-end alignment |
| **Local** (Smith-Waterman) | Best local subsequence alignment |
| **Semi-global** (fitting) | No leading/trailing gap penalties |
| **Affine gap** | Separate gap-open (-5) and gap-extend (-1) costs |

Test data includes human and orangutan mitochondrial sequences.

### Weeks 5 & 6 — Variant Calling + Single-Cell RNA-seq

**Week 5:** Align Illumina and PacBio reads to hg38 chr10, call variants with bcftools, phase with HapCUT2, compare haplotypes, and determine CYP2C19/CYP2C9/CYP2C8 star-alleles via PharmVar.

**Week 6:** End-to-end single-cell RNA-seq pipeline in a self-contained Jupyter notebook — simpleaf/alevin-fry quantification, scanpy preprocessing, and CellTypist cell-type annotation.

## Getting Started

### Prerequisites

- Python 3.11+
- [Codon compiler](https://github.com/exaloop/codon) (weeks 2–4)
- NumPy (`pip install numpy`)

macOS users also need OpenMP for Codon:

```bash
brew install libomp
export DYLD_LIBRARY_PATH="/opt/homebrew/opt/libomp/lib:$DYLD_LIBRARY_PATH"
```

### Running Evaluations

Each ported week has a benchmark script that runs both Python and Codon and prints a timing comparison:

```bash
# Week 3 — Phylogenetic trees
bash week3/evaluate.sh

# Week 4 — Sequence alignment (builds Codon binary automatically)
bash week4/evaluate.sh
```

### Running Tests Individually

```bash
# Week 2 — Python
cd week2/python && pip install -e . && pytest tests/

# Week 2 — Codon
cd week2/codon && for t in tests/*.codon; do codon run "$t"; done

# Week 3 — Python
cd week3/python/test && python3 run_timed_tests.py

# Week 3 — Codon
cd week3/codon && codon run test_phylo.py

# Week 4 — single alignment (Python)
python3 week4/python/main.py run global \
  week4/python/data/MT-human.fa MT_human \
  week4/python/data/MT-orang.fa MT_orang
```

### Week 6 Notebook

The single-cell RNA-seq notebook is fully self-contained and runs in CI. To execute locally:

```bash
cd week6
conda install -c conda-forge -c bioconda simpleaf
pip install scanpy celltypist leidenalg anndata pyroe
jupyter nbconvert --to notebook --execute week6.ipynb --inplace
```

## CI/CD

GitHub Actions runs the week 6 notebook on every push via micromamba with conda-forge and bioconda channels. The executed notebook is archived as a build artifact. See [`.github/workflows/actions.yml`](.github/workflows/actions.yml).

## Performance

Codon typically delivers **2–10x speedup** over CPython for loop-heavy bioinformatics code. Gains come from native compilation and static type inference. NumPy-dominated workloads see smaller improvements since NumPy already delegates to C.

## Key Porting Lessons

| Python / Cython | Codon |
|---|---|
| `Dict[str, Any]` mixed-type dicts | Uniform-type nested dicts only |
| `None` initialization | Must initialize with concrete type |
| `**kwargs` | Explicit keyword parameters |
| `dtype=object` arrays | Typed int/float arrays |
| `.format()`, `%` strings | f-strings only |
| pytest fixtures & parametrize | Plain functions with `assert` |
| `cdef class`, memoryviews | Standard classes, NumPy arrays |

## Technologies

**Languages:** Python, Codon, Cython (source)
**Scientific stack:** NumPy, SciPy, scanpy, anndata
**Bioinformatics tools:** minimap2, samtools, bcftools, HapCUT2, simpleaf, alevin-fry, CellTypist, IGV
**Infrastructure:** GitHub Actions, Jupyter, conda/micromamba
