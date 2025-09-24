Codon versions of Week 1 scripts

This folder contains Codon equivalents of the Python scripts in `../py`.

Contents:
- `main.codon`: Builds a DBG and writes contigs into the specified dataset directory
- `utils.codon`: FASTA reading helpers
- `dbg.codon`: DBG implementation with integer node IDs
- `dbg_kmer_as_key.codon`: DBG implementation with k-mer strings as keys
- `compute_metrics_table.codon`: Computes simple metrics (num contigs, total bases, N50)

Datasets for convenience are under `../codon/genome-assembly/` (`data1`, `data2`, ...).

Prerequisites
- Codon installed (`codon` in your PATH). See installation guide: https://docs.exaloop.io/codon/getting-started

Run: assemble a dataset
Run from `fall25-csc-bioinf/week1/code` and pass a dataset directory.

```bash
cd fall25-csc-bioinf/week1/code
codon run @codon/main.codon codon/genome-assembly/data1
```

This writes `contig.fasta` into the passed dataset directory (e.g., `codon/genome-assembly/data1/contig.fasta`).

Notes:
- `main.codon` accepts one argument: the dataset directory (relative or absolute).
- The de Bruijn graph k is set inside `main.codon` (default `k = 25`).

Run: compute metrics table
Option A (auto-detect `data*`): run from the directory containing the datasets so auto-discovery works.

```bash
cd fall25-csc-bioinf/week1/code/codon/genome-assembly
codon run ../../@codon/compute_metrics_table.codon metrics.csv
```

Option B (explicit datasets): run from anywhere and pass dataset paths explicitly.

```bash
codon run fall25-csc-bioinf/week1/code/@codon/compute_metrics_table.codon metrics.csv \
  fall25-csc-bioinf/week1/code/codon/genome-assembly/data1 \
  fall25-csc-bioinf/week1/code/codon/genome-assembly/data2
```

Build native binaries (optional)
```bash
cd fall25-csc-bioinf/week1/code
codon build -release -o build/main @codon/main.codon
codon build -release -o build/metrics @codon/compute_metrics_table.codon
```

Implementation notes
- Where Python's stdlib is used (e.g., `os.path.join`), we import via `from python import os` and wrap joins with `str(...)`.
- `dbg.codon` and `dbg_kmer_as_key.codon` are libraries used by `main.codon`; they are not standalone entry points.


