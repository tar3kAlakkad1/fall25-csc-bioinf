#!/usr/bin/env python3
import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def read_fasta_lengths(fasta_path):
    """Return list of sequence lengths from a FASTA file."""
    lengths = []
    current_len = 0
    with open(fasta_path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_len > 0:
                    lengths.append(current_len)
                current_len = 0
            else:
                current_len += len(line)
        if current_len > 0:
            lengths.append(current_len)
    return lengths


def calculate_n50(lengths):
    """Compute N50 given a list of contig lengths."""
    if not lengths:
        return 0
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    half = total / 2.0
    acc = 0
    for length in lengths:
        acc += length
        if acc >= half:
            return length
    return 0


def compute_metrics_for_dataset(dataset_dir):
    dataset_dir = Path(dataset_dir)
    contigs = dataset_dir / 'contig.fasta'
    if not contigs.exists():
        raise FileNotFoundError(f"contig.fasta not found in {dataset_dir}")

    ref = dataset_dir / 'long.fasta'
    if not ref.exists():
        raise FileNotFoundError(
            f"No reference provided and {ref} not found. Provide --reference."
        )

    # N50 from contigs
    contig_lengths = read_fasta_lengths(contigs)
    n50 = calculate_n50(contig_lengths)
    num_contigs = len(contig_lengths)
    total_contig_bases = sum(contig_lengths)

    return {
        'dataset': dataset_dir.name,
        'num_contigs': num_contigs,
        'total_contig_bases': total_contig_bases,
        'N50': n50,
    }


def find_datasets(base_dir):
    base = Path(base_dir)
    candidates = []
    for name in os.listdir(base):
        if name.startswith('data'):
            d = base / name
            if (d / 'contig.fasta').exists():
                candidates.append(str(d))
    return sorted(candidates)


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Compute assembly metrics per dataset directory. By default, uses each dataset\'s long.fasta as reference.'
        )
    )
    parser.add_argument('--datasets', nargs='*', default=None, help='Dataset directories (default: auto-detect data*)')
    parser.add_argument('--out', default='metrics.csv', help='Output CSV path (default: metrics.csv)')

    args = parser.parse_args()

    if args.datasets:
        datasets = args.datasets
    else:
        datasets = find_datasets('.')
        if not datasets:
            print('No datasets found (expected directories like data1, data2 with contig.fasta).', file=sys.stderr)
            sys.exit(1)

    rows = []
    for d in datasets:
        print(f'Computing metrics for {d} ...', file=sys.stderr)
        stats = compute_metrics_for_dataset(d)
        rows.append(stats)

    # Write CSV
    fieldnames = [
        'dataset',
        'num_contigs',
        'total_contig_bases',
        'N50',
    ]
    with open(args.out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f'Wrote {args.out} with metrics for {len(rows)} dataset(s).')


if __name__ == '__main__':
    main()


