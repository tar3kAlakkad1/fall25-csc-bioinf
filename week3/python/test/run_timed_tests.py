#!/usr/bin/env python3
"""
Wrapper script to run Python tests and output timing in milliseconds.
"""
import time
import sys
import os

# Add the parent directory to the path to import test modules
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import biotite.sequence.phylo as phylo
from util import data_dir
from os.path import join

def load_upgma_newick():
    """Load the UPGMA newick reference data."""
    with open(join(data_dir("sequence"), "newick_upgma.txt"), "r") as file:
        newick = file.read().strip()
    return newick

def load_distances():
    """Load the distance matrix."""
    return np.loadtxt(join(data_dir("sequence"), "distances.txt"), dtype=int)

def test_distances():
    """Test that all leaves have equal distance to root (UPGMA property)."""
    distances = load_distances()
    tree = phylo.upgma(distances)
    
    dist = tree.root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        assert leaf.distance_to(tree.root) == dist
    # Example topological distances
    assert tree.get_distance(0, 19, True) == 9
    assert tree.get_distance(4, 2, True) == 10

def test_upgma():
    """Compare the results of upgma() with DendroUPGMA."""
    distances = load_distances()
    tree = phylo.upgma(distances)
    upgma_newick = load_upgma_newick()
    
    ref_tree = phylo.Tree.from_newick(upgma_newick)
    # Cannot apply direct tree equality assertion because the distance
    # might not be exactly equal due to floating point rounding errors
    for i in range(len(tree)):
        for j in range(len(tree)):
            # Check for equal distances and equal topologies
            assert abs(tree.get_distance(i, j) - ref_tree.get_distance(i, j)) < 1e-3
            assert tree.get_distance(i, j, topological=True) == ref_tree.get_distance(
                i, j, topological=True
            )

def test_neighbor_joining():
    """Compare the results of neighbor_join() with a known tree."""
    dist = np.array([
        [ 0,  5,  4,  7,  6,  8],
        [ 5,  0,  7, 10,  9, 11],
        [ 4,  7,  0,  7,  6,  8],
        [ 7, 10,  7,  0,  5,  9],
        [ 6,  9,  6,  5,  0,  8],
        [ 8, 11,  8,  9,  8,  0],
    ])  # fmt: skip

    ref_tree = phylo.Tree(
        phylo.TreeNode(
            [
                phylo.TreeNode(
                    [
                        phylo.TreeNode(
                            [
                                phylo.TreeNode(index=0),
                                phylo.TreeNode(index=1),
                            ],
                            [1, 4],
                        ),
                        phylo.TreeNode(index=2),
                    ],
                    [1, 2],
                ),
                phylo.TreeNode(
                    [
                        phylo.TreeNode(index=3),
                        phylo.TreeNode(index=4),
                    ],
                    [3, 2],
                ),
                phylo.TreeNode(index=5),
            ],
            [1, 1, 5],
        )
    )

    test_tree = phylo.neighbor_joining(dist)

    assert test_tree == ref_tree

def format_time(elapsed_seconds):
    """Format time in milliseconds, showing decimals for very fast tests."""
    elapsed_ms = elapsed_seconds * 1000
    rounded_ms = round(elapsed_ms)
    if rounded_ms == 0:
        # Show precise time with 2 decimal places for very fast tests
        return f"{elapsed_ms:.2f}ms"
    else:
        return f"{rounded_ms}ms"

if __name__ == "__main__":
    # Use Python's perf_counter for high-resolution timing
    overall_start = time.perf_counter()
    
    try:
        # Test distances
        start_time = time.perf_counter()
        test_distances()
        print(f"  test_distances: {format_time(time.perf_counter() - start_time)}")
        
        # Test UPGMA
        start_time = time.perf_counter()
        test_upgma()
        print(f"  test_upgma: {format_time(time.perf_counter() - start_time)}")
        
        # Test neighbor joining
        start_time = time.perf_counter()
        test_neighbor_joining()
        print(f"  test_neighbor_joining: {format_time(time.perf_counter() - start_time)}")
        
        # Total time
        print(f"python      {format_time(time.perf_counter() - overall_start)}")
        sys.exit(0)
    except Exception as e:
        end_time = time.perf_counter()
        total_elapsed_ms = round((end_time - overall_start) * 1000)
        
        print(f"python      {total_elapsed_ms}ms", file=sys.stderr)
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

