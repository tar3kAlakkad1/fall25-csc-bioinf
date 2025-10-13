# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

import numpy as np
from tree import Tree, TreeNode
from upgma import upgma
from nj import neighbor_joining
import time

def approx_equal(a: float, b: float, abs_tol: float = 1e-3):
    return abs(a - b) <= abs_tol

def upgma_newick():
    with open("newick_upgma.txt", "r") as file:
        newick = file.read().strip()
    return newick

def distances():
    # Distances are based on the example
    # "Dendrogram of the BLOSUM62 matrix"
    # with the small modification M[i,j] += i+j
    # to reduce ambiguity in the tree construction.
    # Load and convert to float64 for proper processing
    with open("distances.txt", "r") as f:
        lines = f.readlines()
    data = []
    for line in lines:
        values = [float(x) for x in line.strip().split()]
        data.extend(values)
    n = int(len(data) ** 0.5)
    return np.array(data, dtype=np.float64).reshape(n, n)

def tree(distances):
    return upgma(distances)

@test
def test_distances():
    # Tree is created via UPGMA
    # -> The distances to root should be equal for all leaf nodes
    test_distances_data = distances()
    test_tree = tree(test_distances_data)
    
    dist = test_tree.root.distance_to(test_tree.leaves[0])
    for leaf in test_tree.leaves:
        assert leaf.distance_to(test_tree.root) == dist
    # Example topological distances
    assert test_tree.get_distance(0, 19, True) == 9
    assert test_tree.get_distance(4, 2, True) == 10


@test
def test_upgma():
    """
    Compare the results of `upgma()` with DendroUPGMA.
    """
    test_distances_data = distances()
    test_tree = tree(test_distances_data)
    test_upgma_newick = upgma_newick()
    
    ref_tree = Tree.from_newick(test_upgma_newick)
    # Cannot apply direct tree equality assertion because the distance
    # might not be exactly equal due to floating point rounding errors
    for i in range(len(test_tree)):
        for j in range(len(test_tree)):
            # Check for equal distances and equal topologies
            assert approx_equal(test_tree.get_distance(i, j), ref_tree.get_distance(i, j), 1e-3)
            assert test_tree.get_distance(i, j, topological=True) == ref_tree.get_distance(
                i, j, topological=True
            )


@test
def test_neighbor_joining():
    """
    Compare the results of `neighbor_join()` with a known tree.
    """
    dist_data = [
         0.0,  5.0,  4.0,  7.0,  6.0,  8.0,
         5.0,  0.0,  7.0, 10.0,  9.0, 11.0,
         4.0,  7.0,  0.0,  7.0,  6.0,  8.0,
         7.0, 10.0,  7.0,  0.0,  5.0,  9.0,
         6.0,  9.0,  6.0,  5.0,  0.0,  8.0,
         8.0, 11.0,  8.0,  9.0,  8.0,  0.0,
    ]  # fmt: skip
    dist = np.array(dist_data, dtype=np.float64).reshape(6, 6)

    ref_tree = Tree(
        TreeNode(
            [
                TreeNode(
                    [
                        TreeNode(
                            [
                                TreeNode(index=0),
                                TreeNode(index=1),
                            ],
                            [1.0, 4.0],
                        ),
                        TreeNode(index=2),
                    ],
                    [1.0, 2.0],
                ),
                TreeNode(
                    [
                        TreeNode(index=3),
                        TreeNode(index=4),
                    ],
                    [3.0, 2.0],
                ),
                TreeNode(index=5),
            ],
            [1.0, 1.0, 5.0],
        )
    )

    test_tree = neighbor_joining(dist)

    assert test_tree == ref_tree


def format_time(elapsed_seconds: float) -> str:
    """Format time in milliseconds, showing decimals for very fast tests."""
    elapsed_ms = elapsed_seconds * 1000.0
    rounded_ms = round(elapsed_ms)
    if rounded_ms == 0:
        # Show precise time with 2 decimal places for very fast tests
        return f"{elapsed_ms:.2f}ms"
    else:
        return f"{rounded_ms}ms"


if __name__ == "__main__":
    # Use Codon's perf_counter for high-resolution timing
    overall_start = time.perf_counter()
    
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
    print(f"codon       {format_time(time.perf_counter() - overall_start)}")