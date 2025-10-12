# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

from python import numpy as pnp
import numpy.pybridge
import numpy as np
from tree import Tree, TreeNode
from upgma import upgma
from nj import neighbor_joining
from python import sys
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
    distances: np.ndarray[int,2] = pnp.loadtxt("distances.txt", dtype=pnp.int64)
    return distances

def tree(distances):
    return upgma(distances)

def session_timing():
    start_time = time.perf_counter()
    yield
    end_time = time.perf_counter()
    total_elapsed = end_time - start_time
    print(f"\n{'='*60}")
    print(f"TOTAL TIME FOR ALL TESTS: {total_elapsed:.4f} seconds")
    print(f"{'='*60}")

def timing(request, session_timing):
    """Fixture to time each individual test execution"""
    test_name = request.node.name
    start_time = time.perf_counter()
    yield
    end_time = time.perf_counter()
    elapsed = end_time - start_time
    print(f"\n  ⏱️  {test_name} took {elapsed:.4f} seconds")


def test_distances(tree):
    # Tree is created via UPGMA
    # -> The distances to root should be equal for all leaf nodes
    dist = tree.root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        assert leaf.distance_to(tree.root) == dist
    # Example topological distances
    assert tree.get_distance(0, 19, True) == 9
    assert tree.get_distance(4, 2, True) == 10


def test_upgma(tree, upgma_newick):
    """
    Compare the results of `upgma()` with DendroUPGMA.
    """
    labels = [""]
    ref_tree = Tree.from_newick(upgma_newick, labels=labels)
    # Cannot apply direct tree equality assertion because the distance
    # might not be exactly equal due to floating point rounding errors
    for i in range(len(tree)):
        for j in range(len(tree)):
            # Check for equal distances and equal topologies
            assert approx_equal(tree.get_distance(i, j), ref_tree.get_distance(i, j), 1e-3)
            assert tree.get_distance(i, j, topological=True) == ref_tree.get_distance(
                i, j, topological=True
            )


def test_neighbor_joining():
    """
    Compare the results of `neighbor_join()` with a known tree.
    """
    dist = np.array([
        [ 0,  5,  4,  7,  6,  8],
        [ 5,  0,  7, 10,  9, 11],
        [ 4,  7,  0,  7,  6,  8],
        [ 7, 10,  7,  0,  5,  9],
        [ 6,  9,  6,  5,  0,  8],
        [ 8, 11,  8,  9,  8,  0],
    ])  # fmt: skip

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
                            [1, 4],
                        ),
                        TreeNode(index=2),
                    ],
                    [1, 2],
                ),
                TreeNode(
                    [
                        TreeNode(index=3),
                        TreeNode(index=4),
                    ],
                    [3, 2],
                ),
                TreeNode(index=5),
            ],
            [1, 1, 5],
        )
    )

    test_tree = neighbor_joining(dist)

    assert test_tree == ref_tree


def test_suite():
    passed = 0
    failed = 0

    test_distances_data = distances()
    test_tree = tree(test_distances_data)
    test_upgma_newick = upgma_newick()

    test_distances(test_tree)
    print("✓ distances")

    test_upgma(test_tree, test_upgma_newick)
    print("✓ upgma")

    test_neighbor_joining()
    print("✓ neighbor_joining")

if __name__ == "__main__":
    success = test_suite()
    sys.exit(0 if success else 1)