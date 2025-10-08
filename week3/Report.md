# Jotting Notes

# Python steps
## Step 1: Downloaded necessary files. 
Downloaded the necessary files from the provided Github [link](https://github.com/biotite-dev/biotite/tree/v1.4.0/src/biotite/sequence/phylo) from Piazza. 
 

### Downloaded files
#### Source code
1) 1) `__init__.py`: some necessary stuff for managing the Python package.

2) `nj.pyx`: a Cython file that only contains the function `neighbor_joining(np.ndarray distances)`. The function performs "hierarchical clustering using the *neighbor joining* algorithm." In contrast to UPGMA this algorithm does not assume a constant evolution rate. The resulting tree is considered to be unrooted.

3) `tree.pyx`: a Cython file that defines an object-oriented implementation of a `class Tree` with it's relevant functions. As described in the docstrings, a `class Tree` "represents a rooted tree (e.g. alignment guide tree or phylogenetic tree)."

4) `upgma.pyx`: a Cython file that only contains the function `upgma(np.ndarray distances)`. The function performs "hierarchical clustering using the *unweighted pair group method with arithmetic mean* (UPGMA)." This algorithm produces leaf nodes with the same distance to the root node. In the context of evolution this means a constant evolution rate (molecular clock).

#### Test file
According to the deliverable instructions posted by the instructor on [piazza](https://piazza.com/class/mevafycmxgp28j/post/53), we only need the tests `test_distances`, `test_upgma`, and `test_neighbour_joining` from the test file `test_phylo.py`. 

Therefore, I took the following steps to get tests running locally:

1) Download the file `test_phylo.py` from the source code on [Github](https://github.com/biotite-dev/biotite/blob/v1.4.0/tests/sequence/test_phylo.py).

2) Downloaded the necessary data files for test:
    - `distances.txt` ([ref](https://github.com/biotite-dev/biotite/blob/v1.4.0/tests/sequence/data/distances.txt)) used in `test_distances` 
    - and `newick_upgma.txt` ([ref](https://github.com/biotite-dev/biotite/blob/v1.4.0/tests/sequence/data/newick_upgma.txt))

3) Downloaded necessary `util.py` file from the source code on [Github](https://github.com/biotite-dev/biotite/blob/v1.4.0/tests/util.py).

4) Commented out the unwanted tests from `test_phylo.py`.

5) Ran `pyenv virtualenv 3.13 deliverable-3` to create a Python virtual environment.

6) Ran `pip install biotite numpy pytest` to downloaded necessary packages: 
    - `pytest` for tests,
    - `numpy` as it's a dependency for the code,
    - and `biotite` as specified under the "Important" section  on deliverable instructions on Piazza,

7) Ran `pip freeze > requirements.txt` to save the dependencies for later (for CI).

8) Ran `pytest test_phylo.py`