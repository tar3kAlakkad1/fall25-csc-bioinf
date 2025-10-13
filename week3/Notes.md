# Goal
Port `biotite`’s phylo package to Codon. The source code is written in Cython, not pure Python.

# Source code
The source can be found on [GitHub](https://github.com/biotite-dev/biotite/tree/v1.4.0/src/biotite/sequence/phylo).

# Steps
## Getting necessary files
### Source Code
There are four files found under `biotite/src/biotite/sequence/phylo`:

1) `__init__.py`: some necessary stuff for managing the Python package.

2) `nj.pyx`: a Cython file that only contains the function `neighbor_joining(np.ndarray distances)`. The function performs "hierarchical clustering using the *neighbor joining* algorithm." In contrast to UPGMA this algorithm does not assume a constant evolution rate. The resulting tree is considered to be unrooted.

3) `tree.pyx`: a Cython file that defines an object-oriented implementation of a `class Tree` with it's relevant functions. As described in the docstrings, a `class Tree` "represents a rooted tree (e.g. alignment guide tree or phylogenetic tree)."

4) `upgma.pyx`: a Cython file that only contains the function `upgma(np.ndarray distances)`. The function performs "hierarchical clustering using the *unweighted pair group method with arithmetic mean* (UPGMA)." This algorithm produces leaf nodes with the same distance to the root node. In the context of evolution this means a constant evolution rate (molecular clock).

### Tests
According to the deliverable instructions posted by the instructor on [piazza](https://piazza.com/class/mevafycmxgp28j/post/53), we only need the tests `test_distances`, `test_upgma`, and `test_neighbour_joining` from the test file `test_phylo.py`. 

The [test file](https://github.com/biotite-dev/biotite/blob/v1.4.0/tests/sequence/test_phylo.py) can be found on Github as well.

## Setting up Cython/Python
With the code files needed identified, main objective now is to get the Cython/Python tests running.

### Source code path
The Cython/Python src code will be in `fall25-csc-bioinf/week3/python/src` and the test file will be in `fall25-csc-bioinf/week3/python/test`.

### Getting tests running
No clue how to do this at the time of writing this tbh but we will figure it out lol.

## Setting up Codon
Once Cython/Python side of things is working, we will start porting the identified files to Codon (including the test file). 

### Source code path
The Codon src code will be in `fall25-csc-bioinf/week3/codon/src` and the test file will be in `fall25-csc-bioinf/week3/python/test`. 

## Setting up a CI/CD test
### Timing Python and Codon results
Must run Python and Codon tests and output the time of all tests in milliseconds. Use the time module in Python and Codon (do not time this from your evaluate.sh).

### Creating the `evaluate.sh`
Create a `evaluate.sh` bash script to use in `.github/workflows/actions.yml`

### Update `actions.yml`
Use `evaulate.sh` script in `.github/workflows/actions.yml` and change the code found in `.github/workflows/actions.yml` to trigger Week 3's evaluate script. 


