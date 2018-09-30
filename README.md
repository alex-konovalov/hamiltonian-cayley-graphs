[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/alex-konovalov/hamiltonian-cayley-graphs/master)
# hamiltonian-cayley-graphs

## Reproducible experiment for arXiv:1805.00149

This repository provides a working environment needed to reproduce
calculations described in the arXiv preprint "Cayley graphs of order
kp are hamiltonian for k < 48" by Dave Witte Morris and Kirsten Wilk
(<https://arxiv.org/abs/1805.00149>).

It uses the Docker container with the latest public release of GAP with
added LKH - a C implementation of the Lin-Kernighan heuristic for solving
the traveling salesman problem (<http://akira.ruc.dk/~keld/research/LKH/>)
by Keld Helsgaun. The container is maintained in a separate repository 
at <https://github.com/alex-konovalov/gap-docker-lkh>.

The `gap` directory contains the GAP source code and Jupyter notebooks,
and the `log` directory contains output files. The GAP input and output
code from <https://arxiv.org/src/1805.00149v1/anc> was used to initialise
this repository.

To use the code on Binder (<https://mybinder.org/>), perform the following
steps:

1. Click on the "launch binder" badge in this README file on GitHub, or open
<https://mybinder.org/v2/gh/alex-konovalov/hamiltonian-cayley-graphs/master>
in your browser.

2. First "Loading repository: alex-konovalov/hamiltonian-cayley-graphs/master"
will be displayed, followed by a non-interactive preview. When the environment
will be ready, you will see the main Jupyter screen with a list of files.

3. Click on any `.ipynb` file to open it. This will start GAP Jupyter notebook.
You should be able to run the code in the notebook, e.g. by executing one cell
after another, or by running all cells, and perform other actions, e.g. add
new cells to the notebook, or start a new notebook. Your changes will not be
preserved after the window will be closed, but you should be able to download
the notebook in various formats. For further information about Jupyter, see 
<https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/index.html>.

Note that GitHub renders Jupyter notebooks in a human readable form, see e.g.
<https://github.com/alex-konovalov/hamiltonian-cayley-graphs/blob/master/gap/4-5-Valence2.ipynb>. 