# NeighbourNet

The NeighbourNet (NNet) package is currently under development. 

The data and the reproducible R scripts used to generate the figures in the manuscript can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15031726.svg)](https://doi.org/10.5281/zenodo.15031726).

Note that the functions implemented in the final R package will differ from those used in the analysis presented in our paper.

The main difference lies in how NNet prunes the inferred co-expression networks. The functions used in our paper's analysis are provided in [this script](./tests/script.R) and can be directly sourced into R’s global environment. For now, if you wish to use the alternative pruning strategy intended for the R package, please source the R script located [here](./tests/script.new.R).

In comparison, the pruning strategy described in the paper is heuristic and less statistically rigorous, although it yields better results in our numerical evaluation. We conducted an analysis, available [here](./tests/investigate.pruning.md), to compare and understand the differences between the two strategies.

