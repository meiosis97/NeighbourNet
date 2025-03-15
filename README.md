# NeighbourNet

The NeighbourNet (NNet) package is currently under development. Note that the functions implemented in the final R package will differ from those used in the analysis presented in our paper.

The main difference lies in how NNet prunes the inferred co-expression networks. The functions used in our paper's analysis are provided in [this script](./tests/script.R) and can be directly sourced into R’s global environment. For now, if you wish to use the alternative pruning strategy intended for the R package, please source the R script located [here](./tests/script.new.R).

In comparison, the pruning strategy described in the paper is heuristic and less statistically rigorous, although it yields better results in our numerical evaluation. We conducted an analysis, available [here](./tests/investigate.pruning.Rmd), to compare and understand the differences between the two strategies.
