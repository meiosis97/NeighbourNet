# NeighbourNet

The NeighbourNet (NNet) package is currently under development. 

The data and the reproducible R scripts used to generate the figures in the manuscript can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15031726.svg)](https://doi.org/10.5281/zenodo.15031726).

Note that the functions implemented in the final R package will differ from those used in the analysis presented in our paper.

The main difference lies in how NNet prunes the inferred co-expression networks. The functions used in our paper's analysis are provided in [this script](./tests/script.R) and can be directly sourced into R’s global environment. For now, if you wish to use the alternative pruning strategy intended for the R package, please source the R script located [here](./tests/script.new.R).

In comparison, the pruning strategy described in the paper is heuristic and less statistically rigorous, although it yields better results in our numerical evaluation. We conducted an analysis, available [here](./tests/investigate.pruning.md), to compare and understand the differences between the two strategies.

## ── NeighbourNet: rapid test drive ─────────────────────────────────────────────
## 0. House‑keeping ----------------------------------------------------------------
pkgs <- c("Seurat","dplyr","Matrix","ggplot2","ggraph",
          "scatterpie","ggrepel","ggpubr","igraph")
if(any(miss <- !pkgs %in% installed.packages()[,1]))
    install.packages(pkgs[miss], repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

## 1. Load priors + demo Seurat object -------------------------------------------
lapply(c("gene.list","sig.graph","gr.graph","receptor.ppr"),
       \(f) load(file.path("../data", paste0(f,".rda"))))
load("data/luad.rda")                                           # -->  obj

## 2. Minimal preprocessing (QC, PCA, KNN, subsample, local‑var settings) ---------
rt.ppr <- get.ppr()                                             # receptor–target matrix
genes  <- select.gene(obj, min.cells = 10)                      # TF / target lists

obj <- obj              |>                                      # original Seurat obj
       prepare.seurat(genes = genes$genes) |>                   # scaling + PCA
       prepare.graph()              |>                          # KNN construction
       select.cell()                |>                          # subsample if n > 5k
       prepare.reg(predictors = genes$tfs,
                   responses  = genes$targets)                  # local variance etc.

## 3. Run NeighbourNet on first 10 response genes + build meta‑networks ----------
top10 <- head(genes$targets, 10)
obj   <- run.nn.reg(obj, responses = top10, return.p.val = TRUE) |>
         build.meta.network()

## 4. One‑click snapshot: plot meta‑network #1 -----------------------------------
cut   <- mean(apply(obj@misc$mod$meta.network$p.val[,,1], 1, max))
visualise.network(obj, 1, meta.network = TRUE, cutoff = cut,
                  radius = c(.4,.7,.85,1), pie.radius = .04, text.size = 5)

## You now have:
##   obj@misc$mod$effect         # cell‑level TF‑target networks
##   obj@misc$mod$meta.network   # meta‑networks for downstream exploration
##   receptor.activity(obj)      # receptor activity scores (optional)
