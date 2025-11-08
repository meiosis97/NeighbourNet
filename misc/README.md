# NeighbourNet

<p align="center">
<img width="30%" src="https://github.com/meiosis97/NeighbourNet/blob/main/misc/logo/NNet.png" > 
<img width="50%" src="https://github.com/meiosis97/NeighbourNet/blob/main/misc/logo/An_cell_specific_network.png"> 
</p>

|                                               |                                               |
|-----------------------------------------------|-----------------------------------------------|
|<img width="30%" src="https://github.com/meiosis97/NeighbourNet/blob/main/misc/logo/NNet.png" >  | <img width="50%" src="https://github.com/meiosis97/NeighbourNet/blob/main/misc/logo/An_cell_specific_network.png">  |


The **NeighbourNet (NNet)** package is currently under development.

The data and reproducible R scripts used to generate the figures in the
manuscript can be found at  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15031726.svg)](https://doi.org/10.5281/zenodo.15031726).

Note that the functions implemented in the final R package will differ
from those used in the analysis presented in our paper.  
The main difference lies in how NNet prunes the inferred co-expression
networks.  
The functions used in our paper’s analysis are provided in [this
script](https://github.com/meiosis97/NeighbourNet/blob/main/misc/script.R)
and can be directly sourced into R’s global environment.

In comparison, the pruning strategy implemented by the R package is more
statistically rigorous, although it does not necessarily yield better
results in our numerical evaluation.  
We conducted an analysis, available
[here](https://github.com/meiosis97/NeighbourNet/blob/main/misc/investigate.pruning.md),
to compare and understand the differences between the two strategies.

# ⚡ 5-Minute Quick-Start

This mini-walkthrough reproduces the **full vignette** in about 25
lines.  
The full vignette can be found
[here](https://github.com/meiosis97/NeighbourNet/blob/main/misc/vignettes.for.script.md).

------------------------------------------------------------------------

## 1 · Install / load core packages

``` r
remotes::install_github("https://github.com/meiosis97/NeighbourNet")
```

## 2 · Fetch demo dataset (can be found in the data folder)

``` r
# Load the package
require(NeighbourNet)

# Seurat object with ~4k LUAD cells
data(luad)     # loads `obj`
```

## 3 · One-liner preprocessing

``` r
genes  <- select.gene(obj, min.cells = 10) # QC → TF / target lists

obj <- obj |>
  prepare.seurat(genes = genes$genes) |>   # scale + PCA
  prepare.graph() |>                       # 30-NN graph
  select.cell() |>                         # subsample 
  prepare.reg(predictors = genes$tfs,      # local variance scaffolding
              responses  = genes$targets)
```

## 4 · Run NeighbourNet and build meta-networks

``` r
top10 <- head(genes$targets, 10)           # demo: first 10 targets
obj   <- run.nn.reg(obj, responses = top10, return.p.val = TRUE) |>
         build.meta.network() 
```

`obj@misc$mod` now contains:

| slot           | description                                        |
|----------------|----------------------------------------------------|
| `effect`       | (response × predictor × cell) co-expression tensor |
| `p.val`        | matching significance tensor                       |
| `meta.network` | (response × predictor × meta-cell) ensemble        |

## 5.0 Prepare Visualisation

``` r
ctr.genes <- select.central.genes(obj)  
obj <- prepare.visualise(obj, central.genes = ctr.genes)
```

## 5.1 · Snapshot plot (Cell \#1)

``` r
visualise.network(obj, 1, 
                  radius = c(.4,.7,.85,1), pie.radius = .04,
                  text.size = 5)
```

## 5.2 · Snapshot plot (meta-network \#1)

``` r
cut <- mean(apply(obj@misc$mod$meta.network$p.val[,,1], 1, max))
visualise.network(obj, 1, meta.network = TRUE, cutoff = cut,
                  radius = c(.4,.7,.85,1), pie.radius = .04,
                  text.size = 5)
```

## 6 · (Option) Receptor activity per cell

``` r
act  <- receptor.activity(obj)             # matrix: receptor × cell
lrp6 <- act$receptor.act["LRP6", ]
```

Plot on PCA:

``` r
library(ggplot2)
ggplot(Embeddings(obj, "pca"), aes(PC_1, PC_2)) +
  geom_point(alpha = .2, size = .6, colour = "grey80") +
  geom_point(data = Embeddings(obj, "pca")[names(lrp6), ],
             aes(col = as.numeric(lrp6)), size = 1)
```

### ✨ That’s it

You now have a full **NeighbourNet** analysis in under two minutes of
run-time. Tweak gene lists, increase `responses`, or switch
visualisation parameters exactly as in the full vignette.
