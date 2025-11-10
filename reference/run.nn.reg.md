# Run Nearest-Neighbor PC Regression for Gene Co-Expression Analysis

This function performs nearest-neighbor regression on a Seurat object to
measure gene co-expression. It uses principal components (PCs) to smooth
response gene expression and calculates permutation feature importance
as a co-expression measure. The results are stored as a list named
`NNet.mod` in the `misc` slot of the Seurat object for downstream
analysis.

## Usage

``` r
run.nn.reg(
  seurat.obj,
  responses = NULL,
  Y = NULL,
  predictors = NULL,
  t = 3,
  k = NULL,
  remove.self.loops = T,
  f = function(x) 2 * x^2,
  assay = c("effect", "p.val"),
  prune = TRUE,
  cutoff = 0.95,
  return.p.val = FALSE,
  return.smooth = TRUE,
  return.prune = FALSE
)
```

## Arguments

- seurat.obj:

  A `Seurat` object on which
  [`prepare.graph`](https://meiosis97.github.io/NeighbourNet/reference/prepare.graph.md)
  and
  [`prepare.reg`](https://meiosis97.github.io/NeighbourNet/reference/prepare.reg.md)
  have been run.

- responses:

  A character vector of response gene names to return for the network
  ensemble. If `NULL`, uses all response genes from `NNet.setting`.

- Y:

  A custom response matrix (cells x responses). If provided, overrides
  `responses`.

- predictors:

  A character vector of predictor gene names to return for the network
  ensemble. If `NULL`, uses all predictor genes from `NNet.setting`.

- t:

  A numeric value indicating the power to scale the singular values of
  the Laplacian operator. Default is 3.

- k:

  An integer specifying the number of singular vectors to use in the
  Laplacian operator. Default is the number of nearest neighbors.

- remove.self.loops:

  A logical indicating whether to remove self-loops from networks for
  downstream analysis. Default is `TRUE`.

- f:

  A function to map the effect estimation to importance scores for
  downstream analysis. Default is `function(x) 2*x^2` that yields
  permutation feature importamce.

- assay:

  A character string indicating the co-expression measure to use in
  downstream analysis. Options are `"effect"` or `"p.val"`. Default is
  `"effect"`.

- prune:

  A logical indicating whether to prune networks based on p-values for
  downstream analysis. Default is `TRUE`.

- cutoff:

  A numeric value specifying the p-value threshold for pruning for
  downstream analysis.. Default is 0.95.

- return.p.val:

  A logical indicating whether to return p-values in addition to the
  effect tensor. Set to `FALSE` to save memory. Default is `FALSE`.

- return.smooth:

  A logical indicating whether to return the smoothed effect tensor. If
  set to `TRUE`, networks will not be further smoothed in downstream
  analysis. Default is `TRUE`.

- return.prune:

  A logical indicating whether to return the pruned effect tensor. If
  set to `TRUE`, networks will not be further pruned in downstream
  analysis. Default is `FALSE`.

## Value

A `Seurat` object with a list named `NNet.mod` stored in its `misc`
slot. `NNet.mod` is a list containing:

- effect:

  A tensor of effect sizes (responses x predictors x cells).

- p.val:

  A tensor of p-values (responses x predictors x cells) if
  `return.p.val = TRUE`.

- meta.network:

  Reserved for
  [`build.meta.network`](https://meiosis97.github.io/NeighbourNet/reference/build.meta.network.md).

- meta.response:

  Reserved for
  [`build.meta.response`](https://meiosis97.github.io/NeighbourNet/reference/build.meta.network.md).

- mus:

  A named vector of mean log-transformed noise distributions for each
  response.

- sigmas:

  A named vector of standard deviations for log-transformed noise
  distributions for each response.

- subsampled:

  A logical indicating whether cell subsampling was performed.

- smoothed:

  A logical indicating whether smoothed effects are returned.

- pruned:

  A logical indicating whether pruned effects are returned.

- gene.sets:

  A list of selected predictors, responses, and all genes used in the
  analysis. Includes TFs and targets subsets.

- cells:

  A vector of selected cell indices used in the analysis.

- defaults:

  Default settings used for network extraction, such as cutoff and
  pruning function.

- custom.y:

  A logical indicating whether a custom response matrix (`Y`) was
  provided.

- w:

  A list containing the Laplacian operator components (`u` and `vd`).

## Details

This function performs nearest-neighbor PC regression to assess
co-expression relationships between response and predictor genes within
a Seurat object. The co-expression measure is based on permutation
feature importance estimated by local gene variances, enabling scalable
analysis of transcriptional regulation. Parameters `remove.self.loops`,
`f`, `assay`, `prune`, and `cutoff` define default settings for how
networks will be post-processed before downstream analysis, such as
meta-network learning. However, these settings do not immediately affect
the regression analysis performed by this function. These defaults can
be updated later using
[`set.defaults`](https://meiosis97.github.io/NeighbourNet/reference/set.defaults.md).

## See also

[`prepare.graph`](https://meiosis97.github.io/NeighbourNet/reference/prepare.graph.md),
[`prepare.reg`](https://meiosis97.github.io/NeighbourNet/reference/prepare.reg.md),
[`set.defaults`](https://meiosis97.github.io/NeighbourNet/reference/set.defaults.md)

## Examples

``` r
# Assuming `seurat.obj` has been pre-processed with `prepare.reg`:
responses <- Seurat::Misc(seurat.obj, "NNet.setting")$responses %>%
             head(n =5)
#> Error: object 'seurat.obj' not found

# Run PC regression
seurat.obj <- run.nn.reg(
  seurat.obj,
  responses = responses,
)
#> Error: object 'seurat.obj' not found

# Check regression outcome settings
str(Seurat::Misc(seurat.obj, "NNet.mod"))
#> Error: object 'seurat.obj' not found
```
