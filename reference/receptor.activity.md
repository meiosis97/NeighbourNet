# Infer Receptor Activity from NeighbourNet TF–Target Networks

This function infers receptor, transcription factor (TF), and target
activity by combining a receptor–TF prior model (personalised PageRank)
with NeighbourNet TF–target co-expression networks. Activity can be
evaluated at the per-cell level or on meta-networks, and can be
optionally filtered by receptor expression and prior gene regulatory
network (GRN) evidence.

## Usage

``` r
receptor.activity(
  seurat.obj,
  i = NULL,
  meta.network = FALSE,
  cutoff = NULL,
  check.receptor.expression = TRUE,
  check.gr.evidence = TRUE,
  t = 2,
  p = NULL,
  scale.ppr = TRUE,
  scale.network = TRUE,
  as.tfs = c("predictors", "responses"),
  receptors = NULL,
  tfs = NULL,
  targets = NULL,
  receptor.activity = c("cprod", "dist")
)
```

## Arguments

- seurat.obj:

  A `Seurat` object with a `NNet.mod` list stored in the `misc` slot.
  This list is created by
  [`run.nn.reg`](https://meiosis97.github.io/NeighbourNet/reference/run.nn.reg.md).
  To optionally infer receptor activity on meta-networks (set
  `meta.network = TRUE`), the `NNet.mod` should be updated by a
  subsequent run of
  [`build.meta.network`](https://meiosis97.github.io/NeighbourNet/reference/build.meta.network.md).

- i:

  Optional index of cells or meta-network components for which to
  compute receptor activity. If `NULL`, uses `NNet.mod$cells` for
  per-cell networks or all meta-network components when
  `meta.network = TRUE`.

- meta.network:

  A logical indicating whether to compute activity on NeighbourNet
  meta-networks instead of per-cell networks. If `TRUE`, the
  `"meta.network"` assay is used and each component is treated as a
  meta-cell. Default is `FALSE`.

- cutoff:

  A numeric value specifying the p-value threshold to prune
  insignificant co-expression values. If `NULL`, the default value is
  taken from `NNet.mod$defaults`.

- check.receptor.expression:

  A logical indicating whether receptor activity should be down-weighted
  by receptor expression. Default is `TRUE`.

- check.gr.evidence:

  A logical indicating whether TF–target edges should be filtered by
  prior GRN evidence from
  [`get.gr.adj`](https://meiosis97.github.io/NeighbourNet/reference/get.gr.adj.md).
  Default is `TRUE`.

- t:

  A numeric value passed to
  [`get.gr.adj`](https://meiosis97.github.io/NeighbourNet/reference/get.gr.adj.md)
  that controls the depth or order of TF–target adjacency used as prior
  GRN evidence. Default is `2`.

- p:

  A numeric value passed to
  [`get.prior.model`](https://meiosis97.github.io/NeighbourNet/reference/get.prior.model.md)
  specifying the quantile threshold for pruning the PPR matrix. If
  `NULL`, the default threshold stored in
  [`receptor.ppr`](https://meiosis97.github.io/NeighbourNet/reference/receptor.ppr.md)`$ltf`
  is used. Default is `NULL`.

- scale.ppr:

  A logical indicating whether to normalise the receptor–TF prior model
  (PageRank matrix) by receptor-wise norms. Default is `TRUE`.

- scale.network:

  A logical indicating whether to normalise the TF–target co-expression
  network by target-wise norms before computing activity. Default is
  `TRUE`.

- as.tfs:

  A character string indicating whether TFs are encoded as
  `"predictors"` (columns) or `"responses"` (rows) in the NeighbourNet
  networks. Default is `"predictors"`.

- receptors:

  A character vector of receptor genes to include. If `NULL`, all
  receptors present in the prior model returned by
  [`get.prior.model`](https://meiosis97.github.io/NeighbourNet/reference/get.prior.model.md)
  are used.

- tfs:

  A character vector of TFs to include. If `NULL`, TFs are taken from
  the default gene set stored in `NNet.mod$defaults`.

- targets:

  A character vector of target genes to include. If `NULL`, targets are
  taken from the default gene set stored in `NNet.mod$defaults`.

- receptor.activity:

  A character string specifying how receptor–TF and TF–target
  information should be combined. Options are `"cprod"` (cross-product
  aggregration that measures cosine similarity, default) and `"dist"`
  (max-over-path aggregation).

## Value

If multiple cells or meta-network components are requested
(`length(i) > 1` or `i = NULL`), the function returns a list with:

- receptor.act:

  A matrix of receptor activity scores (receptors x cells/components).

- tf.act:

  A matrix of TF activity scores (TFs x cells/components).

- target.act:

  A matrix of target activity scores (targets x cells/components).

If a single cell or component is requested, the function returns a list
with:

- act.mat:

  A receptor–target activity matrix (receptors x targets).

- ppr:

  The receptor–TF prior used after any scaling and expression weighting.

- network:

  The TF–target network used for activity computation.

## Details

This function integrates three hierarchical components to infer receptor
activity: (i) a receptor–TF prior model obtained from
[`get.prior.model`](https://meiosis97.github.io/NeighbourNet/reference/get.prior.model.md),
(ii) TF–target co-expression networks derived from NeighbourNet, and
(iii) optional TF–target adjacency matrices from
[`get.gr.adj`](https://meiosis97.github.io/NeighbourNet/reference/get.gr.adj.md)
as prior gene regulatory network (GRN) evidence. Together, these
components define a receptor–TF–target signalling cascade that
quantifies how upstream receptors potentially influence downstream
target genes through TFs.

For each cell or meta-network component, the receptor–TF prior model
encodes the probability of regulatory transmission from receptors to TF,
while the TF–target network represents local co-expression relationships
inferred by NeighbourNet. Receptor activity is computed by aggregating
these two layers to estimate receptor influence on downstream targets,
mediated by TFs. The resulting matrices describe receptor-, TF-, and
target-level activity for each cell or component, providing an
interpretable summary of upstream signalling dynamics.

Activity can be inferred from per-cell networks (default) or from
meta-networks when `meta.network = TRUE`. The arguments `scale.ppr` and
`scale.network` apply L2-normalisation to the prior model and
co-expression network respectively, reducing the influence of highly
connected nodes and ensuring comparable contribution across genes. The
options `check.receptor.expression` and `check.gr.evidence` introduce
two filters on receptor activity : the first enforces that receptors
must be transcriptionally detected, while the second prunes TF–target
co-expression to those supported by prior GRN knowledge.

The argument `receptor.activity` controls how the two layers are
combined: `"cprod"` performs a cross-product aggregation equivalent to
measuring cosine similarity between receptor–TF and TF–target profiles,
while `"dist"` applies a max-over-path operation that highlights the
strongest possible receptor–target link through any intermediate TF.

## See also

`get.prior.model`, `get.gr.adj`

## Examples

``` r
# Per-cell receptor activity using default settings
act <- receptor.activity(seurat.obj)
#> Error: object 'seurat.obj' not found

# Meta-network receptor activity for the first three components
act <- receptor.activity(
  seurat.obj,
  meta.network = TRUE,
  i = 1:3
)
#> Error: object 'seurat.obj' not found

```
