# Prepare Gene and Receptor Settings for NeighbourNet Visualisation

This function prepares gene sets, hub genes, receptor priors, and
GRN-based evidence required for visualising NeighbourNet co-expression
networks. It selects and orders central genes, identifies
cluster-specific hub genes based on PC loadings, constructs a
receptor–gene prior matrix, and encodes optional prior GRN evidence for
edges. The resulting settings are stored in the `misc` slot of the
`Seurat` object for downstream plotting and visualisation routines.

## Usage

``` r
prepare.visualise(
  seurat.obj,
  n.clu = 4,
  central.genes = NULL,
  check.gr.evidence = TRUE,
  t = 2,
  p = NULL,
  as.g2 = c("predictors", "responses"),
  g1 = NULL,
  g2 = NULL,
  receptors = NULL
)
```

## Arguments

- seurat.obj:

  A `Seurat` object with a `NNet.mod` list stored in the `misc` slot.
  This list is created by
  [`run.nn.reg`](https://meiosis97.github.io/NeighbourNet/reference/run.nn.reg.md).

- n.clu:

  An integer specifying the maximum number of clusters used to group
  genes based on their PC loadings when defining hub genes. Default is
  `4`.

- central.genes:

  An optional list of central genes, typically the output of
  [`select.central.genes`](https://meiosis97.github.io/NeighbourNet/reference/select.central.genes.md),
  containing elements such as `central.responses` and
  `central.predictors`. If provided and `g1` / `g2` are `NULL`, these
  central genes are used to initialise `g1` and `g2`.

- check.gr.evidence:

  A logical indicating whether prior GRN evidence should be used to
  annotate and weight gene–gene relationships. If `TRUE`, TF–target
  adjacency from
  [`get.gr.adj`](https://meiosis97.github.io/NeighbourNet/reference/get.gr.adj.md)
  and a GRN graph `gr.graph` (in the NeighbourNet namespace) are used to
  distinguish supported and unsupported edges. Default is `TRUE`.

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

- as.g2:

  A character string indicating whether `g2` (the second gene layer)
  should correspond to `"predictors"` or `"responses"` in the
  NeighbourNet model. This choice determines which side of the effect
  tensor is treated as the upstream or downstream gene set. Default is
  `"predictors"`.

- g1:

  An optional character vector specifying the first gene layer to
  visualise (e.g. responses if `as.g2 = "predictors"`). If `NULL`, it is
  derived from `central.genes` and intersected with the corresponding
  default gene set in the NeighbourNet model.

- g2:

  An optional character vector specifying the second gene layer to
  visualise (e.g. predictors if `as.g2 = "predictors"`). If `NULL`, it
  is derived from `central.genes` and intersected with the corresponding
  default gene set in the NeighbourNet model.

- receptors:

  An optional character vector of receptors to visualise. If `NULL`, use
  all receptors available in the prior knowledge model obtained from
  [`get.prior.model`](https://meiosis97.github.io/NeighbourNet/reference/get.prior.model.md).

## Value

A `Seurat` object with an additional list named `NNet.visual.setting`
stored in its `misc` slot. This list contains:

- g1:

  The ordered first-layer gene set used for visualisation.

- g2:

  The ordered second-layer gene set used for visualisation.

- clu.g12:

  A `hclust` object describing hierarchical clustering of selected genes
  based on PCs.

- clu.g12:

  A `hclust` object describing hierarchical clustering of `g2` genes
  based on PCs.

- hubs:

  A character vector of hub genes, one per cluster, selected to
  represent major loading patterns.

- receptors:

  The subset of receptors to visualise, if any.

- ppr:

  A receptor–gene prior model (receptors x `g2` set).

- as.g2:

  The role of `g2` in the NeighbourNet model (`"predictors"` or
  `"responses"`).

- evidence:

  A matrix encoding prior GRN evidence for gene–gene relationships
  between `g1` and the `g2` set, used for visual annotation.

## Details

This function organises and annotates genes for visualising NeighbourNet
co-expression networks. Starting from either user-specified gene sets
(`g1`, `g2`) or central genes obtained from
[`select.central.genes`](https://meiosis97.github.io/NeighbourNet/reference/select.central.genes.md),
it defines two gene layers aligned with the response and predictor sides
of the NeighbourNet model, depending on `as.g2`. These genes are
projected into the PC space, and clustered via hierarchical clustering
on their PC correlations. Within each cluster, a representative hub gene
is selected, yielding one hub per cluster that captures dominant
variation patterns across PCs.

In parallel, a receptor–g2 prior model is extracted and stored as `ppr`.
It receptors' prior regulatory potential to the second gene layer and
can be used to overlay upstream receptor influence in downstream
visualisation.

If `check.gr.evidence = TRUE`, GRN-based evidence is extrated from
[`get.gr.adj`](https://meiosis97.github.io/NeighbourNet/reference/get.gr.adj.md)
to construct an `evidence` matrix that annotates co-expression between
g1 and g2 as supported or unsupported by prior regulatory knowledge (and
optionally distinguishes activation and inhibition). If
`check.gr.evidence = FALSE`, all co-expression edges are treated as
stimulatory in the visualisation.

## See also

[`select.central.genes`](https://meiosis97.github.io/NeighbourNet/reference/select.central.genes.md),
[`get.gr.adj`](https://meiosis97.github.io/NeighbourNet/reference/get.gr.adj.md)

## Examples

``` r
# Select central genes from meta-networks
sel <- select.central.genes(seurat.obj)
#> Error: object 'seurat.obj' not found

# Prepare visualisation settings using central genes (predictors as g2)
seurat.obj <- prepare.visualise(
  seurat.obj,
  n.clu = 4,
  as.g2 = "predictors",
  central.genes = sel
)
#> Error: object 'seurat.obj' not found

# Access visual settings
vis.set <- Seurat::Misc(seurat.obj, "NNet.visual.setting")
#> Error: object 'seurat.obj' not found
str(vis.set)
#> Error: object 'vis.set' not found
```
