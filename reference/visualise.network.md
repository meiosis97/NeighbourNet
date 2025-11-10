# Visualise Receptor–TF–Target Pathways on NeighbourNet Cell-Specific Networks.

This function visualises NeighbourNet co-expression networks as a
layered receptor–TF–target signalling diagram. It takes pre-computed
visualisation settings from
[`prepare.visualise`](https://meiosis97.github.io/NeighbourNet/reference/prepare.visualise.md)
and, for a given cell or meta-network component, builds a multi-layer
graph linking upstream receptors to downstream gene modules. Edges are
weighted by co-expression strength and annotated by prior GRN evidence,
and node shapes are represented as cluster-level pie charts.  
  
**Interpretation** When NeighbourNet is used to build TF–target
co-expression networks (for example, targets as responses and TFs as
predictors), this function provides prior knowledge enriched
visualisation of networks as four concentric layers. From the innermost
to the outermost layer, these are:

- `g1` (targets): target genes, placed at the centre of the plot.

- `g2` (TFs): transcription factors surrounding the targets. TFs are
  grouped by expression-based clustering (node colour), and each TF is
  assigned a significance profile summarised as a pie chart (the
  fraction of the pie filled reflects the likelihood of co-expression
  with the `g1` targets).

- `g3` (intermediate signalling nodes): optional intermediate signalling
  molecules placed between TFs and receptors. These nodes represent the
  shortest signalling paths that connect receptors to TFs in a prior
  signalling graph.

- `g4` (receptors): the most active upstream receptors on the outermost
  layer. Each receptor is linked to a TF that is inferred to mediate its
  influence on the targets. Receptors are coloured in accordance with TF
  clusters, and the proportion of colour fill reflects their relative
  activity or expression.

The connections between `g1` and `g2` encode TF–target co-expression
annotated by prior gene-regulatory evidence. Edges represent three
levels of support: no edge (no significant relationship), dashed edges
for significant co-expression only, and solid edges when significant
co-expression is also supported by prior GRN evidence (with different
arrow styles indicating activation versus repression). The `g3` and `g4`
layers summarise upstream signalling pathways that exert strong
influence on the observed TF–target co-expression patterns.

## Usage

``` r
visualise.network(
  seurat.obj,
  i,
  meta.network = FALSE,
  fix.cluster = TRUE,
  hubs = NULL,
  n.clu = 4,
  cutoff = NULL,
  show.pathways = TRUE,
  change.receptors = TRUE,
  receptor.activity = c("cprod", "dist"),
  check.receptor.expression = TRUE,
  scale.ppr = FALSE,
  scale.network = FALSE,
  scale.signifiance = FALSE,
  swap.layers = FALSE,
  k = 2,
  n.per.component = 10,
  radius = NULL,
  pie.radius = 0.05,
  text.size = 4
)
```

## Arguments

- seurat.obj:

  A `Seurat` object with a `NNet.visualise.setting` list stored in the
  `misc` slot. This list is created by
  [`prepare.visualise`](https://meiosis97.github.io/NeighbourNet/reference/prepare.visualise.md).

- i:

  Index of the cell or meta-network component to visualise.

- meta.network:

  A logical indicating whether to visualise a NeighbourNet meta-network
  component instead of a per-cell network. If `TRUE`, the
  `"meta.network"` assay is used. Default is `FALSE`.

- fix.cluster:

  A logical indicating whether to reuse the gene clustering from
  `NNet.visualise.setting` when assigning clusters to the second gene
  layer (`g2`). If `FALSE`, clusters are re-estimated from the current
  network. Default is `TRUE`.

- hubs:

  An optional character vector of hub genes in `g2`. If supplied and
  `fix.cluster = FALSE`, clusters are inferred by assigning each `g2`
  gene to the closest hub. If `NULL`, hubs are not used for clustering.

- n.clu:

  An integer specifying the maximum number of clusters for `g2` genes.
  Default is `4`. The effective number of clusters is truncated when few
  genes are available.

- cutoff:

  A numeric value specifying the p-value threshold used to define
  significant co-expression for edge evidence. If `NULL`, the default is
  taken from `NNet.mod$defaults$cutoff` stored in the `misc` of the
  Seurat object. Values less than 0 are clamped to 0.

- show.pathways:

  A logical indicating whether receptor layers and intermediate
  signalling paths (via a prior signalling graph) should be included in
  the visualisation. If `FALSE`, only the two gene layers (`g1` and
  `g2`) are plotted. Default is `TRUE`.

- change.receptors:

  A logical indicating whether to select a subset of receptors for
  plotting based on activity (e.g., eigen-like centrality on inferred
  receptor–target activity). If `FALSE`, receptors supplied in
  `NNet.visualise.setting` are used as is. Default is `TRUE`.

- receptor.activity:

  A character string specifying how the prior model of receptor–TF
  interaction and TF–target (`g1`–`g2`) co-exression network are
  combined when computing receptor–g1 activity. Options are `"cprod"`
  (cross-product aggregation) and `"dist"` (max-over-path aggregation).
  Default is `"cprod"`.

- check.receptor.expression:

  A logical indicating whether receptor activity should be weighted by
  receptor expression in the selected cell or component. Default is
  `TRUE`.

- scale.ppr:

  A logical indicating whether to normalise the receptor–TF prior model
  (PageRank matrix) by receptor-wise norms. Default is `TRUE`.

- scale.network:

  A logical indicating whether to normalise the co-expression network by
  `g1`-wise norms before computing activity. Default is `TRUE`.

- scale.signifiance:

  A logical indicating whether node and edge weights should be scaled by
  their maximum significance (derived from p-values), rather than raw
  co-expression magnitudes. Default is `FALSE`.

- swap.layers:

  A logical indicating whether to reverse the order of `g1` and `g2` in
  the concentric layout (inner vs outer ring). Default is `FALSE`.

- k:

  An integer controlling the number of singular vectors used when
  selecting active receptors (layer `g4`) based on receptor–target
  activity, analogous to an eigenvector centrality based selection in
  [`select.central.genes`](https://meiosis97.github.io/NeighbourNet/reference/select.central.genes.md)
  Default is `2`.

- n.per.component:

  An integer indicating the number of top receptors per component (or
  singular vector) to consider when selecting candidate receptors.
  Default is `10`.

- radius:

  Optional numeric vector control the radii of concentric layers.

- pie.radius:

  A numeric value giving the base radius of node pies in the scatterpie
  representation. Default is `0.05`.

- text.size:

  A numeric value giving the size of gene labels in the plotted network.
  Default is `4`.

## Value

A `ggplot` object (built with ggraph and ggplot2) representing a layered
receptor–TF–target network. The plot includes:

- concentric layers for `g1` and `g2`, and optionally intermediate
  signalling nodes (`g3`) and receptors (`g4`);

- edge arrows and line types encoding prior evidence and direction;

- node pies summarising cluster assignments and activity contributions.

The underlying graph structure is built internally using igraph, but
only the final `ggplot` object is returned.

## Details

This function visualises a single NeighbourNet network (either a
cell-level network or a meta-network component) using gene sets and
priors prepared by
[`prepare.visualise`](https://meiosis97.github.io/NeighbourNet/reference/prepare.visualise.md).
The first two layers (`g1` and `g2`) correspond to response and
predictor gene sets. Co-expression between `g1` and `g2` is extracted
from the selected NeighbourNet assay (effect, p-value, or meta-network)
via
[`get.network`](https://meiosis97.github.io/NeighbourNet/reference/get.network.md),
with associated p-values used to define an edge-specific significance
measure.

Node weights for `g1` and `g2` are derived from the maximum significance
of their incident edges and can be optionally rescaled by
`scale.signifiance`. Edge weights are then aggregated into cluster-level
pies, showing how strongly each layer contributes to different `g2`
modules. The `evidence` matrix from `NNet.visual.setting` encodes prior
GRN evidence for `g1`–`g2` pairs and is used to annotate edges as
stimulatory or inhibitory where such information is available.

When `show.pathways = TRUE`, the function also incorporates a
receptor–`g2` prior (generated by
[`get.prior.model`](https://meiosis97.github.io/NeighbourNet/reference/get.prior.model.md))
and a signalling graph (generated by
[`get.gr.adj`](https://meiosis97.github.io/NeighbourNet/reference/get.gr.adj.md))
to construct upstream receptor layers. Receptor-to-`g1` activity is
first inferred by combining the receptor–`g2` prior model with the
pruned `g1`–`g2` co-expression network, and weighted by receptor
expression in the selected cell or component
(`check.receptor.expression`). A subset of receptors (`g4`) and
intermediate signalling nodes (`g3`) is then chosen based on activity
and shortest paths in `sig.graph`, and added as outer layers in the
layout.

Finally, the multi-layer graph is laid out using a concentric layout,
with the order of layers controlled by `swap.layers` and `radius`, and
rendered via ggraph. Pies from ggforce::geom_scatterpie summarise
cluster-level contributions for each node, while edge arrows and line
types encode directionality and prior evidence.

## See also

[`prepare.visualise`](https://meiosis97.github.io/NeighbourNet/reference/prepare.visualise.md),
[`receptor.activity`](https://meiosis97.github.io/NeighbourNet/reference/receptor.activity.md)

## Examples

``` r
# Assume NeighbourNet has been fitted and visual settings prepared:
# seurat.obj <- prepare.graph(...)
# seurat.obj <- prepare.reg(seurat.obj, ...)
# seurat.obj <- run.nn.reg(seurat.obj, ...)
# seurat.obj <- build.meta.network(seurat.obj, ...)
# seurat.obj <- prepare.visualise(seurat.obj)

# Visualise the network for a single cell
p <- visualise.network(seurat.obj, i = 1)
#> Error: object 'seurat.obj' not found
print(p)
#> Error: object 'p' not found

# Visualise a meta-network component with receptor pathways
p.meta <- visualise.network(
  seurat.obj,
  i = 1,
  meta.network = TRUE,
)
#> Error: object 'seurat.obj' not found
print(p.meta)
#> Error: object 'p.meta' not found
```
