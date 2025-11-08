# Extract Per-Cell Gene Co-Expression Networks

Retrieves effect or p-value networks (or meta-networks) from a `Seurat`
object prepared by the NeighbourNet workflow. Supports optional pruning,
transformation, self-loop removal, and flexible cell indexing. By
default, argument values fall back to `NNet.mod$defaults` (`defaults`).

## Usage

``` r
get.network(
  seurat.obj,
  i = NULL,
  assay = NULL,
  remove.self.loops = NULL,
  responses = NULL,
  predictors = NULL,
  f = NULL,
  drop = TRUE,
  cutoff = NULL
)
```

## Arguments

- seurat.obj:

  A `Seurat` object with a `NNet.mod` list stored in the `misc` slot.
  This list is created by
  [`run.nn.reg`](https://meiosis97.github.io/NeighbourNet/reference/run.nn.reg.md).

- i:

  Cells for which to extract networks. Can be:

  - `NULL` (default): use the cells stored in `NNet.mod$cells`.

  - Integer indices matching `NNet.setting$all.cells`.

  - A character vector of cell names matching `NNet.setting$all.cells`.

- assay:

  Which network to extract. One of `"effect"`, `"p.val"`, or
  `"meta.network"`. If `NULL`, uses `defaults$assay`.

- remove.self.loops:

  Logical; whether to zero out self-connections. If `NULL`, uses
  `defaults$remove.self.loops`.

- responses:

  A character vector of response genes to extract. If `NULL`, uses
  `defaults$responses`.

- predictors:

  A character vector of predictor genes to extract. If `NULL`, uses
  `defaults$predictors`.

- f:

  A function to transform effect values for downstream use (e.g., an
  importance score). If `NULL`, uses `defaults$f`. Ignored when
  `assay = "p.val"` or `"meta.network"`.

- drop:

  Logical; if `TRUE` (default), drops the one slice arraies to matrices
  for convenience.

- cutoff:

  Numeric p-value threshold in \[0, 1\] for pruning based on edge
  significance stored in `NNet.mod$p.val`. If `NULL`, uses
  `defaults$cutoff`. The `p.val` tensor will be calculated instantly
  when it is not present in the `NNet.mod` list.

## Value

A network for the requested `assay` with dimensions:

- **effect / p.val**: (responses × predictors × cells) array. If
  `drop = TRUE` and any dimension equals 1, a matrix is returned:
  (responses x predictors) (a single cell), (predictors x cells) (a
  single response), or (responses x cells) (a single predictor).

- **meta.network**: the corresponding slice(s) from
  `NNet.mod$meta.network$meta.network[responses, predictors, i]` with
  the same dropping behavior.

## Details

**Assay behavior**

- *effect*: If `NNet.mod$smoothed = TRUE`, returns the stored smoothed
  effects for the requested cells (cells not in `NNet.mod$cells` will
  cause an error). If `NNet.mod$smoothed = FALSE`, effects are projected
  to the requested cells by applying the stored Laplacian operator
  `NNet.mod$w`.

- *p.val*: If `NNet.mod$p.val` exists and all requested cells are in
  `NNet.mod$cells`, returns the stored p-value tensor. Otherwise,
  p-values are computed from the (smoothed or projected) effects using
  per-response null distribution.

- *meta.network*: Returns slices of
  `NNet.mod$meta.network$meta.network`.

**Pruning**

- For *effect* assay, pruning uses p-values. Effects with corresponding
  p-values less than `cutoff` are zeroed.

- For *p.val* assay, entries less than `cutoff` are set to zero.

**Transformation**

- When `assay = "effect"`, the returned tensor/matrix is transformed by
  `f` (default `function(x) 2*x^2`) to yield an importance score. No
  transform is applied for `"p.val"` or `"meta.network"`.

## See also

[`run.nn.reg`](https://meiosis97.github.io/NeighbourNet/reference/run.nn.reg.md),
[`set.defaults`](https://meiosis97.github.io/NeighbourNet/reference/set.defaults.md)

## Examples

``` r
mod <- Seurat::Misc(seurat.obj, "NNet.mod")
#> Error: object 'seurat.obj' not found

# Default extraction (uses mod$defaults): smoothed effects, default cells
net <- get.network(seurat.obj)
#> Error: object 'seurat.obj' not found

# Effects for a specific set of cells (by name) and a gene panel
cells <- head(names(mod$cells), 3)
#> Error: object 'mod' not found
responses <- head(mod$gene.sets$responses$genes, 5)
#> Error: object 'mod' not found
predictors <- head(mod$gene.sets$predictors$genes, 10)
#> Error: object 'mod' not found
net.eff <- get.network(seurat.obj, i = cells,
                       assay = "effect", responses = responses, predictors = predictors)
#> Error: object 'seurat.obj' not found

# Extract P-values
net.p <- get.network(seurat.obj, assay = "p.val", cutoff = 0.9)
#> Error: object 'seurat.obj' not found

# Meta-network slices (3rd component) without dropping dimensions
net.meta <- get.network(seurat.obj, assay = "meta.network", i = 3, drop = FALSE)
#> Error: object 'seurat.obj' not found

# Effect network with pruning + custom transform, keeping self-loops
net.imp <- get.network(seurat.obj,
  assay = "effect", cutoff = 0.95,
  f = function(x) abs(x),
  remove.self.loops = FALSE
)
#> Error: object 'seurat.obj' not found
```
