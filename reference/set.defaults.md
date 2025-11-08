# Set Default Parameters for Network Extraction

This function updates or resets default settings used in downstream
network extraction and analysis within a Seurat object. Defaults include
pruning thresholds, response and predictor gene sets, co-expression
measures, and transformation functions for effect estimates.

## Usage

``` r
set.defaults(seurat.obj, clean.up = FALSE, defaults = list())
```

## Arguments

- seurat.obj:

  A `Seurat` object with a `NNet.mod` list stored in the `misc` slot.
  This list is created by
  [`run.nn.reg`](https://meiosis97.github.io/NeighbourNet/reference/run.nn.reg.md).

- clean.up:

  A logical indicating whether to reset defaults to initial values
  (effect-based assay, quadratic transformation, self-loop removal, and
  pruning at cutoff 0.95). Default is `FALSE`.

- defaults:

  A named list of default parameters to update. Valid entries include:

  - f A function mapping effect estimates to importance scores. Default
    is `function(x) 2*x^2`.

  - remove.self.loops Logical, whether to remove self-loops. Default:
    `TRUE`.

  - assay Co-expression measure to use for downstream analysis. Options:
    `"effect"`, `"p.val"`, or `"meta.network"`.

  - predictors Character vector of predictor genes.

  - responses Character vector of response genes.

  - cutoff Numeric threshold for pruning, between 0 and 1. Default:
    `0.95`.

## Value

A `Seurat` object with updated default settings stored in
`misc$NNet.mod$defaults`.

## Details

The `defaults` list defines how networks are extracted and
post-processed in downstream analyses (e.g., pruning, TF activity
inference, meta-network learning). By setting `clean.up = TRUE`, all
defaults are reset to standardized values. Only valid arguments present
in the current `NNet.mod$defaults` are updated.

## See also

[`run.nn.reg`](https://meiosis97.github.io/NeighbourNet/reference/run.nn.reg.md),
[`prepare.reg`](https://meiosis97.github.io/NeighbourNet/reference/prepare.reg.md)

## Examples

``` r
# Reset defaults to initial values
seurat.obj <- set.defaults(seurat.obj, clean.up = TRUE)
#> Error: object 'seurat.obj' not found

# Update only the cutoff threshold and assay
seurat.obj <- set.defaults(
  seurat.obj,
  defaults = list(cutoff = 0.5, assay = "p.val")
)
#> Error: object 'seurat.obj' not found

# Check updated defaults
Seurat::Misc(seurat.obj, "mod")$defaults
#> Error: object 'seurat.obj' not found
```
