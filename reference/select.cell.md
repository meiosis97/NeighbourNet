# Select Cells for Downstream PC Regression Analysis

This function selects a subset of cells from a `Seurat` object for
downstream analysis based on their principal component (PC) embeddings.
The selection ensures that the subset represents the distribution of the
full dataset while maintaining computational efficiency.

## Usage

``` r
select.cell(seurat.obj, p = 0.1, n = NULL, all = FALSE, ...)
```

## Arguments

- seurat.obj:

  A `Seurat` object that has been processed with
  [`prepare.graph`](https://meiosis97.github.io/NeighbourNet/reference/prepare.graph.md).
  The PC embeddings and graph structure must be stored in the object's
  `NNet.setting`.

- p:

  A numeric value between 0 and 1, specifying the proportion of cells to
  select. If `n` is provided, `p` is ignored. Default is 0.1.

- n:

  An integer specifying the total number of cells to select. If `NULL`,
  the number is calculated as `ceiling(p * total cells)`. Default is
  `NULL`.

- all:

  A logical value indicating whether to select all cells. If `TRUE`, all
  cells are used, and no subset selection is performed. Default is
  `FALSE`.

- ...:

  Additional parameters passed to
  [`kmeans`](https://rdrr.io/r/stats/kmeans.html) for clustering.

## Value

A `Seurat` object with updated `NNet.setting`, where the selected cells
are stored as a vector of indexes named by cell barcodes in the `cells`
field of the `NNet.setting` list in the `misc` slot.

## Details

This function provides flexibility in selecting cells for downstream PC
regression analysis:

**Subset Selection**: If `all = FALSE`, a subset of cells is selected
based on k-means clustering in the PC space. This ensures that the
selected cells represent the overall distribution of the data.

**Balanced Sampling**: The function identifies cluster centers using
k-means and selects the nearest neighbors to achieve a balanced
representation.

**Full Dataset Usage**: If `all = TRUE`, no cells are excluded, and the
entire dataset is used.

By default, the function selects 10 users can provide `n` to select an
exact number of cells or use `all = TRUE` to include all cells.

## See also

[`prepare.seurat`](https://meiosis97.github.io/NeighbourNet/reference/prepare.seurat.md),
[`kmeans`](https://rdrr.io/r/stats/kmeans.html)

## Examples

``` r
# Assuming `seurat.obj` is a Seurat object containing `NNet.setting` in its `misc` slot.

# Select 10% of the cells
seurat.obj <- select.cell(seurat.obj, p = 0.1)
#> Error: object 'seurat.obj' not found

# Select a fixed number of cells
seurat.obj <- select.cell(seurat.obj, n = 500)
#> Error: object 'seurat.obj' not found

# Use all cells for downstream analysis
seurat.obj <- select.cell(seurat.obj, all = TRUE)
#> Error: object 'seurat.obj' not found
```
