# Prepare a Seurat Object with Selected Genes and PCA

This function scales the data and performs Principal Component Analysis
(PCA) on a `Seurat` object using a specified set of genes. A placeholder
for PC regression settings for co-expression inference is created and
stored in the Seurat object's `misc` slot under `NNet.setting`.

## Usage

``` r
prepare.seurat(
  seurat.obj,
  genes,
  npcs = 100,
  truncated = TRUE,
  ScaleData.ctrl = list(),
  RunPCA.ctrl = list()
)
```

## Arguments

- seurat.obj:

  A `Seurat` object containing a `data` layer.

- genes:

  A character vector of gene names to use for scaling and PCA. Only
  genes present in `seurat.obj` will be used.

- npcs:

  An integer specifying the maximum number of principal components to
  compute. Default is `100`.

- truncated:

  A logical value indicating whether to select a subset of significant
  PCs based on their standard deviation. If `TRUE`, only significant PCs
  are used. Default is `TRUE`.

- ScaleData.ctrl:

  A list of additional parameters to pass to
  [`ScaleData`](https://satijalab.org/seurat/reference/ScaleData.html).
  The `features` and `verbose` arguments are set automatically.

- RunPCA.ctrl:

  A list of additional parameters to pass to
  [`RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html). The
  `features`, `npcs`, and `verbose` arguments are set automatically.

## Value

A `Seurat` object with PC regression settings stored in its `misc` slot
as a list named `NNet.setting`. The list includes:

- pcs:

  A matrix of PC embeddings for all cells. Rows correspond to cells, and
  columns correspond to PCs.

- loadings:

  A matrix of PC loadings. Rows correspond to genes, and columns
  correspond to PCs.

- p:

  Placeholder for the KNN affinity matrix.

- nn.idx:

  Placeholder for the indices of nearest neighbors.

- nn.w:

  Placeholder for the weights of nearest neighbors.

- all.cells:

  A character vector of cells used to ran PCA.

- all.genes:

  A character vector of genes used to ran PCA.

- cells:

  Placeholder for selected cells.

- predictors:

  Placeholder for selected predictor genes.

- responses:

  Placeholder for selected response genes.

- genes:

  Placeholder for the set of selected genes.

- lra:

  Placeholder for the low-rank approximation matrix.

- scale.gene:

  Placeholder for global gene scales.

- nn.scale.gene:

  Placeholder for gene-level variance scaling.

- nn.scale.pc:

  Placeholder for PC-level variance scaling.

- n.eff:

  Placeholder for effective neighborhood sizes.

## Details

This function prepares a Seurat object for PC regression analysis. It
filters the provided gene list to include only those present in the
Seurat object, scales the data using
[`ScaleData`](https://satijalab.org/seurat/reference/ScaleData.html),
and performs PCA using
[`RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html). If
`truncated = TRUE`, the number of significant PCs is determined using
their standard deviations. Users can customize the scaling and PCA
processes by providing additional control parameters via
`ScaleData.ctrl` and `RunPCA.ctrl`.

The prepared Seurat object is stored with placeholders for regression
settings to enable downstream co-expression analysis using `NNet`.

## See also

[`ScaleData`](https://satijalab.org/seurat/reference/ScaleData.html),
[`RunPCA`](https://satijalab.org/seurat/reference/RunPCA.html)

## Examples

``` r
# Select genes for PC regression
genes <- select.gene(seurat.obj)$genes
#> Error: object 'seurat.obj' not found

# Prepare the Seurat object with default settings
seurat.obj <- prepare.seurat(seurat.obj, genes, npcs = 30)
#> Error: object 'genes' not found

# Customize scaling and PCA parameters
seurat.obj <- prepare.seurat(seurat.obj, genes, npcs = 50,
                             ScaleData.ctrl = list(do.center = FALSE))
#> Error: object 'genes' not found
```
