# Prepare a Seurat Object for PC Regression Analysis

This function prepares a `Seurat` object for PC regression analysis by
selecting response and predictor genes, as well as the cells to include
in the analysis. It calculates local variances for the selected genes
and PCs based on the KNN graph, which will be used for estimating
permutation feature importance as co-expression measure. The results are
stored in the `NNet.setting` field of the `misc` slot of the Seurat
object.

## Usage

``` r
prepare.reg(
  seurat.obj,
  responses = NULL,
  predictors = NULL,
  cells = NULL,
  check.expressed = FALSE
)
```

## Arguments

- seurat.obj:

  A `Seurat` object on which
  [`prepare.graph`](https://meiosis97.github.io/NeighbourNet/reference/prepare.graph.md)
  has been ran.

- responses:

  A character vector of response gene names. Response gene expression
  will be low-rank approximated. If `NULL`, all genes of PC embedding
  are used. Default is `NULL`.

- predictors:

  A character vector of predictor gene names. Local variance will be
  calculated on the provided predictors. If `NULL`, all genes of PC
  embedding are used. Default is `NULL`.

- cells:

  A vector of cell indices to include. If `NULL`, uses pre-selected
  cells by
  [`select.cell`](https://meiosis97.github.io/NeighbourNet/reference/select.cell.md)
  or all cells if no pre-selection was done. Default is `NULL`.

- check.expressed:

  A logical value indicating whether to prune local variance by if genes
  are expressed in the selected cells. Default is `TRUE`.

## Value

A `Seurat` object with the updated `NNet.setting` stored in the `misc`
slot. The updated settings include:

- responses:

  A character vector of selected response genes.

- predictors:

  A character vector of selected predictor genes.

- genes:

  A character vector of all selected genes (responses and predictors).

- cells:

  A vector of selected cell indices, named by cell barcodes.

- lra:

  A low-rank approximation matrix representing the reconstructed
  expression of response genes based on PCs.

- scale.gene:

  A numeric vector of global gene scales.

- nn.scale.gene:

  A sparse matrix of local variances for the selected genes in each
  selected cell.

- nn.scale.pc:

  A matrix of local variances for the PCs in each selected cell.

- n.eff:

  A numeric vector of effective neighborhood sizes for each selected
  cell.

## Details

This function sets up the data and settings required for performing PC
regression on a Seurat object, with results stored in the `NNet.setting`
field of the `misc` slot. It supports the analysis of gene co-expression
by calculating local variances for genes and PCs, which are essential
for estimating permutation feature importance.

The function allows users to specify response and predictor genes, as
well as subsets of cells, for the analysis. If no responses or
predictors are provided, all genes in the PC loadings are used by
default. Local variances are computed using the KNN graph stored in
`NNet.setting`, enabling downstream co-expression analysis.

## See also

[`prepare.graph`](https://meiosis97.github.io/NeighbourNet/reference/prepare.graph.md),
[`run.nn.reg`](https://meiosis97.github.io/NeighbourNet/reference/run.nn.reg.md)

## Examples

``` r
# Assuming `seurat.obj` is a Seurat object on which `prepare.graph` has been ran.

# seurat.obj <- select.cell(seurat.obj) # Optional if the number of cells is large

# Select genes for PC regression
gene.list <- select.gene(seurat.obj)
#> Error: object 'seurat.obj' not found

# Prepare regression setting, use transcriptional factors as responses and targets as predictors
seurat.obj <- prepare.reg(seurat.obj, responses = gene.list$tfs, predictors = gene.list$targets)
#> Error: object 'seurat.obj' not found

# Check PC regression settings
str(Seurat::Misc(seurat.obj, "NNet.setting"))
#> Error: object 'seurat.obj' not found
```
