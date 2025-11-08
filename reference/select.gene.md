# Select Genes for Analysis

This function filters and selects transcription factors (TFs), target
genes, and background genes from a Seurat object based on their
expression in a minimum number of cells.

## Usage

``` r
select.gene(seurat.obj, tfs = NULL, targets = NULL, bgs = NULL, min.cells = 20)
```

## Arguments

- seurat.obj:

  A `Seurat` object containing gene expression data.

- tfs:

  A vector of transcription factor (TF) gene names. If `NULL`, all TFs
  in the
  [`gene.list`](https://meiosis97.github.io/NeighbourNet/reference/gene.list.md)`$tfs`
  will be considered. Default is `NULL`.

- targets:

  A vector of target gene names. If `NULL`, all target genes in
  [`gene.list`](https://meiosis97.github.io/NeighbourNet/reference/gene.list.md)`$targets`
  will be considered. Default is `NULL`.

- bgs:

  A vector of background gene names. If `NULL`, no additional background
  genes are considered. Default is `NULL`.

- min.cells:

  An integer specifying the minimum number of cells in which a gene must
  be expressed to be included. Default is 20.

## Value

A list containing:

- tfs:

  Filtered TFs that meet the minimum cell expression criteria.

- targets:

  Filtered target genes that meet the minimum cell expression criteria.

- bgs:

  Filtered background genes, excluding any genes in `tfs` or `targets`.

- genes:

  A unique set of all selected genes.

## Details

The function first filters genes in the `seurat.obj` to retain only
those expressed in more than `min.cells` cells. It then intersects the
filtered gene set with provided or default lists of TFs, targets, and
background genes.

## See also

[`gr.graph`](https://meiosis97.github.io/NeighbourNet/reference/pkn.md)

## Examples

``` r
select.gene(seurat.obj)
#> Error: object 'seurat.obj' not found
```
