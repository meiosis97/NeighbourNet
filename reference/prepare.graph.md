# Build a K-Nearest Neighbors (KNN) Graph and Store in Seurat Object

This function constructs a K-Nearest Neighbors (KNN) graph based on
selected principal components (PCs) from a `Seurat` object. The
resulting graph and associated settings are stored in the `misc` slot of
the Seurat object under `NNet.setting`.

## Usage

``` r
prepare.graph(seurat.obj, knn = 30)
```

## Arguments

- seurat.obj:

  A `Seurat` object prepared using
  [`prepare.seurat`](https://meiosis97.github.io/NeighbourNet/reference/prepare.seurat.md),
  which includes PCA results in its `NNet.setting`.

- knn:

  An integer specifying the number of nearest neighbors to use for
  building the KNN graph. Default is 30.

## Value

A `Seurat` object with its `NNet.setting` updated to include:

- p:

  A sparse cell-by-cell affinity matrix.

- nn.idx:

  A cell-by-neighbour matrix where each row contains the indices of the
  nearest neighbors for the corresponding cell.

- nn.w:

  A matrix describing weights of connections between neighbouring cells
  described by `nn.idx`

## Details

The function constructs a KNN graph using the most informative PCs
learned during the execution of
[`prepare.seurat`](https://meiosis97.github.io/NeighbourNet/reference/prepare.seurat.md).
The graph is built with the
[`build.graph`](https://meiosis97.github.io/NeighbourNet/reference/build.graph.md)
function, which computes a sparse affinity matrix, nearest neighbor
indices, and weights for each cell. These results are stored in the
Seurat object's `NNet.setting`, which is used for downstream
co-expression analysis, including PC regression and network inference.

The KNN graph is a critical structure that captures relationships
between cells in the reduced PC space, enabling scalable and
biologically meaningful analyses.

## See also

[`build.graph`](https://meiosis97.github.io/NeighbourNet/reference/build.graph.md),
[`prepare.seurat`](https://meiosis97.github.io/NeighbourNet/reference/prepare.seurat.md),
[`prepare.reg`](https://meiosis97.github.io/NeighbourNet/reference/prepare.reg.md)

## Examples

``` r
# Select genes for PC regression
genes <- select.gene(seurat.obj)$genes
#> Error: object 'seurat.obj' not found

# Scale data and perform PCA
seurat.obj <- prepare.seurat(seurat.obj, genes = genes)
#> Error: object 'genes' not found

# Build the KNN graph
seurat.obj <- prepare.graph(seurat.obj, knn = 30)
#> Error: object 'seurat.obj' not found

# Inspect the settings in NNet.setting
str(Seurat::Misc(seurat.obj, "NNet.setting"))
#> Error: object 'seurat.obj' not found
```
