# Build Graph Based on Principal Components and K-Nearest Neighbors

Constructs an affinity graph from principal components using adaptive
kernel scaling approach adapted from Uniform Manifold Approximation and
Projection (UMAP) (Becht et al. 2019).

## Usage

``` r
build.graph(pcs, knn = 30)
```

## Arguments

- pcs:

  A numeric matrix of principal components (cells x PCs).

- knn:

  Number of nearest neighbors to consider for each cell. Default is 30.

## Value

A list containing:

- p:

  A sparse symmetric affinity matrix representing transition
  probability.

- nn.idx:

  Matrix of nearest neighbor indices for each cell.

- nn.w:

  Matrix of normalized affinity weights for the nearest neighbors.

## References

Becht, E., McInnes, L., Healy, J., Dutertre, C.-A., Kwok, I. W. H., Ng,
L. G., Ginhoux, F., & Newell, E. W. (2019). Dimensionality reduction for
visualizing single-cell data using UMAP. *Nature Biotechnology*, 37(1),
38–44. [doi:10.1038/nbt.4314](https://doi.org/10.1038/nbt.4314)
