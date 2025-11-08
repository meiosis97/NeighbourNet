# Build Meta-Networks (Meta-response) via Network Embeddings (PCA / Non-Negative PCA)

`build.meta.network` aggregates per-cell networks into a low-dimensional
set of meta-networks by network embedding (nPCA on vectorised networks).
Concatenating network vectors yields an edge-by-cell “long-and-tall”
matrix. For efficiency, the edge dimension is first embedded by SVD/PCA,
and then nPCA is applied to the embedded matrix. It accepts either a
`Seurat` object prepared by the NeighbourNet workflow or a precomputed
3D network tensor. Similarly, using `build.meta.response`,
meta-responses representing functional gene modules can be built by nPCA
embedding on vectorised response profiles
(cell\*predictors-by-response). The detailed arguments and usage of
`build.meta.network` is documented, and that of `build.meta.response` is
alike.

## Usage

``` r
build.meta.network(
  seurat.obj = NULL,
  network = NULL,
  k = 100,
  cutoff = NULL,
  big.memory = FALSE,
  scale = TRUE,
  truncated = TRUE,
  n.net = 20,
  non.neg = TRUE,
  max.iter = 1000,
  tol = 1e-10,
  return.p.val = TRUE
)

build.meta.response(
  seurat.obj = NULL,
  network = NULL,
  k = 100,
  cutoff = NULL,
  big.memory = FALSE,
  scale = TRUE,
  truncated = TRUE,
  n.net = 20,
  non.neg = TRUE,
  max.iter = 1000,
  tol = 1e-10,
  return.p.val = TRUE
)
```

## Arguments

- seurat.obj:

  A `Seurat` object with a `NNet.mod` list stored in the `misc` slot.
  This list is created by
  [`run.nn.reg`](https://meiosis97.github.io/NeighbourNet/reference/run.nn.reg.md).
  If provided, `assay` and gene panels are taken from
  `NNet.mod$defaults`. Ignored when `network` is supplied.

- network:

  Optional 3D tensor of networks with dimension order (responses x
  predictors x cells). If supplied, `return.p.val` is forced to `FALSE`.

- k:

  A integer specifying ranks to keep for the “long-and-tall” matrix of
  concatenated network vectors. Limited to
  `min(k, n.response * n.predictor, n.cell)`. Default: `100`.

- cutoff:

  A numeric value specifying the p-value threshold (in `[0, 1]`) passed
  to
  [`get.network`](https://meiosis97.github.io/NeighbourNet/reference/get.network.md)
  when pulling per-cell networks. `NULL` uses defaults from
  `NNet.mod$defaults$cutoff`.

- big.memory:

  A logical indicating whether to load the full network tensor
  internally to reduce repeated reads. Default: `FALSE`.

- scale:

  A logical indicating whether to scale each single-cell network by its
  Frobenius norm across edges before aggregation (per-cell
  normalization). Default: `TRUE`.

- truncated:

  A logical indicating whether to retain only significant edge
  embeddings (spectral-gap heuristic) before running nPCA; otherwise,
  use exactly `k` embeddings. Default: `TRUE`. If `k < 100`, this is set
  to `FALSE`.

- n.net:

  A integer specifying the number of meta-networks to output. Default:
  `20`.

- non.neg:

  A logical indicating whether to perform nPCA (iterative deflation with
  non-negativity constraints) on the edge-embedded matrix, then embeds
  cell networks with nPCA loadings. Otherwise, embeds with standard PCA
  loadings. Default: `TRUE`.

- max.iter:

  A integer specifying the maximum iterations for calculating each nPCA
  component. Default: `1000`.

- tol:

  A numeric value specifying the convergence tolerance for nPCA
  iterations. Default: `1e-10`.

- return.p.val:

  Logical; if `TRUE` and `assay == "effect"`, also returns
  per–meta-network p-values by embedding the p-value tensor
  (`NNet.mod$p.val`) with the nPCA loadings. Forced to `FALSE` when
  `network` is supplied or when `non.neg = FALSE`. Default: `TRUE`.

## Value

If Seurat object is provided, then update it with
`NNet.mod$meta.network` populated, which is, a list containing:

- `meta.network`: Array of size (responses x predictors x `n.net`).

- `p.val`: An Array (same size) if `return.p.val = TRUE`; else `NULL`.

- `pcs`: A matrix of PC scores learnt on per-cell networks; rows =
  cells, columns = components.

- `pca.loadings`: A matrix of PCA loadings learnt on per-cell networks.

- `pca.sd`: A numeric vector of singular values from network covariance.

- `npca.loadings`: A matrix of PC scores learnt on per-cell networks; if
  `non.neg = TRUE`; else `NULL`.

- `scale`, `non.neg`: Echo of the input flags.

- `setting`: Copy of `NNet.mod$defaults` if `seurat.obj` was provided.

## Details

**Edge embedding**: If both `responses` and `predictors` have size \> 1,
a cell-by-cell kernel that represents co-variation in network structure
across cells is computed. Eigen-decomposition is ran to construct a
edge-embedding x cell matrix. If one axis is size 1, perform SVD
directly on the (edges × cells) matrix.

**Rank selection**: If `truncated = TRUE`, a simple spectral-gap rule
over the tail of singular-value differences selects the working rank of
the edge emebdding before nPCA.

**Meta-networks**: Project per-cell network onto the first `n.net`
directions (nPCA or PCA) to obtain meta-network slices.

## See also

[`run.nn.reg`](https://meiosis97.github.io/NeighbourNet/reference/run.nn.reg.md),
[`get.network`](https://meiosis97.github.io/NeighbourNet/reference/get.network.md)

## Examples

``` r
# From a Seurat object (built by run.nn.reg)
seurat.obj <- build.meta.network(seurat.obj, k = 100, n.net = 20, non.neg = TRUE)
#> Error: object 'seurat.obj' not found
```
