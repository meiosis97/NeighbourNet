# Non-Negative PCA with Deflation (nPCA)

Computes a non-negative PCA (nPCA) basis via iterative regression and
deflation on an input matrix `x`.

## Usage

``` r
npca(x, k, tol = 1e-10, max.iter = 1000, return.score = TRUE)
```

## Arguments

- x:

  A numeric matrix of size (p x n). The row dimension will be embedded

- k:

  Integer; number of nPCA components to extract.

- tol:

  Numeric; convergence tolerance for nPCA iterations. Default: `1e-10`.

- max.iter:

  Integer; maximum iterations for calculating each nPCA component.
  Default: `1000`.

- return.score:

  Logical; If `TRUE` return the nPCA projection of `x`.

## Value

A list with:

- `scores`: Matrix (n x k) of non-negative nPCA scores (columns are
  components).

- `loadings`: Matrix (p x k) of non-negative nPCA loadings

- `sd`: Numeric vector (`k`) of nPCA “standard deviations”.

## Details

This function implements non-negative PCA using restricted single-value
decomposition, which is solved using the iterative regression approach
introduced by Sigg and Buhmann (2008) and deflation approach introduced
by Mackey (2008).

## References

Sigg, C. D., & Buhmann, J. M. (2008, July). Expectation-maximization for
sparse and non-negative PCA. *In Proceedings of the 25th international
conference on Machine learning* (pp. 960-967).

Mackey, L. (2008). Deflation methods for sparse PCA. Advances in neural
information processing systems, 21.

## Examples

``` r
x <- matrix(rnorm(1000), nrow - 100)
#> Error in nrow - 100: non-numeric argument to binary operator
res <- nn.pca(x, k = 10, tol = 1e-10, n.iter = 1000)
#> Error in nn.pca(x, k = 10, tol = 1e-10, n.iter = 1000): could not find function "nn.pca"
str(res$loadings); str(res$sd)
#> Error: object 'res' not found
```
