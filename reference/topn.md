# Select indices of the top-n highest values

Utility function returning indices of the top-n elements in a numeric
vector. Ties are broken randomly.

## Usage

``` r
topn(x, n)
```

## Arguments

- x:

  A numeric vector.

- n:

  Number of top elements to select. If `n > length(x)`, all indices are
  returned.

## Value

An integer vector of indices corresponding to the top-n elements in `x`.

## Examples

``` r
x <- c(5, 1, 3, 9, 7)
topn(x, 3)
#> [1] 1 4 5
```
