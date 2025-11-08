# Find Significant Principal Components

This helper function identifies the number of significant principal
components (PCs) based on their standard deviations. It uses a heuristic
approach to determine where the standard deviations significantly
change. This method is based on the work of Linderman et al. (2022).

## Usage

``` r
find.significant.pcs(sd)
```

## Arguments

- sd:

  A numeric vector of standard deviations of PCs.

## References

Linderman, G. C., Zhao, J., Roulis, M., Bielecki, P., Flavell, R. A.,
Nadler, B., & Kluger, Y. (2022). Zero-preserving imputation of
single-cell RNA-seq data. *Nature Communications*, 13(1), 192.
[doi:10.1038/s41467-021-27923-7](https://doi.org/10.1038/s41467-021-27923-7)
