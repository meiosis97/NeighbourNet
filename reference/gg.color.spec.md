# Generate a Color Palette for ggplot Visualizations

This function creates a reversed Spectral color palette from the
RColorBrewer package for use in ggplot or other visualizations.
Optionally, a background color can be added.

## Usage

``` r
gg.color.spec(n = 11, background = FALSE)
```

## Arguments

- n:

  An integer specifying the number of colors to generate. Default is 11.

- background:

  A logical value indicating whether to include a background color
  (`'#D9DDDC'`) as the first element of the palette. Default is `FALSE`.

## Value

A character vector of hexadecimal color codes. If `background = TRUE`,
the first color is a light grey background color (`'#D9DDDC'`), followed
by a reversed Spectral color palette.

## Details

The function generates a reversed Spectral palette using the
[`RColorBrewer::brewer.pal`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html)
function. The palette can be customized by specifying the number of
colors (`n`) and optionally including a background color.

## Examples

``` r
# Default palette with 11 colors, no background color
gg.color.spec()
#>  [1] "#5E4FA2" "#3288BD" "#66C2A5" "#ABDDA4" "#E6F598" "#FFFFBF" "#FEE08B"
#>  [8] "#FDAE61" "#F46D43" "#D53E4F" "#9E0142"

# Palette with 7 colors, including a background color
gg.color.spec(n = 7, background = TRUE)
#> [1] "#D9DDDC" "#3288BD" "#99D594" "#E6F598" "#FFFFBF" "#FEE08B" "#FC8D59"
#> [8] "#D53E4F"
```
