# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Generate a Color Palette for ggplot Visualizations
#'
#' This function creates a reversed Spectral color palette from the \pkg{RColorBrewer} package
#' for use in ggplot or other visualizations. Optionally, a background color can be added.
#'
#' @param n An integer specifying the number of colors to generate. Default is 11.
#' @param background A logical value indicating whether to include a background color
#' (\code{'#D9DDDC'}) as the first element of the palette. Default is \code{FALSE}.
#'
#' @return A character vector of hexadecimal color codes. If \code{background = TRUE},
#' the first color is a light grey background color (\code{'#D9DDDC'}), followed by
#' a reversed Spectral color palette.
#'
#' @details
#' The function generates a reversed Spectral palette using the \code{RColorBrewer::brewer.pal}
#' function. The palette can be customized by specifying the number of colors (\code{n})
#' and optionally including a background color.
#'
#' @examples
#' # Default palette with 11 colors, no background color
#' gg.color.spec()
#'
#' # Palette with 7 colors, including a background color
#' gg.color.spec(n = 7, background = TRUE)
#'
#' # Use the palette in a ggplot
#' library(ggplot2)
#' ggplot(mpg, aes(x = class, fill = class)) +
#'   geom_bar() +
#'   scale_fill_manual(values = gg.color.spec(7))
#'
#' @export
gg.color.spec <- function(n = 11, background = FALSE) {
  # Generate the reversed Spectral palette
  palette <- rev(RColorBrewer::brewer.pal(n = n, name = "Spectral"))

  # Optionally add a background color
  if (background) {
    palette <- c('#D9DDDC', palette)
  }

  return(palette)
}
