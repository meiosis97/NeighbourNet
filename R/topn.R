#' Select indices of the top-n highest values
#'
#' Utility function returning indices of the top-n elements in a numeric vector.
#' Ties are broken randomly.
#'
#' @param x A numeric vector.
#' @param n Number of top elements to select. If \code{n > length(x)}, all indices are returned.
#'
#' @return An integer vector of indices corresponding to the top-n elements in \code{x}.
#'
#' @examples
#' x <- c(5, 1, 3, 9, 7)
#' topn(x, 3)
#'
#' @export
topn <- function(x, n) {
  x <- as.numeric(x)
  n <- min(n, length(x))
  if (anyNA(x)) x[is.na(x)] <- -Inf
  which(rank(-x, ties.method = "random") <= n)
}
