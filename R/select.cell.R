# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Select Cells for Downstream Analysis
#'
#' This function selects a subset of cells from a Seurat object for PC regression
#' based on principal component (PC) embedding.
#'
#' @param seurat.obj A \code{Seurat} object on which \code{\link{prepare.graph}} has been ran.
#' @param p A numeric value between 0 and 1, specifying the proportion of cells to select. Default is 0.1.
#' @param n An integer specifying the total number of cells to select. If \code{NULL},
#' it is calculated based on \code{p}. Default is \code{NULL}.
#' @param all A logical value indicating whether to select all cells. If \code{TRUE},
#' no cells are excluded. Default is \code{FALSE}.
#' @param ... Passed to \code{\link[stats]{kmeans}}
#'
#' @return A \code{Seurat} object with the selected cells stored in the \code{NNet.setting} in Seurat objects' \code{misc} slot.
#'
#' @details
#' If \code{all = TRUE}, all cells are used for downstream analysis. Otherwise,
#' the function uses k-means clustering using \code{\link[stats]{kmeans}} on the PC embedding to identify a subset of cells
#' that represent the data distribution. Nearest neighbors are computed to ensure
#' a balanced selection of cells.
#'
#' @examples
#' # Assuming `seurat.obj` is a Seurat object containing `NNet.setting` in its `misc` slot.
#' seurat.obj <- select.cell(seurat.obj, p = 0.1)
#'
#' # Select a fixed number of cells:
#' seurat.obj <- select.cell(seurat.obj, n = 500)
#'
#' # Select all cells:
#' seurat.obj <- select.cell(seurat.obj, all = TRUE)
#'
#' @seealso \code{\link{prepare.graph}}, \code{\link[stats]{kmeans}}
#'
#' @export
select.cell <- function(seurat.obj, p = 0.1, n = NULL, all = FALSE, ...) {
  # Retrieve the stored settings from the Seurat object
  setting <- Seurat::Misc(seurat.obj, "NNet.setting")

  # Ensure the graph preparation step has been performed
  if (is.null(setting)) stop("Run prepare.graph first.")

  if (all) {
    # Use all cells if `all` is TRUE
    setting["cells"] <- list(cells = NULL)

  } else {
    # Calculate the number of cells to select if not provided
    if (is.null(n)) {
      n <- nrow(setting$pcs)
      n <- ceiling(n * p)
    }

    # Use half the calculated cells as cluster centers
    n <- ceiling(n / 2)

    # Perform k-means clustering to select representative cells
    centers <- stats::kmeans(setting$pcs, centers = n, iter.max = 100, ...)$centers

    # Find the nearest neighbors of the cluster centers
    cells <- RANN::nn2(setting$pcs, centers, k = 2)$nn.idx %>% as.numeric()

    # Remove duplicate cell indices and map them to row names
    cells <- unique(cells)
    names(cells) <- rownames(setting$pcs)[cells]

    # Store the selected cell indices
    setting$cells <- cells
  }

  # Update the Seurat object metadata
  suppressWarnings(
    Seurat::Misc(seurat.obj, "NNet.setting") <- setting
  )

  # Return the updated Seurat object
  return(seurat.obj)
}
