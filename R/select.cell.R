# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Select Cells for Downstream PC Regression Analysis
#'
#' This function selects a subset of cells from a \code{Seurat} object for downstream analysis
#' based on their principal component (PC) embeddings. The selection ensures that the subset
#' represents the distribution of the full dataset while maintaining computational efficiency.
#'
#' @param seurat.obj A \code{Seurat} object that has been processed with \code{\link{prepare.graph}}.
#'   The PC embeddings and graph structure must be stored in the object's \code{NNet.setting}.
#' @param p A numeric value between 0 and 1, specifying the proportion of cells to select.
#'   If \code{n} is provided, \code{p} is ignored. Default is 0.1.
#' @param n An integer specifying the total number of cells to select. If \code{NULL},
#'   the number is calculated as \code{ceiling(p * total cells)}. Default is \code{NULL}.
#' @param all A logical value indicating whether to select all cells. If \code{TRUE},
#'   all cells are used, and no subset selection is performed. Default is \code{FALSE}.
#' @param ... Additional parameters passed to \code{\link[stats]{kmeans}} for clustering.
#'
#' @return A \code{Seurat} object with updated \code{NNet.setting}, where the selected cells
#'   are stored as a vector of indexes named by cell barcodes in the \code{cells} field of the \code{NNet.setting} list in the \code{misc} slot.
#'
#' @details
#' This function provides flexibility in selecting cells for downstream PC regression analysis:
#'
#' - **Subset Selection**: If \code{all = FALSE}, a subset of cells is selected based on k-means clustering
#'   in the PC space. This ensures that the selected cells represent the overall distribution of the data.
#' - **Balanced Sampling**: The function identifies cluster centers using k-means and selects
#'   the nearest neighbors to achieve a balanced representation.
#' - **Full Dataset Usage**: If \code{all = TRUE}, no cells are excluded, and the entire dataset is used.
#'
#' By default, the function selects 10% of the cells (\code{p = 0.1}). For large datasets,
#' users can provide \code{n} to select an exact number of cells or use \code{all = TRUE} to include all cells.
#'
#' @examples
#' # Assuming `seurat.obj` is a Seurat object containing `NNet.setting` in its `misc` slot.
#'
#' # Select 10% of the cells
#' seurat.obj <- select.cell(seurat.obj, p = 0.1)
#'
#' # Select a fixed number of cells
#' seurat.obj <- select.cell(seurat.obj, n = 500)
#'
#' # Use all cells for downstream analysis
#' seurat.obj <- select.cell(seurat.obj, all = TRUE)
#'
#' @seealso \code{\link{prepare.seurat}}, \code{\link[stats]{kmeans}}
#'
#' @export
select.cell <- function(seurat.obj, p = 0.1, n = NULL, all = FALSE, ...) {
  # Retrieve the stored NNet.setting from the Seurat object
  setting <- Seurat::Misc(seurat.obj, "NNet.setting")

  # Ensure the graph preparation step has been performed
  if (is.null(setting)) stop("Run prepare.seurat first.")

  if (all) {
    # Use all cells if `all` is TRUE
    setting$cells <- NULL

  } else {
    # Calculate the number of cells to select if `n` is not provided
    if (is.null(n)) {
      n <- nrow(setting$pcs)
      n <- ceiling(n * p)
    }

    # Use half the calculated number of cells as k-means cluster centers
    n <- ceiling(n / 2)

    # Perform k-means clustering to identify representative cluster centers
    centers <- stats::kmeans(setting$pcs, centers = n, iter.max = 100, ...)$centers

    # Find the nearest neighbors for each cluster center
    cells <- RANN::nn2(setting$pcs, centers, k = 2)$nn.idx %>% as.numeric()

    # Remove duplicate cell indices and assign cell names
    cells <- unique(cells)
    names(cells) <- rownames(setting$pcs)[cells]

    # Store the selected cell indices
    setting$cells <- cells
  }

  # Update the Seurat object's NNet.setting
  suppressWarnings(
    Seurat::Misc(seurat.obj, "NNet.setting") <- setting
  )

  # Return the updated Seurat object
  return(seurat.obj)
}
