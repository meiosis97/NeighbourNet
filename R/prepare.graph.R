# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Build a K-Nearest Neighbors (KNN) Graph and Store in Seurat Object
#'
#' This function constructs a K-Nearest Neighbors (KNN) graph based on selected
#' principal components (PCs) from a \code{Seurat} object. The resulting graph
#' and associated settings are stored in the \code{misc} slot of the Seurat object
#' under \code{NNet.setting}.
#'
#' @param seurat.obj A \code{Seurat} object prepared using \code{\link{prepare.seurat}},
#'   which includes PCA results in its \code{NNet.setting}.
#' @param knn An integer specifying the number of nearest neighbors to use for building the KNN graph. Default is 30.
#'
#' @return A \code{Seurat} object with its \code{NNet.setting} updated to include:
#' \item{p}{A sparse cell-by-cell affinity matrix.}
#' \item{nn.idx}{A cell-by-neighbour matrix where each row contains the indices of the nearest neighbors for the corresponding cell.}
#' \item{nn.w}{A matrix describing weights of connections between neighbouring cells described by \code{nn.idx}}
#'
#' @details
#' The function constructs a KNN graph using the most informative PCs
#' learned during the execution of \code{\link{prepare.seurat}}. The graph is built with the
#' \code{\link{build.graph}} function, which computes a sparse affinity matrix, nearest neighbor indices,
#' and weights for each cell. These results are stored in the Seurat object's \code{NNet.setting},
#' which is used for downstream co-expression analysis, including PC regression and network inference.
#'
#' The KNN graph is a critical structure that captures relationships between cells in the reduced
#' PC space, enabling scalable and biologically meaningful analyses.
#'
#' @examples
#' # Select genes for PC regression
#' genes <- select.gene(seurat.obj)$genes
#'
#' # Scale data and perform PCA
#' seurat.obj <- prepare.seurat(seurat.obj, genes = genes)
#'
#' # Build the KNN graph
#' seurat.obj <- prepare.graph(seurat.obj, knn = 30)
#'
#' # Inspect the settings in NNet.setting
#' str(Seurat::Misc(seurat.obj, "NNet.setting"))
#'
#' @seealso
#' \code{\link{build.graph}}, \code{\link{prepare.seurat}}, \code{\link{prepare.reg}}
#'
#' @export
prepare.graph <- function(seurat.obj, knn = 30) {
  # Ensure Seurat object has been prepared with PCs
  setting <- Seurat::Misc(seurat.obj, "NNet.setting")
  if (is.null(setting) || is.null(setting$pcs)) {
    stop("Run prepare.seurat before prepare.graph.")
  }

  # Extract PCs from settings
  pcs <- setting$pcs

  # Build the KNN graph using the selected PCs
  message("Building KNN graph...")
  graph.result <- NeighbourNet::build.graph(pcs, knn = knn)

  # Update the NNet.setting with graph-related structures
  setting$p <- graph.result$p  # Affinity matrix
  setting$nn.idx <- graph.result$nn.idx  # Nearest neighbor indices
  setting$nn.w <- graph.result$nn.w  # Nearest neighbor weights

  # Store updated settings back in the Seurat object
  suppressWarnings(
    Seurat::Misc(seurat.obj, "NNet.setting") <- setting
  )

  # Return the updated Seurat object
  return(seurat.obj)
}
