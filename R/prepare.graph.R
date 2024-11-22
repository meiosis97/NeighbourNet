# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Build KNN Graph and Store Settings in Seurat Object
#'
#' This function builds a k-nearest neighbors (KNN) graph based on selected
#' principal components (PCs) from a Seurat object and stores the resulting
#' graph and related PC regression settings in the object's metadata.
#'
#' @param seurat.obj A \code{Seurat} object containing PCA results.
#' @param knn An integer specifying the number of nearest neighbors for the KNN graph. Default is 30.
#' @param truncated A logical value indicating whether to select a subset of important PCs
#' based on their standard deviation. If \code{TRUE}, only significant PCs are used. Default is \code{TRUE}.
#'
#' @return A \code{Seurat} object with the KNN graph and related PC regression settings stored in its \code{misc} slot,
#' as a list named \code{NNet.setting}. Refer to \code{\link{prepare.reg}} for the detailed description of the regression settings.
#'
#' @details
#' The function first selects the most informative principal components (PCs) based
#' on their standard deviations if \code{truncated = TRUE}. A KNN graph is then
#' constructed using \code{\link{build.graph}} with these PCs. The PC regression settings,
#' including PCA results and neighborhood indices, are stored in the Seurat object's \code{misc} slot.
#'
#' @examples
#' # Select genes for PC regression
#' genes <- select.gene(seurat.obj)$genes
#'
#' # Scale data and perform PCA
#' seurat.obj <- prepare.seurat(seurat.obj, genes = genes)
#'
#' # Prepare the KNN graph
#' seurat.obj <- prepare.graph(seurat.obj)
#'
#' # Check PC regression settings
#' str(Seurat::Misc(seurat.obj, "NNet.setting"))
#'
#' @seealso \code{\link{build.graph}}, \code{\link{prepare.reg}}
#'
#' @export
prepare.graph <- function(seurat.obj, knn = 30, truncated = TRUE) {
  # Extract standard deviations of principal components
  sd <- Seurat::Reductions(seurat.obj, "pca")@stdev

  # Determine the number of significant PCs
  npcs <- if (truncated) {
    find.significant.pcs(sd)
  } else {
    length(sd)  # Use all PCs if no truncation
  }

  # Extract embeddings and loadings for the selected PCs
  pcs <- Seurat::Embeddings(seurat.obj, "pca")[, 1:npcs]
  loadings <- Seurat::Reductions(seurat.obj, "pca")@feature.loadings[, 1:npcs]

  # Calculate sparsity
  genes <- rownames(loadings)
  sparsity <- mean(SeuratObject::LayerData(seurat.obj, "counts")[genes,] == 0)

  # Build the KNN graph using the selected PCs
  message("Building KNN graph...")
  graph.result <- NeighbourNet::build.graph(pcs, knn = knn)

  # Create a settings list to store in the Seurat object
  setting <- list(
    pcs = pcs,
    loadings = loadings,
    p = graph.result$p,
    nn.idx = graph.result$nn.idx,
    nn.w = graph.result$nn.w,
    sparsity = sparsity,
    cells = NULL,
    predictors = NULL,
    responses = NULL,
    genes = NULL,
    lra = NULL,
    nn.scale.gene = NULL,
    nn.scale.pc = NULL,
    n.eff = NULL
  )

  # Set a class
  class(setting) <- "NNet.setting"

  # Store settings in Seurat object metadata
  suppressWarnings(
    Seurat::Misc(seurat.obj, "NNet.setting") <- setting
  )

  # Return the updated Seurat object
  return(seurat.obj)
}

#' Find Significant Principal Components
#'
#' This helper function identifies the number of significant principal components (PCs)
#' based on their standard deviations. It uses a heuristic approach to determine
#' where the standard deviations significantly change. This method is based on the work of Linderman et al. (2022).
#'
#' @references
#' Linderman, G. C., Zhao, J., Roulis, M., Bielecki, P., Flavell, R. A., Nadler, B., & Kluger, Y. (2022).
#' Zero-preserving imputation of single-cell RNA-seq data.
#' \emph{Nature Communications}, 13(1), 192. \doi{10.1038/s41467-021-27923-7}
#'
#' @param sd A numeric vector of standard deviations of PCs.
find.significant.pcs <- function(sd) {
  k <- length(sd)
  s <- abs(diff(sd))
  mu <- mean(s[(0.8 * k):(k - 1)])
  sigma <- sd(s[(0.8 * k):(k - 1)])
  sk <- mu + 6 * sigma
  npcs <- max(which(s > sk)) + 1

  # Handle case where no significant PCs are identified
  if (length(npcs) == 0) {
    warning("No significant components identified. Using all components.")
    return(k)
  }

  return(npcs)
}
