# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Prepare a Seurat Object for PC Regression Analysis
#'
#' This function prepares a \code{Seurat} object for PC regression analysis by selecting response and predictor genes,
#' as well as the cells to include in the analysis. It calculates local variances for the selected genes and PCs
#' based on the KNN graph, which will be used for estimating permutation feature importance as co-expression measure.
#' The results are stored in the \code{NNet.setting} field of  the \code{misc} slot of the Seurat object.
#'
#'
#' @param seurat.obj A \code{Seurat} object on which \code{\link{prepare.graph}} has been ran.
#' @param responses A character vector of response gene names. Response gene expression will be low-rank approximated.
#' If \code{NULL}, all genes of PC embedding are used. Default is \code{NULL}.
#' @param predictors A character vector of predictor gene names. Local variance will be calculated on the provided predictors.
#' If \code{NULL}, all genes of PC embedding are used. Default is \code{NULL}.
#' @param cells A vector of cell indices to include. If \code{NULL}, uses pre-selected cells by \code{\link{select.cell}} or all cells if no
#' pre-selection was done. Default is \code{NULL}.
#' @param check.expressed A logical value indicating whether to prune local variance by if genes are expressed in the selected cells. Default is \code{TRUE}.
#'
#' @return A \code{Seurat} object with the updated \code{NNet.setting} stored in the \code{misc} slot.
#'   The updated settings include:
#'   \item{responses}{A character vector of selected response genes.}
#'   \item{predictors}{A character vector of selected predictor genes.}
#'   \item{genes}{A character vector of all selected genes (responses and predictors).}
#'   \item{cells}{A vector of selected cell indices, named by cell barcodes.}
#'   \item{lra}{A low-rank approximation matrix representing the reconstructed expression of response genes based on PCs.}
#'   \item{nn.scale.gene}{A sparse matrix of local variances for the selected genes in each selected cell.}
#'   \item{nn.scale.pc}{A matrix of local variances for the PCs in each selected cell.}
#'   \item{n.eff}{A numeric vector of effective neighborhood sizes for each selected cell.}
#'
#' @details
#' This function sets up the data and settings required for performing PC regression
#' on a Seurat object, with results stored in the \code{NNet.setting} field of the \code{misc} slot.
#' It supports the analysis of gene co-expression by calculating local variances
#' for genes and PCs, which are essential for estimating permutation feature importance.
#'
#' The function allows users to specify response and predictor genes, as well as subsets of cells, for the analysis.
#' If no responses or predictors are provided, all genes in the PC loadings are used by default.
#' Local variances are computed using the KNN graph stored in \code{NNet.setting}, enabling downstream co-expression analysis.
#'
#' @examples
#' # Assuming `seurat.obj` is a Seurat object on which `prepare.graph` has been ran.
#'
#' # seurat.obj <- select.cell(seurat.obj) # Optional if the number of cells is large
#'
#' # Select genes for PC regression
#' gene.list <- select.gene(seurat.obj)
#'
#' # Prepare regression setting, use transcriptional factors as responses and targets as predictors
#' seurat.obj <- prepare.reg(seurat.obj, responses = gene.list$tfs, predictors = gene.list$targets)
#'
#' # Check PC regression settings
#' str(Seurat::Misc(seurat.obj, "NNet.setting"))
#'
#' @seealso \code{\link{prepare.graph}}, \code{\link{run.nn.reg}}
#'
#' @export
prepare.reg <- function(seurat.obj, responses = NULL, predictors = NULL, cells = NULL, check.expressed = TRUE) {

  # Retrieve the stored settings from the Seurat object
  setting <- Seurat::Misc(seurat.obj, "NNet.setting")

  # Ensure that prepare.seurat has been run
  if (is.null(setting)) stop("Run prepare.seurat first, then prepare.graph and prepare.reg.")

  # Ensure that prepare.graph has been run
  if (is.null(setting$p)) stop("Run prepare.graph first and then prepare.reg.")

  # Select response genes
  # If responses are not provided, use all available genes in the PCA loadings.
  if (is.null(responses)) {
    responses <- rownames(setting$loadings)
  } else {
    responses <- intersect(responses, rownames(setting$loadings))
  }

  # Select predictor genes
  # If predictors are not provided, use all available genes in the PCA loadings.
  if (is.null(predictors)) {
    predictors <- rownames(setting$loadings)
  } else {
    predictors <- intersect(predictors, rownames(setting$loadings))
  }

  # Combine responses and predictors into a unique set of genes
  genes <- unique(c(responses, predictors))

  # Select cells for analysis
  if (is.null(cells)) {
    # Use pre-selected cells if available, otherwise use all cells
    cells <- setting$cells
    if (is.null(cells)) {
      cells <- 1:nrow(setting$pcs)
      names(cells) <- rownames(setting$pcs)
    }
  } else {
    # Map the provided cell indices to their cell names
    names(cells) <- rownames(setting$pcs)[cells]
    # Store selected cells in settings
    setting$cells <- cells
  }

  # Prepare data structures for regression
  pcs <- setting$pcs[cells, ]  # Select PCs for the selected cells
  setting$responses <- responses
  setting$predictors <- predictors
  setting$genes <- genes

  n.cell <- length(cells)
  n.gene <- length(genes)
  n.pc <- ncol(pcs)

  # Low-rank approximation of selected responses using PCs
  setting$lra <- tcrossprod(setting$pcs, setting$loadings[responses, ])

  # Initialize matrices for local variances
  message("Calculating local variance.")
  nn.scale.gene <- matrix(0, nrow = n.gene, ncol = n.cell,
                          dimnames = list(genes, names(cells)))
  nn.scale.pc <- matrix(0, nrow = n.pc, ncol = n.cell,
                        dimnames = dimnames(pcs) %>% rev)

  # Effective neighborhood size for variance calculation
  n.eff <- apply(setting$nn.w[cells, ], 1, function(x) (sum(x)^2) / sum(x^2))
  names(n.eff) <- names(cells)

  # Retrieve scaled data for selected genes and cells
  scale.data <- SeuratObject::LayerData(seurat.obj, layer = "scale.data")[genes, names(cells)]

  # Loop through each selected cell to calculate local variances
  for (i in 1:n.cell) {
    idx <- setting$nn.idx[cells[i], ]  # Nearest neighbors for the current cell
    w <- setting$nn.w[cells[i], ]      # Weights for the neighbors

    # Local gene scales
    w.mean <- as.numeric(scale.data[, idx] %*% w)  # Weighted mean for genes
    res <- scale.data[, idx] - w.mean              # Residuals
    nn.scale.gene[, i] <- as.numeric(res^2 %*% w) * n.eff[i] / (n.eff[i] - 1)

    # Local PC scales
    w.mean <- as.numeric(w %*% setting$pcs[idx, ]) # Weighted mean for PCs
    res <- t(setting$pcs[idx, ]) - w.mean          # Residuals for PCs
    nn.scale.pc[, i] <- as.numeric(res^2 %*% w) * n.eff[i] / (n.eff[i] - 1)
  }

  # Optionally filter for expressed genes to exclude zero counts
  if (check.expressed) {
    nn.scale.gene <- nn.scale.gene *
      (SeuratObject::LayerData(seurat.obj, layer = "counts")[genes, names(cells)] != 0)
  }

  # Store calculated variances and effective neighborhood size
  setting$nn.scale.gene <- Matrix::Matrix(sqrt(nn.scale.gene))  # Gene-level variances
  setting$nn.scale.pc <- sqrt(nn.scale.pc)                      # PC-level variances
  setting$n.eff <- n.eff                                        # Effective neighborhood sizes

  # Update the Seurat object with the new settings
  suppressWarnings(
    Seurat::Misc(seurat.obj, "NNet.setting") <- setting
  )

  return(seurat.obj)  # Return the updated Seurat object
}
