# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Prepare a Seurat Object with Selected Genes and PCA
#'
#' This function scales the data and performs Principal Component Analysis (PCA) on a \code{Seurat} object
#' using a specified set of genes. A placeholder for PC regression settings for co-expression inference
#' is created and stored in the Seurat object's \code{misc} slot under \code{NNet.setting}.
#'
#' @param seurat.obj A \code{Seurat} object containing a \code{data} layer.
#' @param genes A character vector of gene names to use for scaling and PCA.
#'   Only genes present in \code{seurat.obj} will be used.
#' @param npcs An integer specifying the maximum number of principal components to compute. Default is 100.
#' @param truncated A logical value indicating whether to select a subset of significant PCs
#'   based on their standard deviation. If \code{TRUE}, only significant PCs are used. Default is \code{TRUE}.
#' @param ScaleData.ctrl A list of additional parameters to pass to \code{\link[Seurat]{ScaleData}}.
#'   The \code{features} and \code{verbose} arguments are set automatically.
#' @param RunPCA.ctrl A list of additional parameters to pass to \code{\link[Seurat]{RunPCA}}.
#'   The \code{features}, \code{npcs}, and \code{verbose} arguments are set automatically.
#'
#' @return A \code{Seurat} object with PC regression settings stored in its \code{misc} slot
#'   as a list named \code{NNet.setting}. The list includes:
#'   \item{pcs}{A matrix of PC embeddings for all cells. Rows correspond to cells, and columns correspond to PCs.}
#'   \item{loadings}{A matrix of PC loadings. Rows correspond to genes, and columns correspond to PCs.}
#'   \item{p}{Placeholder for the KNN affinity matrix.}
#'   \item{nn.idx}{Placeholder for the indices of nearest neighbors.}
#'   \item{nn.w}{Placeholder for the weights of nearest neighbors.}
#'   \item{cells}{Placeholder for selected cells.}
#'   \item{predictors}{Placeholder for selected predictor genes.}
#'   \item{responses}{Placeholder for selected response genes.}
#'   \item{genes}{Placeholder for the set of selected genes.}
#'   \item{lra}{Placeholder for the low-rank approximation matrix.}
#'   \item{nn.scale.gene}{Placeholder for gene-level variance scaling.}
#'   \item{nn.scale.pc}{Placeholder for PC-level variance scaling.}
#'   \item{n.eff}{Placeholder for effective neighborhood sizes.}
#'
#' @details
#' This function prepares a Seurat object for PC regression analysis. It filters the provided gene list
#' to include only those present in the Seurat object, scales the data using \code{\link[Seurat]{ScaleData}},
#' and performs PCA using \code{\link[Seurat]{RunPCA}}. If \code{truncated = TRUE}, the number of significant PCs
#' is determined using their standard deviations. Users can customize the scaling and PCA processes by
#' providing additional control parameters via \code{ScaleData.ctrl} and \code{RunPCA.ctrl}.
#'
#' The prepared Seurat object is stored with placeholders for regression settings to enable downstream
#' co-expression analysis using \code{NNet}.
#'
#' @examples
#' # Select genes for PC regression
#' genes <- select.gene(seurat.obj)$genes
#'
#' # Prepare the Seurat object with default settings
#' seurat.obj <- prepare.seurat(seurat.obj, genes, npcs = 30)
#'
#' # Customize scaling and PCA parameters
#' seurat.obj <- prepare.seurat(seurat.obj, genes, npcs = 50,
#'                              ScaleData.ctrl = list(do.center = FALSE))
#'
#' @seealso \code{\link[Seurat]{ScaleData}}, \code{\link[Seurat]{RunPCA}}
#'
#' @export
prepare.seurat <- function(seurat.obj, genes, npcs = 100, truncated = TRUE, ScaleData.ctrl = list(), RunPCA.ctrl = list()) {

  # Ensure `genes` is valid and present in `seurat.obj`
  genes <- intersect(genes, rownames(seurat.obj))
  if (length(genes) == 0) {
    stop("No valid genes found in the Seurat object. Ensure `genes` is a valid subset of rownames(seurat.obj).")
  }

  # Validate `npcs`
  if (!is.numeric(npcs) || npcs <= 0) {
    stop("`npcs` must be a positive numeric value.")
  }

  message("Running Seurat scaling and PCA...")

  # Configure and run ScaleData
  ScaleData.ctrl$features <- genes
  ScaleData.ctrl$verbose <- FALSE
  ScaleData.ctrl$do.scale <- TRUE
  args <- paste(names(ScaleData.ctrl), "=ScaleData.ctrl$",
                names(ScaleData.ctrl), sep = "", collapse = ",")
  args <- paste("Seurat::ScaleData(object = seurat.obj,", args, ")", sep = "")
  seurat.obj <- eval(parse(text = args))

  # Configure and run RunPCA
  RunPCA.ctrl$npcs <- npcs
  RunPCA.ctrl$features <- genes
  RunPCA.ctrl$verbose <- FALSE
  args <- paste(names(RunPCA.ctrl), "=RunPCA.ctrl$",
                names(RunPCA.ctrl), sep = "", collapse = ",")
  args <- paste("Seurat:: RunPCA(object = seurat.obj,", args, ")", sep = "")
  seurat.obj <- eval(parse(text = args))

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

  # Create a settings list to store in the Seurat object
  setting <- list(
    pcs = pcs,
    loadings = loadings,
    p = NULL,
    nn.idx = NULL,
    nn.w = NULL,
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

