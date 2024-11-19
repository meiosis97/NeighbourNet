# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Prepare a Seurat Object with Selected Genes and PCA.
#'
#' This function scales and performs PCA on a Seurat object using a specified set of genes
#' and number of principal components.
#'
#' @param seurat.obj A \code{Seurat} object containing \code{data} layer.
#' @param genes A character vector of gene names to use for scaling and PCA.
#'   Only genes present in \code{seurat.obj} will be used.
#' @param npcs An integer specifying the number of principal components to compute. Default is 100.
#' @param ScaleData.ctrl A list of additional parameters to pass to \code{\link[Seurat]{ScaleData}}.
#'   The \code{features} and \code{verbose} arguments will be set automatically.
#' @param RunPCA.ctrl A list of additional parameters to pass to \code{\link[Seurat]{RunPCA}}.
#'   The \code{features}, \code{npcs}, and \code{verbose} arguments will be set automatically.
#'
#' @return A \code{Seurat} object with scaled data and PCA results added.
#'
#' @details
#' The function first filters the provided gene list to include only those present
#' in the Seurat object. It then scales the data using \code{\link[Seurat]{ScaleData}} and
#' performs PCA using \code{\link[Seurat]{RunPCA}}. Users can customize the scaling and PCA
#' processes by providing additional control parameters through \code{ScaleData.ctrl}
#' and \code{RunPCA.ctrl}. See documentation of \pkg{Seurat} for details.
#'
#' @examples
#' # Select genes for PC regression
#' genes <- select.gene(seurat.obj)$genes
#'
#' # Prepare the Seurat object with default settings
#' seurat.obj <- prepare.seurat(seurat.obj, genes, npcs = 30)
#'
#' # Customize the scaling and PCA
#' seurat.obj <- prepare.seurat(seurat.obj, genes, npcs = 50,
#'                              ScaleData.ctrl = list(do.center = FALSE))
#'
#'
#' @seealso \code{\link[Seurat]{ScaleData}}, \code{\link[Seurat]{RunPCA}}
#'
#' @export
prepare.seurat <- function(seurat.obj, genes, npcs = 100, ScaleData.ctrl = list(), RunPCA.ctrl = list()) {

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

  return(seurat.obj)
}
