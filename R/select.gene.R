# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Select Genes for Analysis
#'
#' This function filters and selects transcription factors (TFs), target genes,
#' and background genes from a Seurat object based on their expression in a
#' minimum number of cells.
#'
#' @param seurat.obj A \code{Seurat} object containing gene expression data.
#' @param tfs A vector of transcription factor (TF) gene names. If \code{NULL},
#' all TFs in the \code{\link{gene.list}$tfs} will be considered. Default is \code{NULL}.
#' @param targets A vector of target gene names. If \code{NULL}, all target
#' genes in \code{\link{gene.list}$targets} will be considered. Default is \code{NULL}.
#' @param bgs A vector of background gene names. If \code{NULL}, no additional
#' background genes are considered. Default is \code{NULL}.
#' @param min.cells An integer specifying the minimum number of cells in which
#' a gene must be expressed to be included. Default is 20.
#'
#' @return A list containing:
#' \item{tfs}{Filtered TFs that meet the minimum cell expression criteria.}
#' \item{targets}{Filtered target genes that meet the minimum cell expression criteria.}
#' \item{bgs}{Filtered background genes, excluding any genes in \code{tfs} or \code{targets}.}
#' \item{genes}{A unique set of all selected genes.}
#'
#' @details
#' The function first filters genes in the \code{seurat.obj} to retain only those
#' expressed in more than \code{min.cells} cells. It then intersects the filtered
#' gene set with provided or default lists of TFs, targets, and background genes.
#'
#' @examples
#' selected.genes <- select.gene(seurat.obj)
#'
#' @seealso \code{\link{gr.graph}}
#'
#' @export
select.gene <- function(seurat.obj,
                        tfs = NULL,
                        targets = NULL,
                        bgs = NULL,
                        min.cells = 20){
  # Load the default gene list from NeighbourNet
  gene.list <- NeighbourNet::gene.list

  # Select TFs: default to all TFs in gene.list if tfs is NULL
  tfs <- if (is.null(tfs)) unique(gene.list$tfs) else intersect(gene.list$tfs, tfs)

  # Select targets: default to all targets in gene.list if targets is NULL
  targets <- if (is.null(targets)) gene.list$targets else intersect(gene.list$targets, targets)

  # Filter genes based on expression in at least min.cells
  genes <- rownames(seurat.obj)
  n.cell.expressed <- rowSums(SeuratObject::LayerData(obj, "counts")[genes, ] > 0)
  genes <- genes[n.cell.expressed > min.cells]

  # Filter TFs, targets, and background genes from the expressed genes
  tfs <- intersect(genes, tfs)
  targets <- intersect(genes, targets)

  # Handle background genes if provided
  if (!is.null(bgs)) {
    bgs <- intersect(genes, bgs)
    bgs <- bgs[!bgs %in% c(tfs, targets)]  # Exclude TFs and targets from background genes
  } else {
    bgs <- character(0)
  }

  # Combine all selected genes
  genes <- unique(c(tfs, targets, bgs))

  # Return the filtered gene lists
  list(tfs = tfs,
       targets = targets,
       bgs = bgs,
       genes = genes
      )
}
