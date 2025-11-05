# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Gene List from Integrated Prior Knowledge Network
#'
#' A list of gene sets derived from an integrated prior knowledge network used
#' in the NeighbourNet (NNet) framework. This list includes transcription factors (TFs),
#' target genes, receptors, and ligands.
#'
#' @format A list with four entries:
#' \describe{
#'   \item{\code{tfs}}{A character vector of transcription factor (TF) gene names. These are regulators that control gene expression.}
#'   \item{\code{targets}}{A character vector of target gene names. These genes are regulated by TFs.}
#'   \item{\code{receptors}}{A character vector of receptor gene names. These genes encode proteins that receive extracellular signals.}
#'   \item{\code{ligands}}{A character vector of ligand gene names. These genes encode signaling molecules that bind to receptors.}
#' }
#'
#' @details
#' This list is crucial for gene regulatory network (GRN) inference and upstream
#' signaling pathway (USP) analysis in NNet.
#'
#' @source The integrated prior knowledge network used in NNet analysis.
#'
"gene.list"
