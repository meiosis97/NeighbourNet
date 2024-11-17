# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Generate a Pruned Prior Model from Personalized PageRank Results
#'
#' This function retrieves and prunes the personalized PageRank (PPR) matrix
#' to generate a prior model of receptor regulatory potential on target genes.
#' The pruning is based on a quantile threshold applied to each receptor's regulatory
#' potential values.
#'
#' @param p A numeric value specifying the quantile threshold for pruning the
#' PPR matrix. If \code{NULL}, the default threshold stored in \code{\link{receptor.ppr}$ltf}
#' is used.
#'
#' @return A sparse matrix of the same dimensions as the PPR matrix, where
#' receptor-target regulatory potentials below the specified quantile threshold
#' are set to zero.
#'
#' @details
#' The function applies element-wise pruning to the PPR matrix, \code{\link{receptor.ppr}$ppr},
#' based on the quantile threshold provided by the user or the default learned
#' threshold \code{\link{receptor.ppr}$ltf}. For each receptor (row),
#' regulatory potentials below the specified quantile are zeroed out, retaining
#' only high-confidence interactions.
#'
#' This pruned prior model can be used in downstream analyses to study regulatory
#' interactions with higher confidence.
#'
#' @examples
#' # Generate a pruned prior model using the default threshold
#' pruned_model <- get.prior.model()
#'
#' # Generate a pruned prior model using a custom threshold (e.g., 0.75 quantile)
#' pruned_model <- get.prior.model(p = 0.75)
#'
#' # Check the dimensions of the pruned model
#' dim(pruned_model)
#'
#' @seealso \code{\link{receptor.ppr}}
#'
#' @export
get.prior.model <- function(p = NULL){

  # Load the personalized PageRank (PPR) matrix from the NeighbourNet package
  pruned.ppr <- NeighbourNet::receptor.ppr$ppr

  # Use the default learned quantile threshold if no custom threshold is provided
  if(is.null(p)) p <- NeighbourNet::receptor.ppr$ltf

  # Calculate the cutoff value for each row (receptor) in the PPR matrix
  # The cutoff is the specified quantile of each receptor's regulatory potential values
  cutoff <- apply(pruned.ppr, 1, quantile, p)

  # Prune the PPR matrix: set values below the quantile threshold to zero
  pruned.ppr <- pruned.ppr * (pruned.ppr > cutoff)

  # Convert the pruned matrix to a sparse matrix for efficient storage and computation
  Matrix::Matrix(pruned.ppr)
}
