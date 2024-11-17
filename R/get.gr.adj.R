# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Create an Adjacency Matrix of Gene Regulation with Propagation
#'
#' This function generates an adjacency matrix from a gene regulatory network
#' (\code{\link{gr.graph}}) with propagated regulatory effects. The propagation
#' iteratively multiplies the adjacency matrix to capture higher-order regulatory
#' interactions.
#'
#' @param t An integer specifying the number of propagation steps. Default is 2.
#' If \code{t = 1}, the function returns the direct adjacency matrix without propagation.
#'
#' @return A numeric adjacency matrix representing the propagated gene regulatory
#' interactions. Each entry in the matrix reflects the weighted regulatory influence
#' of a transcription factor (TF) on a target gene, adjusted by consensus stimulation
#' and inhibition.
#'
#' @details
#' The function starts by modifying the weights of \code{\link{gr.graph}} edges to reflect
#' the regulatory effects. The edge weights are adjusted by multiplying the original
#' weights with the difference between the indicator attributes \code{consensus_stimulation} and
#' \code{consensus_inhibition}.
#'
#' The resulting weighted adjacency matrix is then propagated by matrix
#' multiplication for \code{t} steps to capture indirect regulatory interactions.
#' Positive weights suggest stimulation effect, whereas negative weights indicate
#' inhibition effect
#'
#' @examples
#' # Generate an adjacency matrix with the default 2-step propagation
#' adj_matrix <- get.gr.adj()
#'
#' # Generate an adjacency matrix with 3-step propagation
#' adj_matrix <- get.gr.adj(t = 3)
#'
#' @seealso \code{\link{gr.graph}}
#'
#' @export
get.gr.adj <- function(t = 2) {

  # Load the gene regulatory graph
  gr.graph <- NeighbourNet::gr.graph

  # Adjust edge weights using indicators of consensus stimulation and inhibition
  # Positive weight: stimulation
  # Negative weight: inhibition
  weight <- igraph::edge.attributes(gr.graph)$weight
  stimulation <- igraph::edge.attributes(gr.graph)$consensus_stimulation
  inhibition <- igraph::edge.attributes(gr.graph)$consensus_inhibition
  igraph::edge.attributes(gr.graph)$weight <- weight * (stimulation - inhibition)

  # Create the adjacency matrix
  adj <- igraph::as_adjacency_matrix(gr.graph, attr = "weight")

  # Use matrix exponentiation for propagation
  if (t > 1) {
    adj <- adj %^% t
  }

  adj
}

# Helper function for matrix exponentiation
`%^%` <- function(mat, power) {
  if (power == 1) return(mat) else power <- power - 1
  for(i in 1:power) mat <- mat %*% mat
  mat
}
