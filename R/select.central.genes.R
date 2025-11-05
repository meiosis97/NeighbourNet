# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Select Central Genes from Meta-Networks via Eigenvector Centrality
#'
#' Identifies "central" response and predictor genes from NeighbourNet
#' meta-networks by decomposing each meta-network into principal components
#' (via singular value decomposition, SVD) and selecting genes with the
#' highest absolute loadings on leading singular vectors. These loadings
#' correspond to each gene’s **eigenvector centrality** within the meta-network,
#' capturing genes that are most influential in the overall co-expression
#' structure. The function optionally operates directly on a supplied
#' meta-network tensor or retrieves it from a \code{Seurat} object prepared by
#' the NeighbourNet workflow.
#'
#' @param seurat.obj A \code{Seurat} object with a \code{NNet.mod} list stored in the
#'   \code{misc} slot. This list is created by  \code{\link{run.nn.reg}}. 
#'   If provided, \code{assay} and gene panels are taken from
#'   \code{NNet.mod$defaults}. Ignored when \code{network} is supplied.
#' @param network Optional 3D tensor of networks with dimension order
#'   (responses x predictors x cells).
#' @param n.net Integer; number of leading meta-networks from which
#'  central genes are extracted
#' @param k Integer; number of leading eigenvectors (components) to extract from
#'   each meta-network when scoring central genes.
#' @param n.per.component Integer; number of top-scoring genes to select per
#'   eigenvector and per meta-network (default 4).
#' @param keep.responses Logical; controls whether responses are also selected.
#'   \itemize{
#'     \item If \code{FALSE} (default), both responses and predictors are
#'       selected based on eigenvector centrality scores.
#'     \item If \code{TRUE}, all responses are retained (no selection) and only
#'       predictors are selected as central genes.
#'   }
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{central.responses}: Character vector of selected central
#'       response genes. When \code{keep.responses = TRUE}, this may contain
#'       all responses.
#'     \item \code{central.predictors}: Character vector of selected central
#'       predictor genes.
#'     \item \code{response.module}: Factor giving the meta-network membership
#'       label (\code{"M1"}, \code{"M2"}, ...) for each unique
#'       \code{central.responses}.
#'       Set to \code{NULL} when all responses are kept without assignment.
#'     \item \code{predictor.module}: Factor giving the meta-network membership
#'       label for each unique \code{central.predictors}. Set to \code{NULL} when all predictors are
#'       kept without assignment.
#'   }
#'
#' @details
#' \strong{Eigenvector centrality interpretation}
#' \itemize{
#'   \item Each meta-network is treated as a weighted bipartite graph between
#'         response and predictor genes.
#'   \item The singular vectors from SVD represent the principal axes of
#'         variation in this bipartite graph. Genes with high absolute loadings
#'         on leading singular vectors have strong co-expression links with
#'         other highly connected genes. This is analogous to high
#'         \emph{eigenvector centrality} in network theory.
#'   \item Selecting top-loading genes across singular vectors therefore yields
#'         a set of “hub-like” genes that summarise the
#'         co-expression pattern of each meta-network.
#' }
#'
#' \strong{Selection strategy}
#' \itemize{
#'   \item For each meta-network, a truncated SVD is applied to the
#'         weighted adjancency matrix. 
#'   \item Absolute singular vector loadings are used as eigenvector centrality
#'         scores.
#'   \item Genes with the highest loadings per component (up to
#'         \code{n.per.component}) are selected.
#'   \item Previously selected genes are zeroed out before selecting additional
#'         components, ensuring diverse representation across meta-networks.
#' }
#'
#' @examples
#' # From a Seurat object with NNet meta-networks
#' sel <- select.central.genes(seurat.obj,
#'                             n.net = 5,
#'                             k = 2,
#'                             n.per.component = 3)
#'
#' # View top eigen-central genes
#' head(sel$central.predictors)
#'
#' @seealso \code{\link{build.meta.network}}
#'
#' @export

select.central.genes <- function(seurat.obj = NULL,
                                 network = NULL,
                                 n.net = NULL,
                                 k = 1,
                                 n.per.component = 4,
                                 keep.responses = FALSE) {

  # 1. Retrieve meta-network tensor
  if (is.null(network)) {
    network <- NeighbourNet::get.network(seurat.obj, assay = "meta.network", drop = FALSE)
  }

  # Number of meta-networks to use
  total.net <- dim(network)[3]
  if (is.null(n.net)) {
    n.net <- total.net
  } else {
    n.net <- min(total.net, n.net)
  }
  network <- network[, , 1:n.net, drop = FALSE]

  # Extract dimnames / sizes
  responses  <- dimnames(network)[[1]]
  predictors <- dimnames(network)[[2]]
  n.response  <- length(responses)
  n.predictor <- length(predictors)

  # 2. Helper: indices of top n loadings (ties broken at random)
  which.max.n <- function(x, n = n.per.component){
    which(rank(-x, ties.method = "random" ) <= n)
  }

  # Storage for selected genes
  central.predictors <- character()
  central.responses  <- character()

  # Pre-allocate module labels (one label per candidate position)
  predictor.module <- rep(paste0("M", 1:n.net),
                          each = k * n.per.component)
  response.module  <- predictor.module

  # Effective k: cannot exceed either gene dimension
  k.use <- min(k, n.predictor, n.response)

  # 3. Main selection

  # Case A: many responses, keep all responses and select predictors only
  if (n.response > n.predictor && keep.responses) {
    central.predictors <- predictors
    central.responses  <- responses
    predictor.module <- NULL
    response.module  <- NULL

  # Case B: select both responses and predictors
  } else if (!keep.responses) {

    # General SVD-based case: at least 2 responses and 2 predictors
    if (min(n.predictor, n.response) != 1) {
      for (i in 1:n.net) {
        mat <- network[, , i]

        # Truncated SVD: base svd for small k, RSpectra for larger k
        if (k.use <= 3) {
          svd.mod <- suppressWarnings(svd(mat, nu = k.use, nv = k.use))
        } else {
          svd.mod <- suppressWarnings(RSpectra::svds(mat, k = k.use))
        }

        # Absolute loadings as importance scores
        v <- abs(svd.mod$v)  # predictors
        u <- abs(svd.mod$u)  # responses
        rownames(v) <- predictors
        rownames(u) <- responses

        # Zero out already-selected genes to encourage diversity
        if (length(central.predictors)) v[central.predictors, ] <- 0
        if (length(central.responses))  u[central.responses, ]  <- 0

        # Extract top genes per component
        new.preds <- apply(v, 2, function(col) predictors[which.max.n(col)])
        new.resps <- apply(u, 2, function(col) responses[which.max.n(col)])

        central.predictors <- c(central.predictors, as.vector(new.preds))
        central.responses  <- c(central.responses,  as.vector(new.resps))
      }

    # Special case: only one predictor, select responses per component
    } else if (n.predictor == 1) {
      mat <- matrix(abs(network),
                    ncol = n.net,
                    dimnames = dimnames(network)[c(1, 3)])  # responses × components
      if (length(central.responses)) mat[central.responses, ] <- 0

      new.resps <- apply(mat, 2, function(col) responses[which.max.n(col)])
      central.responses  <- as.vector(new.resps)
      central.predictors <- predictors
      predictor.module   <- NULL

    # Special case: only one response, select predictors per component
    } else if (n.response == 1) {
      mat <- matrix(abs(network),
                    ncol = n.net,
                    dimnames = dimnames(network)[c(2, 3)])  # predictors × components
      if (length(central.predictors)) mat[central.predictors, ] <- 0

      new.preds <- apply(mat, 2, function(col) predictors[which.max.n(col)])
      central.predictors <- as.vector(new.preds)
      central.responses  <- responses
      response.module    <- NULL
    }

  # Case C: keep.responses = TRUE, but not covered by the first branch
  } else {
    # Keep all responses, select predictors only
    central.responses <- responses
    response.module   <- NULL

    if (min(n.predictor, n.response) != 1) {
      for (i in 1:n.net) {
        mat <- network[, , i]

        if (k.use <= 3) {
          svd.mod <- suppressWarnings(svd(mat, nu = k.use, nv = k.use))
        } else {
          svd.mod <- suppressWarnings(RSpectra::svds(mat, k = k.use))
        }

        v <- abs(svd.mod$v)
        rownames(v) <- predictors
        if (length(central.predictors)) v[central.predictors, ] <- 0

        new.preds <- apply(v, 2, function(col) predictors[which.max.n(col)])
        central.predictors <- c(central.predictors, as.vector(new.preds))
      }

    } else if (n.predictor == 1) {
      # All predictors kept; no module assignment needed
      central.predictors <- predictors
      predictor.module   <- NULL

    } else if (n.response == 1) {
      mat <- matrix(abs(network),
                    ncol = n.net,
                    dimnames = dimnames(network)[c(2, 3)])  # predictors × components
      if (length(central.predictors)) mat[central.predictors, ] <- 0

      new.preds <- apply(mat, 2, function(col) predictors[which.max.n(col)])
      central.predictors <- as.vector(new.preds)
    }
  }

  # 4. Deduplicate and align module labels
  central.predictors <- unique(central.predictors)
  central.responses  <- unique(central.responses)

  if (!is.null(predictor.module)) {
    predictor.module <- predictor.module[seq_along(central.predictors)]
    predictor.module <- factor(predictor.module,
                               levels = paste0("M", 1:n.net))
  }

  if (!is.null(response.module)) {
    response.module <- response.module[seq_along(central.responses)]
    response.module <- factor(response.module,
                              levels = paste0("M", 1:n.net))
  }

  # 5. Output
  list(
    central.responses  = central.responses,
    central.predictors = central.predictors,
    response.module    = response.module,
    predictor.module   = predictor.module
  )
}
