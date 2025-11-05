# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Build Meta-Networks (Meta-response) via Network Embeddings (PCA / Non-Negative PCA)
#'
#' \code{build.meta.network} aggregates per-cell networks into a low-dimensional set of
#' meta-networks by network embedding (nPCA on vectorised networks). Concatenating
#' network vectors yields an edge-by-cell “long-and-tall” matrix. For efficiency,
#' the edge dimension is first embedded by SVD/PCA, and then nPCA is applied to
#' the embedded matrix. It accepts either a \code{Seurat} object prepared by the
#' NeighbourNet workflow or a precomputed 3D network tensor. Similarly, using
#' \code{build.meta.response}, meta-responses representing functional gene modules can
#' be built by nPCA embedding on vectorised response profiles (cell*predictors-by-response). 
#' The detailed arguments and usage of \code{build.meta.network} is documented, and that of
#' \code{build.meta.response} is alike. 
#'
#' @param seurat.obj A \code{Seurat} object with a \code{NNet.mod} list stored in the
#'   \code{misc} slot. This list is created by  \code{\link{run.nn.reg}}. 
#'   If provided, \code{assay} and gene panels are taken from
#'   \code{NNet.mod$defaults}. Ignored when \code{network} is supplied.
#' @param network Optional 3D tensor of networks with dimension order
#'   (responses x predictors x cells).
#'   If supplied, \code{return.p.val} is forced to \code{FALSE}.
#' @param k Integer; Ranks to keep for the “long-and-tall” matrix of concatenated network vectors.
#'   Limited to \code{min(k, n.response * n.predictor, n.cell)}. Default: \code{100}.
#' @param cutoff Numeric p-value threshold (in \code{[0, 1]}) passed to
#'   \code{\link{get.network}} when pulling per-cell networks. \code{NULL} uses
#'   defaults from \code{NNet.mod$defaults$cutoff}.
#' @param big.memory Logical; if \code{TRUE}, loads the full network tensor internally via
#'   \code{get.network(..., drop = FALSE, cutoff = cutoff)} to reduce repeated reads.
#'   Default: \code{FALSE}.
#' @param scale Logical; if \code{TRUE}, scales each single-cell network by its Frobenius
#'   norm across edges before aggregation (per-cell normalization). Default: \code{TRUE}.
#' @param truncated Logical; if \code{TRUE}, retain only significant edge embeddings
#'   (spectral-gap heuristic) before running nPCA; otherwise, use exactly \code{k}
#'   embeddings. Default: \code{TRUE}. If \code{k < 100}, this is set to \code{FALSE}.
#' @param n.net Integer; number of meta-networks to output. Default: \code{20}.
#' @param non.neg Logical; if \code{TRUE}, performs nPCA (iterative deflation with
#'   non-negativity constraints) on the edge-embedded matrix, then embeds cell networks
#'   with nPCA loadings. Otherwise, embeds with standard PCA loadings. Default: \code{TRUE}.
#' @param max.iter Integer; maximum iterations for calculating each nPCA component.
#'   Default: \code{1000}.
#' @param tol Numeric; convergence tolerance for nPCA iterations. Default: \code{1e-10}.
#' @param return.p.val Logical; if \code{TRUE} and \code{assay == "effect"}, also returns
#'   per–meta-network p-values by embedding the p-value tensor (\code{NNet.mod$p.val})
#'   with the nPCA loadings. Forced to \code{FALSE} when \code{network} is supplied or
#'   when \code{non.neg = FALSE}. Default: \code{TRUE}.
#'
#' @return If Seurat object is provided, then update it with
#'   \code{NNet.mod$meta.network} populated, which is, a list containing:
#'   \itemize{
#'     \item \code{meta.network}: Array of size
#'       (responses x predictors x \code{n.net}).
#'     \item \code{p.val}: An Array (same size) if \code{return.p.val = TRUE}; else \code{NULL}.
#'     \item \code{pcs}: A matrix of PC scores learnt on per-cell networks; rows = cells, columns = components.
#'     \item \code{pca.loadings}: A matrix of PCA loadings learnt on per-cell networks.
#'     \item \code{pca.sd}: A numeric vector of singular values from network covariance.
#'     \item \code{npca.loadings}: A matrix of PC scores learnt on per-cell networks; 
#'                                 if \code{non.neg = TRUE}; else \code{NULL}.
#'     \item \code{scale}, \code{non.neg}: Echo of the input flags.
#'     \item \code{setting}: Copy of \code{NNet.mod$defaults} if \code{seurat.obj} was provided.
#'   }
#'
#' @details
#' \strong{Edge embedding}: If both \code{responses} and \code{predictors} have size > 1, a
#' cell-by-cell kernel that represents co-variation in network structure across cells is computed. 
#' Eigen-decomposition is ran to construct a edge-embedding x cell matrix.
#' If one axis is size 1, perform SVD directly on the (edges × cells) matrix.
#'
#' \strong{Rank selection}: If \code{truncated = TRUE}, a simple spectral-gap rule over the
#' tail of singular-value differences selects the working rank of the edge emebdding before nPCA.
#'
#' \strong{Meta-networks}: Project per-cell network onto the first
#' \code{n.net} directions (nPCA or PCA) to obtain meta-network slices.
#'
#' @examples
#' # From a Seurat object (built by run.nn.reg)
#' seurat.obj <- build.meta.network(seurat.obj, k = 100, n.net = 20, non.neg = TRUE)
#'
#' @seealso \code{\link{run.nn.reg}}, \code{\link{get.network}}
#'
#' @export
build.meta.network <- function(
  seurat.obj = NULL,
  network = NULL,
  k = 100,
  cutoff = NULL,
  big.memory = FALSE,
  scale = TRUE,
  truncated = TRUE,
  n.net = 20,
  non.neg = TRUE,
  max.iter = 1000,
  tol = 1e-10,
  return.p.val = TRUE
) {
  # Resolve inputs: Seurat-driven vs. direct tensor; support NNet.mod and mod
  if (!is.null(seurat.obj)) {
    # Model extraction
    mod <- Seurat::Misc(seurat.obj, "NNet.mod")
    if (is.null(mod)) {
      stop("No model found in misc. Run prepare.reg / run.nn.reg first.")
    }

    assay <- mod$defaults$assay
    assay <- match.arg(assay, choices = c("effect", "p.val"))
    cells <- names(mod$cells)
    predictors <- mod$defaults$predictors
    responses  <- mod$defaults$responses

    # Optionally materialize full tensor to avoid repeated get.network calls
    if (big.memory) {
      network <- NeighbourNet::get.network(seurat.obj, drop = FALSE, cutoff = cutoff)
    }
  } else if (!is.null(network)) {
    assay <- NULL
    cells <- 1:dim(network)[[3]]
    names(cells) <- dimnames(network)[[3]]
    predictors <- dimnames(network)[[2]]
    responses  <- dimnames(network)[[1]]
    return.p.val <- FALSE
  } else {
    stop("Neither a Seurat object nor a network tensor is provided.")
  }

  n.cell      <- length(cells)
  n.response  <- length(responses)
  n.predictor <- length(predictors)

  # Per-cell scaling (Frobenius norm) if requested
  if (scale) scales <- rep(0, n.cell)

  # Number of eigenvectors to extract
  k <- min(k, n.predictor * n.response, n.cell)
  if (k < 100) truncated <- FALSE

  # Build edge embedding: covariance (multi-edge) or SVD (single-edge axis)
  message("Now construct the covariance matrix.")

  if (min(n.predictor, n.response) != 1) {
    # Retrieved from @https://www.dummies.com/article/technology/programming-web-design/r/how-to-generate-your-own-error-messages-in-r-175112/
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                   total = n.response,
                                   complete = "=",   # Completion bar character
                                   incomplete = "-", # Incomplete bar character
                                   current = ">",    # Current bar character
                                   clear = FALSE,    # If TRUE, clears the bar when finish
                                   width = 100)      # Width of the progress bar

    # Multiple predictors and responses: accumulate K = sum_r (A_r' A_r)
    K <- matrix(0, nrow = n.cell, ncol = n.cell,
                dimnames = list(cells, cells))

    for (r in responses) {
      pb$tick()
      A <- if (is.null(network)) {
        NeighbourNet::get.network(seurat.obj, responses = r, cutoff = cutoff)
      } else {
        network[r, , ]
      }
      if (scale) scales <- scales + colSums(A^2)
      A <- Matrix::Matrix(A)
      K <- K + crossprod(A)
    }

    if (scale) {
      scales <- sqrt(scales)
      scales[scales == 0] <- 1
      K <- K / outer(scales, scales)
    }

    # Top-k eigendecomposition of K
    message("Eigen decomposition.")
    eigs <- if(k < n.cell){
      RSpectra::eigs(K, k = k)
    } else {
      as.matrix(K) %>% eigen()
    }
    v <- eigs$vectors
    d <- sqrt(eigs$values)
    vd <- v %*% diag(d)
    rownames(v)  <- rownames(vd) <- cells
    colnames(v)  <- colnames(vd) <- paste("component", 1:k, sep = "_")

  } else {
    # One axis has size 1: SVD on (edges × cells)
    A <- if (is.null(network)) {
      NeighbourNet::get.network(seurat.obj, drop = FALSE, cutoff = cutoff)
    } else {
      network
    }
    # Coerce to edges × cells
    if (n.predictor == 1) {
      A <- matrix(A, ncol = n.cell, dimnames = dimnames(A)[c(1, 3)])
    }
    if (n.response == 1) {
      A <- matrix(A, ncol = n.cell, dimnames = dimnames(A)[c(2, 3)])
    }

    if (scale) {
      scales <- sqrt(colSums(A^2))
      scales[scales == 0] <- 1
      A <- sweep(A, 2, scales, "/")
    }

    message("Single value decomposition.")
    if (min(dim(A)) <= 3) {
      svd.mod <- svd(A, nu = 0, nv = k)
    } else {
      svd.mod <- RSpectra::svds(A, k = k)
    }
    v <- svd.mod$v
    d <- svd.mod$d[1:k]
    vd <- v %*% diag(d)
    rownames(v)  <- rownames(vd) <- cells
    colnames(v)  <- colnames(vd) <- paste("component", 1:k, sep = "_")
  }

  # Rank selection via simple spectral-gap on tail differences (if truncated)
  if (truncated) {
    rank <- k
  } else {
    rank <- find.significant.pcs(d) # find.significant.pcs is defined in prepare.seurat
    vd <- vd[, 1:rank, drop = FALSE]
    v  <-  v[, 1:rank,  drop = FALSE]
  }

  # nPCA (optional) and meta-network projection
  message("Non-negative PCA.")
  n.net <- min(n.net, abs(rank - 1))

  meta.network <- array(
    dim = c(n.response, n.predictor, n.net),
    dimnames = list(responses,
                    predictors,
                    paste("component", 1:n.net, sep = "_")
                   )
  )

  # p-values can be returned only with Seurat context + effect assay + non.neg
  if (is.null(assay) || !non.neg)  return.p.val <- FALSE
  
  p.val <- if (return.p.val) meta.network else NULL

  if (non.neg) {
    # Iterative nPCA with deflation in cell space
    npca.res <- NeighbourNet::npca(vd, n.net, tol, max.iter, return.score = FALSE)
    npca.loadings <- npca.res$loadings
    npca.sd <- npca.res$sd

    # Project each response’s networks
    for (r in responses) {
      # Extract network
      A <- if (is.null(network)){
        NeighbourNet::get.network(seurat.obj, responses = r, cutoff = cutoff)
      } else {
        network[r, , ]
      }

      # Calculate p-values
      if (return.p.val) {
        P <- if (assay == "effect") {
          NeighbourNet::get.network(seurat.obj, responses = r, assay = "p.val", cutoff = cutoff)
        } else {
          A
        } 
        p.val[r, , ] <- P %*% npca.loadings
      }

      # Projection
      if (scale) A <- sweep(A, 2, scales, "/")
      meta.network[r, , ] <- A %*% npca.loadings
    }

    # Scale meta-network p-value, assume at least one significant connection
    if (return.p.val) {
      for (i in 1:n.net) p.val[, , i] <- p.val[, , i] / max(p.val[, , i])
    }

  } else {
    # Standard PCA projection
    npca.loadings <- NULL
    npca.sd <- NULL
    for (r in responses) {
      A <- if (is.null(network)) {
        NeighbourNet::get.network(seurat.obj, responses = r, cutoff = cutoff)
      } else {
        network[r, , ]
      }

      # Projection
      if (scale) A <- sweep(A, 2, scales, "/")
      meta.network[r, , ] <- A %*% v[, 1:n.net, drop = FALSE]
    }
  }

  # Package outputs and persist if a Seurat object was provided
  out <- list(
    meta.network  = meta.network,
    p.val         = p.val,
    pcs           = vd,
    pca.loadings  = v,
    pca.sd        = d,
    npca.loadings = npca.loadings,
    npca.sd       = npca.sd,
    scale         = scale,
    non.neg       = non.neg,
    setting       = if (!is.null(seurat.obj)) mod$defaults else NULL
  )

  if (!is.null(seurat.obj)) {
    mod$meta.network <- out
    suppressWarnings(Seurat::Misc(seurat.obj, "NNet.mod") <- mod)
    return(seurat.obj)
  } else {
    return(out)
  }
}

#' @export
#' @rdname build.meta.network
build.meta.response <- function(
  seurat.obj = NULL,
  network = NULL,
  k = 100,
  cutoff = NULL,
  big.memory = FALSE,
  scale = TRUE,
  truncated = TRUE,
  n.net = 20,
  non.neg = TRUE,
  max.iter = 1000,
  tol = 1e-10,
  return.p.val = TRUE
) {
  # Resolve inputs: Seurat-driven vs. direct tensor (NNet.mod convention)
  if (!is.null(seurat.obj)) {
    # Model extraction
    mod <- Seurat::Misc(seurat.obj, "NNet.mod")
    if (is.null(mod)) stop("No model found in misc. Run prepare.reg / run.nn.reg first.")

    assay <- mod$defaults$assay
    assay <- match.arg(assay, choices = c("effect", "p.val"))
    cells <- mod$cells                     # integer indices named by cell barcodes
    predictors <- mod$defaults$predictors
    responses  <- mod$defaults$responses

    # Optionally materialize full tensor to avoid repeated get.network calls
    if (big.memory) {
      network <- NeighbourNet::get.network(seurat.obj, drop = FALSE, cutoff = cutoff)
    }
  } else if (!is.null(network)) {
    assay <- NULL
    cells <- 1:dim(network)[[3]]
    names(cells) <- dimnames(network)[[3]]
    predictors <- dimnames(network)[[2]]
    responses  <- dimnames(network)[[1]]
    return.p.val <- FALSE
  } else {
    stop("Neither a Seurat object nor a network tensor is provided.")
  }

  n.cell      <- length(cells)
  n.response  <- length(responses)
  n.predictor <- length(predictors)

  # Per-response scaling factors if requested
  if (scale) scales <- rep(0, n.response)

  # Effective k along the response dimension
  k <- min(k, n.predictor * n.cell, n.response)
  if (k < 100) truncated <- FALSE

  # Build response embedding: covariance (multi-edge) or SVD (single-edge axis)
  message("Now construct the covariance matrix.")

  if (min(n.predictor, n.response) != 1) {
    # Retrieved from @https://www.dummies.com/article/technology/programming-web-design/r/how-to-generate-your-own-error-messages-in-r-175112/
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                   total = n.cell,
                                   complete = "=",   # Completion bar character
                                   incomplete = "-", # Incomplete bar character
                                   current = ">",    # Current bar character
                                   clear = FALSE,    # If TRUE, clears the bar when finish
                                   width = 100)      # Width of the progress bar

    # Accumulate K = sum_c (A_c A_c^T), A_c is (responses x predictors)
    K <- matrix(0, nrow = n.response, ncol = n.response,
                dimnames = list(responses, responses))

    for (i in 1:n.cell) {
      pb$tick()
      A <- if (is.null(network)) {
        NeighbourNet::get.network(seurat.obj, i = cells[i], cutoff = cutoff)
      } else {
        network[, , i]
      }
      if (scale) scales <- scales + rowSums(A^2)
      A <- Matrix::Matrix(A)
      K <- K + tcrossprod(A)
    }

    if (scale) {
      scales <- sqrt(scales)
      scales[scales == 0] <- 1
      K <- K / outer(scales, scales)
    }

    # Top-k eigendecomposition of K
    message("Eigen decomposition.")
    eigs <- if(k < n.response){
      RSpectra::eigs(K, k = k)
    } else {
      as.matrix(K) %>% eigen()
    }
    v <- eigs$vectors
    d <- sqrt(eigs$values)
    vd <- v %*% diag(d)
    rownames(v)  <- rownames(vd) <- responses
    colnames(v)  <- colnames(vd) <- paste("component", 1:k, sep = "_")

  } else {
    # One axis has size 1: SVD on (responses × edges) matrix
    A <- if (is.null(network)) {
      NeighbourNet::get.network(seurat.obj, drop = FALSE, cutoff = cutoff)
    } else {
      network
    }

    # Coerce to responses × edges
    if (n.predictor == 1) {
      # responses × cells (when predictors = 1)
      A <- matrix(A, ncol = n.response, dimnames = dimnames(A)[c(3, 1)])
    }
    if (n.cell == 1) {
      # responses × predictors (when single cell)
      A <- matrix(A, ncol = n.response, dimnames = dimnames(A)[c(2, 1)])
    }

    if (scale) {
      scales <- sqrt(colSums(A^2))
      scales[scales == 0] <- 1
      A <- sweep(A, 2, scales, "/")
    }

    message("Single value decomposition.")
    if (min(dim(A)) <= 3) {
      svd.mod <- svd(A, nu = 0, nv = k)
    } else {
      svd.mod <- RSpectra::svds(A, k = k)
    }
    v <- svd.mod$v
    d <- svd.mod$d[1:k]
    vd <- v %*% diag(d)
    rownames(v)  <- rownames(vd) <- responses
    colnames(v)  <- colnames(vd) <- paste("component", 1:k, sep = "_")
  }
  
  # Rank selection via simple spectral-gap on tail differences (if truncated)
  if (truncated) {
    rank <- k
  } else {
    # find.significant.pcs is defined in prepare.seurat
    rank <- find.significant.pcs(d)
    vd <- vd[, 1:rank, drop = FALSE]
    v  <-  v[, 1:rank,  drop = FALSE]
  }

  # nPCA (optional) and meta-network projection
  message("Non-negative PCA.")
  n.net <- min(n.net, abs(rank - 1))

  meta.response <- array(
    dim = c(n.net, n.predictor, n.cell),
    dimnames = list(paste("component", 1:n.net, sep = "_"),
    predictors,
    names(cells)
    )
  )

  # p-values can be returned only with Seurat context + effect assay + non.neg
  if (is.null(assay) || !non.neg)  return.p.val <- FALSE
  
  p.val <- if (return.p.val) meta.response else NULL
  
  if (non.neg) {
    # Iterative nPCA with deflation in cell space
    npca.res <- NeighbourNet::npca(vd, n.net, tol, max.iter, return.score = FALSE)
    npca.loadings <- npca.res$loadings
    npca.sd <- npca.res$sd

    # Project each cell’s networks (component x predictors x cell)
    for (i in 1:n.cell) {
      A <- if (is.null(network)) {
        NeighbourNet::get.network(seurat.obj, i = cells[i], cutoff = cutoff)
      } else {
        network[, , i]
      }

      # Optional p-values projection
      if (return.p.val) {
        P <- if (assay == "effect") {
          NeighbourNet::get.network(seurat.obj, i = cells[i], assay = "p.val", cutoff = cutoff)
        } else {
          A
        }
        p.val[, , i] <- crossprod(npca.loadings, P)
      }

      # Projection
      if (scale) A <- A / scales
      meta.response[, , i] <- crossprod(npca.loadings, A)
    }

  } else {
    # Standard PCA projection with response loadings
    npca.loadings <- NULL
    npca.sd <- NULL
    for (i in 1:n.cell) {
      A <- if (is.null(network)) {
        NeighbourNet::get.network(seurat.obj, i = cells[i], cutoff = cutoff)
      } else {
        network[, , i]
      }

      # Projection
      if (scale) A <- A / scales
      meta.response[, , i] <- crossprod(v[, 1:n.net, drop = FALSE], A)
    }
  }

  # Package outputs and persist if a Seurat object was provided
  out <- list(
    meta.response  = meta.response,
    p.val         = if (isTRUE(return.p.val)) p.val else NULL,
    pcs           = vd,
    pca.loadings  = v,
    pca.sd        = d,
    npca.loadings = npca.loadings,
    npca.sd       = npca.sd,
    scale         = scale,
    non.neg       = non.neg,
    setting       = if (!is.null(seurat.obj)) mod$defaults else NULL
  )

  if (!is.null(seurat.obj)) {
    mod$meta.response <- out
    suppressWarnings(Seurat::Misc(seurat.obj, "NNet.mod") <- mod)
    return(seurat.obj)
  } else {
    return(out)
  }
}
