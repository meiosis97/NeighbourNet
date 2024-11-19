# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Build Graph Based on Principal Components and K-Nearest Neighbors
#'
#' Constructs an affinity graph from principal components using adaptive kernel scaling approach
#' adapted from Uniform Manifold Approximation and Projection (UMAP) (Becht et al. 2019).
#'
#' @param pcs A numeric matrix of principal components (cells x PCs).
#' @param knn Number of nearest neighbors to consider for each cell. Default is 30.
#'
#' @return A list containing:
#' \item{p}{A sparse symmetric affinity matrix representing transition probability.}
#' \item{nn.idx}{Matrix of nearest neighbor indices for each cell.}
#' \item{nn.w}{Matrix of normalized affinity weights for the nearest neighbors.}
#'
#' @references
#' Becht, E., McInnes, L., Healy, J., Dutertre, C.-A., Kwok, I. W. H., Ng, L. G., Ginhoux, F., & Newell, E. W. (2019).
#' Dimensionality reduction for visualizing single-cell data using UMAP.
#' \emph{Nature Biotechnology}, 37(1), 38--44. \doi{10.1038/nbt.4314}
#'
#' @export
build.graph <- function(pcs, knn = 30){

  n <- nrow(pcs)  # Number of cells
  nn2.result <- RANN::nn2(pcs, k = knn)  # Find k-nearest neighbors
  a <- 2 * log2(knn)  # Scaling parameter for sigma
  sigma <- numeric(n)  # Initialize vector to store sigma values

  # Scale distances for better numerical stability
  scaled.dists <- nn2.result$nn.dists
  scaled.dists <- scaled.dists - scaled.dists[,2]  # Normalize distances
  scaled.dists[,1] <- 0  # Zero out self-distances

  # Compute sigma for each cell
  for (i in 1:n) {
    dk <- scaled.dists[i, 1:knn]  # Distances to k-nearest neighbors
    sigma[i] <- FindSigma(dk, a, knn)  # Optimal sigma for cell i
  }

  # Calculate affinity weights
  nn.aff <- exp(-scaled.dists / sigma)
  nn.idx <- nn2.result$nn.idx  # Indices of k-nearest neighbors
  nn.w <- nn.aff / rowSums(nn.aff)  # Normalize affinity weights

  # Create a sparse affinity matrix
  i <- rep(1:n, times = knn)  # Row indices
  j <- as.numeric(nn.idx)  # Column indices
  x <- as.numeric(nn.aff)  # Affinity values
  aff <- Matrix::sparseMatrix(i = i, j = j, x = x)  # Sparse matrix

  # Symmetrise the affinity matrix using the probabilistic t-norm method
  aff <- aff + t(aff) - aff * t(aff)

  # Normalize the affinity matrix (I - Laplacian)
  d <- rowSums(aff)  # Degree vector
  p <- aff
  p <- 0.5 * (sweep(p, 2, sqrt(d), "/") / sqrt(d))  # Symmetric normalization
  diag(p) <- diag(p) + 0.5  # Adjust diagonal entries

  # Return the results
  return(list(p = p, nn.idx = nn.idx, nn.w = nn.w))
}

# Helper function to find the optimal sigma for adaptive kernel scaling
# dk: distance vector to nearest neighbors
# a: scaling factor
# knn: number of nearest neighbors
FindSigma <- function(dk, a, knn) {
  lower <- 0
  upper <- Inf
  cur <- dk[knn]  # Initial guess for sigma based on the kth nearest neighbor distance
  while (TRUE) {
    psum <- sum(exp(-dk / cur))  # Compute sum of exponentiated distances
    if (psum > a) {
      upper <- cur
      cur <- (lower + cur) / 2  # Narrow the range
    } else if (psum < a) {
      if (is.infinite(upper)) {
        lower <- cur
        cur <- 2 * cur  # Expand the range if upper bound is infinite
      } else {
        lower <- cur
        cur <- (upper + cur) / 2  # Narrow the range
      }
    }
    if (abs(psum - a) < 1e-5) break  # Stop when the sum is close enough to target
  }
  cur  # Return optimal sigma
}

