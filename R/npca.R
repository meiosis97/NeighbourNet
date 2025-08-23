# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Non-Negative PCA with Deflation (nPCA)
#'
#' Computes a non-negative PCA (nPCA) basis via iterative regression and
#' deflation on an input matrix \code{x}. 
#' 
#' @param x A numeric matrix of size (p x n). The row dimension will be embedded
#' @param k Integer; number of nPCA components to extract.
#' @param max.iter Integer; maximum iterations for calculating each nPCA component.
#'   Default: \code{1000}.
#' @param tol Numeric; convergence tolerance for nPCA iterations. Default: \code{1e-10}.
#' @param return.score Logical; If \code{TRUE} return the nPCA projection of \code{x}.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{scores}: Matrix (n x k) of non-negative nPCA scores
#'           (columns are components).
#'     \item \code{loadings}: Matrix (p x k) of non-negative nPCA loadings
#'     \item \code{sd}: Numeric vector (\code{k}) of nPCA “standard deviations”.
#'   }
#'
#' @details
#' This function implements non-negative PCA using restricted single-value decomposition, 
#' which is solved using the iterative regression approach introduced by Sigg and Buhmann (2008) and deflation
#' approach introduced by Mackey (2008).
#' 
#' @references 
#' Sigg, C. D., & Buhmann, J. M. (2008, July). Expectation-maximization for sparse and non-negative PCA.
#' \emph{In Proceedings of the 25th international conference on Machine learning} (pp. 960-967).
#' 
#' Mackey, L. (2008). Deflation methods for sparse PCA. Advances in neural information processing systems, 21.
#'
#' @examples
#' x <- matrix(rnorm(1000), nrow - 100)
#' res <- nn.pca(x, k = 10, tol = 1e-10, n.iter = 1000)
#' str(res$loadings); str(res$sd)
#'
#' @export
npca <- function(x, k, tol = 1e-10, max.iter = 1000, return.score = TRUE) {
  # Input checks and setup
  if (!is.matrix(x)) x <- as.matrix(x)
  p <- nrow(x)
  n <- ncol(x)
  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("k must be a positive integer.")
  }

  # Metric for normalization (identity deflation metric, as in your code)
  B <- diag(p)

  # Outputs
  loadings <- matrix(0, nrow = p, ncol = k,
                     dimnames = list(rownames(x), paste0("component_", 1:k)))
  sd  <- numeric(k)

  # Track total “variance” as in your previous implementation
  tot.var <- sum(x^2) / n

  # Regression and deflation
  for (i in 1:k) {
    # Non-negative init on unit sphere
    w <- runif(p)
    w <- w / sqrt(sum(w^2))
    for (j in 1:max.iter) {
      w.old <- w
      u <- crossprod(x, w)
      w <- x %*% u
      w <- w / sum(w)
      w[w < 0] <- 0
      w <- w / sqrt(sum(w^2))
      if (sum((w - w.old)^2) < tol) break
    }

    # Normalize w in metric B, then deflate
    w <- w / as.numeric(sqrt(crossprod(w, B) %*% w))
    q <- B %*% w
    x <- x - q %*% (t(q) %*% x)
    B <- B - B %*% q %*% t(q)

    # Store loadings and component “sd” (matching your original formula)
    loadings[, i] <- w / sqrt(sum(w^2))
    var.diff <- tot.var - sum(x^2) / n
    sd[i] <- sqrt(max(var.diff, 0))
    tot.var <- tot.var - var.diff
  }
  
  
  scores <- if(return.score) {
    tcrossprod(x, loadings) 
  } else {
    NULL
  }

  list(scores = scores, loadings = loadings, sd = sd)
}