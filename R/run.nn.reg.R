#' Run Nearest-Neighbor PC Regression for Gene Co-Expression Analysis
#'
#' This function performs nearest-neighbor regression on a Seurat object to measure
#' gene co-expression. It uses principal components (PCs) to smooth response
#' gene expression and calculates permutation feature importance as a co-expression measure.
#' The results are stored as a list named \code{NNet.mod} in the \code{misc}
#' slot of the Seurat object for downstream analysis.
#'
#' @param seurat.obj A \code{Seurat} object on which \code{\link{prepare.graph}} and \code{\link{prepare.reg}} have been run.
#' @param responses A character vector of response gene names to return for the network ensemble. If \code{NULL}, uses all response genes from \code{NNet.setting}.
#' @param Y A custom response matrix (cells x responses). If provided, overrides \code{responses}.
#' @param predictors A character vector of predictor gene names to return for the network ensemble. If \code{NULL}, uses all predictor genes from \code{NNet.setting}.
#' @param t A numeric value indicating the power to scale the singular values of the Laplacian operator. Default is 3.
#' @param k An integer specifying the number of singular vectors to use in the Laplacian operator. Default is the number of nearest neighbors.
#' @param remove.self.loops A logical indicating whether to remove self-loops from networks for downstream analysis. Default is \code{TRUE}.
#' @param f A function to map the effect estimation to importance scores for downstream analysis. Default is \code{function(x) 2*x^2} that yields permutation feature importamce.
#' @param assay A character string indicating the co-expression measure to use in downstream analysis. Options are \code{"effect"} or \code{"p.val"}. Default is \code{"effect"}.
#' @param prune A logical indicating whether to prune networks based on p-values for downstream analysis. Default is \code{TRUE}.
#' @param cutoff A numeric value specifying the p-value threshold for pruning for downstream analysis.. Default is 0.5.
#' @param return.p.val A logical indicating whether to return p-values in addition to the effect tensor. Set to \code{FALSE} to save memory. Default is \code{FALSE}.
#' @param return.smooth A logical indicating whether to return the smoothed effect tensor.
#' If set to \code{TRUE}, networks will not be further smoothed in downstream analysis. Default is \code{TRUE}.
#' @param return.prune A logical indicating whether to return the pruned effect tensor.
#'  If set to \code{TRUE}, networks will not be further pruned in downstream analysis. Default is \code{FALSE}.
#'
#' @return A \code{Seurat} object with a list named \code{NNet.mod} stored in its \code{misc} slot.
#' \code{NNet.mod} is a list containing:
#' \item{effect}{A tensor of effect sizes (responses x predictors x cells).}
#' \item{p.val}{A tensor of p-values (responses x predictors x cells) if \code{return.p.val = TRUE}.}
#' \item{meta.network}{Reserved for future use.}
#' \item{mus}{A named vector of mean log-transformed noise distributions for each response.}
#' \item{sigmas}{A named vector of standard deviations for log-transformed noise distributions for each response.}
#' \item{subsampled}{A logical indicating whether cell subsampling was performed.}
#' \item{smoothed}{A logical indicating whether smoothed effects are returned.}
#' \item{pruned}{A logical indicating whether pruned effects are returned.}
#' \item{gene.sets}{A list of selected predictors, responses, and all genes used in the analysis. Includes TFs and targets subsets.}
#' \item{cells}{A vector of selected cell indices used in the analysis.}
#' \item{defaults}{Default settings used for network extraction, such as cutoff and pruning function.}
#' \item{custom.y}{A logical indicating whether a custom response matrix (\code{Y}) was provided.}
#' \item{w}{A list containing the Laplacian operator components (\code{u} and \code{vd}).}
#'
#' @details
#' This function performs nearest-neighbor PC regression to assess co-expression relationships between response
#' and predictor genes within a Seurat object. The co-expression measure is based on
#' permutation feature importance estimated by local gene variances,
#' enabling scalable analysis of transcriptional regulation. Parameters \code{remove.self.loops},
#' \code{f}, \code{assay}, \code{prune}, and \code{cutoff} define default settings for how networks will be
#' post-processed before downstream analysis, such as meta-network learning. However, these settings do not
#' immediately affect the regression analysis performed by this function. These defaults can be updated later
#' using \code{\link{set.defaults}}.
#'
#' @examples
#' # Assuming `seurat.obj` has been pre-processed with `prepare.reg`:
#' responses <- Seurat::Misc(seurat.obj, "NNet.setting")$responses %>%
#'              head(n =5)
#'
#' # Run PC regression
#' seurat.obj <- run.nn.reg(
#'   seurat.obj,
#'   responses = responses,
#' )
#'
#' # Check regression outcome settings
#' str(Seurat::Misc(seurat.obj, "NNet.mod"))
#'
#' @seealso \code{\link{prepare.graph}}, \code{\link{prepare.reg}}, \code{\link{set.defaults}}
#'
#' @export
run.nn.reg <- function(seurat.obj, responses = NULL, Y = NULL,
                       predictors = NULL, t = 3, k = NULL,
                       remove.self.loops = T, f = function(x) 2*x^2, assay = c("effect", "p.val"),
                       prune = TRUE, cutoff = 0.5,
                       return.p.val = FALSE, return.smooth = TRUE, return.prune = FALSE) {
  # Retrieve the stored settings from the Seurat object
  setting <- Seurat::Misc(seurat.obj, "NNet.setting")

  # Ensure that prepare.seurat has been run
  if (is.null(setting)) stop("Run prepare.seurat first, then prepare.graph and prepare.reg.")

  # Ensure that prepare.graph has been run
  if (is.null(setting$p)) stop("Run prepare.graph first and then prepare.reg.")

  # Ensure that prepare.reg has been run
  if (is.null(setting$nn.scale.gene)) stop("Run prepare.reg first.")

  # Match the assay argument to valid options
  assay <- match.arg(assay)

  # Extract PC embedding (X matrix)
  X <- setting$pcs

  # Determine if subsampling of cells has been applied
  subsampled <- !is.null(setting$cells)
  if (!subsampled) {
    cells <- 1:nrow(X)  # Use all cells if no subsampling
    names(cells) <- rownames(X)
  } else {
    cells <- setting$cells  # Use preselected cells
    names(cells) <- rownames(X)[cells]
  }

  # Extract response matrix (Y) or use a custom one
  custom.y <- !is.null(Y)
  if (!custom.y) {
    # Select response genes with non-zero local variances
    responses <- intersect(setting$responses, responses)
    responses <- responses[rowSums(setting$nn.scale.gene[responses, names(cells), drop = F]) > 0]
    Y <- setting$lra[, responses, drop = F]  # Low-rank approximated response matrix
  } else {
    Y <- as.matrix(Y)
    responses <- colnames(Y)
    # Assign default names if the responses are unnamed
    if (is.null(responses)) responses <- colnames(Y) <- paste("Y", 1:ncol(Y), sep = "")
    message("Custom response is provided, will not prune or return p-value.")
    prune <- FALSE
    return.p.val <- FALSE
  }

  # Select predictors or use all available predictors
  if (is.null(predictors)) {
    predictors <- setting$predictors
  } else {
    predictors <- intersect(predictors, setting$genes)
  }

  # Determine the gene set to be used in the analysis
  if (!custom.y) {
    genes <- unique(c(responses, predictors))
  } else {
    genes <- predictors
  }

  # Normalize PC loadings
  loadings <- setting$loadings[genes, , drop = F]
  loading.scale <- sqrt(rowSums(loadings^2))
  loadings <- loadings / loading.scale

  # Initialize variables for dimensions
  n.cell <- length(cells)
  n.gene <- length(genes)
  n.response <- length(responses)
  n.predictor <- length(predictors)
  n.pc <- ncol(X)
  if (is.null(k)) k <- ncol(setting$nn.idx)

  # Adjust logical flags for smooth and pruned effects
  if (!return.smooth) return.prune <- FALSE
  if (return.prune) prune <- FALSE

  # Display messages about the analysis settings
  message(ifelse(return.smooth, "Return smoothed effect, can only generate networks for sampled cells.", "Return raw effect."))
  message(ifelse(return.prune, "Return pruned effect.", "Return unpruned effect."))
  message(ifelse(return.p.val, "Return p-value.", "Will not return p-value."))
  message(ifelse(assay == "effect", "Downstream analysis will be performed on the effect tensor.", "Downstream analysis will be performed on the p-val tensor."))
  message(ifelse(prune, "Downstream analysis will perform network pruning.", "Downstream analysis will not perform network pruning."))

  # Retrieve transcription factors (TFs) and targets for reporting
  gene.list <- NeighbourNet::gene.list
  tfs.in.responses <- responses[responses %in% gene.list$tfs]
  tfs.in.predictors <- predictors[predictors %in% gene.list$tfs]
  targets.in.responses <- responses[responses %in% gene.list$targets]
  targets.in.predictors <- predictors[predictors %in% gene.list$targets]

  # Initialize a matrix for local variances of response genes
  nn.scale.y <- matrix(0, nrow = n.response, ncol = n.cell, dimnames = dimnames(Y[cells, ]) %>% rev)

  # Calculate local variances for responses
  for (i in 1:n.cell) {
    j <- cells[i]
    idx <- setting$nn.idx[j, ]
    w <- setting$nn.w[j, ]

    # Compute residuals for response genes
    w.mean <- as.numeric(w %*% Y[idx, ])
    res <- t(Y[idx, ]) - w.mean
    nn.scale.y[, i] <- as.numeric(res^2 %*% w) * setting$n.eff[i] / (setting$n.eff[i] - 1)
  }
  nn.scale.y <- sqrt(nn.scale.y)

  # Build the Laplacian operator for smoothing
  message("Build the Laplacian operator.")
  svds.p <- RSpectra::svds(A = setting$p, k = k)  # Singular value decomposition
  u <- svds.p$u
  vd <- sweep(svds.p$v, 2, svds.p$d^t, "*")  # Scale singular vectors
  rownames(u) <- rownames(vd) <- rownames(X)
  d <- rowSums(tcrossprod(u, vd[cells, ]))  # Normalize scaling factor
  u <- u / d

  # Initialize tensors for storing effects and p-values
  effect.tensor <- array(dim = c(n.response, n.gene, n.cell), dimnames = list(responses, genes, names(cells)))
  p.val.tensor <- if (return.p.val) effect.tensor else NULL

  # Initialize containers for noise distribution statistics
  mus <- c()
  sigmas <- c()

  # Perform regression for each response gene
  message("Now regress.")

  # Retrieved from @https://www.dummies.com/article/technology/programming-web-design/r/how-to-generate-your-own-error-messages-in-r-175112/
  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                   total = n.response,
                                   complete = "=",   # Completion bar character
                                   incomplete = "-", # Incomplete bar character
                                   current = ">",    # Current bar character
                                   clear = FALSE,    # If TRUE, clears the bar when finish
                                   width = 100)      # Width of the progress bar

  for (i in 1:n.response) {
    pb$tick()

    # Initialize a matrix to store regression coefficients
    b <- matrix(0, nrow = n.cell, ncol = n.pc)
    for (j in 1:n.cell) {
      idx <- setting$nn.idx[cells[j], ]
      w <- setting$nn.w[cells[j], ]

      # Perform local scaling of data
      w.mean <- as.numeric(w %*% Y[idx, i])
      y <- scale(Y[idx, i], center = w.mean, scale = nn.scale.y[i, j])
      w.mean <- as.numeric(w %*% X[idx, ])
      x <- scale(X[idx, ], center = w.mean, scale = setting$nn.scale.pc[, j])

      # Perform local regression
      if (!is.na(attr(y, "scaled:scale"))) {
        lambda <- 5
        qr.mod <- qr(rbind(x * sqrt(w), diag(sqrt(lambda), ncol(x))))
        v <- qr.qty(qr.mod, c(y * sqrt(w), rep(0, ncol(x))))
        b[j, ] <- backsolve(qr.R(qr.mod), v) / attr(x, "scaled:scale") * attr(y, "scaled:scale")
      } else {
        b[j, ] <- rep(0, ncol(x))
      }
    }

    # Compute noise distribution for pruning
    rand.loadings <- replicate(100, rnorm(n.pc)) %>% t
    rand.loadings <- rand.loadings/sqrt(rowSums(rand.loadings^2))
    noise <- tcrossprod(rand.loadings, b) %>% as.matrix
    noise <- tcrossprod(noise %*% vd[cells, ], u[cells, ])
    noise <- log(abs(noise))
    mus[i] <- mean(noise)
    sigmas[i] <- sd(noise)

    # Transform regression coefficients to effects
    b <- tcrossprod(loadings, b) %>% as.matrix

    # Calculate effect
    effect <- (b * setting$nn.scale.gene[genes,names(cells),drop=F]) %>%
      as.matrix

    if (return.smooth | return.p.val) effect.hat <- tcrossprod(effect %*% vd[cells, ], u[cells, ])

    # Compute p-values if required
    if (return.p.val | return.prune) {
      p.val <- pnorm(log(abs(effect.hat)), mus[i], sigmas[i])
      if(return.p.val) p.val.tensor[i, , ] <- p.val
    }

    # Store smoothed or raw effects in the tensor
    if (return.smooth) {
      if (return.prune) effect.hat[p.val < cutoff] <- 0
      effect.tensor[i, , ] <- effect.hat
    } else {
      effect.tensor[i, , ] <- effect
    }

  }
  names(mus) <- names(sigmas) <- responses

  mod <- list(
    effect = effect.tensor, p.val = p.val.tensor, meta.network = NULL,
    mus = mus, sigmas = sigmas, subsampled = subsampled,
    smoothed = return.smooth, pruned = return.prune,
    gene.sets = list(
      predictors = list(genes = predictors, tfs = tfs.in.predictors, targets = targets.in.predictors),
      responses = list(genes = responses, tfs = tfs.in.responses, targets = targets.in.responses),
      genes = genes
    ),
    cells = cells,
    defaults = list(f = f, remove.self.loops = remove.self.loops, assay = assay, predictors = predictors, responses = responses, cutoff = cutoff, prune = prune),
    custom.y = custom.y, w = list(u = u, vd = vd)
  )

  class(mod) <- "NNet.mod"

  # Store results in the Seurat object
  suppressWarnings(
    Seurat::Misc(seurat.obj, "NNet.mod") <- mod
  )

  return(seurat.obj)  # Return the updated Seurat object
}
