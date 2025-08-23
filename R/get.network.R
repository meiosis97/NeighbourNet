# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Extract Per-Cell Gene Co-Expression Networks
#'
#' Retrieves effect or p-value networks (or meta-networks) from a \code{Seurat}
#' object prepared by the NeighbourNet workflow. Supports optional pruning,
#' transformation, self-loop removal, and flexible cell indexing. By default,
#' argument values fall back to \code{NNet.mod$defaults} (\code{defaults}).
#'
#' @param seurat.obj A \code{Seurat} object with a \code{NNet.mod} list stored in the \code{misc} slot.
#'                   This list is created by \code{\link{prepare.reg}} and updated by
#'                   \code{\link{run.nn.reg}}.
#' @param i Cells for which to extract networks. Can be:
#'   \itemize{
#'     \item \code{NULL} (default): use the cells stored in \code{NNet.mod$cells}.
#'     \item Integer indices matching \code{NNet.setting$all.cells}.
#'     \item A character vector of cell names matching \code{NNet.setting$all.cells}.
#'   }
#' @param assay Which network to extract. One of \code{"effect"}, \code{"p.val"},
#'   or \code{"meta.network"}. If \code{NULL}, uses \code{defaults$assay}.
#' @param remove.self.loops Logical; whether to zero out self-connections. If \code{NULL},
#'   uses \code{defaults$remove.self.loops}.
#' @param responses A character vector of response genes to extract. If \code{NULL},
#'   uses \code{defaults$responses}. 
#' @param predictors A character vector of predictor genes to extract. If \code{NULL},
#'   uses \code{defaults$predictors}. 
#' @param f A function to transform effect values for downstream use (e.g.,
#'   an importance score). If \code{NULL}, uses \code{defaults$f}. Ignored
#'   when \code{assay = "p.val"} or \code{"meta.network"}.
#' @param drop Logical; if \code{TRUE} (default), drops the one slice arraies to
#'   matrices for convenience. 
#' @param cutoff Numeric p-value threshold in [0, 1] for pruning based on 
#'   edge significance stored in \code{NNet.mod$p.val}.
#'   If \code{NULL}, uses \code{defaults$cutoff}. The \code{p.val} tensor will 
#'   be calculated instantly when it is not present in the \code{NNet.mod} list. 
#'
#' @return A network for the requested \code{assay} with dimensions:
#'   \itemize{
#'     \item \strong{effect / p.val}: (responses × predictors × cells) array.
#'       If \code{drop = TRUE} and any dimension equals 1, a matrix is returned:
#'       (responses x predictors) (a single cell),
#'       (predictors x cells) (a single response),
#'       or  (responses x cells) (a single predictor).
#'     \item \strong{meta.network}: the corresponding slice(s) from
#'       \code{NNet.mod$meta.network$meta.network[responses, predictors, i]} with the
#'       same dropping behavior.
#'   }
#'
#' @details
#' \strong{Assay behavior}
#' \itemize{
#'   \item \emph{effect}: If \code{NNet.mod$smoothed = TRUE}, returns the stored smoothed
#'         effects for the requested cells (cells not in \code{NNet.mod$cells} will cause an error). If \code{NNet.mod$smoothed = FALSE}, effects
#'         are projected to the requested cells by applying the stored Laplacian
#'         operator \code{NNet.mod$w}.
#'   \item \emph{p.val}: If \code{NNet.mod$p.val} exists and all requested cells are in
#'         \code{NNet.mod$cells}, returns the stored p-value tensor. Otherwise, p-values are
#'         computed from the (smoothed or projected) effects using per-response
#'         null distribution.
#'   \item \emph{meta.network}: Returns slices of \code{NNet.mod$meta.network$meta.network}.
#' }
#'
#' \strong{Pruning}
#' \itemize{
#'   \item For \emph{effect} assay, pruning uses p-values. 
#'   Effects with corresponding p-values less than \code{cutoff} are zeroed.
#'   \item For \emph{p.val} assay, entries less than \code{cutoff} are set to zero.
#' }
#'
#' \strong{Transformation}
#' \itemize{
#'   \item When \code{assay = "effect"}, the returned tensor/matrix is transformed
#'         by \code{f} (default often \code{function(x) 2*x^2}) to yield an
#'         importance score. No transform is applied for \code{"p.val"} or
#'         \code{"meta.network"}.
#' }
#'
#' @examples
#' mod <- Seurat::Misc(seurat.obj, "NNet.mod")
#'
#' # Default extraction (uses mod$defaults): smoothed effects, default cells
#' net <- get.network(seurat.obj)
#'
#' # Effects for a specific set of cells (by name) and a gene panel
#' cells <- head(names(mod$cells), 3)
#' responses <- head(mod$gene.sets$responses$genes, 5)
#' predictors <- head(mod$gene.sets$predictors$genes, 10)
#' net.eff <- get.network(seurat.obj, i = cells,
#'                        assay = "effect", responses = responses, predictors = predictors)
#'
#' # Extract P-values
#' net.p <- get.network(seurat.obj, assay = "p.val", cutoff = 0.9)
#'
#' # Meta-network slices (3rd component) without dropping dimensions
#' net.meta <- get.network(seurat.obj, assay = "meta.network", i = 3, drop = FALSE)
#'
#' # Effect network with pruning + custom transform, keeping self-loops
#' net.imp <- get.network(seurat.obj,
#'   assay = "effect", cutoff = 0.95,
#'   f = function(x) abs(x),
#'   remove.self.loops = FALSE
#' )
#'
#' @seealso \code{\link{run.nn.reg}}, \code{\link{set.defaults}}
#'
#' @export
get.network <- function(seurat.obj,
                        i = NULL,
                        assay = NULL,
                        remove.self.loops = NULL,
                        responses = NULL,
                        predictors = NULL,
                        f = NULL,
                        drop = TRUE,
                        cutoff = NULL
                        ) {
  mod <- Seurat::Misc(seurat.obj, "NNet.mod")
  if (is.null(mod)) {
    stop("No model found in misc. Run prepare.reg / run.nn.reg first.")
  }

  # Resolve defaults / validate inputs
  # Assay
  if (is.null(assay)) {
    assay <- mod$defaults$assay
  } else {
    assay <- match.arg(assay, choices = c("effect", "p.val", "meta.network"))
  }

  # Cutoff
  cutoff <- mod$defaults$cutoff
  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff)) {
    stop("`cutoff` must be a numeric scalar between 0 and 1.")
  }
  # Gentle clamping with a message (preserves your previous behavior for <0)
  if (cutoff < 0) {
    warning("`cutoff` < 0 detected; clamping to 0.")
    mod$defaults$cutoff <- 0
  } else if (cutoff > 1) {
    warning("`cutoff` > 1 detected; clamping to 1.")
    mod$defaults$cutoff <- 1
  }

  # Validate `remove.self.loops`
  if (is.null(remove.self.loops)) {
    remove.self.loops <- mod$defaults$remove.self.loops
  } 
  if (!is.logical(remove.self.loops) || length(remove.self.loops) != 1L) {
    stop("remove.self.loops must be a single logical (TRUE/FALSE).")
  }

  # Validate `f`
  if (is.null(f)) {
    f <- mod$defaults$f
  }
  if (!is.function(f)) {
    stop("`f` must be a function (e.g., default: function(x) 2*x^2).")
  }

  # Gene panels
  if (is.null(responses)) {
    responses <- mod$defaults$responses
  } else {
    responses <- intersect(responses, mod$gene.sets$responses$genes)
  }

  if (is.null(predictors)) {
    predictors <- mod$defaults$predictors
  } else {
    predictors <- intersect(predictors, mod$gene.sets$genes)
  }

  w <- mod$w
  smoothed <- mod$smoothed
  all.cells <- rownames(w$u)

  # Resolve cell indices
  if (assay != "meta.network") {
    if (is.null(i)) {
      # Use stored cells
      i <- mod$cells
    } else if (is.numeric(i)) {
      # Name numerical indeces with cell names
      names(i) <- all.cells[i]
    } else if (is.character(i)) {
      # Character names → indices via all cell names
      i.tmp <- which(all.cells %in% i)
      names(i.tmp) <- all.cells[i.tmp]
      i <- i.tmp[i]
    } else {
      stop("Invalid index for `i`. Provide NULL, numeric indices, or character cell names.")
    }

    # Can we impute to cells outside mod$cells?
    need_impute <- !all(i %in% mod$cells)
    if (smoothed && need_impute) {
      stop("Unable to impute cell networks: effect tensor is already smoothed.")
    }
  } else {
    # meta.network slices
    if (is.null(i)) {
      i <- 1:dim(mod$meta.network$meta.network)[[3]]
      names(i) <- dimnames(mod$meta.network$meta.network)[[3]]
    } else {
      names(i) <- dimnames(mod$meta.network$meta.network)[[3]][i]
    }
    need_impute <- FALSE
  }

  # Dims
  n.cell <- length(i)
  n.response <- length(responses)
  n.predictor <- length(predictors)

  # Build / fetch the network tensor
  if (assay == "effect") {

    if (smoothed) {
      # Direct subset of smoothed effects
      network <- mod$effect[responses, predictors, names(i), drop = FALSE]

    } else {
      # Project raw effect tensor to requested cells via Laplacian operator
      network <- array(
        0, dim = c(n.response, n.predictor, n.cell),
        dimnames = list(responses, predictors, names(i))
      )
      for (r in responses) {
        network[r, , ] <- tcrossprod(
          mod$effect[r, predictors, ] %*% w$vd[mod$cells, ],
          w$u[i, , drop = FALSE]
        )
      }
    }

  } else if (assay == "p.val") {

    if (!need_impute && !is.null(mod$p.val)) {
      # Stored p-values available for these cells
      network <- mod$p.val[responses, predictors, names(i), drop = FALSE]

    } else if (smoothed) {
      # Compute p-values from smoothed effects
      network <- array(
        0, dim = c(n.response, n.predictor, n.cell),
        dimnames = list(responses, predictors, names(i))
      )
      for (r in responses) {
        network[r, , ] <- pnorm(
          log(abs(mod$effect[r, predictors, names(i)])),
          mean = mod$mus[r], sd = mod$sigmas[r]
        )
      }

    } else {
      # Project effects, then compute p-values
      network <- array(
        0, dim = c(n.response, n.predictor, n.cell),
        dimnames = list(responses, predictors, names(i))
      )
      for (r in responses) {
        projected <- tcrossprod(
          mod$effect[r, predictors, ] %*% w$vd[mod$cells, ],
          w$u[i, , drop = FALSE]
        )
        network[r, , ] <- pnorm(
          log(abs(projected)),
          mean = mod$mus[r], sd = mod$sigmas[r]
        )
      }
    }

  } else { # assay == "meta.network"
    network <- mod$meta.network$meta.network[responses, predictors, i, drop = FALSE]
  }

  # Pruning
  if (cutoff > 0) { 
    if (assay == "effect") {
      # Per-response pruning using p-values
      for (r in responses) {
        if (need_impute || is.null(mod$p.val)) {
          p.val <- network[r, predictors, ] %>%
                    abs %>% 
                    log %>% 
                    pnorm(., mod$mus[r], mod$sigmas[r])
        } else {
          p.val <- mod$p.val[r, predictors, names(i)]
        }
        # Keep edges with p > cutoff; zero-out edges with p <= cutoff
        network[r, , ] <- network[r, , ] * (p.val > cutoff)
      }
    } else if (assay == "p.val") {
      network[network < cutoff] <- 0
    }
    # meta.network has no pruning branch here by design
  }

  # Transform (effects only)
  if (assay == "effect") {
    network <- f(network)
  }

  # Remove self-loops
  if (remove.self.loops) {
    diag_genes <- intersect(responses, predictors)
    if (length(diag_genes)) {
      for (r in diag_genes) network[r, r, ] <- 0
    }
  }

  #Drop length-1 dimensions for convenience
  if (drop) {
    if (n.cell == 1) {
      network <- matrix(network, nrow = n.response,
                        dimnames = list(responses, predictors))
    } else if (n.response == 1) {
      network <- matrix(network, nrow = n.predictor,
                        dimnames = list(predictors, names(i)))
    } else if (n.predictor == 1) {
      network <- matrix(network, nrow = n.response,
                        dimnames = list(responses, names(i)))
    }
  }

  network
}
