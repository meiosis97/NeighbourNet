# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Infer Receptor Activity from NeighbourNet TF–Target Networks
#'
#' This function infers receptor, transcription factor (TF), and target activity
#' by combining a receptor–TF prior model (personalised PageRank) with NeighbourNet
#' TF–target co-expression networks. Activity can be evaluated at the per-cell level
#' or on meta-networks, and can be optionally filtered by receptor expression and 
#' prior gene regulatory network (GRN) evidence.
#'
#' @param seurat.obj A \code{Seurat} object with a \code{NNet.mod} list stored in the
#'   \code{misc} slot. This list is created by \code{\link{run.nn.reg}}.
#'   To optionally infer receptor activity on meta-networks (set \code{meta.network = TRUE}), the \code{NNet.mod}
#'   should be updated by a subsequent run of \code{\link{build.meta.network}}.
#' @param i Optional index of cells or meta-network components for which to
#'   compute receptor activity. If \code{NULL}, uses \code{NNet.mod$cells} for
#'   per-cell networks or all meta-network components when
#'   \code{meta.network = TRUE}.
#' @param meta.network A logical indicating whether to compute activity on
#'   NeighbourNet meta-networks instead of per-cell networks. If \code{TRUE},
#'   the \code{"meta.network"} assay is used and each component is treated as a
#'   meta-cell. Default is \code{FALSE}.
#' @param cutoff A numeric value specifying the p-value threshold to prune insignificant
#' co-expression values. If \code{NULL}, the default value is taken from \code{NNet.mod$defaults}.
#' @param check.receptor.expression A logical indicating whether receptor
#'   activity should be down-weighted by receptor expression. Default is
#'   \code{TRUE}.
#' @param check.gr.evidence A logical indicating whether TF–target edges should
#'   be filtered by prior GRN evidence from
#'   \code{\link{get.gr.adj}}. Default is \code{TRUE}.
#' @param t A numeric value passed to \code{\link{get.gr.adj}} that
#'   controls the depth or order of TF–target adjacency used as prior GRN
#'   evidence. Default is \code{2}.
#' @param p A numeric value passed to \code{\link{get.prior.model}} specifying the quantile threshold
#' for pruning the PPR matrix. If \code{NULL}, the default threshold stored in \code{\link{receptor.ppr}$ltf}
#' is used. Default is \code{NULL}.
#' @param scale.ppr A logical indicating whether to normalise the receptor–TF
#'   prior model (PageRank matrix) by receptor-wise norms. Default is \code{TRUE}.
#' @param scale.network A logical indicating whether to normalise the TF–target 
#'   co-expression network by target-wise norms before computing activity. Default is
#'   \code{TRUE}.
#' @param as.tfs A character string indicating whether TFs are encoded as
#'   \code{"predictors"} (columns) or \code{"responses"} (rows) in the
#'   NeighbourNet networks. Default is \code{"predictors"}.
#' @param receptors A character vector of receptor genes to include. If
#'   \code{NULL}, all receptors present in the prior model returned by
#'   \code{\link{get.prior.model}} are used.
#' @param tfs A character vector of TFs to include. If \code{NULL},
#'   TFs are taken from the default gene set stored in \code{NNet.mod$defaults}.
#' @param targets A character vector of target genes to include. If \code{NULL},
#'   targets are taken from the default gene set stored in \code{NNet.mod$defaults}.
#' @param receptor.activity A character string specifying how receptor–TF and
#'   TF–target information should be combined. Options are \code{"cprod"}
#'   (cross-product aggregration that measures cosine similarity, default)
#'   and \code{"dist"} (max-over-path aggregation).
#'
#' @return
#' If multiple cells or meta-network components are requested (\code{length(i) > 1}
#' or \code{i = NULL}), the function returns a list with:
#' \item{receptor.act}{A matrix of receptor activity scores
#'                     (receptors x cells/components).}
#' \item{tf.act}{A matrix of TF activity scores (TFs x cells/components).}
#' \item{target.act}{A matrix of target activity scores
#'                   (targets x cells/components).}
#'
#' If a single cell or component is requested, the function returns a list with:
#' \item{act.mat}{A receptor–target activity matrix (receptors x targets).}
#' \item{ppr}{The receptor–TF prior used after any scaling and expression
#'            weighting.}
#' \item{network}{The TF–target network used for activity computation.}
#'
#' @details
#' This function integrates three hierarchical components to infer receptor activity:
#' (i) a receptor–TF prior model obtained from \code{\link{get.prior.model}},
#' (ii) TF–target co-expression networks derived from NeighbourNet, and
#' (iii) optional TF–target adjacency matrices from \code{\link{get.gr.adj}} as
#' prior gene regulatory network (GRN) evidence. Together, these components
#' define a receptor–TF–target signalling cascade that quantifies how upstream
#' receptors potentially influence downstream target genes through TFs.
#'
#' For each cell or meta-network component, the receptor–TF prior model encodes
#' the probability of regulatory transmission from receptors to TF, 
#' while the TF–target network represents local co-expression
#' relationships inferred by NeighbourNet. Receptor activity is computed by
#' aggregating these two layers to estimate receptor influence on downstream
#' targets, mediated by TFs. The resulting matrices describe receptor-, TF-, and
#' target-level activity for each cell or component, providing an interpretable
#' summary of upstream signalling dynamics.
#'
#' Activity can be inferred from per-cell networks (default) or from
#' meta-networks when \code{meta.network = TRUE}. The arguments \code{scale.ppr} and
#' \code{scale.network} apply L2-normalisation to the prior model and
#' co-expression network respectively, reducing the influence of highly connected
#' nodes and ensuring comparable contribution across genes. The options
#' \code{check.receptor.expression} and \code{check.gr.evidence} introduce two
#' filters on receptor activity : the first enforces that receptors must be transcriptionally
#' detected, while the second prunes TF–target co-expression to those 
#' supported by prior GRN knowledge.
#'
#' The argument \code{receptor.activity} controls how the two layers are
#' combined: \code{"cprod"} performs a cross-product aggregation equivalent to
#' measuring cosine similarity between receptor–TF and TF–target profiles,
#' while \code{"dist"} applies a max-over-path operation that highlights the
#' strongest possible receptor–target link through any intermediate TF.
#'
#' @examples
#' # Per-cell receptor activity using default settings
#' act <- receptor.activity(seurat.obj)
#'
#' # Meta-network receptor activity for the first three components
#' act <- receptor.activity(
#'   seurat.obj,
#'   meta.network = TRUE,
#'   i = 1:3
#' )
#' 
#'
#' @seealso  \code{get.prior.model}, \code{get.gr.adj}
#'
#' @export
receptor.activity <- function(seurat.obj,
                              i = NULL,
                              meta.network = FALSE,
                              cutoff = NULL,
                              check.receptor.expression = TRUE,
                              check.gr.evidence = TRUE,
                              t = 2,
                              p = NULL,
                              scale.ppr = TRUE,
                              scale.network = TRUE,
                              as.tfs = c("predictors", "responses"),
                              receptors = NULL,
                              tfs = NULL,
                              targets = NULL,
                              receptor.activity = c("cprod", "dist")) {
  # 1. Retrieve settings and model objects 
  setting <- Seurat::Misc(seurat.obj, "NNet.setting")
  mod <- Seurat::Misc(seurat.obj, "NNet.mod")
  if (is.null(mod)) {
    stop("No model found in misc. Run prepare.reg / run.nn.reg first.")
  }

  # 2. Resolve arguments and defaults
  as.tfs <- match.arg(as.tfs)
  receptor.activity <- match.arg(receptor.activity)

  assay <- mod$defaults$assay
  if (meta.network) {
    assay <- "meta.network"
  }
  f <- mod$defaults$f

  # Cutoff: fall back to defaults and clamp lower bound
  if (is.null(cutoff)) {
    cutoff <- mod$defaults$cutoff
  } else if (!is.numeric(cutoff)) {
    stop("cutoff must be numeric.")
  }
  if (cutoff < 0) {
    cutoff <- 0
  }

  # Resolve cells or components
  if (is.null(i)) {
    if (assay != "meta.network") {
      cells <- mod$cells
    } else {
      n.net <- dim(mod$meta.network$meta.network)[[3]]
      cells <- 1:n.net
      names(cells) <- paste0("component_", cells)
    }
  } else if (length(i) > 1) {
    if (assay != "meta.network") {
      cells <- i
      names(cells) <- setting$all.cells[i]
    } else {
      cells <- i
      names(cells) <- paste0("component_", cells)
    }
  } else {
    cells <- i
  }
  n.cells <- length(cells)

  # 3. Resolve TFs and targets from defaults and gene.list
  if (as.tfs == "predictors") {
    if (is.null(tfs)) {
      tfs <- mod$gene.sets$predictors$tfs
    }
    tfs <- intersect(tfs, mod$defaults$predictors)
    if (is.null(targets)) {
      targets <- mod$defaults$responses
    }
    targets <- intersect(targets, mod$defaults$responses)
  } else {
    if (is.null(tfs)) {
      tfs <- mod$gene.sets$responses$tfs
    }
    tfs <- intersect(tfs, mod$defaults$responses)
    if (is.null(targets)) {
      targets <- mod$defaults$predictors
    }
    targets <- intersect(targets, mod$defaults$predictors)
  }

  # Restrict to genes recognised by prior knowledge (gene.list)
  targets <- intersect(targets, NeighbourNet::gene.list$targets)
  tfs <- intersect(tfs, NeighbourNet::gene.list$tfs)
  n.tfs <- length(tfs)
  n.targets <- length(targets)

  # 4. Build receptor–TF PageRank matrix (ppr)
  ppr.full <- NeighbourNet::get.prior.model(p = p)
  receptors.full <- rownames(ppr.full)
  mediators <- intersect(colnames(ppr.full), tfs)

  ppr <- Matrix::Matrix(
    0,
    nrow = length(receptors.full),
    ncol = n.tfs,
    dimnames = list(receptors.full, tfs)
  )

  if (length(mediators)) {
    ppr[, mediators] <- ppr.full[, mediators, drop = FALSE]
  }

  # Restrict to requested receptors (if any)
  if (is.null(receptors)) {
    receptors <- receptors.full
  } else {
    receptors <- intersect(receptors.full, receptors)
  }
  ppr <- ppr[receptors, , drop = FALSE]
  n.receptors <- length(receptors)

  # Optional Laplacian-like scaling on receptor–TF PageRank
  if (scale.ppr) {
    receptor.scale <- sqrt(rowSums(ppr^2))
    ppr <- ppr / receptor.scale
    ppr[is.na(ppr)] <- 0
  }

  # 5. Receptor expression check and propagation 
  if (check.receptor.expression) {
    counts <- SeuratObject::LayerData(seurat.obj, "counts")
    all.cells <- setting$all.cells
    all.genes <- setting$all.genes
    detected.receptors <- intersect(all.genes, receptors)

    # Binary detection matrix: receptors × pca.cells
    expr.mask <- Matrix::Matrix(
      0,
      nrow = n.receptors,
      ncol = length(all.cells),
      dimnames = list(receptors, all.cells)
    )
    if (length(detected.receptors) > 0) {
      expr.mask[detected.receptors, ] <- counts[detected.receptors, all.cells, drop = FALSE] != 0
    }

    # Network propagation using NeighbourNet Laplacian
    if (assay != "meta.network") {
      w <- tcrossprod(mod$w$u[cells, , drop = FALSE], mod$w$vd)
      expr.mask <- tcrossprod(expr.mask, w)
    } else {
      w <- mod$w$vd %*%
        crossprod(mod$w$u[mod$cells, , drop = FALSE],
                  mod$meta.network$npca.loadings[, cells, drop = FALSE])
      w <- sweep(w, 2, colSums(w), "/")
      expr.mask <- expr.mask %*% w
    }

    # Clamp to [0, 1]
    expr.mask[expr.mask > 1] <- 1
    expr.mask[expr.mask < 0] <- 0

  }

  # 6. Prior GRN evidence (TF–target adjacency)
  if (check.gr.evidence) {
    adj <- t(NeighbourNet::get.gr.adj(t = t))
  }

  # 7. Pre-allocate outputs for multi-cell/component case
  if (n.cells > 1) {
    receptor.act <- data.frame(
      matrix(0, nrow = n.receptors, ncol = n.cells),
      row.names = receptors
    )
    colnames(receptor.act) <- names(cells)

    tf.act <- data.frame(
      matrix(0, nrow = n.tfs, ncol = n.cells),
      row.names = tfs
    )
    colnames(tf.act) <- names(cells)

    target.act <- data.frame(
      matrix(0, nrow = n.targets, ncol = n.cells),
      row.names = targets
    )
    colnames(target.act) <- names(cells)

    ppr.colsums <- colSums(ppr)

    # TF 
    message("Now infer receptor activity.")

    # Retrieved from @https://www.dummies.com/article/technology/programming-web-design/r/how-to-generate-your-own-error-messages-in-r-175112/
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                    total = n.cells,
                                    complete = "=",   # Completion bar character
                                    incomplete = "-", # Incomplete bar character
                                    current = ">",    # Current bar character
                                    clear = FALSE,    # If TRUE, clears the bar when finish
                                    width = 100)      # Width of the progress bar


  }

  # 8. Loop over cells/components and compute activity 
  for (i in seq_along(cells)) {
    # 8.1 Extract TF–target network for this cell/component
    network <- NeighbourNet::get.network(
      seurat.obj,
      i = cells[i],
      drop = TRUE,
      assay = assay,
      cutoff = cutoff
    )
    responses <- rownames(network)
    predictors <- colnames(network)

    # Meta-network-specific p-value pruning
    if (assay == "meta.network" && cutoff > 0) {
      if (is.null(mod$meta.network$p.val)) {
        stop("Meta-network p-val not found.")
      }
      prob.network <- mod$meta.network$p.val[
        mod$defaults$responses,
        mod$defaults$predictors,
        cells[i],
        drop = FALSE
      ]
      network[prob.network <= cutoff] <- 0
    }

    # Orient network so that TFs lie in columns
    if (as.tfs == "responses") {
      network <- t(network)
    }

    # Restrict to targets and TFs and optionally enforce GRN evidence
    network <- network[targets, tfs, drop = FALSE]
    if (check.gr.evidence) {
      network <- network * (adj[targets, tfs, drop = FALSE] != 0)
    }

    # 8.2 Scale TF–target network if requested
    if (scale.network) {
      target.scale <- sqrt(rowSums(network^2))
      network <- network / target.scale
      network[is.na(network)] <- 0
    }

    # 8.3 Combine receptor–TF and TF–target scores
    if (n.cells > 1) {
      pb$tick()

      # Multi-cell/component case
      if (receptor.activity == "cprod") {
        # Column sums of TF–target network
        network.colsums <- colSums(network)

        ppr.i <-ppr
        if (check.receptor.expression) {
          ppr.i <- ppr * expr.mask[, i]
        }
        ppr.colsums <- colSums(ppr.i)

        receptor.act[, i] <- as.numeric(ppr.i %*% network.colsums)
        tf.act[, i] <- as.numeric(ppr.colsums * network.colsums)
        target.act[, i] <- as.numeric(network %*% ppr.colsums)

      } else {
        # "dist" metric: max-over-paths style scoring
        network.colmax <- apply(network, 2, max)

        ppr.i <-ppr
        if (check.receptor.expression) {
          ppr.i <- ppr * expr.mask[, i]
        }
        ppr.colsums <- apply(ppr.i, 2, max)

        # Receptor activity: max over TFs of (receptor–TF × TF strength)
        receptor.act[, i] <- apply(
          sweep(ppr.i, 2, network.colmax, "*"),
          1,
          max
        )
        tf.act[, i] <- as.numeric(network.colmax * ppr.colmax)
        target.act[, i] <- apply(
          sweep(network, 2, ppr.colmax, "*"),
          1,
          max
        )

      }
    } else {
      if (check.receptor.expression) {
        ppr <- ppr * as.numeric(expr.mask)
      }

      # Single cell/component case
      if (receptor.activity == "cprod") {
        act.mat <- tcrossprod(ppr, network)
      } else {
        act.mat <- matrix(
          0,
          nrow = n.receptors,
          ncol = n.targets,
          dimnames = list(receptors, targets)
        )
        for (i in targets) {
          act.mat[, i] <- apply(
            sweep(ppr, 2, network[i, ], "*"),
            1,
            max
          )
        }
      }
    }
  }

  # 9. Return results
  if (n.cells > 1) {
    return(invisible(list(
      receptor.act = receptor.act,
      tf.act = tf.act,
      target.act = target.act
    )))
  } else {
    return(invisible(list(
      act.mat = act.mat,
      ppr = ppr,
      network = network
    )))
  }
}
