# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Prepare Gene and Receptor Settings for NeighbourNet Visualisation
#'
#' This function prepares gene sets, hub genes, receptor priors, and GRN-based
#' evidence required for visualising NeighbourNet co-expression networks. It
#' selects and orders central genes, identifies cluster-specific hub genes
#' based on PC loadings, constructs a receptor–gene prior matrix, and encodes
#' optional prior GRN evidence for edges. The resulting settings are stored in
#' the \code{misc} slot of the \code{Seurat} object for downstream plotting
#' and visualisation routines.
#'
#' @param seurat.obj A \code{Seurat} object with a \code{NNet.mod} list stored in the
#'   \code{misc} slot. This list is created by \code{\link{run.nn.reg}}. 
#' @param n.clu An integer specifying the maximum number of clusters used to
#'   group genes based on their PC loadings when defining hub genes. Default is
#'   \code{4}.
#' @param central.genes An optional list of central genes, typically the output
#'   of \code{\link{select.central.genes}}, containing elements such as
#'   \code{central.responses} and \code{central.predictors}. If provided and
#'   \code{g1} / \code{g2} are \code{NULL}, these central genes are used to
#'   initialise \code{g1} and \code{g2}.
#' @param check.gr.evidence A logical indicating whether prior GRN evidence
#'   should be used to annotate and weight gene–gene relationships. If
#'   \code{TRUE}, TF–target adjacency from \code{\link{get.gr.adj}} and a GRN
#'   graph \code{gr.graph} (in the NeighbourNet namespace) are used to
#'   distinguish supported and unsupported edges. Default is \code{TRUE}.
#' @param t A numeric value passed to \code{\link{get.gr.adj}} that controls
#'   the depth or order of TF–target adjacency used as prior GRN evidence.
#'   Default is \code{2}.
#' @param p A numeric value passed to \code{\link{get.prior.model}} specifying the quantile threshold
#'   for pruning the PPR matrix. If \code{NULL}, the default threshold stored in \code{\link{receptor.ppr}$ltf}
#'   is used. Default is \code{NULL}.
#' @param as.g2 A character string indicating whether \code{g2} (the second
#'   gene layer) should correspond to \code{"predictors"} or \code{"responses"}
#'   in the NeighbourNet model. This choice determines which side of the
#'   effect tensor is treated as the upstream or downstream gene set. Default
#'   is \code{"predictors"}.
#' @param g1 An optional character vector specifying the first gene layer to
#'   visualise (e.g. responses if \code{as.g2 = "predictors"}). If \code{NULL},
#'   it is derived from \code{central.genes} and intersected with the
#'   corresponding default gene set in the NeighbourNet model.
#' @param g2 An optional character vector specifying the second gene layer to
#'   visualise (e.g. predictors if \code{as.g2 = "predictors"}). If \code{NULL},
#'   it is derived from \code{central.genes} and intersected with the
#'   corresponding default gene set in the NeighbourNet model.
#' @param receptors An optional character vector of receptors to visualise.
#'   If \code{NULL}, use all receptors available in the prior knowledge model
#'   obtained from \code{\link{get.prior.model}}.
#'
#' @return
#' A \code{Seurat} object with an additional list named \code{NNet.visual.setting}
#' stored in its \code{misc} slot. This list contains:
#' \item{g1}{The ordered first-layer gene set used for visualisation.}
#' \item{g2}{The ordered second-layer gene set used for visualisation.}
#' \item{clu.g12}{A \code{hclust} object describing hierarchical clustering of
#'                selected genes based on PCs.}
#' \item{clu.g12}{A \code{hclust} object describing hierarchical clustering of
#'                \code{g2} genes based on PCs.}
#' \item{hubs}{A character vector of hub genes, one per cluster, selected to
#'             represent major loading patterns.}
#' \item{receptors}{The subset of receptors to visualise, if any.}
#' \item{ppr}{A receptor–gene prior model (receptors x \code{g2}
#'                set).}
#' \item{as.g2}{The role of \code{g2} in the NeighbourNet model
#'              (\code{"predictors"} or \code{"responses"}).}
#' \item{evidence}{A matrix encoding prior GRN evidence for gene–gene
#'                 relationships between \code{g1} and the \code{g2}
#'                 set, used for visual annotation.}
#'
#' @details
#' This function organises and annotates genes for visualising NeighbourNet
#' co-expression networks. Starting from either user-specified gene sets
#' (\code{g1}, \code{g2}) or central genes obtained from
#' \code{\link{select.central.genes}}, it defines two gene layers aligned with
#' the response and predictor sides of the NeighbourNet model, depending on
#' \code{as.g2}. These genes are projected into the PC space, and clustered via
#' hierarchical clustering on their PC correlations. Within each cluster,
#' a representative hub gene is selected, yielding one hub per cluster that captures
#' dominant variation patterns across PCs.
#'
#' In parallel, a receptor–g2 prior model is extracted and stored as \code{ppr}. It receptors'
#' prior regulatory potential to the second gene layer and can be used to overlay
#' upstream receptor influence in downstream visualisation.
#' 
#' If \code{check.gr.evidence = TRUE}, GRN-based evidence is extrated from \code{\link{get.gr.adj}} 
#' to construct an \code{evidence} matrix that annotates co-expression between g1 and g2 as supported or
#' unsupported by prior regulatory knowledge (and optionally distinguishes
#' activation and inhibition). If \code{check.gr.evidence = FALSE}, all
#' co-expression edges are treated as stimulatory in the visualisation. 
#'
#' @examples
#' # Select central genes from meta-networks
#' sel <- select.central.genes(seurat.obj)
#'
#' # Prepare visualisation settings using central genes (predictors as g2)
#' seurat.obj <- prepare.visualise(
#'   seurat.obj,
#'   n.clu = 4,
#'   as.g2 = "predictors",
#'   central.genes = sel
#' )
#'
#' # Access visual settings
#' vis.set <- Seurat::Misc(seurat.obj, "NNet.visual.setting")
#' str(vis.set)
#'
#' @seealso \code{\link{select.central.genes}}, \code{\link{get.gr.adj}}
#'
#' @export
prepare.visualise <- function(seurat.obj,
                              n.clu = 4,
                              central.genes = NULL,
                              check.gr.evidence = TRUE,
                              t = 2,
                              p = NULL,
                              as.g2 = c("predictors", "responses"),
                              g1 = NULL,
                              g2 = NULL,
                              receptors = NULL) {
  # Retrieve settings and model objects 
  setting <- Seurat::Misc(seurat.obj, "NNet.setting")
  mod <- Seurat::Misc(seurat.obj, "NNet.mod")
  if (is.null(mod)) {
    stop("No model found in misc. Run prepare.reg / run.nn.reg first.")
  }
  
  as.g2 <- match.arg(as.g2)
  assay <- mod$defaults$assay
  
  if(as.g2 == "predictors"){
    if(is.null(g1)) g1 <- central.genes$central.responses
    if(is.null(g2)) g2 <- central.genes$central.predictors
    g1 <- intersect(g1, mod$defaults$responses)
    g2 <- intersect(g2, mod$defaults$predictors)
    g2.full <- mod$defaults$predictors
  }else{
    if(is.null(g1)) g1 <- central.genes$central.predictors
    if(is.null(g2)) g2 <- central.genes$central.responses
    g1 <- intersect(g1, mod$defaults$predictors)
    g2 <- intersect(g2, mod$defaults$responses)
    g2.full <- mod$defaults$responses
  }
  
  # Set the first two layers
  if(mod$custom.y){
    g12 <- g2
  }else{
    g12 <- unique(c(g1, g2))
  }
  d <- apply(setting$pcs, 2, sd)
  pcs <- setting$loadings[g12, ] %>% sweep(., 2, d, "*")
  cor.mat <- cor(t(pcs))
  clu.g12 <- hclust(cor.mat %>% dist) # Will be used to order layouts.
  clu.g2 <- hclust(cor.mat[g2,] %>% dist)  # Will be used to identify hub genes.
  
  # Set hub genes
  if(length(g2) > 1){
    n.clu <- min(n.clu, length(g2))
    hubs <- cbind(data.frame(cluster = cutree(clu.g2, k = n.clu), genes = g2),
                  data.frame(pcs[g2,])
                 ) %>%
            plyr::ddply(., "cluster", function(x){
                      genes <- x[,2]
                      x <- x[,-c(1,2)]
                      v <- svd(as.matrix(x), nu=0, nv=1)$v[,1,drop=FALSE]
                      genes[cor(t(x), v) %>% which.max()]
                   }
                  )
    hubs <- hubs[,2]
  }else{
    hubs <- g2
  }
  
  g12 <- g12[clu.g12$order]
  if(!mod$custom.y) g1 <- g12[g12 %in% g1]
  g2 <- g12[g12 %in% g2]
  
  # Create a receptor to g2 personalized page rank score matrix.
  ppr.full <- NeighbourNet::get.prior.model(p = p)
  receptors.full <- rownames(ppr.full)
  targets.full <- colnames(ppr.full)
  ppr <- Matrix::Matrix(0,
                        nrow = length(receptors.full),
                        ncol = length(g2.full),
                        dimnames = list(receptors.full, g2.full)
                       )
  ppr[,intersect(g2.full, targets.full)] <- ppr.full[,intersect(g2.full, targets.full)]
  
  # Set receptors
  if(!is.null(receptors)){
    receptors <- intersect(receptors, receptors.full)
  }
  
  # Build up evidences
  if(check.gr.evidence){
    evidence <- matrix(1, 
                       nrow = length(g1),
                       ncol = length(g2.full),
                       dimnames = list(g1,g2.full)
                      )
    adj <- t(NeighbourNet::get.gr.adj(t = t))
    GRN.genes <- names(igraph::V(gr.graph))
    known.g1 <- g1[g1 %in% GRN.genes]
    known.g2 <- g2.full[g2.full %in% GRN.genes]
    
    evidence[known.g1,known.g2] <- as.matrix(
        evidence[known.g1,known.g2] +
        (adj[known.g1, known.g2] > 0) + 
        2*(adj[known.g1, known.g2] < 0)
    )
  }else{
    message("check.gr.evidence is set FLASE, all co-expression will be considerred as stimulation.")
    evidence <- matrix(2, 
                       nrow = length(g1),
                       ncol = length(g2.full),
                       dimnames = list(g1,g2.full)
                      )
  }

  visual.setting <- list(g1 = g1,
                         g2 = g2,
                         clu.g12 = clu.g12,
                         clu.g2  = clu.g2,
                         hubs = hubs,
                         receptors = receptors,
                         ppr = ppr,
                         as.g2 = as.g2,
                         evidence = evidence
                        )

  suppressWarnings(
    Seurat::Misc(seurat.obj, "NNet.visual.setting") <- visual.setting
  )
  seurat.obj
}