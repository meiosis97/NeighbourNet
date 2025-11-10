# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Visualise Receptor–TF–Target Pathways on NeighbourNet Cell-Specific Networks.
#'
#' This function visualises NeighbourNet co-expression networks as a layered
#' receptor–TF–target signalling diagram. It takes pre-computed visualisation
#' settings from \code{\link{prepare.visualise}} and, for a given cell or
#' meta-network component, builds a multi-layer graph linking upstream receptors
#' to downstream gene modules. Edges are weighted by co-expression strength and
#' annotated by prior GRN evidence, and node shapes are represented as
#' cluster-level pie charts.
#' \cr
#' \cr
#' \strong{Interpretation}
#' When NeighbourNet is used to build TF–target co-expression networks
#' (for example, targets as responses and TFs as predictors), this function
#' provides prior knowledge enriched visualisation of networks as four concentric
#' layers. From the innermost to the outermost layer, these are:
#' \itemize{
#'   \item \code{g1} (targets): target genes, placed at the centre of the plot.
#'   \item \code{g2} (TFs): transcription factors surrounding the targets.
#'         TFs are grouped by expression-based clustering (node colour), and
#'         each TF is assigned a significance profile summarised as a pie chart
#'         (the fraction of the pie filled reflects the likelihood of
#'         co-expression with the \code{g1} targets).
#'   \item \code{g3} (intermediate signalling nodes): optional intermediate
#'         signalling molecules placed between TFs and receptors. These nodes
#'         represent the shortest signalling paths that connect receptors to
#'         TFs in a prior signalling graph.
#'   \item \code{g4} (receptors): the most active upstream receptors on the
#'         outermost layer. Each receptor is linked to a TF that is inferred to
#'         mediate its influence on the targets. Receptors are coloured in
#'         accordance with TF clusters, and the proportion of colour fill
#'         reflects their relative activity or expression.
#' }
#' The connections between \code{g1} and \code{g2} encode TF–target
#' co-expression annotated by prior gene-regulatory evidence. Edges represent
#' three levels of support: no edge (no significant relationship), dashed edges
#' for significant co-expression only, and solid edges when significant
#' co-expression is also supported by prior GRN evidence (with different
#' arrow styles indicating activation versus repression). The \code{g3} and
#' \code{g4} layers summarise upstream signalling pathways that exert strong
#' influence on the observed TF–target co-expression patterns.

#'
#' @param seurat.obj A \code{Seurat} object with a \code{NNet.visualise.setting}
#'   list stored in the \code{misc} slot. This list is created by \code{\link{prepare.visualise}}. 
#' @param i Index of the cell or meta-network component to visualise. 
#' @param meta.network A logical indicating whether to visualise a NeighbourNet
#'   meta-network component instead of a per-cell network. If \code{TRUE}, the
#'   \code{"meta.network"} assay is used. Default is \code{FALSE}.
#' @param fix.cluster A logical indicating whether to reuse the gene clustering
#'   from \code{NNet.visualise.setting} when assigning clusters to the second
#'   gene layer (\code{g2}). If \code{FALSE}, clusters are re-estimated from the
#'   current network. Default is \code{TRUE}.
#' @param hubs An optional character vector of hub genes in \code{g2}. If
#'   supplied and \code{fix.cluster = FALSE}, clusters are inferred by assigning
#'   each \code{g2} gene to the closest hub. If \code{NULL}, hubs are not used
#'   for clustering.
#' @param n.clu An integer specifying the maximum number of clusters for
#'   \code{g2} genes. Default is \code{4}. The effective number of clusters is
#'   truncated when few genes are available.
#' @param cutoff A numeric value specifying the p-value threshold used to define
#'   significant co-expression for edge evidence. If \code{NULL}, the default is
#'   taken from \code{NNet.mod$defaults$cutoff} stored in the \code{misc} of the Seurat object.
#'   Values less than 0 are clamped to 0.
#' @param show.pathways A logical indicating whether receptor layers and
#'   intermediate signalling paths (via a prior signalling graph) should be
#'   included in the visualisation. If \code{FALSE}, only the two gene layers
#'   (\code{g1} and \code{g2}) are plotted. Default is \code{TRUE}.
#' @param change.receptors A logical indicating whether to select a subset of
#'   receptors for plotting based on activity (e.g., eigen-like centrality on
#'   inferred receptor–target activity). If \code{FALSE}, receptors supplied in
#'   \code{NNet.visualise.setting} are used as is. Default is \code{TRUE}.
#' @param receptor.activity A character string specifying how the prior model of 
#'   receptor–TF interaction and TF–target (\code{g1}–\code{g2}) co-exression network are
#'   combined when computing receptor–g1 activity. Options are \code{"cprod"}
#'   (cross-product aggregation) and \code{"dist"} (max-over-path aggregation).
#'   Default is \code{"cprod"}.
#' @param check.receptor.expression A logical indicating whether receptor
#'   activity should be weighted by receptor expression in the selected cell or
#'   component. Default is \code{TRUE}.
#' @param scale.ppr A logical indicating whether to normalise the receptor–TF
#'   prior model (PageRank matrix) by receptor-wise norms. Default is \code{TRUE}.
#' @param scale.network A logical indicating whether to normalise the 
#'   co-expression network by \code{g1}-wise norms before computing activity. Default is
#'   \code{TRUE}.
#' @param scale.signifiance A logical indicating whether node and edge weights
#'   should be scaled by their maximum significance (derived from p-values),
#'   rather than raw co-expression magnitudes. Default is \code{FALSE}.
#' @param swap.layers A logical indicating whether to reverse the order of
#'   \code{g1} and \code{g2} in the concentric layout (inner vs outer ring).
#'   Default is \code{FALSE}.
#' @param k An integer controlling the number of singular vectors used when
#'   selecting active receptors (layer \code{g4}) based on receptor–target
#'   activity, analogous to an eigenvector centrality based selection in 
#'   \code{\link{select.central.genes}} Default is \code{2}.
#' @param n.per.component An integer indicating the number of top receptors per
#'   component (or singular vector) to consider when selecting candidate
#'   receptors. Default is \code{10}.
#' @param radius Optional numeric vector control the radii of concentric layers.
#' @param pie.radius A numeric value giving the base radius of node pies in
#'   the scatterpie representation. Default is \code{0.05}.
#' @param text.size A numeric value giving the size of gene labels in the
#'   plotted network. Default is \code{4}.
#'
#' @return A \code{ggplot} object (built with \pkg{ggraph} and
#'   \pkg{ggplot2}) representing a layered receptor–TF–target network. The plot
#'   includes:
#' \itemize{
#'   \item concentric layers for \code{g1} and \code{g2}, and optionally
#'         intermediate signalling nodes (\code{g3}) and receptors (\code{g4});
#'   \item edge arrows and line types encoding prior evidence and direction;
#'   \item node pies summarising cluster assignments and activity contributions.
#' }
#' The underlying graph structure is built internally using \pkg{igraph}, but
#' only the final \code{ggplot} object is returned.
#'
#' @details 
#' This function visualises a single NeighbourNet network (either a cell-level
#' network or a meta-network component) using gene sets and priors prepared by
#' \code{\link{prepare.visualise}}. The first two layers (\code{g1} and
#' \code{g2}) correspond to response and predictor gene sets. Co-expression between
#' \code{g1} and \code{g2} is extracted from the selected NeighbourNet assay
#' (effect, p-value, or meta-network) via \code{\link{get.network}}, with
#' associated p-values used to define an edge-specific significance measure.
#'
#' Node weights for \code{g1} and \code{g2} are derived from the maximum
#' significance of their incident edges and can be optionally rescaled by
#' \code{scale.signifiance}. Edge weights are then aggregated into cluster-level
#' pies, showing how strongly each layer contributes to different \code{g2}
#' modules. The \code{evidence} matrix from \code{NNet.visual.setting} encodes prior GRN evidence
#' for \code{g1}–\code{g2} pairs and is used to annotate edges as stimulatory
#' or inhibitory where such information is available.
#'
#' When \code{show.pathways = TRUE}, the function also incorporates a
#' receptor–\code{g2} prior (generated by \code{\link{get.prior.model}}) and a signalling graph
#' (generated by \code{\link{get.gr.adj}}) to construct upstream receptor layers. 
#' Receptor-to-\code{g1} activity is first inferred by combining the receptor–\code{g2} prior model with the
#' pruned \code{g1}–\code{g2} co-expression network, and
#' weighted by receptor expression in the selected cell or component
#' (\code{check.receptor.expression}). A subset of receptors (\code{g4}) and
#' intermediate signalling nodes (\code{g3}) is then chosen based on activity
#' and shortest paths in \code{sig.graph}, and added as outer layers in the
#' layout.
#'
#' Finally, the multi-layer graph is laid out using a concentric layout, with
#' the order of layers controlled by \code{swap.layers} and \code{radius}, and
#' rendered via \pkg{ggraph}. Pies from \pkg{ggforce::geom_scatterpie} summarise
#' cluster-level contributions for each node, while edge arrows and line types
#' encode directionality and prior evidence. 
#' 
#' 
#'
#' @examples
#' # Assume NeighbourNet has been fitted and visual settings prepared:
#' # seurat.obj <- prepare.graph(...)
#' # seurat.obj <- prepare.reg(seurat.obj, ...)
#' # seurat.obj <- run.nn.reg(seurat.obj, ...)
#' # seurat.obj <- build.meta.network(seurat.obj, ...)
#' # seurat.obj <- prepare.visualise(seurat.obj)
#'
#' # Visualise the network for a single cell
#' p <- visualise.network(seurat.obj, i = 1)
#' print(p)
#'
#' # Visualise a meta-network component with receptor pathways
#' p.meta <- visualise.network(
#'   seurat.obj,
#'   i = 1,
#'   meta.network = TRUE,
#' )
#' print(p.meta)
#'
#' @seealso \code{\link{prepare.visualise}}, \code{\link{receptor.activity}}
#'
#' @export
visualise.network <- function(seurat.obj,
                              i,
                              meta.network = FALSE,
                              fix.cluster = TRUE,
                              hubs = NULL,
                              n.clu = 4,
                              cutoff = NULL,
                              show.pathways = TRUE,
                              change.receptors = TRUE,
                              receptor.activity = c("cprod", "dist"),
                              check.receptor.expression = TRUE,
                              scale.ppr = FALSE,
                              scale.network = FALSE,
                              scale.signifiance = FALSE,
                              swap.layers = FALSE,
                              k = 2,
                              n.per.component = 10,
                              radius  = NULL,
                              pie.radius = 0.05,
                              text.size = 4){
  # 1. Retrieve settings and model objects 
  setting <- Seurat::Misc(seurat.obj, "NNet.setting")
  mod <- Seurat::Misc(seurat.obj, "NNet.mod")
  visual.setting <- Seurat::Misc(seurat.obj, "NNet.visual.setting")
  if (is.null(mod)) {
    stop("No model found in misc. Run prepare.reg / run.nn.reg first, and then prepare.visualise.")
  }
  if(is.null(visual.setting)) stop("Run prepare.visualise first.")
  
  if(is.null(cutoff)){
    cutoff <- mod$defaults$cutoff
  }else if(!is.numeric(cutoff)){
    stop("cutoff must be numerical and between 0 and 1.")
  }
  if(cutoff < 0) cutoff <- 0
  
  g1 <- visual.setting$g1
  g2 <- visual.setting$g2
  as.g2 <- visual.setting$as.g2
  evidence <- visual.setting$evidence
  ppr <- visual.setting$ppr
  receptors <- visual.setting$receptors
  receptors.full <- rownames(ppr)
  assay <- if(meta.network) "meta.network" else mod$defaults$assay

  # 2. ######################## network processing ########################
  # Get the network
  network <- get.network(seurat.obj, i = i, drop = TRUE, assay = assay, cutoff = 0)
  responses <- rownames(network)
  predictors <- colnames(network)
  
  # Get the probability network
  if(assay == "effect"){
    prob.network <- get.network(seurat.obj, i = i, drop = TRUE, assay = "p.val", cutoff = 0)
  }else if(assay == "p.val"){
    prob.network <-  network
  }else{
    if(is.null(mod$meta.network$p.val)) stop("Meta-network p-val not found.")
    prob.network <- mod$meta.network$p.val[responses,predictors,i]
    prob.network <- matrix(prob.network,
                           nrow = length(responses),
                           ncol = length(predictors),
                           dimnames = list(responses, predictors))
  }
  
  if(as.g2 == "responses"){
    network <- t(network)
    prob.network <- t(prob.network)
  }
  
  g1 <- intersect(g1, rownames(network))
  g2 <- intersect(g2, colnames(network))

  # Allow duplicated gene names in the first and second layer
  names(g1) <- paste("g1", g1, sep = ".")
  names(g2) <- paste("g2", g2, sep = ".")
  
  # g1 g2 node weight
  g1.weight <- apply(prob.network, 1, max)
  g2.weight <- apply(prob.network, 2, max)
  
  # Scale the weight
  if(scale.signifiance){
    g1.weight <- g1.weight/max(g1.weight)
    g2.weight <- g2.weight/max(g2.weight)
    network <- network/max(network)
  }else if(assay != "p.val"){
    network <- network/max(network)*max(prob.network)
  }
  
  g1.weight <- g1.weight[g1]
  g2.weight <- g2.weight[g2]
  network <- network[g1,g2,drop=FALSE]
  evidence <- evidence[g1,g2,drop=FALSE]
  prob.network <- prob.network[g1,g2,drop=FALSE]
  dimnames(prob.network) <- dimnames(network) <- list(names(g1), names(g2))
  names(g1.weight) <- names(g1)
  names(g2.weight) <- names(g2)
  
  # Check evidence for the g2 to g1 network.
  evidence <- (prob.network > cutoff) * evidence + 1
  
  # 3. ################# Hierarchical clustering g2 genes #################
  if(length(g2)>1){
    n.clu <- min(n.clu, length(g2))
    if(fix.cluster){
      clu.g2 <- cutree(visual.setting$clu.g2, n.clu)[g2]
      clu.g2 <- as.factor(clu.g2) %>% as.numeric()
      n.clu <- max(clu.g2)
      names(clu.g2) <- names(g2)
    }else if(is.null(hubs)){
      if(length(g1)>1){
        clu.g2 <- hclust(cor(network) %>% dist)
      }else{
        clu.g2 <- hclust(network %>% dist)
      }
      clu.g2 <- cutree(clu.g2,n.clu)[names(g2)]
    }else{
      hubs <- intersect(hubs, g2)
      names(hubs) <- paste("g2", hubs, sep = ".")
      if(length(g1)>1){
        clu.g2 <- cor(network[,names(hubs),drop=FALSE], network)
        clu.g2 <- apply(clu.g2, 2, which.max)
      }else{
        clu.g2 <-  sapply(names(hubs), function(x) sqrt(colSums(network-network[,x])^2))
        clu.g2 <- apply(clu.g2, 1, which.min)
      }
      n.clu <- length(hubs)
    }
  }else{
    clu.g2 <- 1
  }
  
  names(clu.g2) <- names(g2)
  
  # 4. ################# Receptor to target activity #################
  if(show.pathways){
    pruned.network <- network * (evidence > 2)
    ppr <- ppr[,g2,drop=FALSE]
    
    if(scale.network){
      g1.scale <- sqrt(rowSums(pruned.network^2))
      pruned.network <- pruned.network/g1.scale
      pruned.network[is.na(pruned.network)] <- 0
    }
    
    if(scale.ppr){
      receptor.scale <- sqrt(rowSums(ppr^2))
      ppr <- ppr / receptor.scale
    }

    ppr <- ppr/max(ppr)
    ppr[is.na(ppr)] <- 0

    # Whether to weight receptors by their expression.
    if(check.receptor.expression){
      counts <- SeuratObject::LayerData(seurat.obj, "counts")
      all.cells <- setting$all.cells
      all.genes <- rownames(counts)
      detected.receptors <- intersect(all.genes, receptors.full)
      
      # Check receptor expression
      expr.mask <- Matrix::Matrix(0,
                           nrow = length(receptors.full),
                           ncol = length(all.cells), 
                           dimnames = list(receptors.full, all.cells))
      expr.mask[detected.receptors, ] <- counts[detected.receptors, all.cells] != 0
      
      # Network propagation
      if(assay != "meta.network"){
        w <- tcrossprod(mod$w$u[i,,drop=FALSE], mod$w$vd)
        expr.mask <- tcrossprod(expr.mask, w)
      }else{
        w <-  mod$w$vd %*%
          crossprod(mod$w$u[mod$cells,], mod$meta.network$npca.loadings[,i])
        w <-  w/sum(w)
        expr.mask <- expr.mask %*% w
      }
      expr.mask[expr.mask > 1] <- 1
      expr.mask[expr.mask < 0] <- 0
      ppr <- ppr * as.numeric(expr.mask)
    }
    
    receptor.activity <- match.arg(receptor.activity)
    if(receptor.activity == "cprod"){
      act <- tcrossprod(ppr, pruned.network)
    }else{
      act <- matrix(0, nrow = length(receptors.full), ncol = length(g1),
                    dimnames = list(receptors.full, g1))
      for(j in g1){
        act[,j] <- sweep(ppr, 2, pruned.network[j,], "*") %>% apply(., 1, max)
      }
    }
    
    # Check maximum activity
    max.act <- apply(act, 2, function(x) max(x))
    
    k.min <- min(length(g1), nrow(ppr), sum(colSums(act) > 0))
    k <- min(k, k.min)
    if(k == 0) show.pathways <- FALSE
    max.act <- max.act[colSums(act)>0]
    act <- act[,colSums(act)>0,drop=FALSE]

    # Find g4
    if(k == 0){
      show.pathways <- FALSE
    }else{
      if(change.receptors){
        # Eigen centrality.
        if(k.min >= 3){
          svd.mod <- RSpectra::svds(act, k = k) %>% suppressWarnings()
          g4 <- receptors.full[apply(abs(svd.mod$u), 2, NeighbourNet::topn, n = n.per.component) %>%
                               as.numeric
                              ]
        }else{
          g4 <- receptors.full[apply(act, 2, NeighbourNet::topn, n = n.per.component) %>%
                               as.numeric
                              ]
        }
        g4 <- unique(g4)
        
      }else{
        g4 <- unique(receptors)
      }
        
      # Receptor weight
      receptor.weight <- sweep(act, 2, max.act, "/") %>%
                        apply(., 1, max)
      active.receptors <- names(receptor.weight)[receptor.weight>0]
      g4 <- intersect(g4, active.receptors)
      if(length(g4) == 0) show.pathways <- FALSE
    }
  }
  
  # 5. ################# Evidence processing   #################
  # Check evidence for the receptor to TF network.
  check.evidence <- function(graph){
    from <- graph[,1]
    to <- graph[,2]
    consensus_stimulation <- igraph::E(NeighbourNet::sig.graph)$consensus_stimulation
    consensus_inhibition <- igraph::E(NeighbourNet::sig.graph)$consensus_inhibition
    
    evidence <- rep(2, nrow(graph))
    e.ids <- igraph::get.edge.ids(NeighbourNet::sig.graph,
                                  data.frame(from, to) %>% t)
    
    evidence <- evidence + consensus_stimulation[e.ids] + 2*consensus_inhibition[e.ids]
    evidence[evidence>4] <- 2
    evidence
  }
  
  # 6. Build network
  ################# g2 to g1 layer   #################
  network.df <- reshape2::melt(t(network))
  colnames(network.df) <- c("from","to","weight")
  network.df[,1:2] <- apply(network.df[,1:2], 2, as.character)
  network.df$cluster <- 1 + clu.g2[network.df[,1]]
  network.df$evidence <- reshape2::melt(t(evidence))[,3]
  network.df$weight[network.df$weight < 0.05] <- 0.05
  graph <- igraph::graph_from_data_frame(network.df,directed = T)
  
  if(show.pathways){
    ################# g4 to g3 layer  #################
    # Check activated TFs
    activeTF <- (colSums(pruned.network) >0) & (colSums(ppr) >0)

    # Create g3 to receptor network
    network1 <- sweep(ppr, 2, g2.weight * activeTF,"*")
    network1 <- network1[active.receptors,,drop=FALSE]

    # Add extra g4
    activeTF <- g2[activeTF]
    closest.receptors <- active.receptors[apply(network1[,activeTF,drop=FALSE], 2, which.max)]
    
    # Get closest TFs
    closestTF <- g2[apply(network1[g4,,drop=FALSE], 1, which.max)]
    
    # Create g3 to g4 network
    g4 <- c(g4, closest.receptors)
    closestTF <- c(closestTF, activeTF)

    # Find shortest paths between TFs and their closest receptors.
    paths <- c()
    for (i in seq(length(g4))){
      r <- g4[i]
      sp <- igraph::shortest_paths(NeighbourNet::sig.graph,
                                   r,
                                   closestTF[i],
                                   weights = 1/igraph::E(NeighbourNet::sig.graph),
                                   output = c("both"))
      paths <- c(paths , sp$vpath)
    }
    names(paths) <- closestTF
    paths <- paths[order(sapply(closestTF, function(x) which(g2==x)[1]))]
    closestTF <- names(paths)
    pass1 <- sapply(paths, function(x) names(x)[2])
    pass2 <- sapply(paths, function(x) names(x)[length(x)-1])
    g3 <- sapply(paths, function(x) names(x)[-c(1,length(x))], simplify = FALSE)
    g3 <- sapply(g3, function(x) paste(x, collapse = "->"))
    g4 <- sapply(paths, function(x) names(x)[1])
    if(sum(g3 == "") > 0){
      g30 <- paste("NULL", seq(sum(g3 == "")), sep = "")
      g3[g3 == ""] <- g30
    }else{
      g30 <- NULL
    }
    g4.weight <- receptor.weight[g4]
    network1 <- network1[g4,,drop=FALSE]

    names(g3) <- paste("g3", g3, sep = ".")
    names(g4) <- paste("g4", g4, sep = ".")
    if(!is.null(g30)) names(g30) <- paste("g3", g30, sep = ".")
    names(pass1) <- paste("g3", pass1, sep = "." )
    names(pass1)[g3 %in% g30] <- paste("g2", pass1, sep = "." )[g3 %in% g30]
    names(pass2) <- paste("g3", pass2, sep = "." )
    names(pass2)[g3 %in% g30] <- paste("g4", pass2, sep = "." )[g3 %in% g30]
    
    network.df <- data.frame(from = g4, to = pass1, weight =  network1[cbind(g4, closestTF)], cluster = 1)
    network.df$evidence <- check.evidence(network.df)
    network.df$from <- names(g4)
    network.df$to <- names(pass1)
    network.df[!g3%in%g30,2] <-  names(g3)[!g3%in%g30]
    network.df$weight[network.df$weight < 0.05] <- 0.05
    
    graph <- rbind(igraph::as_data_frame(graph),
                   network.df
                  ) %>% 
             igraph::graph_from_data_frame()
    
    ################# g3 to g2 layer  #################
    network.df <- data.frame(from = pass2,
                             to = closestTF,
                             weight = network.df$weight,
                             cluster = 1)
    network.df$evidence <- check.evidence(network.df)
    network.df$to <- paste("g2", closestTF, sep = ".")
    network.df$from <- names(g3)
    network.df <- network.df[!g3%in%g30,]
    network.df$weight[network.df$weight < 0.05] <- 0.05
    graph <- rbind(igraph::as_data_frame(graph),
                   network.df
                  ) %>%
             igraph::graph_from_data_frame()
    
    ################# Clean up  #################
    rm(network.df)
    g3 <- g3[!duplicated(g3)]
    g3.weight <- rep(0, length(g3))
    g4 <- g4[!duplicated(g4)]
    g4.weight <- g4.weight[g4]
    if(!is.null(g30)){
        graph <- graph %>% igraph::add_vertices(nv = length(g30), name = names(g30))
    } 
    graph <- igraph::permute.vertices(graph,
                                      permutation = 
                                      sapply(
                                             names(igraph::V(graph)),
                                             function(x){
                                               which(
                                                 c(names(g1), names(g2), names(g3), names(g4)) == x
                                               )
                                             }
                                            )
                                     )
  }
  
  ################# Create graph annotation  #################
  # Vertex cluster
  igraph::V(graph)[names(g1)]$type <- 1
  igraph::V(graph)[names(g2)]$type <- 2
  if(show.pathways){
    igraph::V(graph)[names(g3)]$type <- 3
    igraph::V(graph)[names(g4)]$type <- 4
    if(!is.null(g30)) igraph::V(graph)[names(g30)]$type <- 5
  }
  igraph::V(graph)$cluster <- 1
  igraph::V(graph)[names(g2)]$cluster <- clu.g2[names(g2)]+1
  
  # Edge evidence
  igraph::E(graph)$arrow <- "-"
  igraph::E(graph)[igraph::E(graph)$evidence == 3]$arrow <- ">"
  igraph::E(graph)[igraph::E(graph)$evidence == 4]$arrow <- ">"
  
  # Layout
  if(show.pathways){
    if(swap.layers){
      l <-  layout.concentric(graph, concentric = list(names(g4), names(g3), names(g2), names(g1)))
    }else{
      l <-  layout.concentric(graph, concentric = list(names(g1), names(g2), names(g3), names(g4)), radius = radius)
    }
  }else{
    if(swap.layers){
      l <-  layout.concentric(graph, concentric = list(names(g2), names(g1)))
    }else{
      l <-  layout.concentric(graph, concentric = list(names(g1), names(g2)))
    }
  }
  rownames(l) <- names(igraph::V(graph))
  colnames(l) <- c("x","y")
  
  ################# Create pie chart  #################
  # Pie for g1 and g2
  pie.values <- t(network) %>% as.data.frame()
  pie.values$cluster <- clu.g2
  pie.values <- plyr::ddply(pie.values, "cluster",colMeans) %>%
                t %>%
                as.data.frame()
  pie.values <- pie.values[-nrow(pie.values),,drop=FALSE]
  pie.values <- pie.values/rowSums(pie.values)*g1.weight
  pie.values <- rbind(pie.values,
                      matrix(0, nrow = length(g2), ncol = ncol(pie.values),
                             dimnames = list(names(g2), colnames(pie.values))))
  for(i in names(clu.g2)) pie.values[i,clu.g2[i]] <- g2.weight[i]
  pie.values <- cbind(data.frame(N=1-rowSums(pie.values)), pie.values)
  
  if(show.pathways){
    # Pie for g3
    pie.values1 <- matrix(0,
                          ncol = length(unique(clu.g2)),
                          nrow = length(g3)
                          ) %>% 
                          as.data.frame()
    rownames(pie.values1) <- names(g3)
    pie.values1 <- cbind(
                         data.frame(N=1-rowSums(pie.values1)),
                         pie.values1
                        )
    
    # Pie for g4
    pie.values2 <- t(network1) %>% as.data.frame()
    pie.values2$cluster <- clu.g2
    pie.values2 <- plyr::ddply(pie.values2, "cluster", colMeans) %>%
                   t %>%
                   as.data.frame()
    pie.values2 <- pie.values2[-nrow(pie.values2),,drop=FALSE]
    pie.values2 <- pie.values2/rowSums(pie.values2)*g4.weight
    pie.values2 <- cbind(data.frame(N=1-rowSums(pie.values2)), pie.values2)
    rownames(pie.values2) <- paste("g4",rownames(pie.values2),sep = ".")
    
    # Combine Pies
    pie.values <- rbind(pie.values, pie.values1, pie.values2)
    pie.values <- pie.values[names(igraph::V(graph)),,drop=FALSE]
    pie.values <- cbind(data.frame(cluster = igraph::V(graph)$cluster),pie.values)
    pie.values <- cbind(data.frame(scale = 1),pie.values)
    pie.values[names(g3),]$scale <- 0
    pie.values <- cbind(l,pie.values)
    
  }else{
    pie.values <- pie.values[names(igraph::V(graph)),,drop=FALSE]
    pie.values <- cbind(data.frame(cluster = igraph::V(graph)$cluster),pie.values)
    pie.values <- cbind(data.frame(scale = 1),pie.values)
    pie.values <- cbind(l,pie.values)
    
  }
  
  ################# Aes  #################
  col <- RColorBrewer::brewer.pal("Set1", n = 9)
  col1 <- c("white",col) %>% ggplot2::alpha(.,0.8)
  col2 <-  c("grey",col)
  label.col <- c("darkblue","darkblue","grey2","darkblue","darkblue")
  vertex.label <- names(igraph::V(graph))
  if(show.pathways){
    vertex.label <- c(g1,g2,g3,g4)[vertex.label]
    vertex.label[vertex.label %in% g30] <- ""
  }else{
    vertex.label <- c(g1,g2)[vertex.label]
  }
  
  lty <- c("blank","dashed","solid", "solid")
  arrow.angle <- c(0,0,15,90)
  segement_data <- make_segements(graph, l)
  
  p <- ggraph::ggraph(graph, layout = l)  +
    
    geom_segment(data = segement_data,aes(x=x,y=y,xend=xend,yend=yend),
                 arrow = grid::arrow(type = "closed",
                                     angle = arrow.angle[igraph::E(graph)$evidence],
                                     length = unit(0.15, "inches")),
                 color = col2[igraph::E(graph)$cluster],
                 linewidth = igraph::E(graph)$weight*2,
                 linetype = lty[igraph::E(graph)$evidence])+
    
    scatterpie::geom_scatterpie(
      data = pie.values,
      cols = colnames(pie.values)[-c(1,2,3,4)],
      aes(x = x,y = y, col= as.character(cluster), r = pie.radius * scale),
      alpha = 0.5,
    )+
    
    ggraph::geom_node_text(aes(label = vertex.label), col = label.col[igraph::V(graph)$type], size = text.size)+
    
    scale_fill_manual(values = col1)+
    scale_color_manual(values = col2)+
    theme_classic()+ theme(axis.line=element_blank(),axis.text.x=element_blank(),
                           axis.text.y=element_blank(),axis.ticks=element_blank(),
                           axis.title.x=element_blank(),
                           axis.title.y=element_blank(),legend.position="none",
                           panel.border=element_blank(),panel.grid.major=element_blank(),
                           panel.grid.minor=element_blank(),plot.background=element_blank())
  p
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper functions, some are acquired from rTRM package
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_segements <- function(g,
                           l,
                           shorten.start = 0.07,
                           shorten.end = 0.07){
    es <- igraph::get.edgelist(g)
    data <- data.frame(matrix(0 , nrow = length(igraph::E(g)), ncol = 4))
    colnames(data) <- c("x","y","xend", "yend")
    for(i in 1:nrow(data)){
        data[i,1:2] <- l[es[i,1],]
        data[i,3:4] <- l[es[i,2],]
        
    }
    data$dx = data$xend - data$x
    data$dy = data$yend - data$y
    data$dist = sqrt( data$dx^2 + data$dy^2 )
    data$px = data$dx/data$dist
    data$py = data$dy/data$dist

    data$x = data$x + data$px * shorten.start
    data$y = data$y + data$py * shorten.start
    data$xend = data$xend - data$px * shorten.end
    data$yend = data$yend - data$py * shorten.end
    data[,1:4]
}

layout.concentric = function (g,
                              concentric = NULL,
                              radius = NULL,
                              order.by){
  if(is.null(concentric))
    concentric = list(igraph::V(g)$name)
  
  all_c = unlist(concentric, use.names = FALSE)
  
  if (!.checkValid(all_c))
    stop("Duplicated nodes in layers!")
  
  if (!.checkValid(radius))
    stop("Duplicated radius in layers!")
  
  all_n = igraph::V(g)$name
  sel_other = all_n[ ! all_n %in% all_c ]
  
  if(length(sel_other) > 0)
    concentric[[length(concentric)+1]] = sel_other
  
  if(is.null(radius)) {
    radius = seq(0, 1, 1/(length(concentric)))
    if(length(concentric[[1]]) == 1)
      radius = radius[-length(radius)]
    else
      radius = radius[-1]
  }
  
  if( ! missing(order.by) )
    order.values = lapply(order.by, function(b) igraph::get.vertex.attribute(g, b))
  
  res = matrix(NA, nrow = length(all_n), ncol = 2)
  for(k in 1:length(concentric)) {
    r = radius[k]
    l = concentric[[k]]
    
    i = which(igraph::V(g)$name %in% l) - 1
    i_o = i
    if (!missing(order.by)) {
      ob = lapply(order.values, function(v) v[i + 1])
      ord = do.call(order, ob)
      i_o = i_o[ord]
    }
    res[i_o+1, ] = .getCoordinates(i_o, r)
    
  }
  res
}

.getCoordinates = function(x, r) {
  l = length(x)
  d = 360/l
  c1 = seq(0, 360, d)
  c1 = c1[1:(length(c1)-1)]
  tmp = t(sapply(c1, function(cc) c(cos(cc*pi/180)*r, sin(cc*pi/180)*r)))
  rownames(tmp) = x
  tmp
}

.checkValid = function(x) {
  if(any(table(x) > 1)) FALSE else TRUE
}
