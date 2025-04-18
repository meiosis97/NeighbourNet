require(Seurat)
require(dplyr)
require(Matrix)
require(ggplot2)
require(ggraph)
require(scatterpie)
require(ggrepel)
require(ggpubr)

gg.color.spec <- function(n = 11, background = F) {
  if(background){
    c('#D9DDDC',rev(RColorBrewer::brewer.pal(n = n,name = "Spectral")))
  }else{
    rev(RColorBrewer::brewer.pal(n = n,name = "Spectral"))
  }
}
get.ppr <- function(p = NULL){
  rt.ppr <- receptor.ppr$ppr
  if(is.null(p)) p <- receptor.ppr$ltf
  cutoff <- apply(rt.ppr, 1, quantile, p)
  rt.ppr <- rt.ppr * (rt.ppr > cutoff)
  Matrix::Matrix(rt.ppr)
}
get.gr.adj <- function(t = 2){
  tmp <- gr.graph
  igraph::edge.attributes(tmp)$weight <- igraph::edge.attributes(tmp)$weight *
    (igraph::edge.attributes(tmp)$consensus_stimulation - igraph::edge.attributes(tmp)$consensus_inhibition)
  adj <- igraph::as_adjacency_matrix(tmp, attr = "weight")
  if(t > 1){
    for(i in 1:(t-1)) adj <- adj %*% adj
  }
  adj
}

build.graph <- function(pcs, knn = 30){

  FindSigma <- function(dk, a, knn) {
    lower <- 0
    upper <- Inf
    cur <- dk[knn]
    while (T) {
      psum <- sum(exp(-dk / cur))
      if (psum > a) {
        upper <- cur
        cur <- (lower + cur) / 2
      } else if (psum < a) {
        if (is.infinite(upper)) {
          lower <- cur
          cur <- 2 * cur
        } else {
          lower <- cur
          cur <- (upper + cur) / 2
        }
      }
      if (abs(psum - a) < 1e-5) break
    }
    cur
  }

  n <- nrow(pcs)
  nn2.result <- RANN::nn2(pcs, k = knn)
  a <- 2*log2(knn)
  sigma <- c()

  scaled.dists <- nn2.result$nn.dists
  scaled.dists <- scaled.dists - scaled.dists[,2]
  scaled.dists[,1] <- 0

  for (i in 1:n) {
    dk <- scaled.dists[i, 1:knn]
    sigma[i] <- FindSigma(dk, a, knn)
  }

  nn.aff <- exp(-scaled.dists/sigma)
  nn.idx <- nn2.result$nn.idx
  nn.w <- nn.aff/rowSums(nn.aff)

  i <- rep(1:n, times = knn)
  j <- as.numeric(nn.idx)
  x <- as.numeric(nn.aff)
  aff <- Matrix::sparseMatrix(i = i,j = j,x = x)

  aff <- aff + t(aff) - aff * t(aff)

  d <- rowSums(aff)
  p <- aff
  p <- 0.5*(sweep(p, 2, sqrt(d), "/") / sqrt(d))
  diag(p) <-  diag(p) + 0.5

  return(list(p = p, nn.idx = nn.idx, nn.w = nn.w))
}
select.gene <- function(seurat.obj,
                        tfs = NULL,
                        targets = NULL,
                        bg = NULL, min.cells = 20){
  if(is.null(tfs)){
    tfs <- unique(gene.list$tfs)
  }else{
    tfs <- intersect(gene.list$tfs, tfs)
  }

  if(is.null(targets)){
    targets <- gene.list$targets
  }else{
    targets <- intersect(gene.list$targets, targets)
  }

  genes <- rownames(seurat.obj)
  n.cell.expressed <- rowSums(seurat.obj@assays$RNA$counts[genes,]>0)
  genes <- genes[n.cell.expressed > min.cells]

  tfs <- intersect(genes, tfs)
  targets <- intersect(genes, targets)
  bg <- intersect(genes, bg)
  bg <- bg[!bg %in% c(tfs, targets)]
  genes <- unique(c(tfs, targets, bg))

  list(tfs=tfs, targets = targets, bg = bg, genes = genes)

}

prepare.seurat <- function(seurat.obj, genes = VariableFeatures(seurat.obj), k = 100){
  genes <- intersect(genes, rownames(seurat.obj))
  message("Run Seurat.")
  seurat.obj <- Seurat::ScaleData(seurat.obj, features = genes, verbose = F, scale.max = Inf, do.scale = T)
  seurat.obj <- Seurat::RunPCA(seurat.obj, npcs = k, features = genes, verbose = F)
  seurat.obj
}
prepare.graph <- function(seurat.obj,
                          knn = 30, truncated = TRUE){
  sd <- Seurat::Reductions(seurat.obj, "pca")@stdev
  k <- length(sd)

  # Select important pcs
  if(truncated){
    s <- abs(diff(sd))
    mu <- mean(s[(0.8*k):k-1])
    sigma <- sd(s[(0.8*k):k-1])
    sk <- mu + 6*sigma
    npcs <- max(which(s > sk))+1
    if(length(npcs)==0){
      warning(npcs, " components could be insufficient to approximate the data.")
      npcs <- k
    }
  }else{
    npcs <- k
  }

  pcs <- Seurat::Embeddings(seurat.obj, "pca")[,1:npcs]
  loadings <- Seurat::Reductions(seurat.obj, "pca")@feature.loadings[,1:npcs]
  sd <- sd[1:npcs]

  # Build graph
  message("Now building knn graph.")

  graph.result <- build.graph(pcs, k = knn)

  setting <- list(pcs = pcs, loadings = loadings,
                  predictors = NULL, responses = NULL, cells = NULL,
                  p = graph.result$p, nn.idx = graph.result$nn.idx, nn.w = graph.result$nn.w,
                  lra = NULL, scale.gene = NULL,
                  nn.scale.gene = NULL, nn.scale.pc = NULL, n.eff = NULL)

  suppressWarnings(
    Seurat::Misc(seurat.obj, "setting") <- setting
  )
  seurat.obj
}
select.cell <- function(seurat.obj, p = 0.1, n = NULL, all = FALSE){
  setting <- Seurat::Misc(seurat.obj, "setting")

  if(is.null(setting)) stop("Run prepare.graph first.")

  if(all){
    setting$cells <- NULL
  }else{
    if(is.null(n)){
      n <- ceiling(nrow(setting$pcs) * p)
    }
    centers <- stats::kmeans(setting$pcs, n, iter.max = 100)$centers
    ids <- RANN::nn2(setting$pcs, centers, k = 2)$nn.idx %>% as.numeric()
    ids <- unique(ids)
    names(ids) <- rownames(setting$pcs)[ids]
    setting$cells <- ids
  }

  suppressWarnings(
    Seurat::Misc(seurat.obj, "setting") <- setting
  )

  seurat.obj
}
prepare.reg <- function(seurat.obj, responses = NULL, predictors = NULL, cells = NULL, check.expressed = TRUE){
  setting <- Seurat::Misc(seurat.obj, "setting")

  if(is.null(setting)) stop("Run prepare.graph first.")

  if(is.null(responses)){
    responses <- rownames(setting$loadings)
  }else{
    responses <- intersect(responses, rownames(setting$loadings))
  }

  if(is.null(predictors)){
    predictors <- rownames(setting$loadings)
  }else{
    predictors <- intersect(predictors, rownames(setting$loadings))
  }
  genes <- c(responses, predictors) %>% unique()
  if(is.null(cells)){
    cells <- setting$cells
    if(is.null(cells)){
      cells <- 1:nrow(setting$pcs)
      names(cells) <- rownames(setting$pcs)
    }
  }else{
    names(cells) <- rownames(setting$pcs)[cells]
    setting$cells <- cells
  }
  pcs <- setting$pcs[cells,]
  n.gene <- length(genes)
  n.cell <- length(cells)
  n.pc <- ncol(pcs)
  setting$responses <- responses
  setting$predictors <- predictors
  setting$genes <- genes

  # Low rank approximation.
  setting$lra <- tcrossprod(setting$pcs, setting$loadings[responses,])

  # Local variances.
  message("Calculating local variance.")

  nn.scale.gene <-  matrix(0, nrow = n.gene, ncol = n.cell,
                           dimnames = list(genes, rownames(pcs)))
  nn.scale.pc <- matrix(0, nrow = n.pc, ncol = n.cell,
                        dimnames = dimnames(pcs) %>% rev)
  n.eff <- apply(setting$nn.w[cells,], 1, function(x) (sum(x)^2)/sum(x^2))
  names(n.eff) <- rownames(pcs)

  scale.data <- SeuratObject::LayerData(seurat.obj, layer = "scale.data")[genes,rownames(setting$pcs)]
  setting$scale.gene <- SeuratObject::LayerData(seurat.obj, layer = "data")[genes,rownames(setting$pcs)] %>%
    apply(., 1, sd)
  if(check.expressed) expression <- SeuratObject::LayerData(seurat.obj, layer = "counts")[genes,rownames(setting$pcs)] != 0

  for(i in 1:n.cell){
    idx <- setting$nn.idx[cells[i],]
    w <- setting$nn.w[cells[i],]

    # Local gene scales
    w.mean <- as.numeric(scale.data[,idx] %*% w)
    res <- (scale.data[,idx] - w.mean)
    nn.scale.gene[,i] <- as.numeric(res^2 %*% w)*n.eff[i]/(n.eff[i]-1)
    if(check.expressed) nn.scale.gene[,i] <- nn.scale.gene[,i] * as.numeric(expression[,idx] %*% w)

    # Local pc scales
    w.mean <- as.numeric(w %*% setting$pcs[idx,])
    res <- t(setting$pcs[idx,]) - w.mean
    nn.scale.pc[,i] <- as.numeric(res^2 %*% w)*n.eff[i]/(n.eff[i]-1)
  }

  setting$nn.scale.gene <- Matrix::Matrix(sqrt(nn.scale.gene))
  setting$nn.scale.pc <- sqrt(nn.scale.pc)
  setting$n.eff <- n.eff

  suppressWarnings(
    Seurat::Misc(seurat.obj, "setting") <- setting
  )

  seurat.obj
}
run.nn.reg <- function(seurat.obj, responses = NULL, Y = NULL,
                       predictors = NULL, t = 3, k = NULL,
                       remove.self.loops = T, f = function(x) 2*x^2, assay = c("effect", "p.val"),
                       prune = TRUE, cutoff = 0.5,
                       return.p.val = FALSE, return.smooth = TRUE, return.prune = FALSE){
  setting <- Seurat::Misc(seurat.obj, "setting")

  if(is.null(setting)) stop("Run prepare.graph first, and then prepare.reg.")

  assay <- match.arg(assay)

  # Extract X
  X <- setting$pcs

  # Extract Cells
  subsampled <- !is.null(setting$cells)
  if(!subsampled){
    cells <- 1:nrow(X)
    names(cells) <- rownames(X)
  }else{
    cells <- setting$cells
  }

  # Extract Y
  custom.y <- !is.null(Y)
  if(!custom.y){
    # Check responses.
    responses <- intersect(setting$responses, responses)
    responses <- responses[rowSums(setting$nn.scale.gene[responses,names(cells),drop=F]) > 0]
    Y <- setting$lra[,responses, drop = F]
  }else{
    Y <- as.matrix(Y)
    responses <- colnames(Y)
    if(is.null(responses)) responses <- colnames(Y) <- paste("Y", 1:ncol(Y), sep = "")
  }

  # Select regression coefficients of genes to be kept.
  if(is.null(predictors)){
    predictors <- setting$predictors
  }else{
    predictors <- intersect(predictors, setting$genes)
  }

  if(!custom.y){
    genes <- unique(c(responses, predictors))
  }else{
    genes <- predictors
  }

  # Extract loadings
  loadings <- setting$loadings[genes,, drop = F]
  loading.scale <- rowSums(loadings^2) %>% sqrt
  loadings <- loadings/loading.scale

  n.cell <- length(cells)
  n.gene <- length(genes)
  n.response <- length(responses)
  n.predictor <- length(predictors)
  n.pc <- ncol(X)
  if(is.null(k)) k <- ncol(setting$nn.idx)

  if(!return.smooth) return.prune <- FALSE
  if(return.prune) prune <- FALSE

  # Generate message
  if(return.smooth){
    message("Return smoothed effect, can only generate networks for sampled cells.")
  }else{
    message("Return raw effect.")
  }
  if(return.prune){
    message("Return pruned effect.")
  }else{
    message("Return unpruned effect.")
  }
  if(return.p.val){
    message("Return p-value.")
  }else{
    message("Will not return p-value.")
  }
  if(assay == "effect"){
    message("By default, downstream analysis will be performed on the effect tensor.")
  }else{
    message("By default, downstream analysis will be performed on the p-val tensor.")
  }
  if(prune){
    message("By default, downstream analysis will perform network prunning.")
  }else{
    message("By default, downstream analysis will not perform network prunning")
  }

  # Get TF information
  tfs.in.responses <- responses[responses%in%gene.list$tfs]
  tfs.in.predictors <- predictors[predictors%in%gene.list$tfs]
  targets.in.responses <- responses[responses%in%gene.list$targets]
  targets.in.predictors <- predictors[predictors%in%gene.list$targets]

  # Local variances of y
  nn.scale.y <- matrix(0, nrow = n.response,
                       ncol = n.cell, dimnames = dimnames(Y[cells,]) %>% rev)

  for(i in 1:n.cell){
    j <- cells[i]
    idx <- setting$nn.idx[j,]
    w <- setting$nn.w[j,]

    # Local response scales
    w.mean <- as.numeric(w %*% Y[idx,])
    res <- t(Y[idx,]) - w.mean
    nn.scale.y[,i] <- as.numeric(res^2 %*% w)*setting$n.eff[i]/(setting$n.eff[i]-1)
  }
  nn.scale.y <- sqrt(nn.scale.y)

  # Build the Laplacian operator.
  message("Build the Laplacian operator.")
  svds.p <-RSpectra::svds(A = setting$p, k = k)
  u <- svds.p$u
  vd <- sweep(svds.p$v, 2, svds.p$d^t, "*")
  rownames(u) <- rownames(vd) <- rownames(X)
  d <- tcrossprod(u, vd[cells,]) %>% rowSums()
  u <- u/d

  effect.tensor <-
    array(dim = c(n.response, n.gene, n.cell),
          dimnames = list(responses, genes, rownames(X)[cells]))
  p.val.tensor <- if(return.p.val) effect.tensor else NULL

  mus <- c()
  sigmas <- c()

  ############ Regression starts here ############
  message("Now regress.")
  for(j in 1:n.response){
    r <- responses[j]
    message(r)

    b <- matrix(0, nrow = n.cell, ncol = n.pc)
    for(k in 1:n.cell){
      i <- cells[k]
      idx <- setting$nn.idx[i,]
      w <- setting$nn.w[i,]

      # Local data scaling
      n <- length(idx)
      w.mean <- as.numeric(w %*% Y[idx,j])
      y <- Y[idx,j] %>% scale(center = w.mean, scale = nn.scale.y[j,k])
      w.mean <- as.numeric(w %*% X[idx,])
      x <- X[idx,] %>% scale(center = w.mean, scale = setting$nn.scale.pc[,k])
      y.scale <- attr(y,"scaled:scale")
      x.scale <- attr(x,"scaled:scale")

      # Local regression
      if(y.scale){
        lambda <- 5
        qr.mod <- rbind(x *sqrt(w),diag(sqrt(lambda), ncol(x))) %>% qr
        v <- qr.qty(qr.mod, c(y *sqrt(w), rep(0, ncol(x))))
        b[k,] <- backsolve(qr.R(qr.mod), v)/x.scale*y.scale
      }else{
        b[k,] <- rep(0, ncol(x))
      }

    }

    # Dot product
    b <- tcrossprod(loadings, b) %>% as.matrix

    # Account for response variances inflated by LRA.
    if(!custom.y){
      y.factor <- setting$nn.scale.gene[r,names(cells)]/nn.scale.y[j,]*setting$scale.gene[r]
    }else{
      y.factor <- 1
    }

    # Calculate effect
    effect <- b * setting$nn.scale.gene[genes,names(cells),drop=F] %>%
      sweep(2, y.factor , "*") %>% as.matrix
    if(return.smooth|return.p.val) effect.hat <-  tcrossprod(effect %*% vd[cells,], u[cells,])

    # Calculate p-value
    if(!custom.y){
      noise <- effect[j,]
      noise <- replicate(100, sample(noise))
    }else{
      ref <- apply(abs(effect), 1, max) %>% which.max
      noise <- effect[ref,]
      noise <- replicate(100, sample(noise))
    }
    noise <- u[cells,] %*% crossprod(vd[cells,], noise)
    noise <- log(abs(noise))
    mu <- mean(noise)
    sigma <- sd(noise)
    mus[j] <- mu
    sigmas[j]  <- sigma
    p.val <- if(return.p.val) pnorm(log(abs(effect.hat)), mu, sigma) else NULL

    if(return.smooth){
      if(return.prune){
        if(is.null(p.val)) p.val <- pnorm(log(abs(effect.hat)), mu, sigma)
        effect.hat[p.val < cutoff] <- 0
      }
      effect.tensor[j,,] <- effect.hat
    }else{
      effect.tensor[j,,] <- effect
    }

    if(return.p.val) p.val.tensor[j,,] <- p.val

  }
  names(mus) <- responses
  names(sigmas) <- responses

  # Default controls for extracting sub networks.
  defaults <- list(f = f, remove.self.loops = remove.self.loops,
                   assay = assay, predictors = predictors, responses = responses,
                   cutoff = cutoff, prune = prune)
  if(!custom.y){
    gene.sets <- list(predictors = list(genes = predictors, tfs = tfs.in.predictors, targets = targets.in.predictors),
                      responses = list(genes = responses, tfs = tfs.in.responses, targets = targets.in.responses),
                      genes = genes)
  }else{
    gene.sets <- list(predictors = list(genes = predictors, tfs = tfs.in.predictors, targets = targets.in.predictors),
                      responses = list(genes = responses, tfs = NULL, targets = NULL),
                      genes = genes)
  }


  suppressWarnings(
    Seurat::Misc(seurat.obj, "mod") <- list(effect = effect.tensor, p.val = p.val.tensor,
                                            meta.network = NULL, mus = mus, sigmas = sigmas,
                                            subsampled = subsampled,
                                            smoothed = return.smooth, pruned = return.prune,
                                            gene.sets = gene.sets, cells = cells,
                                            defaults = defaults, custom.y = custom.y, w = list(u = u, vd = vd),  return.smooth = return.smooth)
  )

  seurat.obj
}
set.defaults <- function(seurat.obj, clean.up = FALSE, defaults = list()){
  mod <- Seurat::Misc(seurat.obj, "mod")

  if(is.null(mod)) stop("Run prepare.reg first.")

  if(clean.up){
    mod$defaults <- list(f = function(x) 2*x^2,
                         remove.self.loops = T,
                         assay = "effect",
                         predictors = mod$gene.sets$predictors$genes,
                         responses = mod$gene.sets$responses$genes,
                         cutoff = 0.5, prune = TRUE)
  }

  # Filter valid arguments
  valid.args <- intersect(names(defaults), names(mod$defaults))
  if(length(valid.args)) defaults <- defaults[valid.args]

  mod$defaults[valid.args] <- defaults[valid.args]

  # Filter valid responses and predictors
  mod$defaults$predictors <-
    intersect(mod$defaults$predictors, mod$gene.sets$genes)
  mod$defaults$responses <-
    intersect(mod$defaults$responses, mod$gene.sets$responses$genes)

  # Check assay input
  mod$defaults$assay <- match.arg(mod$defaults$assay, choices = c("effect", "p.val", "meta.network"))

  # Check remove.self.loops input
  if(!is.logical(mod$defaults$remove.self.loops)) stop("remove.self.loops must be logical.")

  # Check prune input
  if(!is.logical(mod$defaults$prune)) stop("prune must be logical.")

  # Check f input
  if(!is.function(mod$defaults$f)) stop("f must be a function.")

  # Check cutoff input
  cutoff <- mod$defaults$cutoff
  if(!is.numeric(cutoff)) stop("cutoff must be numerical and between 0 and 1.")

  suppressWarnings(
    Seurat::Misc(seurat.obj, "mod") <- mod
  )
  seurat.obj
}
get.network <- function(seurat.obj, i = NULL, assay = NULL, remove.self.loops = NULL,
                        responses = NULL, predictors = NULL, f = NULL, drop = T,
                        prune = NULL, cutoff = NULL){
  mod <- Seurat::Misc(seurat.obj, "mod")
  if(is.null(mod)) stop("Run prepare.reg first.")

  if(is.null(assay)){
    assay <- mod$defaults$assay
  }else{
    assay <- match.arg(assay, choices = c("effect", "p.val", "meta.network"))
  }

  if(is.null(cutoff)){
    cutoff <- mod$defaults$cutoff
  }else if(!is.numeric(cutoff)){
    stop("cutoff must be numerical and between 0 and 1.")
  }

  if(is.null(remove.self.loops)){
    remove.self.loops <- mod$defaults$remove.self.loops
  }else if(!is.logical(remove.self.loops)){
    stop("remove.self.loops must be logical.")
  }

  if(is.null(prune)){
    prune <- mod$defaults$prune
  }else if(!is.logical(prune)){
    stop("prune must be logical.")
  }

  if(prune & mod$pruned & assay != "meta.network"){
    if(assay == "effect"){
      warning("The effect ensemble is already pruned. ",
              "Use set.defaults to set pruned to FALSE to enable further pruning")
      prune <- FALSE
    }else if(!is.null(mod$p.val)){
      warning("The effect ensemble is already pruned. ",
              "Use set.defaults to set pruned to FALSE to enable further pruning")
      prune <- FALSE
    }
  }

  if(is.null(responses)){
    responses <- mod$defaults$responses
  }else{
    responses <- intersect(responses, mod$gene.sets$responses$genes)
  }

  if(is.null(predictors)){
    predictors <- mod$defaults$predictors
  }else{
    predictors <- intersect(predictors, mod$gene.sets$genes)
  }

  if(is.null(f)){
    f <- mod$defaults$f
  }else if(!is.function(f)){
    stop("f must be a function.")
  }

  w <- mod$w
  smoothed <- mod$smoothed
  if(assay != "meta.network"){
    if(is.null(i)){
      i <- mod$cells
    }else{
      names(i) <- rownames(w$u)[i]
    }
    need.impute <- !all(i %in% mod$cells)
    if(smoothed & need.impute) stop("Unable to impute cell networks due to the effect network is already smoothed.")
  }else if(is.null(i)){
    i <- 1:dim(mod$meta.network$meta.network)[[3]]
    names(i) <- dimnames(mod$meta.network$meta.network)[[3]]
  }else{
    names(i) <- dimnames(mod$meta.network$meta.network)[[3]][i]
  }
  n.cell <- length(i)
  n.response <- length(responses)
  n.predictor <- length(predictors)

  # Get the network
  if(assay=="effect"){

    if(smoothed){
      network <- mod$effect[responses,predictors,names(i),drop=F]
    }else{
      network <- array(0, dim = c(n.response, n.predictor, n.cell),
                       dimnames = list(responses, predictors, names(i)))
      for(r in responses) network[r,,] <- tcrossprod(mod$effect[r,predictors,] %*% w$vd[mod$cells,], w$u[i,,drop=F])
    }

  }else if(assay == "p.val"){

    if(!need.impute & !is.null(mod$p.val)){
      network <- mod$p.val[responses,predictors,names(i),drop=F]

    }else if(smoothed){
      network <- array(0, dim = c(n.response, n.predictor, n.cell),
                       dimnames = list(responses, predictors, names(i)))

      for(r in responses){
        network[r,,] <- pnorm(log(abs(mod$effect[r,predictors,names(i)])), mod$mus[r], mod$sigmas[r])
      }

    }else{
      network <- array(0, dim = c(n.response, n.predictor, n.cell),
                       dimnames = list(responses, predictors, names(i)))
      for(r in responses){
        network[r,,] <- tcrossprod(mod$effect[r,predictors,] %*% w$vd[mod$cells,], w$u[i,,drop=F])
        network[r,,] <- pnorm(log(abs(network[r,,])), mod$mus[r], mod$sigmas[r])
      }
    }

  }else{
    network <- mod$meta.network$meta.network[responses,predictors,i,drop=F]
  }

  # Prune the effect network
  if(prune){
    if(assay == "effect"){
      # Prune by each response
      for(r in responses){
        if(need.impute | is.null(mod$p.val)){
          p.val <- pnorm(log(abs(network[r,predictors,])), mod$mus[r], mod$sigmas[r])
        }else{
          p.val <- mod$p.val[r,predictors,names(i)]
        }
        network[r,,] <- network[r,,] * (p.val > cutoff)
      }
    }else if(assay == "p.val"){
      network[network < cutoff] <- 0
    }
  }

  # Transform
  if(assay == "effect") network <- f(network)

  # Remove self-loop
  if(remove.self.loops){
    for(r in intersect(responses,predictors)) network[r,r,] <- 0
  }

  if(drop){
    if(n.cell == 1){
      network <- matrix(network, nrow = n.response, dimnames = list(responses, predictors))
    }else if(n.response == 1){
      network <- matrix(network, nrow = n.predictor, dimnames = list(predictors, names(i)))
    }else if(n.predictor == 1){
      network <- matrix(network, nrow = n.response, dimnames = list(responses, names(i)))
    }
  }
  network
}
build.meta.network <- function(seurat.obj = NULL, network = NULL, k = 100,
                               big.memory = FALSE,
                               scale = TRUE,
                               truncated =TRUE, n.net = 20, non.neg = T,
                               max.iter = 1000, tol = 1e-10, return.p.val = TRUE){

  if(!is.null(seurat.obj)){
    mod <- Seurat::Misc(seurat.obj, "mod")
    if(is.null(mod)) stop("Run run.nn.reg first.")
    assay <- mod$defaults$assay
    assay <- match.arg(assay, choices = c("effect", "p.val"))
    cells <- names(mod$cells)
    predictors <- mod$defaults$predictors
    responses <- mod$defaults$responses
    if(big.memory) network <- get.network(seurat.obj, drop = F)
  }else if(!is.null(network)){
    assay <- NULL
    cells <- dimnames(network)[[3]]
    predictors <- dimnames(network)[[2]]
    responses <- dimnames(network)[[1]]
  }else{
    stop("Neither a Seurat object nor a network tensor is provided.")
  }

  n.cell <- length(cells)
  n.response <- length(responses)
  n.predictor <- length(predictors)

  # Scale the size of each single cell network.
  if(scale) scales <- rep(0, n.cell)


  # Number of eigen vectors to extract
  k <- min(k, n.predictor*n.response, n.cell)
  if(k < 100) truncated = FALSE

  # Now construct the covariance matrix.
  message("Now construct the covariance matrix.")

  if(min(n.predictor, n.response) != 1){ # If there are more than one response or predictor.
    K <- matrix(0, nrow = n.cell, ncol = n.cell)
    colnames(K) <- rownames(K) <- cells

    for(r in responses){
      message(r)
      A <- if(is.null(network)) get.network(seurat.obj, responses = r) else network[r,,]
      # Scale the network
      if(scale) scales <- scales + colSums(A^2)
      A <- Matrix::Matrix(A)
      K <- K + crossprod(A)
    }

    if(scale){
      scales <- sqrt(scales)
      scales[scales==0] <- 1
      K <- K / outer(scales, scales)
    }

    # Covariance eigen decomposition to get embedding on edges.
    message("Eigen decomposition.")
    eigs <- RSpectra::eigs(K, k = k)
    v <- eigs$vectors
    d <- sqrt(eigs$values)
    vd <- v %*% diag(d)
    rownames(v) <-  rownames(vd) <- cells
    colnames(v) <- colnames(vd) <- paste("component", 1:k, sep ='_')

  }else{
    A <- if(is.null(network)) get.network(seurat.obj, drop = F) else network

    if(n.predictor == 1){
      A <- matrix(A, ncol = n.cell, dimnames = dimnames(A)[c(1,3)])
    }
    if(n.response == 1){
      A <- matrix(A, ncol = n.cell, dimnames = dimnames(A)[c(2,3)])
    }

    if(scale){
      scales <- sqrt(colSums(A^2))
      scales[scales==0] <- 1
      A <- sweep(A, 2, scales, "/")
    }

    # SVD to get embedding on edges.
    message("Single value decomposition.")
    if(min(dim(A)) <= 3){
      svd.mod <- svd(A, nu = 0, nv = k)
    }else{
      svd.mod <- RSpectra::svds(A, k = k)
    }
    v <- svd.mod$v
    d <- svd.mod$d[1:k]
    vd <- v %*% diag(d)
    rownames(v) <-  rownames(vd) <- cells
    colnames(v) <- colnames(vd) <- paste("component", 1:k, sep ='_')
  }

  # Select the number of components.
  if(!truncated){
    rank <- k
  }else{
    s <- abs(diff(d))
    mu <- mean(s[(0.8*k):k-1])
    sigma <- sd(s[(0.8*k):k-1])
    sk <- mu + 6*sigma
    rank <- max(which(s > sk))+1
    if(length(rank)==0){
      warning(rank, " single values could be insufficient to approximate the data.")
      rank <- k
    }
    vd <- vd[,1:rank]
    v <- v[,1:rank]
  }

  # Select the number of components.
  message("Non-negative PCA.")
  n.net <- min(n.net, abs(rank-1))
  meta.network <-
    array(dim = c(n.response, n.predictor, n.net),
          dimnames = list(responses, predictors, paste("component",1:n.net,sep = "_")))
  if(is.null(assay) | !non.neg){
    return.p.val <- FALSE
  }
  p.val <- if(return.p.val) meta.network else NULL


  if(non.neg){
    x <- vd
    B <- diag(n.cell)
    npca.loadings <- matrix(0, ncol = n.net, nrow = n.cell)
    rownames(npca.loadings) <- rownames(x)
    colnames(npca.loadings) <- paste("component", 1:n.net, sep = "_")
    tot.var <- sum(x^2)/rank
    npca.sd <- c()

    for(j in 1:n.net){
      w <- runif(n.cell)
      w <- w/sqrt(sum(w^2))
      for(i in 1:max.iter){
        w.old <- w
        u <- crossprod(x,w)
        w <- x %*% u
        w <- w/sum(w)
        w[w<0] <- 0
        w <- w/sqrt(sum(w^2))
        if(sum(w-w.old)^2 < tol) break
      }
      w <- w/as.numeric(sqrt(crossprod(w, B) %*% w))

      q <- B %*% w
      x <- x - q %*% (t(q) %*% x)
      B <- B - B %*% q %*% t(q)
      npca.loadings[,j] <- w/sqrt(sum(w^2))
      var.diff <- tot.var - sum(x^2)/rank
      npca.sd[j] <- sqrt(var.diff)
      tot.var <- tot.var - var.diff
    }

    for(r in responses){
      A <- if(is.null(network)) get.network(seurat.obj, responses = r) else network[r,,]
      if(return.p.val){
        B <- if(assay == "effect") get.network(seurat.obj, responses = r, assay = "p.val") else A
        p.val[r,,] <- B %*% npca.loadings
      }
      if(scale) A <- sweep(A, 2, scales, "/")

      meta.network[r,,] <- A %*% npca.loadings
    }

    if(return.p.val){
      for(i in 1:n.net) p.val[,,i] <- p.val[,,i]/max(p.val[,,i])
    }

  }else{
    npca.loadings <- NULL
    npca.sd <- NULL
    for(r in responses){
      A <- if(is.null(network)) get.network(seurat.obj, responses = r) else network[r,,]
      if(scale) A <- sweep(A, 2, scales, "/")

      meta.network[r,,] <- A %*% v[,1:n.net]
    }

  }

  meta.network <- list(meta.network = meta.network, p.val = p.val,
                       pcs = vd, pca.loadings = v, pca.sd = d,
                       npca.loadings = npca.loadings, npca.sd = npca.sd,
                       scale = scale, non.neg = non.neg,
                       setting = mod$defaults)

  if(!is.null(seurat.obj)){
    mod$meta.network <- meta.network
    suppressWarnings(
      Seurat::Misc(seurat.obj, "mod") <- mod
    )
    return(seurat.obj)
  }else{
    return(meta.network)
  }

}
build.meta.response <- function(seurat.obj = NULL, network = NULL, k = 100,
                                big.memory = FALSE,
                                scale = TRUE,
                                truncated =TRUE, n.net = 20, non.neg = TRUE,
                                max.iter = 1000, tol = 1e-10, return.p.val = TRUE){

  if(!is.null(seurat.obj)){
    mod <- Seurat::Misc(seurat.obj, "mod")
    if(is.null(mod)) stop("Run run.nn.reg first.")
    assay <- mod$defaults$assay
    assay <- match.arg(assay, choices = c("effect", "p.val"))
    cells <- mod$cells
    predictors <- mod$defaults$predictors
    responses <- mod$defaults$responses
    if(big.memory) network <- get.network(seurat.obj, drop = F)
  }else if(!is.null(network)){
    assay <- NULL
    cells <- 1:dim(network)[[3]]
    names(cells) <- dimnames(network)[[3]]
    predictors <- dimnames(network)[[2]]
    responses <- dimnames(network)[[1]]
  }else{
    stop("Neither a Seurat object nor a network tensor is provided.")
  }

  n.cell <- length(cells)
  n.response <- length(responses)
  n.predictor <- length(predictors)

  # Scale the size of each single cell network.
  if(scale) scales <- rep(0, n.response)

  # Number of eigen vectors to extract
  k <- min(k, n.predictor*n.cell, n.response)
  if(k < 100) truncated = FALSE

  # Now construct the covariance matrix.
  message("Now construct the covariance matrix.")

  if(min(n.predictor, n.response) != 1){ # If there are more than one response or predictor.
    K <- matrix(0, nrow = n.response, ncol = n.response)
    colnames(K) <- rownames(K) <- responses

    for(i in 1:n.cell){
      message(i)
      A <- if(is.null(network)) get.network(seurat.obj, i = cells[i]) else network[,,i]
      # Scale the network
      if(scale) scales <- scales + rowSums(A^2)
      A <- Matrix::Matrix(A)
      K <- K + tcrossprod(A)
    }
    if(scale){
      scales <- sqrt(scales)
      scales[scales==0] <- 1
      K <- K / outer(scales, scales)
    }

    # Covariance eigen decomposition to get embedding on edges.
    message("Eigen decomposition.")
    eigs <- RSpectra::eigs(K, k = k)
    v <- eigs$vectors
    d <- sqrt(eigs$values)
    vd <- v %*% diag(d)
    rownames(v) <-  rownames(vd) <- responses
    colnames(v) <- colnames(vd) <- paste("component", 1:k, sep ='_')

  }else{
    A <- if(is.null(network)) get.network(seurat.obj, drop = F) else network

    if(n.predictor == 1){
      A <- matrix(A, ncol = n.response, dimnames = dimnames(A)[c(3,1)])
    }
    if(n.cell == 1){
      A <- matrix(A, ncol = n.response, dimnames = dimnames(A)[c(2,1)])
    }

    if(scale){
      scales <- sqrt(colSums(A^2))
      scales[scales==0] <- 1
      A <- sweep(A, 2, scales, "/")
    }

    # SVD to get embedding on edges.
    message("Single value decomposition.")
    if(min(dim(A)) <= 3){
      svd.mod <- svd(A, nu = 0, nv = k)
    }else{
      svd.mod <- RSpectra::svds(A, k = k)
    }
    v <- svd.mod$v
    d <- svd.mod$d[1:k]
    vd <- v %*% diag(d)
    rownames(v) <-  rownames(vd) <- responses
    colnames(v) <- colnames(vd) <- paste("component", 1:k, sep ='_')
  }

  # Select the number of components.
  if(!truncated){
    rank <- k
  }else{
    s <- abs(diff(d))
    mu <- mean(s[(0.8*k):k-1])
    sigma <- sd(s[(0.8*k):k-1])
    sk <- mu + 6*sigma
    rank <- max(which(s > sk))+1
    if(length(rank)==0){
      warning(rank, " single values could be insufficient to approximate the data.")
      rank <- k
    }
    vd <- vd[,1:rank]
    v <- v[,1:rank]
  }

  # Select the number of components.
  message("Non-negative PCA.")
  n.net <- min(n.net, abs(rank-1))
  meta.network <-
    array(dim = c(n.net, n.predictor, n.cell),
          dimnames = list(paste("component",1:n.net,sep = "_"), predictors, names(cells)))
  if(is.null(assay) | !non.neg){
    return.p.val <- FALSE
  }
  p.val <- if(return.p.val) meta.network else NULL

  if(non.neg){
    x <- vd
    B <- diag(n.response)
    npca.loadings <- matrix(0, ncol = n.net, nrow = n.response)
    rownames(npca.loadings) <- rownames(x)
    colnames(npca.loadings) <- paste("component", 1:n.net, sep = "_")
    tot.var <- sum(x^2)/rank
    npca.sd <- c()

    for(j in 1:n.net){
      w <- runif(n.response)
      w <- w/sqrt(sum(w^2))
      for(i in 1:max.iter){
        w.old <- w
        u <- crossprod(x,w)
        w <- x %*% u
        w <- w/sum(w)
        w[w<0] <- 0
        w <- w/sqrt(sum(w^2))
        if(sum(w-w.old)^2 < tol) break
      }
      w <- w/as.numeric(sqrt(crossprod(w, B) %*% w))

      q <- B %*% w
      x <- x - q %*% (t(q) %*% x)
      B <- B - B %*% q %*% t(q)
      npca.loadings[,j] <- w/sqrt(sum(w^2))
      var.diff <- tot.var - sum(x^2)/rank
      npca.sd[j] <- sqrt(var.diff)
      tot.var <- tot.var - var.diff
    }

    for(i in 1:n.cell){
      A <- if(is.null(network)) get.network(seurat.obj, i = cells[i]) else network[,,i]
      if(return.p.val){
        B <- if(assay == "effect") get.network(seurat.obj, i = cells[i], assay = "p.val") else A
        p.val[,,i] <- crossprod(npca.loadings, B)
      }
      if(scale) A <- A/scales
      meta.network[,,i] <- crossprod(npca.loadings, A)
    }

  }else{
    npca.loadings <- NULL
    npca.sd <- NULL
    for(i in 1:n.cell){
      A <- if(is.null(network)) get.network(seurat.obj, i = cells[i]) else network[,,i]
      if(scale) A <- A/scales
      meta.network[,,i] <- crossprod(v[,1:n.net], A)
    }

  }

  meta.network <- list(meta.network = meta.network, p.val = p.val,
                       pcs = vd, pca.loadings = v, pca.sd = d,
                       npca.loadings = npca.loadings, npca.sd = npca.sd,
                       scale = scale, non.neg = non.neg,
                       setting = mod$defaults)

  return(meta.network)
}

select.central.genes <- function(seurat.obj = NULL, network = NULL, n.net = NULL,
                                 k = 1, n.per.component = 4, keep.responses = FALSE){
  if(!is.null(seurat.obj)){
    network <- get.network(seurat.obj, assay =  "meta.network", drop = F)
  }

  # Number of meta-networks to consider
  if(is.null(n.net)){
    n.net <- dim(network)[3]
  }else{
    n.net <- min(dim(network)[3], n.net)
  }
  network <- network[,,1:n.net,drop=F]

  responses <- dimnames(network)[[1]]
  predictors <- dimnames(network)[[2]]
  n.response <- length(responses)
  n.predictor <- length(predictors)

  central.predictors <- c()
  central.responses <- c()

  # Number of components to extract from each meta-network.
  k <- min(k, n.predictor, n.response)
  predictor.module <- rep(paste("M", 1:n.net, sep = ""), each = (k*n.per.component))
  response.module <- rep(paste("M", 1:n.net, sep = ""), each = (k*n.per.component))

  which.max.n <- function(x, n = n.per.component){
    which(rank(-x, ties.method = "random" ) <= n)
  }

  if(n.response > n.predictor & keep.responses){
    central.predictors <- predictors
    central.responses <- responses
    predictor.module <- NULL
    response.module <- NULL

  }else if(!keep.responses){
    if(min(n.predictor, n.response) != 1){ # If there are more than one response or predictor.
      for(i in 1:n.net){
        if(k <= 3){
          svd.mod <- svd(network[,,i], nu = k, nv = k) %>% suppressWarnings()
        }else{
          svd.mod <- RSpectra::svds(networks[,,i], k = k) %>% suppressWarnings()
        }
        v <- abs(svd.mod$v)
        u <- abs(svd.mod$u)
        rownames(v) <- predictors
        rownames(u) <- responses
        v[central.predictors,] <- 0
        u[central.responses,] <- 0

        # Extract genes with the highest scores on eigenvectors.
        central.predictors <- c(central.predictors,
                                predictors[apply(v, 2, which.max.n) %>% as.numeric])
        central.responses <- c(central.responses,
                               responses[apply(u, 2, which.max.n) %>% as.numeric])
      }
    }else if(n.predictor == 1){
      network <- matrix(abs(network), ncol = n.net, dimnames = dimnames(network)[c(1,3)])
      network[central.responses,] <- 0
      for(i in 1:n.net){
        central.responses <- c(central.responses, responses[which.max.n(network[,i]) %>% as.numeric])
      }
      central.responses <- responses[apply(network, 2, which.max.n) %>% as.numeric]
      central.predictors <- predictors
      predictor.module <- NULL

    }else if(n.response == 1){
      network <- matrix(abs(network), ncol = n.net, dimnames = dimnames(network)[c(2,3)])
      network[central.predictors,] <- 0
      for(i in 1:n.net){
        central.predictors <- c(central.predictors, predictors[which.max.n(network[,i]) %>% as.numeric])
      }
      central.responses <- responses
      response.module <- NULL
    }

  }else{
    central.responses <- responses
    response.module <- NULL
    if(min(n.predictor, n.response) != 1){ # If there are more than one response or predictor.
      for(i in 1:n.net){
        if(k <= 3){
          svd.mod <- svd(network[,,i], nu = k, nv = k) %>% suppressWarnings()
        }else{
          svd.mod <- RSpectra::svds(network[,,i], k = k) %>% suppressWarnings()
        }
        v <- abs(svd.mod$v)
        rownames(v) <- predictors
        v[central.predictors,] <- 0

        # Extract genes with the highest scores on eigenvectors.
        central.predictors <- c(central.predictors,
                                predictors[apply(v, 2, which.max.n) %>% as.numeric])
      }
    }else if(n.predictor == 1){
      central.predictors <- predictors
      predictor.module <- NULL

    }else if(n.response == 1){
      network <- matrix(abs(network), ncol = n.net, dimnames = dimnames(network)[c(2,3)])
      network[central.predictors,] <- 0
      for(i in 1:n.net){
        central.predictors <- c(central.predictors, predictors[which.max.n(network[,i]) %>% as.numeric])
      }
    }
  }

  if(!is.null(predictor.module)){
    predictor.module <- predictor.module[!duplicated(central.predictors)]
    predictor.module <- factor(predictor.module, levels = paste("M", 1:n.net, sep = ""))
  }
  if(!is.null(response.module)){
    response.module <- response.module[!duplicated(central.responses)]
    response.module <- factor(response.module, levels = paste("M", 1:n.net, sep = ""))
  }
  central.predictors <- central.predictors[!duplicated(central.predictors)]
  central.responses <- central.responses[!duplicated(central.responses)]

  list(central.responses = central.responses, central.predictors = central.predictors,
       response.module = response.module, predictor.module = predictor.module)
}


receptor.activity <- function(seurat.obj, i = NULL, meta.network = FALSE, t = 2,
                              prune = NULL,
                              receptors = NULL,
                              cutoff = NULL, check.receptor.expression = TRUE,
                              scale.receptor.activity = F,
                              scale.target.activity = F,
                              tfs = NULL, as.tfs = c("predictors", "responses"),
                              targets = NULL, receptor.activity = c("cprod", "dist")){
  setting <- Seurat::Misc(seurat.obj, "setting")
  if(is.null(setting)) stop("Run prepare.graph first, and then prepare.reg.")
  mod <- Seurat::Misc(seurat.obj, "mod")
  if(is.null(mod)) stop("Run run.nn.reg first.")

  as.tfs <- match.arg(as.tfs)
  assay <- mod$defaults$assay
  if(meta.network) assay <- "meta.network"
  f <- mod$defaults$f
  receptor.activity <- match.arg(receptor.activity)

  if(is.null(i)){
    if(assay != "meta.network"){
      cells <- mod$cells
    }else{
      cells <- 1:dim(mod$meta.network)[[3]]
      names(cells) <- paste("component", cells, sep = "_")
    }
  }else if(length(i)>1){
    if(assay != "meta.network"){
      cells <- i
      names(cells) <- rownames(setting$pcs)[i]
    }else{
      cells <- i
      names(cells) <- paste("component", cells, sep = "_")
    }
  }else{
    cells <- i
  }

  if(as.tfs == "predictors"){
    if(is.null(tfs)) tfs <- mod$gene.sets$predictors$tfs
    tfs <- intersect(tfs, mod$defaults$predictors)
    if(is.null(targets)) targets <- mod$defaults$responses
    targets <- intersect(targets, mod$defaults$responses)
  }else{
    if(is.null(tfs)) tfs <- mod$gene.sets$responses$tfs
    tfs <- intersect(tfs, mod$defaults$responses)
    if(is.null(targets)) targets <- mod$defaults$predictors
    targets <- intersect(targets, mod$defaults$predictors)
  }

  targets <- intersect(targets, names(igraph::V(gr.graph)))
  tfs <- intersect(tfs, names(igraph::V(gr.graph)))

  # Create a receptor to g2 personalized page rank score matrix.
  rtfs.ppr <- Matrix::Matrix(0, nrow = nrow(rt.ppr),
                             ncol = length(tfs), dimnames = list(rownames(rt.ppr), tfs))
  rtfs.ppr[,intersect(tfs, colnames(rt.ppr))] <-
    rt.ppr[,intersect(tfs, colnames(rt.ppr))]

  if(is.null(receptors)){
    receptors <- rownames(rtfs.ppr)
  }else{
    receptors <- intersect(rownames(rtfs.ppr), receptors)
    rtfs.ppr <- rtfs.ppr[receptors,,drop=F]
  }

  if(scale.receptor.activity){
    # Check receptor activity.
    receptor.scale <- rowSums(rtfs.ppr)
    tf.scale <- colSums(rtfs.ppr)
    rtfs.ppr <- rtfs.ppr/receptor.scale^0.5
    rtfs.ppr <- sweep(rtfs.ppr, 2, tf.scale^0.5, "/")
    rtfs.ppr[is.na(rtfs.ppr)] <- 0
  }

  if(check.receptor.expression){
    # Check receptor expression
    rw <- Matrix::Matrix(0, nrow = nrow(rtfs.ppr),
                         ncol = nrow(setting$pcs), dimnames = list(receptors, rownames(setting$pcs)))
    rw[intersect(rownames(seurat.obj), receptors), ] <-
      SeuratObject::LayerData(seurat.obj, "counts")[intersect(rownames(seurat.obj), receptors),rownames(setting$pcs)] != 0

    # Network propagation
    if(assay != "meta.network"){
      w <- tcrossprod(mod$w$u[cells,,drop=F], mod$w$vd)
      rw <- tcrossprod(rw, w)
    }else{
      w <- mod$w$vd %*%
        crossprod(mod$w$u[mod$cells,], mod$meta.network$npca.loadings[,cells])
      w <- sweep(w, 2, colSums(w), "/")
      rw <- rw %*% w
    }
    rw[rw > 1] <- 1
    rw[rw < 0] <- 0
    if(length(cells) == 1) rtfs.ppr <- rtfs.ppr * as.numeric(rw)
  }

  # Get the prior knowledge network
  adj <- get.gr.adj(t = t) %>% t

  # Prepare output
  if(length(cells) > 1){
    receptor.act <- data.frame(matrix(nrow = nrow(rtfs.ppr), ncol = length(cells),
                                      dimnames = list(receptors, names(cells))))
    tf.act <- data.frame(matrix(nrow = length(tfs), ncol = length(cells),
                                dimnames = list(tfs, names(cells))))
    target.act <- data.frame(matrix(nrow = length(targets), ncol = length(cells),
                                    dimnames = list(targets, names(cells))))
    rtfs.ppr.colsums <- colSums(rtfs.ppr)
  }

  # Get the network
  for(j in 1:length(cells)){
    if(is.null(prune)) prune <- mod$defaults$prune
    network <- get.network(seurat.obj, i = cells[j], drop = TRUE, assay = assay, prune = prune, cutoff = cutoff)
    responses <- rownames(network)
    predictors <- colnames(network)

    # Get the probability network
    if(assay == "meta.network" & prune){
      if(is.null(mod$meta.network$p.val)) stop("Meta-network p-val not found.")
      prob.network <- mod$meta.network$p.val[mod$defaults$responses,mod$defaults$predictors,cells[j]]
      network[prob.network < cutoff] <- 0
    }

    if(as.tfs == "responses") network <- t(network)

    network <- network[targets,tfs,drop=F] *
      (adj[targets, tfs,drop=F] != 0)

    if(scale.target.activity){
      # Scale target activity.
      target.scale <- rowSums(network)
      tf.scale <- colSums(network)
      network <- network/target.scale^0.5
      network <- sweep(network, 2, tf.scale^0.5, "/")
      network[is.na(network)] <- 0
    }

    if(length(cells) > 1){
      if(receptor.activity == "cprod"){
        network.colsums <- colSums(network)
        if(check.receptor.expression) rtfs.ppr.colsums <- colSums(rtfs.ppr * rw[,j])
        receptor.act[,j] <- as.numeric(rtfs.ppr %*% network.colsums)
        tf.act[,j] <- as.numeric(rtfs.ppr.colsums*network.colsums)
        target.act[,j] <- as.numeric(network%*%rtfs.ppr.colsums)
        if(check.receptor.expression) receptor.act[,j] <- receptor.act[,j] * rw[,j]
      }else{
        network.colmax <- apply(network, 2, max)
        if(check.receptor.expression) rtfs.ppr.colmax <- apply(rtfs.ppr * rw[,j], 2, max)
        receptor.act[,j] <- sweep(rtfs.ppr, 2, network.colmax, "*") %>% apply(., 1, max)
        tf.act[,j] <- as.numeric(network.colmax*rtfs.ppr.colmax)
        target.act[,j] <- sweep(network, 2, rtfs.ppr.colmax, "*") %>% apply(., 1, max)
        if(check.receptor.expression) receptor.act[,j] <- receptor.act[,j] * rw[,j]
      }
    }else{
      if(receptor.activity == "cprod"){
        rt <- tcrossprod(rtfs.ppr, network)
      }else{
        rt <- matrix(0, nrow = nrow(rtfs.ppr), ncol = length(targets),
                     dimnames = list(receptors, targets))
        for(j in targets){
          rt[,j] <- sweep(rtfs.ppr, 2, network[j,], "*") %>% apply(., 1, max)
        }
      }
    }
  }
  if(length(cells) > 1){
    return(list(receptor.act = receptor.act, tf.act = tf.act, target.act = target.act) %>% invisible)
  }else{
    return(list(rt = rt, rtfs.ppr = rtfs.ppr, network = network) %>% invisible)
  }
}


prepare.visualise <- function(seurat.obj, n.clu = 4, as.g2 = c("predictors", "responses"),
                              central.genes = NULL, g1 = NULL, g2 = NULL, receptors = NULL, t = 2){
  setting <- Seurat::Misc(seurat.obj, "setting")
  if(is.null(setting)) stop("Run prepare.graph first, and then prepare.reg.")
  mod <- Seurat::Misc(seurat.obj, "mod")
  if(is.null(mod)) stop("Run run.nn.reg first.")

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
  loadings <- setting$loadings[g12, ] %>% sweep(., 2, d, "*")
  clu.g12 <- hclust(cor(t(loadings)) %>% dist) #Will be used to order layouts.

  # Set hub genes
  if(length(g2) > 1){
    n.clu <- min(n.clu, length(g2))
    hubs <- cbind(data.frame(cluster = cutree(clu.g12, k = n.clu), genes = g12),
                  data.frame(loadings)) %>%
      plyr::ddply(., "cluster", function(x){
        genes <- x[,2]
        x <- x[,-c(1,2)]
        v <- svd(as.matrix(x), nu=0, nv=1)$v[,1,drop=F]
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
  rg2.ppr <- Matrix::Matrix(0, nrow = nrow(rt.ppr),
                            ncol = length(g2.full), dimnames = list(rownames(rt.ppr), g2.full))
  rg2.ppr[,intersect(g2.full, colnames(rt.ppr))] <-
    rt.ppr[,intersect(g2.full, colnames(rt.ppr))]

  # Set receptors
  if(!is.null(receptors)){
    receptors <- intersect(receptors, rownames(rg2.ppr))
  }

  # Build up evidences
  evidence <- matrix(1, nrow = length(g1), ncol = length(g2.full), dimnames = list(g1,g2.full))
  adj <- get.gr.adj(t = t) %>% t
  known.g1 <- g1[g1%in%names(igraph::V(gr.graph))]
  known.g2 <- g2.full[g2.full%in%names(igraph::V(gr.graph))]
  evidence[known.g1,known.g2] <- as.matrix(
    evidence[known.g1,known.g2] +
      (adj[known.g1, known.g2] > 0) + 2*(adj[known.g1, known.g2] < 0)
  )
  suppressWarnings(
    Seurat::Misc(seurat.obj, "visual.setting") <-
      list(g1 = g1, g2 = g2, clu.g12 = clu.g12, hubs = hubs, receptors = receptors, rg2.ppr = rg2.ppr, as.g2 = as.g2,
           evidence = evidence)
  )
  seurat.obj
}
visualise.network <- function(seurat.obj, i, meta.network = FALSE,
                              fix.cluster = TRUE,
                              hubs = NULL, n.clu = 4,
                              cutoff = 0.5,
                              show.pathways = TRUE,
                              change.receptors = T,
                              receptor.activity = c("cprod", "dist"),
                              check.receptor.expression = TRUE,
                              scale.receptor.activity = F,
                              scale.g1.activity = F,
                              scale = FALSE,
                              swap.layers = FALSE,
                              k = 2,
                              n.per.component = 10, radius  = NULL, pie.radius = 0.05, text.size = 4){
  setting <- Seurat::Misc(seurat.obj, "setting")
  if(is.null(setting)) stop("Run prepare.graph first, and then prepare.reg.")
  mod <- Seurat::Misc(seurat.obj, "mod")
  if(is.null(mod)) stop("Run run.nn.reg first, and then prepare.visualise.")
  visual.setting <- Seurat::Misc(seurat.obj, "visual.setting")
  if(is.null(visual.setting)) stop("Run prepare.visualise first.")

  g1 <- visual.setting$g1
  g2 <- visual.setting$g2
  receptors <- visual.setting$receptors
  as.g2 <- visual.setting$as.g2
  rg2.ppr <- visual.setting$rg2.ppr
  evidence <- visual.setting$evidence
  g2.full <- colnames(rg2.ppr)
  assay <- if(meta.network) "meta.network" else mod$defaults$assay
  remove.self.loops <- mod$defaults$remove.self.loops

  #################### network processing ########################
  # Get the network
  network <- get.network(seurat.obj, i = i, drop = TRUE, assay = assay, prune = FALSE)
  responses <- rownames(network)
  predictors <- colnames(network)

  # Get the probability network
  if(assay == "effect"){
    prob.network <- get.network(seurat.obj, i = i, drop = TRUE, assay = "p.val", prune = FALSE)
  }else if(assay == "p.val"){
    prob.network <-  network
  }else{
    if(is.null(mod$meta.network$p.val)) stop("Meta-network p-val not found.")
    prob.network <- mod$meta.network$p.val[responses,predictors,i]
    prob.network <- matrix(prob.network, nrow = length(responses), ncol = length(predictors),
                           dimnames = list(responses, predictors))
  }

  if(as.g2 == "responses"){
    network <- t(network)
    prob.network <- t(prob.network)
  }

  g1 <- intersect(g1, rownames(network))
  g2 <- intersect(g2, colnames(network))
  if(mod$custom.y){
    common.node <- NULL
  }else{
    common.node <- intersect(g1,g2)
  }

  # Remove self-loop
  for(r in common.node){
    network[r,r] <- 0
    prob.network[r,r] <- 0
  }

  # Allow duplicated gene names in the first and second layer
  names(g1) <- paste("g1", g1, sep = ".")
  names(g2) <- paste("g2", g2, sep = ".")

  # g1 g2 node weight
  g1.weight <- apply(prob.network, 1, max)
  g2.weight <- apply(prob.network, 2, max)

  # Scale the weight
  if(scale){
    g1.weight <- g1.weight/max(g1.weight)
    g2.weight <- g2.weight/max(g2.weight)
    network <- network/max(network)
  }else if(assay != "p.val"){
    network <- network/max(network)*max(prob.network)
  }

  ################# Hierarchical clustering g2 genes   #################
  if(length(g2)>1){
    n.clu <- min(n.clu, length(g2))
    if(fix.cluster){
      clu.g2 <- cutree(visual.setting$clu.g12, n.clu)[g2]
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
      names(hubs) <- paste("g2", hubs, sep = ".")
      if(length(g1)>1){
        clu.g2 <- cor(network[,names(hubs)], network[,g2,drop=F])
        clu.g2 <- apply(clu.g2, 2, which.max)
      }else{
        clu.g2 <-  sapply(names(hubs), function(x) sqrt(colSums(network[,g2,drop=F]-network[,x])^2))
        clu.g2 <- apply(clu.g2, 1, which.min)
      }
      n.clu <- length(hubs)
    }
  }else{
    clu.g2 <- 1
  }

  ################# Receptor to target activity   #################
  if(show.pathways){
    pruned.network <- network[g1,g2,drop=F] *
      (evidence[g1, g2, drop=F] > 1) *
      (prob.network[g1, g2, drop=F] > cutoff)
    rg2.ppr <- rg2.ppr[,g2,drop=F]

    if(scale.g1.activity){
      # g1.scale <- sqrt(rowSums(pruned.network^2))
      g1.scale <- rowSums(pruned.network)
      g2.scale <- colSums(pruned.network)
      pruned.network <- pruned.network/g1.scale^0.5
      pruned.network <- sweep(pruned.network, 2, g2.scale^0.5, "/")
      pruned.network[is.na(pruned.network)] <- 0

    }

    if(scale.receptor.activity){
      # Scale receptor activity
      receptor.scale <- rowSums(rg2.ppr)
      g2.scale <- colSums(rg2.ppr)
      rg2.ppr <- rg2.ppr/receptor.scale^0.5
      rg2.ppr <- sweep(rg2.ppr, 2, g2.scale^0.5, "/")
      rg2.ppr[is.na(rg2.ppr)] <- 0
    }

    receptor.activity <- match.arg(receptor.activity)
    if(receptor.activity == "cprod"){
      rg1 <- tcrossprod(rg2.ppr, pruned.network)
    }else{
      rg1 <- matrix(0, nrow = nrow(rg2.ppr), ncol = length(g1),
                    dimnames = list(rownames(rg2.ppr), g1))
      for(j in g1){
        rg1[,j] <- sweep(rg2.ppr, 2, pruned.network[j,], "*") %>% apply(., 1, max)
      }
    }

    # Check maximum activity
    max.act <- apply(rg1, 2, function(x) max(x))

    # Whether to weight receptors by their expression.
    if(check.receptor.expression){
      # Check receptor expression
      rw <- Matrix::Matrix(0, nrow = nrow(rg1),
                           ncol = nrow(setting$pcs), dimnames = list(rownames(rg1), rownames(setting$pcs)))
      rw[intersect(rownames(seurat.obj), rownames(rg1)),rownames(setting$pcs)] <-
        SeuratObject::LayerData(seurat.obj, "data")[intersect(rownames(seurat.obj), rownames(rg1)),rownames(setting$pcs)] != 0

      # Network propagation
      if(assay != "meta.network"){
        w <- tcrossprod(mod$w$u[i,,drop=F], mod$w$vd)
        rw <- tcrossprod(rw, w)
      }else{
        w <-  mod$w$vd %*%
          crossprod(mod$w$u[mod$cells,], mod$meta.network$npca.loadings[,i])
        w <-  w/sum(w)
        rw <- rw %*% w
      }
      rw[rw > 1] <- 1
      rw[rw < 0] <- 0
      rg1 <- rg1 * as.numeric(rw)
    }

    k.min <- min(length(g1), nrow(rg2.ppr), sum(colSums(rg1) > 0))
    k <- min(k, k.min)
    if(k == 0) show.pathways <- FALSE
    max.act <- max.act[colSums(rg1)>0]
    rg1 <- rg1[,colSums(rg1)>0,drop=F]
  }

  g1.weight <- g1.weight[g1]
  g2.weight <- g2.weight[g2]
  rg2.ppr <- rg2.ppr[,g2,drop=F]
  network <- network[g1,g2,drop=F]
  evidence <- evidence[,g2,drop=F]
  prob.network <- prob.network[g1,g2,drop=F]
  dimnames(prob.network) <- dimnames(network) <- list(names(g1), names(g2))
  names(g1.weight) <- names(g1)
  names(g2.weight) <- names(clu.g2) <- names(g2)

  ################# Evidence processing   #################
  # Check evidence for the receptor to TF network.
  check.evidence <- function(graph){
    from <- graph[,1]
    to <- graph[,2]
    consensus_stimulation <- igraph::E(sig.graph)$consensus_stimulation
    consensus_inhibition <- igraph::E(sig.graph)$consensus_inhibition

    evidence <- rep(2, nrow(graph))
    e.ids <- igraph::get.edge.ids(sig.graph,
                                  data.frame(from, to) %>% t)

    evidence <- evidence + consensus_stimulation[e.ids] + 2*consensus_inhibition[e.ids]
    evidence[evidence>4] <- 2
    evidence
  }

  # Check evidence for the g2 to g1 network.
  if(!is.null(cutoff)){
    evidence <- (prob.network > cutoff) * evidence + 1
  }else{
    evidence <- evidence+1
  }

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
    # Find g4
    if(change.receptors){
      which.max.n <- function(x, n = n.per.component){
        which(rank(-x, ties.method = "random" ) <= n)
      }
      # Eigen centrality.
      if(k.min >= 3){
        svd.mod <- RSpectra::svds(rg1, k = k) %>% suppressWarnings()
        g4 <- rownames(rg2.ppr)[apply(abs(svd.mod$u), 2, which.max.n) %>% as.numeric]
      }else{
        g4 <- rownames(rg2.ppr)[apply(rg1,
                                      2, which.max.n) %>% as.numeric]
      }
      g4 <- unique(g4)

    }else{
      g4 <- unique(receptors)
    }
    g4.weight <- sweep(rg1[g4,,drop=F], 2, max.act, "/") %>%
      apply(., 1, max)
    g4 <- g4[rowSums(rg2.ppr[g4,g2, drop = F])>0]
    g4 <- g4[g4.weight[g4]>0]
    g4.weight <- g4.weight[g4]

    # Create g3 to g4 network
    rg2.ppr <- rg2.ppr/apply(rg2.ppr, 1, max)
    if(check.receptor.expression) rg2.ppr <- rg2.ppr * as.numeric(rw)
    network1 <- sweep(rg2.ppr[g4,,drop=F], 2, g2.weight * (colSums(pruned.network[,g2]) >0),"*")
    network1 <- network1[rowSums(network1) > 0,]
    g4 <- rownames(network1)
    g4.weight <- g4.weight[g4]

    closestTF <- g2[apply(network1, 1, which.max)]

    # Find shortest paths between TFs and their closest receptors.
    paths <- c()
    for (i in seq(length(g4))){
      receptor <- g4[i]
      sp <- igraph::shortest_paths(sig.graph, receptor, closestTF[i], weights = 1/igraph::E(sig.graph),
                                   output = c("both"))
      paths <- c(paths , sp$vpath)
    }
    names(paths) <- closestTF
    paths <- paths[order(sapply(closestTF, function(x) which(g2==x)[1]))]
    closestTF <- names(paths)
    pass1 <- sapply(paths, function(x) names(x)[2])
    pass2 <- sapply(paths, function(x) names(x)[length(x)-1])
    g3 <- sapply(paths, function(x) names(x)[-c(1,length(x))])
    g3 <- sapply(g3, function(x) paste(x, collapse = "->"))
    g4 <- sapply(paths, function(x) names(x)[1])
    if(sum(g3 == "") > 0){
      g30 <- paste("NULL", seq(sum(g3 == "")), sep = "")
      g3[g3 == ""] <- g30
    }else{
      g30 <- NULL
    }
    g4.weight <- g4.weight[g4]
    network1 <- network1[g4,,drop=F]

    names(g3) <- paste("g3", g3, sep = ".")
    names(g4) <- paste("g4", g4, sep = ".")
    if(!is.null(g30)) names(g30) <- paste("g3", g30, sep = ".")
    names(pass1) <- paste("g3", pass1, sep = "." )
    names(pass1)[g3 %in% g30] <- paste("g2", pass1, sep = "." )[g3 %in% g30]
    names(pass2) <- paste("g3", pass2, sep = "." )
    names(pass2)[g3 %in% g30] <- paste("g4", pass2, sep = "." )[g3 %in% g30]

    network.df <- data.frame(from = g4, to = pass1, weight =  g4.weight, cluster = 1)
    network.df$evidence <- check.evidence(network.df)
    network.df$from <- names(g4)
    network.df$to <- names(pass1)
    network.df[!g3%in%g30,2] <-  names(g3)[!g3%in%g30]
    network.df$weight[network.df$weight < 0.05] <- 0.05

    graph <- rbind(igraph::as_data_frame(graph),network.df) %>% igraph::graph_from_data_frame()

    ################# g3 to g2 layer  #################
    network.df <- data.frame(from = pass2, to = closestTF,
                             weight = network.df$weight, cluster = 1)
    network.df$evidence <- check.evidence(network.df)
    network.df$to <- paste("g2", closestTF, sep = ".")
    network.df$from <- names(g3)
    network.df <- network.df[!g3%in%g30,]
    network.df$weight[network.df$weight < 0.05] <- 0.05
    graph <- rbind(igraph::as_data_frame(graph),network.df) %>% igraph::graph_from_data_frame()

    ################# Clean up  #################
    rm(network.df)
    g3 <- g3[!duplicated(g3)]
    g3.weight <- rep(0, length(g3))
    if(!is.null(g30)) graph <- graph %>% igraph::add_vertices(nv = length(g30), name = names(g30))
    graph <- igraph::permute.vertices(graph, permutation = sapply(names(igraph::V(graph)),
                                                                  function(x) which(c(names(g1), names(g2), names(g3), names(g4)) == x)))
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
  pie.values <- plyr::ddply(pie.values, "cluster",colMeans) %>% t %>% as.data.frame()
  pie.values <- pie.values[-nrow(pie.values),,drop=F]
  pie.values <- pie.values/rowSums(pie.values)*g1.weight
  pie.values <- rbind(pie.values,
                      matrix(0, nrow = length(g2), ncol = ncol(pie.values),
                             dimnames = list(names(g2), colnames(pie.values))))
  for(i in names(clu.g2)) pie.values[i,clu.g2[i]] <- g2.weight[i]
  pie.values <- cbind(data.frame(N=1-rowSums(pie.values)), pie.values)

  if(show.pathways){
    # Pie for g3
    pie.values1 <- matrix(0, ncol = length(unique(clu.g2)), nrow = length(g3)) %>% as.data.frame()
    rownames(pie.values1) <- names(g3)
    pie.values1 <- cbind(data.frame(N=1-rowSums(pie.values1)), pie.values1)

    # Pie for g4
    pie.values2 <- t(network1) %>% as.data.frame()
    pie.values2$cluster <- clu.g2
    pie.values2 <- plyr::ddply(pie.values2, "cluster", colMeans) %>% t %>% as.data.frame()
    pie.values2 <- pie.values2[-nrow(pie.values2),,drop=F]
    pie.values2 <- pie.values2/rowSums(pie.values2)*g4.weight
    pie.values2 <- cbind(data.frame(N=1-rowSums(pie.values2)), pie.values2)
    rownames(pie.values2) <- paste("g4",rownames(pie.values2),sep = ".")

    # Combine Pies
    pie.values <- rbind(pie.values, pie.values1, pie.values2)
    pie.values <- pie.values[names(igraph::V(graph)),,drop=F]
    pie.values <- cbind(data.frame(cluster = igraph::V(graph)$cluster),pie.values)
    pie.values <- cbind(data.frame(scale = 1),pie.values)
    pie.values[names(g3),]$scale <- 0
    pie.values <- cbind(l,pie.values)

  }else{
    pie.values <- pie.values[names(igraph::V(graph)),,drop=F]
    pie.values <- cbind(data.frame(cluster = igraph::V(graph)$cluster),pie.values)
    pie.values <- cbind(data.frame(scale = 1),pie.values)
    pie.values <- cbind(l,pie.values)

  }

  ################# Aes  #################
  col <- RColorBrewer::brewer.pal("Set1", n = 9)
  col1 <- c("white",col) %>% alpha(.,0.8)
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
  make_segements <- function(graph, shorten.start = 0.07, shorten.end = 0.07){
    es <- igraph::get.edgelist(graph)
    data <- data.frame(matrix(0 , nrow = length(igraph::E(graph)), ncol = 4))
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
  segement_data <- make_segements(graph)

  p <- ggraph(graph, layout = l)  +

    geom_segment(data = segement_data,aes(x=x,y=y,xend=xend,yend=yend),
                 arrow = grid::arrow(type = "closed",
                                     angle = arrow.angle[igraph::E(graph)$evidence],
                                     length = unit(0.15, "inches")),
                 color = col2[igraph::E(graph)$cluster],
                 linewidth = igraph::E(graph)$weight*2,
                 linetype = lty[igraph::E(graph)$evidence])+

    geom_scatterpie(
      data = pie.values,
      cols = colnames(pie.values)[-c(1,2,3,4)],
      aes(x = x,y = y, col= as.character(cluster), r = pie.radius * scale),
      alpha = 0.5,
    )+

    geom_node_text(aes(label = vertex.label), col = label.col[igraph::V(graph)$type], size = text.size)+

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


plot.receptor.expression <- function(act){
  rt <- act$rt
  rtfs.ppr <- act$rtfs.ppr
  network <- act$network
  targets <- rownames(network)
  tfs <- colnames(act$rtfs.ppr)
  receptors <- rownames(rt)

  # PLS plot
  pls.mod <- svd(rt)

  lx <- pls.mod$u[,1:2]
  lx <- lx/max(sqrt(rowSums(lx^2)))
  ly <- pls.mod$v[,1:2]
  ly <- ly/max(sqrt(rowSums(ly^2)))

  sx <- crossprod(rtfs.ppr, lx)
  sx <- sx/max(sqrt(rowSums(sx^2)))
  sy <- crossprod(network, ly)
  sy <- sy/max(sqrt(rowSums(sy^2)))

  lx <- data.frame(lx)
  ly <- data.frame(ly)
  sx <- data.frame(sx)
  sy <- data.frame(sy)
  receptors <- rownames(rt)
  rownames(lx) <- receptors
  rownames(ly) <- targets
  rownames(sx) <- rownames(sy) <- tfs


  p1 <- ggplot() +
    geom_point(data = sx, aes(X1,X2)) +
    geom_point(data = lx, aes(X1,X2), col = "blue") +
    geom_point(data = ly, aes(X1,X2), col = "red")+
    geom_text_repel(data = sx, aes(X1,X2, label = tfs), col = "black") +
    geom_text_repel(data = lx, aes(X1,X2, label = receptors), col = "blue") +
    geom_text_repel(data = ly, aes(X1,X2, label = targets), col = "red") +
    ylim(c(-1,1))+ xlim(c(-1,1)) + xlab("Component 1") + ylab("Component2 ")+
    theme_classic()+
    theme(text = element_text(size = 15, face = "bold"))

  p2 <- ggplot() +
    geom_point(data = sy, aes(X1,X2, col = "TF")) +
    geom_point(data = lx, aes(X1,X2, col = "Receptor")) +
    geom_point(data = ly, aes(X1,X2, col = "Target"))+
    geom_text_repel(data = sy, aes(X1,X2, label = tfs), col = "black") +
    geom_text_repel(data = lx, aes(X1,X2, label = receptors), col = "blue") +
    geom_text_repel(data = ly, aes(X1,X2, label = targets), col = "red") +
    ylim(c(-1,1))+ xlim(c(-1,1)) + xlab("Component 1") + ylab("Component2 ")+
    scale_color_manual("",values = c(TF = "black", Receptor = "blue", Target = "red")) +
    theme_classic()+
    theme(text = element_text(size = 15, face = "bold"))

  print(p1+p2)

  #if(do.pls) print(p1+p2)

  # lx <- lx[order(-abs(lx[,k])),]
  # lx <- head(lx, n.top)
  # ly <- ly[order(-abs(ly[,k])),]
  # ly <- head(ly, n.top)
  # sx <- sx[order(-abs(sx[,k])),]
  # sx <- head(sx, n.top)
  # sy <- sy[order(-abs(sy[,k])),]
  # sy <- head(sy, n.top)
  #
  # p3 <- ggplot() + geom_col(aes(x=lx[,k], y = factor(rownames(lx), levels = rev(rownames(lx))))) + ylab("") + xlab("Receptor score")
  # p4 <- ggplot() + geom_col(aes(x=ly[,k], y = factor(rownames(ly), levels = rev(rownames(ly))))) + ylab("") + xlab("Target score")
  # p5 <- ggplot() + geom_col(aes(x=sx[,k], y = factor(rownames(sx), levels = rev(rownames(sx))))) + ylab("") + xlab("TF score (Receptor-TF potential)")
  # p6 <- ggplot() + geom_col(aes(x=sy[,k], y = factor(rownames(sy), levels = rev(rownames(sy))))) + ylab("") + xlab("TF score (Target-TF co-expression)")
  # p <- p3 + p4 + p5 + p6
  # p <- p & theme_classic() & theme(text = element_text(face = "bold", size = 15),
  #                                  legend.position = "none",
  #                                  axis.line = element_blank())
  # print(p)

  # mat1 <- rtfs.ppr[rownames(lx), rownames(sy)]
  # mat1[mat1==0] <- NA
  # mat2 <- network[rownames(ly), rownames(sy)]
  # mat2[mat2==0] <- NA
  # pheatmap::pheatmap(mat1, cluster_cols = F, cluster_rows = F,
  #                    color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Reds"))(100), na_col = "white")
  # pheatmap::pheatmap(mat2, cluster_cols = F, cluster_rows = F,
  #                    color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Reds"))(100), na_col = "white")


}


# Acquired from rTRM
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

layout.concentric = function (g, concentric = NULL, radius = NULL, order.by)
{
  if(is.null(concentric))
    concentric = list(V(g)$name)

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
    order.values = lapply(order.by, function(b) get.vertex.attribute(g, b))

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


