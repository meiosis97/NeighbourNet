require(devtools)
document()
build()
install()


require(Seurat)
require(NeighbourNet)
require(ggplot2)
require(dplyr)
require(Matrix)
load("tests/data/luad.rda")
obj <- Seurat::FindVariableFeatures(obj)
genes <- Seurat::VariableFeatures(obj)

get.prior.model() %>% str

get.gr.adj() %>% str

gene.list <- select.gene(obj)

obj <- prepare.seurat(obj, genes = genes)

obj <- prepare.graph(obj)

obj <- select.cell(obj)

obj <- prepare.reg(obj, responses = genes)

obj <- run.nn.reg(obj, responses = genes[1:5], return.p.val = T)

obj <- set.defaults(obj)

obj <- set.defaults(obj, defaults = c(cutoff = 2))

obj <- set.defaults(obj, clean.up = TRUE)

str(Seurat::Misc(obj, "NNet.mod"))

str(Seurat::Misc(obj, "NNet.setting"))

get.network(obj) %>% str

get.network(obj, i = 1)

get.network(obj, i = 2) %>% str

get.network(obj, i = c(2,4)) %>% str

get.network(obj, i = 2,
                       assay = "effect", responses = genes[2:3], predictors = sample(genes, 10)) %>% str

get.network(obj, assay = "p.val", cutoff = 0.9) %>% str

obj <- build.meta.network(obj)

obj <- build.meta.response(obj)

cells <- Misc(obj)$NNet.mod$cells

select.central.genes(obj) %>% str

select.central.genes(obj, k = 2) %>% str

select.central.genes(obj, k = 2, keep.responses = TRUE) %>% str

select.central.genes(obj, k = 2, n.net = 1) %>% str

receptor.activity(obj) %>% str

receptor.activity(obj, meta.network = TRUE) %>% str

receptor.activity(obj, i = 2) %>% str

ctr.genes <- select.central.genes(obj) 

obj <- prepare.visualise(obj,  central.genes = ctr.genes) 

Seurat::Misc(obj, "NNet.visual.setting") %>% str

obj <- prepare.visualise(obj,  central.genes = ctr.genes, as.g2 = "responses") 

Seurat::Misc(obj, "NNet.visual.setting") %>% str

###################  Check results ###################
pcs <- Seurat::Misc(obj, "NNet.setting")$pcs[cells,]
lra <- Seurat::Misc(obj, "NNet.setting")$lra[cells,]
obj <- Seurat::RunUMAP(obj, dims = 1:ncol(pcs))
umap <- data.frame(Seurat::Embeddings(obj, "umap"))[cells,]

gene <- genes[1]
predictor <- rowSums(2*obj@misc$NNet.mod$effect[gene, , ]^2) %>%
  sort(decreasing = T) %>% names %>% dplyr::nth(1)

Seurat::FeaturePlot(obj, features = c(gene, predictor), reduction = "umap")
ggplot() +
  geom_point(data = umap, aes(umap_1, umap_2, col =
                                2*obj@misc$NNet.mod$effect[gene,predictor,]^2))+
  scale_color_gradientn(colours = gg.color.spec(11))

ggplot() +
  geom_point(data = umap, aes(umap_1, umap_2, col =
                                obj@misc$NNet.mod$p.val[gene,predictor,])) +
  scale_color_gradientn(colours = gg.color.spec(11), limits = c(0,1))

ggplot() + geom_point(aes(lra[,predictor],lra[,gene], col =
                            2*obj@misc$NNet.mod$effect[gene,predictor,]^2))+
  scale_color_gradientn(colours = gg.color.spec(11))

ggplot() + geom_point(aes(lra[,predictor],lra[,gene], col =
                            obj@misc$NNet.mod$p.val[gene,predictor,]))+
  scale_color_gradientn(colours = gg.color.spec(11), limits = c(0,1))

ggplot() +
  geom_point(data = umap, aes(umap_1, umap_2, col =
                                obj@misc$NNet.mod$p.val[gene,predictor,] > 0.95))

                                ggplot() +
geom_point(data = umap, aes(umap_1, umap_2, col =
                                obj@misc$NNet.mod$meta.network$npca.loadings[,1]))


###################  Check results2 ###################
# Build the null data matrix
perm.data <- obj@assays$RNA$data[genes, ] %>%
  apply(.,1,sample) %>% t
rownames(perm.data) <- paste("NULL", rownames(perm.data), sep = "-")
n.umi <- colSums(SeuratObject::LayerData(obj, "counts", features = genes))
n.umi <- sample(log(n.umi))
perm.data <- rbind(perm.data, n.umi)

new.obj <- Seurat::CreateSeuratObject(counts = rbind(obj@assays$RNA$data[genes,], perm.data), project = "celline",
                              meta.data = data.frame(obj@meta.data))
SeuratObject::LayerData(new.obj, "data") <- SeuratObject::LayerData(new.obj, "counts")
null.genes <- rownames(perm.data)
expand.genes <- rownames(new.obj)

new.obj <- prepare.seurat(new.obj, genes = expand.genes)
new.obj <- prepare.graph(new.obj)
new.obj <- prepare.reg(new.obj, check.expressed = T)

pcs <- Seurat::Misc(new.obj, "NNet.setting")$pcs
lra <- Seurat::Misc(new.obj, "NNet.setting")$lra
new.obj <- Seurat::RunUMAP(new.obj, dims = 1:ncol(pcs))
umap <- data.frame(Seurat::Embeddings(new.obj, "umap"))

gene <- null.genes[1]
new.obj <- run.nn.reg(new.obj, responses = gene, return.p.val = T)
predictor <-  gene
mu <- new.obj@misc$NNet.mod$mus
sigma <- new.obj@misc$NNet.mod$sigmas

ggplot() + geom_jitter(aes(x= "umi", y= log(abs(new.obj@misc$NNet.mod$effect[gene,"n.umi",])))) +
  geom_jitter(aes(x= "predictor", y= log(abs(new.obj@misc$NNet.mod$effect[gene,predictor,])))) +
  geom_jitter(aes(x= "null", y= log(abs(new.obj@misc$NNet.mod$effect[gene,paste("NULL",predictor, sep = "-"),]))))+
  geom_errorbar(aes(x = "predictor", ymin = mu-1.64*sigma, ymax =mu + 1.64*sigma), color = "darkred",size = 3)+
  geom_point(aes(x = "predictor", y=mu), color = "darkred",size = 5)

rowMeans(log(abs(new.obj@misc$NNet.mod$effect)[1,,]) > mu + 1.64*sigma) %>% sort(decreasing = T) %>%
  plot

log(abs(new.obj@misc$NNet.mod$effect)[1,,])

###################  Check results3 ###################
n.umi <- colSums(SeuratObject::LayerData(obj, "counts", features = genes))
perm.data <- replicate(10, rnorm(ncol(obj))) %>% t #replicate(10, sample(log(n.umi))) %>% t
hk <- paste("HOUSEKEEPING", 1:10, sep="")
rownames(perm.data) <-  hk

new.obj <- Seurat::CreateSeuratObject(counts = rbind(obj@assays$RNA$data[genes,], perm.data), project = "celline",
                                      meta.data = data.frame(obj@meta.data))
SeuratObject::LayerData(new.obj, "data") <- SeuratObject::LayerData(new.obj, "counts")
expand.genes <- rownames(new.obj)

new.obj <- prepare.seurat(new.obj, genes = expand.genes)
new.obj <- prepare.graph(new.obj)
new.obj <- prepare.reg(new.obj)

pcs <- Seurat::Misc(new.obj, "NNet.setting")$pcs
lra <- Seurat::Misc(new.obj, "NNet.setting")$lra
new.obj <- Seurat::RunUMAP(new.obj, dims = 1:ncol(pcs))
umap <- data.frame(Seurat::Embeddings(new.obj, "umap"))

gene <- genes[2]
new.obj <- run.nn.reg(new.obj, responses = gene, return.p.val = T)
predictor <-  gene

null.effect <- as.numeric(new.obj@misc$NNet.mod$effect[gene,hk,])
mu <- mean(log(abs(null.effect)))
sigma <- sd(log(abs(null.effect)))
ggplot() + geom_jitter(aes(x= "umi", y= log(abs(null.effect)))) +
  geom_jitter(aes(x= "predictor", y= log(abs(new.obj@misc$NNet.mod$effect[gene,predictor,])),
             col = as.numeric(SeuratObject::LayerData(obj, "data", features = gene))))+
  geom_errorbar(aes(x = "umi", ymin = mu-1*sigma, ymax =mu + 1*sigma), color = "darkred",size = 3)+
  geom_point(aes(x = "umi", y=mu), color = "darkred",size = 5)+
  scale_color_gradientn(colours = gg.color.spec(11))


###################  Check results4 ###################
source("tests/20241210.R")
load("tests/data/luad.rda")
obj <- Seurat::FindVariableFeatures(obj)
genes <- Seurat::VariableFeatures(obj)
gene.list <-  NeighbourNet::gene.list

n.umi <- colSums(SeuratObject::LayerData(obj, "counts", features = genes))
perm.data <- replicate(10, rnorm(ncol(obj))) %>% t #replicate(10, sample(log(n.umi))) %>% t
hk <- paste("HOUSEKEEPING", 1:10, sep="")
rownames(perm.data) <-  hk

new.obj <- Seurat::CreateSeuratObject(counts = rbind(obj@assays$RNA$data[genes,], perm.data), project = "celline",
                                      meta.data = data.frame(obj@meta.data))
SeuratObject::LayerData(new.obj, "data") <- SeuratObject::LayerData(new.obj, "counts")
expand.genes <- rownames(new.obj)

new.obj <- prepare.seurat(new.obj, genes = expand.genes)
new.obj <- prepare.graph(new.obj)
new.obj <- prepare.reg(new.obj)

pcs <- Seurat::Misc(new.obj, "setting")$pcs
lra <- Seurat::Misc(new.obj, "setting")$lra

gene <- "GAPDH"#genes[1]
new.obj <- run.nn.reg(new.obj, responses = gene, return.p.val = T)
predictor <-  gene

null.effect <- as.numeric(new.obj@misc$mod$effect[gene,hk,])
mu <- new.obj@misc$mod$mus
sigma <- new.obj@misc$mod$sigmas
mu1 <- mean(log(abs(null.effect)))
sigma1 <- sd(log(abs(null.effect)))
ggplot() + geom_jitter(aes(x= "umi", y= log(abs(null.effect)))) +
  geom_jitter(aes(x= "predictor", y= log(abs(new.obj@misc$mod$effect[gene,predictor,])),
                  col = as.numeric(SeuratObject::LayerData(obj, "data", features = gene)))) +
  geom_errorbar(aes(x = "umi", ymin = mu-0.64*sigma, ymax =mu + 0.64*sigma), color = "darkred",size = 3)+
  geom_errorbar(aes(x = "umi", ymin = mu1-1.65*sigma1, ymax =mu1 + 1.65*sigma1), color = "darkblue",size = 3)+
  geom_point(aes(x = "umi", y=mu), color = "darkred",size = 5)+
  geom_point(aes(x = "umi", y=mu1), color = "darkblue",size = 5)+
  scale_color_gradientn("",colours = gg.color.spec(11))

print(mu)
print(sigma)
FeaturePlot(obj,gene)



###################  Check results4 ###################
special.run.nn.reg <- function(seurat.obj, responses = NULL, Y = NULL,
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
    r <- responses[i]

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


    # Transform regression coefficients to effects
    b <- tcrossprod(loadings, b) %>% as.matrix

    # Calculate effect
    effect <- (b[genes,] * setting$nn.scale.gene[genes,names(cells),drop=F]) %>%
      as.matrix

    # Calculate p-value
    if(!custom.y){
      noise <- effect[i,]
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
    mus[i] <- mu
    sigmas[i]  <- sigma

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
special.run.nn.reg <- function(seurat.obj, responses = NULL, Y = NULL,
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
    effect <- (b* setting$nn.scale.gene[genes,names(cells),drop=F]) %>%
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
    custom.y = custom.y, w = list(u = u, vd = vd), noise = noise
  )

  class(mod) <- "NNet.mod"

  # Store results in the Seurat object
  suppressWarnings(
    Seurat::Misc(seurat.obj, "NNet.mod") <- mod
  )

  return(seurat.obj)  # Return the updated Seurat object
}


perm.data <- obj@assays$RNA$data[genes, ] %>%
  apply(.,1,sample) %>% t
rownames(perm.data) <- paste("NULL", rownames(perm.data), sep = "-")
n.umi <- colSums(SeuratObject::LayerData(obj, "counts", features = genes))
n.umi <- sample(log(n.umi))
perm.data <- rbind(perm.data, n.umi)

new.obj <- Seurat::CreateSeuratObject(counts = rbind(obj@assays$RNA$data[genes,], perm.data), project = "celline",
                                      meta.data = data.frame(obj@meta.data))
SeuratObject::LayerData(new.obj, "data") <- SeuratObject::LayerData(new.obj, "counts")
null.genes <- rownames(perm.data)
expand.genes <- rownames(new.obj)

new.obj <- prepare.seurat(new.obj, genes = expand.genes)
new.obj <- prepare.graph(new.obj)
new.obj <- prepare.reg(new.obj)

response.list <- tapply(genes, cut(1:length(genes), breaks = 200), list)

# Calculate power and FDR for one set of responses
one.iter <- function(responses){
  predictors <- paste("NULL",responses,sep="-")
  #new.obj <- run.nn.reg(new.obj, responses, return.p.val = T, predictors = predictors)
  new.obj <- special.run.nn.reg(new.obj, responses, return.p.val = T, predictors = predictors)

  # Calculate power and FDR for one response
  evaluation <- function(p){
    null.p <- paste("NULL", p, sep = "-")
    p.val <- new.obj@misc$NNet.mod$p.val[p,p,]
    null.p.val <- new.obj@misc$NNet.mod$p.val[p,null.p,]
    expressed <- new.obj@assays$RNA$data[p,] != 0
    plot.df <- data.frame(expressed = expressed, p.val = p.val, null.p.val = null.p.val)
    plot.df %>% group_by(expressed) %>%
      dplyr::summarise(power = mean(p.val > 0.9), fdr =  mean(null.p.val > 0.9)) %>%
      dplyr::select(-1)%>%
      unlist(use.names = T)
  }

  sapply(responses, evaluation)
}

# Test
one.iter(response.list[[10]])
noise <- special.run.nn.reg(new.obj, "GAPDH", return.p.val = T)@misc$NNet.mod$noise
hist(noise, breaks = 50)
ggplot() + geom_histogram(aes(noise))
