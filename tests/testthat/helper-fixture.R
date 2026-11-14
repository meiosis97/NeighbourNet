# Small simulated data set run through the public NeighbourNet pipeline, for
# testing visualise.network(). Built once per test session and cached.
#
# Genes are chosen so that the prior knowledge is exercised: TF -> target
# pairs come from gr.graph (some with stimulation / inhibition evidence),
# target expression follows its prior TFs, and receptors from gene.list are
# expressed so that receptor -> TF pathways can be drawn.

.fixture.env <- new.env()

nnet_fixture <- function() {
  if (!is.null(.fixture.env$obj)) return(.fixture.env$obj)
  set.seed(42)

  gl <- NeighbourNet::gene.list
  pri <- igraph::as_data_frame(NeighbourNet::gr.graph, "edges")
  pri <- pri[pri$from %in% gl$tfs & pri$to %in% gl$targets &
             !pri$to %in% gl$tfs & pri$from != pri$to, ]
  # TFs with several prior targets, each contributing up to 4 targets
  tf.n <- sort(table(pri$from), decreasing = TRUE)
  tfs <- names(tf.n)[seq(1, 60, by = 6)]                  # 10 TFs with priors
  pairs <- do.call(rbind, lapply(tfs, function(tf) utils::head(pri[pri$from == tf, ], 4)))
  targets <- unique(pairs$to)
  other.tfs <- utils::head(setdiff(gl$tfs, c(tfs, targets)), 15)
  receptors <- utils::head(setdiff(unique(gl$receptors), c(tfs, other.tfs, targets)), 25)

  n <- 250
  act <- matrix(stats::rnorm(n * length(tfs)), n, dimnames = list(NULL, tfs))
  # smooth latent trajectory so local neighbourhoods carry co-expression
  t <- sort(stats::runif(n))
  act <- act * 0.5 + outer(t, seq_along(tfs), function(t, k) sin(2 * pi * (t + k / 10)))
  mu <- function(x) exp(1 + 0.8 * scale(x))
  expr <- cbind(
    act,
    sapply(targets, function(g) {
      p <- pairs[pairs$to == g, ]
      sgn <- ifelse(p$consensus_inhibition == 1, -1, 1)
      drop(act[, p$from, drop = FALSE] %*% sgn) + stats::rnorm(n, sd = 0.5)
    }),
    matrix(stats::rnorm(n * length(other.tfs)), n, dimnames = list(NULL, other.tfs)),
    matrix(stats::rnorm(n * length(receptors)), n, dimnames = list(NULL, receptors)) +
      act[, rep_len(seq_along(tfs), length(receptors))]
  )
  counts <- matrix(stats::rpois(length(expr), mu(expr)), nrow(expr),
                   dimnames = list(NULL, colnames(expr)))
  counts <- t(counts)
  colnames(counts) <- paste0("cell", seq_len(n))

  obj <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE)))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  # run a pipeline step without its progress output
  quiet <- function(expr) {
    out <- NULL
    utils::capture.output(out <- suppressWarnings(suppressMessages(expr)))
    out
  }
  genes <- NeighbourNet::select.gene(obj, min.cells = 10)
  obj <- quiet(NeighbourNet::prepare.seurat(obj, genes = genes$genes, npcs = 20))
  obj <- quiet(NeighbourNet::prepare.graph(obj, knn = 20))
  obj <- quiet(NeighbourNet::prepare.reg(obj, predictors = genes$tfs,
                                         responses = genes$targets))
  obj <- quiet(NeighbourNet::run.nn.reg(obj, responses = genes$targets,
                                        return.p.val = TRUE))
  obj <- quiet(NeighbourNet::build.meta.network(obj, n.net = 3))
  cg <- quiet(NeighbourNet::select.central.genes(obj, n.net = 2, k = 2,
                                                 n.per.component = 4))
  obj <- quiet(NeighbourNet::prepare.visualise(obj, central.genes = cg, n.clu = 3))
  .fixture.env$obj <- obj
  obj
}

# Mean of each target's strongest meta-network significance: the cutoff used
# in the package article.
meta_cutoff <- function(obj, i = 1) {
  p <- Seurat::Misc(obj, "NNet.mod")$meta.network$p.val[, , i]
  mean(apply(p, 1, max), na.rm = TRUE)
}
