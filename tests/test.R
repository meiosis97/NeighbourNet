require(NeighbourNet)
load("tests/data/luad.rda")

get.prior.model()

get.gr.adj()

gene.list <- select.gene(obj)

obj <- prepare.seurat(obj, genes = gene.list$genes)

obj <- prepare.graph(obj)

obj <- select.cell(obj, all = TRUE)

obj <- prepare.reg(obj, responses = gene.list$tfs, predictors = gene.list$targets)

str(obj@misc$NNet.setting)

