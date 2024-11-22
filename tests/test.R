require(NeighbourNet)
require(ggplot2)
require(dplyr)
require(Matrix)
load("tests/data/luad.rda")
obj <- Seurat::FindVariableFeatures(obj)
genes <- Seurat::VariableFeatures(obj)

get.prior.model()

get.gr.adj()

gene.list <- select.gene(obj)

obj <- prepare.seurat(obj, genes = genes)
n.cell <- ncol(obj)
perm.data <- replicate(10, rnorm(n.cell)) %>% t
rownames(perm.data) <- paste("BACKGROUND", 1:10)
SeuratObject::LayerData(obj, "scale.data") <-
  rbind(SeuratObject::LayerData(obj, "scale.data"), perm.data)

test.function


obj <- prepare.graph(obj)

obj <- select.cell(obj, all = TRUE)

obj <- prepare.reg(obj, responses = genes)

obj <- run.nn.reg(obj, responses = genes[1:5], return.p.val = T)

str(Seurat::Misc(obj, "NNet.mod"))

str(Seurat::Misc(obj, "NNet.setting"))



###################  Check results ###################
pcs <- Seurat::Misc(obj, "NNet.setting")$pcs
lra <- Seurat::Misc(obj, "NNet.setting")$lra
obj <- Seurat::RunUMAP(obj, dims = 1:ncol(pcs))
umap <- data.frame(Seurat::Embeddings(obj, "umap"))

gene <- genes[4]
predictor <- rowSums(2*obj@misc$NNet.mod$effect[gene, , ]^2) %>%
  sort(decreasing = T) %>% names %>% dplyr::nth(1000)

Seurat::FeaturePlot(obj, features = c(gene, predictor))
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
new.obj <- prepare.reg(new.obj)

pcs <- Seurat::Misc(new.obj, "NNet.setting")$pcs
lra <- Seurat::Misc(new.obj, "NNet.setting")$lra
new.obj <- Seurat::RunUMAP(new.obj, dims = 1:ncol(pcs))
umap <- data.frame(Seurat::Embeddings(new.obj, "umap"))

gene <- "GAPDH"#genes[1]
new.obj <- run.nn.reg(new.obj, responses = gene, return.p.val = T)
predictor <-  gene

ggplot() + geom_jitter(aes(x= "umi", y= log(abs(new.obj@misc$NNet.mod$effect[gene,"n.umi",])))) +
  geom_jitter(aes(x= "predictor", y= log(abs(new.obj@misc$NNet.mod$effect[gene,predictor,])))) +
  geom_jitter(aes(x= "null", y= log(abs(new.obj@misc$NNet.mod$effect[gene,paste("NULL",predictor, sep = "-"),]))))


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


