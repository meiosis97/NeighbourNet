## Introduction

This vignette demonstrates the usage of the NeighbourNet by sourcing the
provided R script into the global environment.

## Install dependencies

``` r
uninstalled.packages <- c("Seurat", "dplyr", "Matrix", "ggplot2", "ggraph", "scatterpie", "ggrepel", "ggpubr", "igraph")
uninstalled.packages <- uninstalled.packages[!uninstalled.packages%in%installed.packages()]
if(length(uninstalled.packages)){
  install.packages(uninstalled.packages)
}
```

## Source and read data

Then, load in the pre-computed prior knowledge networks that are
available at our
[github](https://github.com/meiosis97/NeighbourNet/tree/main/data).

``` r
load("../data/gene.list.rda")
load("../data/sig.graph.rda")
load("../data/gr.graph.rda")
load("../data/receptor.ppr.rda")
```

Alternatively, you can also download the prior knowledge networks
complied by NicheNet, which contain more but less confident interactions

``` r
# Not run
gene.list <-  readRDS(url("https://zenodo.org/records/15240485/files/gene.list.rds"))
sig.graph <-  readRDS(url("https://zenodo.org/records/15240485/files/sig.graph.rds"))
gr.graph <-  readRDS(url("https://zenodo.org/records/15240485/files/gr.graph.rds"))
receptor.ppr <-  readRDS(url("https://zenodo.org/records/15240485/files/receptor.ppr.rds"))
```

Read the demonstrating data, which should be a Seurat object.

``` r
load("data/luad.rda")
```

## Preprocess

Prepare receptor prior regulatory potential matrix on targets.

``` r
rt.ppr <- get.ppr()
```

Quality control to filter out lowly expressed genes. Remaining genes are
grouped by their annotations as transcriptional factors (TF), targets,
receptors, and ligands.

``` r
genes <- select.gene(obj, min.cells = 10)
```

Run PCA using genes selected. Co-expression will be measured within
these genes.

``` r
obj <- prepare.seurat(obj, genes = genes$genes)
```

    ## Run Seurat.

Construct a KNN graph

``` r
obj <- prepare.graph(obj)
```

    ## Now building knn graph.

Optionally, for data containing many cells (\>5000), run analysis on a
subset of representative cells.

``` r
obj <- select.cell(obj)
```

Finally, calculate local variance.

``` r
# Local variance will be calculated for predictor genes and response genes
# Reponse genes will be low rank approximated
# We want to measure co-expression between TFs and targets
obj <- prepare.reg(obj, predictors = genes$tfs, responses = genes$targets)
```

    ## Calculating local variance.

Now, the prerequisite settings of NeighbourNet regression are stored in
the `misc` slot of the Seurat object.

``` r
Misc(obj, "setting") %>% str
```

    ## List of 14
    ##  $ pcs          : num [1:3918, 1:49] 4.15 1.61 -18.67 4.71 1.72 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:3918] "Lib90_00000" "Lib90_00001" "Lib90_00002" "Lib90_00003" ...
    ##   .. ..$ : chr [1:49] "PC_1" "PC_2" "PC_3" "PC_4" ...
    ##  $ loadings     : num [1:3995, 1:49] 0.00286 -0.00183 0.03787 -0.00224 -0.00568 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:3995] "URI1" "TADA2A" "NUPR1" "MYRF" ...
    ##   .. ..$ : chr [1:49] "PC_1" "PC_2" "PC_3" "PC_4" ...
    ##  $ predictors   : chr [1:640] "URI1" "TADA2A" "NUPR1" "MYRF" ...
    ##  $ responses    : chr [1:3851] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##  $ cells        : Named num [1:773] 3716 2610 2635 3043 3430 ...
    ##   ..- attr(*, "names")= chr [1:773] "Lib90_03759" "Lib90_02611" "Lib90_02637" "Lib90_03052" ...
    ##  $ p            :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##   .. ..@ i       : int [1:175602] 0 1 4 5 19 22 38 49 61 64 ...
    ##   .. ..@ p       : int [1:3919] 0 54 111 147 255 339 376 408 465 501 ...
    ##   .. ..@ Dim     : int [1:2] 3918 3918
    ##   .. ..@ Dimnames:List of 2
    ##   .. .. ..$ : NULL
    ##   .. .. ..$ : NULL
    ##   .. ..@ x       : num [1:175602] 0.52867 0.0107 0.00334 0.02509 0.00484 ...
    ##   .. ..@ factors : list()
    ##  $ nn.idx       : int [1:3918, 1:30] 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ nn.w         : num [1:3918, 1:30] 0.102 0.102 0.102 0.102 0.102 ...
    ##  $ lra          : num [1:3918, 1:3851] -0.2297 -0.0836 0.0749 -0.2167 -0.2838 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:3918] "Lib90_00000" "Lib90_00001" "Lib90_00002" "Lib90_00003" ...
    ##   .. ..$ : chr [1:3851] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##  $ scale.gene   : Named num [1:3995] 0.1135 0.2385 0.5815 0.3113 0.0858 ...
    ##   ..- attr(*, "names")= chr [1:3995] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##  $ nn.scale.gene:Formal class 'dgeMatrix' [package "Matrix"] with 4 slots
    ##   .. ..@ Dim     : int [1:2] 3995 773
    ##   .. ..@ Dimnames:List of 2
    ##   .. .. ..$ : chr [1:3995] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##   .. .. ..$ : chr [1:773] "Lib90_03759" "Lib90_02611" "Lib90_02637" "Lib90_03052" ...
    ##   .. ..@ x       : num [1:3088135] 0.15 0.222 0.644 0.337 0 ...
    ##   .. ..@ factors : list()
    ##  $ nn.scale.pc  : num [1:49, 1:773] 1.8 2.06 1.91 4.27 4.11 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:49] "PC_1" "PC_2" "PC_3" "PC_4" ...
    ##   .. ..$ : chr [1:773] "Lib90_03759" "Lib90_02611" "Lib90_02637" "Lib90_03052" ...
    ##  $ n.eff        : Named num [1:773] 17.2 21.5 22.2 22.7 21.4 ...
    ##   ..- attr(*, "names")= chr [1:773] "Lib90_03759" "Lib90_02611" "Lib90_02637" "Lib90_03052" ...
    ##  $ genes        : chr [1:3995] "COL27A1" "POLR3D" "SMC4" "URI1" ...

Get cells on which co-expression will be measured.

``` r
cells <- Misc(obj, "setting")$cells
```

## NeighbourNet regression

For the demonstration purpose, we only use ten response genes to
construct co-expression networks.

``` r
(responses <- genes$targets[1:10])
```

    ##  [1] "COL27A1" "POLR3D"  "SMC4"    "URI1"    "BAAT"    "TPM1"    "GAA"    
    ##  [8] "DMPK"    "SV2A"    "EGFR"

NeighbourNet regression.

``` r
obj <- run.nn.reg(obj, responses = responses, return.p.val = T)
```

    ## Return smoothed effect, can only generate networks for sampled cells.

    ## Return unpruned effect.

    ## Return p-value.

    ## By default, downstream analysis will be performed on the effect tensor.

    ## By default, downstream analysis will perform network prunning.

    ## Build the Laplacian operator.

    ## Now regress.

    ## COL27A1

    ## POLR3D

    ## SMC4

    ## URI1

    ## BAAT

    ## TPM1

    ## GAA

    ## DMPK

    ## SV2A

    ## EGFR

The resulting network ensemble can be found also in the `misc` slot
under the name `mod`.

``` r
# effect: the (response x predictor x cell) ensemble of co-expression networks.
# p.val: the ensemble of co-expression significance corresponds to `effect`
Misc(obj, "mod") %>% str
```

    ## List of 14
    ##  $ effect       : num [1:10, 1:649, 1:773] 1.50e-03 -5.00e-04 4.81e-03 1.46e-03 1.43e-08 ...
    ##   ..- attr(*, "dimnames")=List of 3
    ##   .. ..$ : chr [1:10] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##   .. ..$ : chr [1:649] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##   .. ..$ : chr [1:773] "Lib90_03759" "Lib90_02611" "Lib90_02637" "Lib90_03052" ...
    ##  $ p.val        : num [1:10, 1:649, 1:773] 0.9459 0.3937 0.8968 0.6125 0.0017 ...
    ##   ..- attr(*, "dimnames")=List of 3
    ##   .. ..$ : chr [1:10] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##   .. ..$ : chr [1:649] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##   .. ..$ : chr [1:773] "Lib90_03759" "Lib90_02611" "Lib90_02637" "Lib90_03052" ...
    ##  $ meta.network : NULL
    ##  $ mus          : Named num [1:10] -8.68 -7.28 -6.76 -6.85 -11.32 ...
    ##   ..- attr(*, "names")= chr [1:10] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##  $ sigmas       : Named num [1:10] 1.36 1.18 1.13 1.14 2.3 ...
    ##   ..- attr(*, "names")= chr [1:10] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##  $ subsampled   : logi TRUE
    ##  $ smoothed     : logi TRUE
    ##  $ pruned       : logi FALSE
    ##  $ gene.sets    :List of 3
    ##   ..$ predictors:List of 3
    ##   .. ..$ genes  : chr [1:640] "URI1" "TADA2A" "NUPR1" "MYRF" ...
    ##   .. ..$ tfs    : chr [1:640] "URI1" "TADA2A" "NUPR1" "MYRF" ...
    ##   .. ..$ targets: chr [1:496] "URI1" "NUPR1" "TSG101" "HIC1" ...
    ##   ..$ responses :List of 3
    ##   .. ..$ genes  : chr [1:10] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##   .. ..$ tfs    : chr "URI1"
    ##   .. ..$ targets: chr [1:10] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##   ..$ genes     : chr [1:649] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##  $ cells        : Named num [1:773] 3716 2610 2635 3043 3430 ...
    ##   ..- attr(*, "names")= chr [1:773] "Lib90_03759" "Lib90_02611" "Lib90_02637" "Lib90_03052" ...
    ##  $ defaults     :List of 7
    ##   ..$ f                :function (x)  
    ##   ..$ remove.self.loops: logi TRUE
    ##   ..$ assay            : chr "effect"
    ##   ..$ predictors       : chr [1:640] "URI1" "TADA2A" "NUPR1" "MYRF" ...
    ##   ..$ responses        : chr [1:10] "COL27A1" "POLR3D" "SMC4" "URI1" ...
    ##   ..$ cutoff           : num 0.95
    ##   ..$ prune            : logi TRUE
    ##  $ custom.y     : logi FALSE
    ##  $ w            :List of 2
    ##   ..$ u : num [1:3918, 1:30] -0.0609 -0.0677 -0.0707 -0.0753 -0.0748 ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   .. .. ..$ : chr [1:3918] "Lib90_00000" "Lib90_00001" "Lib90_00002" "Lib90_00003" ...
    ##   .. .. ..$ : NULL
    ##   ..$ vd: num [1:3918, 1:30] -0.0166 -0.0191 -0.016 -0.025 -0.022 ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   .. .. ..$ : chr [1:3918] "Lib90_00000" "Lib90_00001" "Lib90_00002" "Lib90_00003" ...
    ##   .. .. ..$ : NULL
    ##  $ return.smooth: logi TRUE

Build meta-networks.

``` r
obj <- build.meta.network(obj)
```

    ## Now construct the covariance matrix.

    ## COL27A1

    ## POLR3D

    ## SMC4

    ## URI1

    ## BAAT

    ## TPM1

    ## GAA

    ## DMPK

    ## SV2A

    ## EGFR

    ## Eigen decomposition.

    ## Non-negative PCA.

Meta-networks are also stored as a (response x predictor x meta-cell)
network ensemble.

``` r
Misc(obj, "mod")$meta.network$memeta.network %>% dim
```

## Visualisation

First, we identify the most connected genes to visualise.

``` r
# All responses will be visualised by `keep.responses = T`
# `n.net` controls the number of meta-networks to examine
# For each meta-network, eigenvector centrality of `k` sub-networks will be calculated
# On each eigenvector, `n.per.component` central genes will be selected to visualise
central.genes  <- select.central.genes(obj, keep.responses = T, n.per.component = 4, n.net = 5, k = 2)
```

Prepare visualisation parameters.

``` r
# TFs will be clustered into `n.clu` groups represented by different colours
obj <- prepare.visualise(obj, central.genes = central.genes, n.clu = 4)
```

Visualise the network of cell `i`.

``` r
i <- cells[10]
visualise.network(obj,i, cutoff = 0.95,
                  radius = c(0.4,0.7,0.85,1),pie.radius = 0.04, text.size = 5)
```

    ## Warning: `get.edge.ids()` was deprecated in igraph 2.1.0.
    ## ℹ Please use `get_edge_ids()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: `permute.vertices()` was deprecated in igraph 2.0.0.
    ## ℹ Please use `permute()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: `get.edgelist()` was deprecated in igraph 2.0.0.
    ## ℹ Please use `as_edgelist()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](vignettes.for.script_files/figure-markdown_github/unnamed-chunk-21-1.png)

Visualise meta-network `i`.

``` r
i <- 1
cutoff <- apply(obj@misc$mod$meta.network$p.val[,,i],1,max) %>% sort(decreasing = T) %>% mean
visualise.network(obj,i, meta.network = T, cutoff = cutoff,
                  radius = c(0.4,0.7,0.85,1),pie.radius = 0.04, text.size = 5)
```

    ## Warning: Removed 5 rows containing non-finite outside the scale range
    ## (`stat_pie()`).

![](vignettes.for.script_files/figure-markdown_github/unnamed-chunk-22-1.png)

Visualise which cells represent meta-network `i`.

``` r
i <- 1
weights <- obj@misc$mod$meta.network$npca.loadings[,i]
ggplot() + geom_point(data = Embeddings(obj,"pca"), aes(PC_1,PC_2)) + 
  geom_point(data = Embeddings(obj,"pca")[cells,], aes(PC_1,PC_2, col = weights))
```

![](vignettes.for.script_files/figure-markdown_github/unnamed-chunk-23-1.png)
\## Receptor activity

Calculate receptor activity for each cell.

``` r
act <- receptor.activity(obj)
```

Visualise the activity distribution of NOTCH2, and identify TFs with
highly correlated activities.

``` r
notch2.act <- act$receptor.act["NOTCH2",] %>% as.numeric()
ggplot() + geom_point(data = Embeddings(obj,"pca"), aes(PC_1,PC_2)) + 
  geom_point(data = Embeddings(obj,"pca")[cells,], aes(PC_1,PC_2, col = notch2.act))
```

![](vignettes.for.script_files/figure-markdown_github/unnamed-chunk-25-1.png)

``` r
cor(notch2.act, t(act$tf.act)) %>% drop %>% sort 
```

    ## Warning in cor(notch2.act, t(act$tf.act)): the standard deviation is zero

    ##        SMAD3      TGFB1I1        FOXA2         E2F1         SFPQ        MYBL2 
    ## -0.176631011 -0.151598307 -0.142649139 -0.141363295 -0.115478842 -0.110618144 
    ##          ID2         HES6        NR1H2        EPAS1          RB1          YY1 
    ## -0.109707241 -0.107903916 -0.106747436 -0.099787304 -0.098437222 -0.094429587 
    ##        ESRRA        SNAI1          REL       SREBF1        KDM5B       STAT5B 
    ## -0.094287161 -0.090010246 -0.080805812 -0.075681691 -0.074109859 -0.072018809 
    ##        MECOM         MAFB          VDR        GATA2         KAT5         DLX2 
    ## -0.070412510 -0.067046197 -0.066752438 -0.065536371 -0.064278329 -0.063631553 
    ##        EWSR1        NR3C1          AHR        CLOCK         MSX2        PPARG 
    ## -0.061142385 -0.060905000 -0.057856843 -0.056576992 -0.055576199 -0.053856641 
    ##        PRDM2         RXRA         USF2        KDM5C          FOS          SP3 
    ## -0.048586971 -0.047296889 -0.045290140 -0.040405432 -0.039984728 -0.036666775 
    ##         E2F4         MTA1         TP53         JUND         BCL3         PBX1 
    ## -0.036461158 -0.036002748 -0.031346192 -0.029533063 -0.027844945 -0.025636753 
    ##        NFKB1        NFIL3        HDAC9      SMARCC1          MAZ        NCOA3 
    ## -0.020725119 -0.017500058 -0.008626192 -0.008620578 -0.005065852 -0.000958519 
    ##          JUN        CREB1          MYC        FOXA1         ELF3         RARA 
    ##  0.008252310  0.029564825  0.047879965  0.064549498  0.065614890  0.065621548 
    ##         HSF1         EGR1         CUX1        FOXM1         IRF2         PURA 
    ##  0.065623569  0.066819047  0.111378831  0.137606275  0.139241206  0.211715862 
    ##         YBX3       BCLAF1         ELK3        MECP2        KMT2A        DNMT1 
    ##  0.220361815  0.373525646  0.393305993  0.407475118  0.407569410  0.444402653 
    ##         NRG1        MEF2C         ENO1       RBFOX2          AIP         NONO 
    ##  0.444696974  0.505091913  0.513123663  0.546218782  0.584553719  0.601727195 
    ##         NFYC        HDAC1       NFE2L2         FHL2        FOSL1         CTCF 
    ##  0.606599318  0.614943742  0.634727759  0.643698695  0.660142413  0.679738940 
    ##        CTBP1        STAT1        MYBL1      SMARCA4         MED1         ETS1 
    ##  0.742598770  0.753939140  0.755040006  0.776785212  0.795550179  0.806757520 
    ##        HIF1A         MITF        TFDP1         MUC1        DMAP1        KLF10 
    ##  0.837509883  0.842854012  0.862571851  0.910622094  0.912416875  0.915882072 
    ##      LRRFIP1       HOXA10        NR2F2         SOX2         RBPJ 
    ##  0.919312798  0.945587426  0.984347510  0.985741533  0.988629735

We identified potential E2F1 mediated NOTCH inhibition of URI1.
