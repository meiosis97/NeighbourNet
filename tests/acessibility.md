## 1 · One‑liner preprocessing

```r
require(NeighbourNet)
require(Seurat)
rt.ppr <- get.ppr()                        # receptor‑target prior matrix
genes  <- select.gene(obj, min.cells = 10) # QC → TF / target lists

# Obj: A Seurat Object
obj <- obj |>
  prepare.seurat(genes = genes$genes) |>   # scale + PCA
  prepare.graph() |>                       # 30‑NN graph
  select.cell() |>                         # subsample 
  prepare.reg(predictors = genes$tfs,      # local variance scaffolding
              responses  = genes$targets)
```

## 2 · NNet regression
```r
top10 <- head(genes$targets, 10)           # demo: first 10 targets
obj   <- run.nn.reg(obj, responses = top10, return.p.val = TRUE) |>
         build.meta.network() |>
         select.central.genes() |>
         prepare.visualise()
```

## 3 · Snapshot plot (Cell #1)
```r
visualise.network(obj, 1)
```

## 4 · Snapshot plot (meta‑network #1)
```r
visualise.network(obj, 1, meta.network = TRUE)
```

## 5. Receptor activity
```r
act  <- receptor.activity(obj)             # matrix: receptor × cell
```
