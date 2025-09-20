## 1 · One‑liner preprocessing

```r
require(NeighbourNet)
rt.ppr <- get.ppr()                        # receptor‑target prior matrix
genes  <- select.gene(obj, min.cells = 10) # QC → TF / target lists

obj <- obj |>
  prepare.seurat(genes = genes$genes) |>   # scale + PCA
  prepare.graph() |>                       # 30‑NN graph
  select.cell() |>                         # subsample 
  prepare.reg(predictors = genes$tfs,      # local variance scaffolding
              responses  = genes$targets)
```

