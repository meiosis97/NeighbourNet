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
visualise.network(obj, 1, 
                  radius = c(.4,.7,.85,1), pie.radius = .04,
                  text.size = 5)
```

## 4 · Snapshot plot (meta‑network #1)
```r
cut <- mean(apply(obj@misc$mod$meta.network$p.val[,,1], 1, max))
visualise.network(obj, 1, meta.network = TRUE, cutoff = cut,
                  radius = c(.4,.7,.85,1), pie.radius = .04,
                  text.size = 5)
```

## 5. Receptor activity
```r
act  <- receptor.activity(obj)             # matrix: receptor × cell
```
