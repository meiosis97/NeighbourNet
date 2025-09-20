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

```mermaid
flowchart TD

    A[1 · One-liner preprocessing<br/><br/>require(NeighbourNet)<br/>rt.ppr <- get.ppr()<br/>genes <- select.gene(obj, min.cells = 10)<br/>obj <- obj |>
      prepare.seurat(genes = genes$genes) |>
      prepare.graph() |>
      select.cell() |>
      prepare.reg(predictors = genes$tfs, responses = genes$targets)] 

    B[2 · NNet regression<br/><br/>top10 <- head(genes$targets, 10)<br/>obj <- run.nn.reg(obj, responses = top10, return.p.val = TRUE) |>
      build.meta.network() |>
      select.central.genes() |>
      prepare.visualise()]

    C[3 · Snapshot plot (Cell #1)<br/><br/>visualise.network(obj, 1)]

    D[4 · Snapshot plot (meta-network #1)<br/><br/>visualise.network(obj, 1, meta.network = TRUE)]

    E[5 · Receptor activity<br/><br/>act <- receptor.activity(obj)]

    A --> B --> C --> D --> E
```
