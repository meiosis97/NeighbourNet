# Gene List from Integrated Prior Knowledge Network

A list of gene sets derived from an integrated prior knowledge network
used in the NeighbourNet (NNet) framework. This list includes
transcription factors (TFs), target genes, receptors, and ligands.

## Usage

``` r
gene.list
```

## Format

A list with four entries:

- `tfs`:

  A character vector of transcription factor (TF) gene names. These are
  regulators that control gene expression.

- `targets`:

  A character vector of target gene names. These genes are regulated by
  TFs.

- `receptors`:

  A character vector of receptor gene names. These genes encode proteins
  that receive extracellular signals.

- `ligands`:

  A character vector of ligand gene names. These genes encode signaling
  molecules that bind to receptors.

## Source

The integrated prior knowledge network used in NNet analysis.

## Details

This list is crucial for gene regulatory network (GRN) inference and
upstream signaling pathway (USP) analysis in NNet.
