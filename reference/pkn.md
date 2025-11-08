# Integrated Weighted Directed Networks for Prior Knowledge in Gene Regulation and Signaling

Two `igraph` objects representing integrated weighted directed networks:
`gr.graph` describes transcription factor (TF) and target gene
regulatory relationships, and `sig.graph` describes intra-cellular
signaling interactions.

## Usage

``` r
gr.graph

sig.graph
```

## Format

Two `igraph` objects with directed and weighted edges, containing the
following edge attributes:

- `weight`:

  A numeric value representing the strength of the interaction.

- `consensus_stimulation`:

  A logical attribute indicating whether the interaction is consistently
  recorded as stimulating across multiple databases in OmniPath.

- `consensus_inhibition`:

  A logical attribute indicating whether the interaction is consistently
  recorded as inhibiting across multiple databases in OmniPath.

An object of class `igraph` of length 9867.

## Source

OmniPath database via the OmniPathR package.

## Details

These prior knowledge networks (PKN) are generated using the OmniPath R
package and provide high-confidence regulatory and signaling
interactions by integrating data from multiple databases. They are
crucial for gene regulatory network (GRN) inference and upstream
signaling pathway (USP) analysis in NNet.

## Examples

``` r
# Summary of the gene regulatory network
summary(gr.graph)
#> IGRAPH 63bbc63 DNW- 6855 41206 -- 
#> + attr: name (v/c), weight (e/n), consensus_direction (e/n),
#> | consensus_stimulation (e/n), consensus_inhibition (e/n)

# Summary of the signaling network
summary(sig.graph)
#> IGRAPH 2d729a0 DNW- 9867 79242 -- 
#> + attr: name (v/c), weight (e/n), consensus_direction (e/n),
#> | consensus_stimulation (e/n), consensus_inhibition (e/n)

# Access edge attributes in gr.graph
edge_attr(gr.graph)
#> Error in edge_attr(gr.graph): could not find function "edge_attr"

# Access edge attributes in sig.graph
edge_attr(sig.graph)
#> Error in edge_attr(sig.graph): could not find function "edge_attr"
```
