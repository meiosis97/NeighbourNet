# NeighbourNet 1.1.0

## New features

* `visualise.network()` gains `interactive = TRUE`, which returns an
  interactive network (a visNetwork htmlwidget) instead of a ggplot. It draws
  the same genes, edges, clusters and layout as the static plot, and adds:
  - hover tooltips with each gene's layer, cluster, edge significance and the
    partners drawn, and each edge's weight, significance and prior evidence;
  - click a gene, or pick it from the menu, to highlight its partners; filter
    by layer; zoom and pan;
  - distinct shapes per layer: `g2` genes as diamonds filled with their
    cluster colour, `g1` genes as circles, receptors as triangles and
    intermediate signalling steps as small labelled squares;
  - plain lines for co-expression only, with an arrow (activation) or a bar
    (inhibition) when prior knowledge supports the edge;
  - greyed-out genes for those with no edge above `cutoff`;
  - a legend below the network.

  It needs the visNetwork, htmlwidgets and htmltools packages, which are now
  suggested dependencies. The static plot (`interactive = FALSE`, the
  default) is unchanged.

## Internal changes

* `visualise.network()` is split into a computation step and a drawing step,
  so that the static and interactive plots share one computation.
* First unit tests (testthat), covering `visualise.network()` on a small
  simulated data set.
