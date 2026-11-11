skip_if_no_widgets <- function() {
  skip_if_not_installed("visNetwork")
  skip_if_not_installed("htmlwidgets")
  skip_if_not_installed("htmltools")
}

# Static plot and interactive widget for the same arguments.
both_views <- function(obj, ...) {
  p <- suppressWarnings(visualise.network(obj, ...))
  w <- suppressWarnings(visualise.network(obj, ..., interactive = TRUE))
  list(p = p, w = w, graph = attr(p$data, "graph"))
}

drawn_edges <- function(graph) {
  e <- igraph::as_data_frame(graph, "edges")
  e[e$evidence >= 2, ]
}

test_that("the default is the static ggplot", {
  obj <- nnet_fixture()
  p <- suppressWarnings(visualise.network(obj, 1, meta.network = TRUE,
                                          cutoff = meta_cutoff(obj)))
  expect_s3_class(p, "ggplot")
})

test_that("interactive view draws the same genes and edges as the static plot", {
  skip_if_no_widgets()
  obj <- nnet_fixture()
  cases <- list(
    meta  = list(i = 1, meta.network = TRUE, cutoff = meta_cutoff(obj)),
    nopath = list(i = 1, meta.network = TRUE, cutoff = meta_cutoff(obj),
                  show.pathways = FALSE),
    cell  = list(i = 1, meta.network = FALSE, cutoff = 0.5)
  )
  for (nm in names(cases)) {
    v <- do.call(both_views, c(list(obj), cases[[nm]]))
    expect_s3_class(v$w, "visNetwork")
    expect_s3_class(v$w, "htmlwidget")

    types <- igraph::V(v$graph)$type
    expect_setequal(v$w$x$nodes$id, igraph::V(v$graph)$name[types != 5])  # spacers dropped

    e <- drawn_edges(v$graph)
    expect_setequal(paste(v$w$x$edges$from, v$w$x$edges$to), paste(e$from, e$to))
  }
})

test_that("node shapes follow the layers, including receptor pathways", {
  skip_if_no_widgets()
  obj <- nnet_fixture()
  v <- both_views(obj, 1, meta.network = TRUE, cutoff = meta_cutoff(obj))
  nodes <- v$w$x$nodes
  # the fixture is built so that meta-network 1 has receptor pathways
  expect_true(all(c("TF", "target", "receptor", "pathway step") %in% nodes$group))

  shape <- tapply(nodes$shape, nodes$group, unique)
  expect_equal(shape[["TF"]], "diamond")
  expect_equal(shape[["target"]], "dot")
  expect_equal(shape[["receptor"]], "triangle")
  expect_equal(shape[["pathway step"]], "square")

  # TFs (the g2 layer) carry their cluster colour; targets are neutral
  drawn <- c(v$w$x$edges$from, v$w$x$edges$to)
  tf <- nodes$group == "TF" & nodes$id %in% drawn
  expect_false(any(nodes$color.background[tf] %in% c("#FFFFFF", "#F2F2F2")))
  expect_true(all(nodes$color.background[nodes$group == "target" & nodes$id %in% drawn] == "#FFFFFF"))
})

test_that("edges mark prior activation with an arrow and inhibition with a bar", {
  skip_if_no_widgets()
  obj <- nnet_fixture()
  v <- both_views(obj, 1, meta.network = TRUE, cutoff = meta_cutoff(obj))
  e <- drawn_edges(v$graph)
  ed <- v$w$x$edges
  key <- function(from, to) paste(from, to)
  ev <- e$evidence[match(key(ed$from, ed$to), key(e$from, e$to))]

  expect_true(all(c(2, 3, 4) %in% ev))                       # fixture covers all three
  expect_true(all(!ed$arrows.to.enabled[ev == 2]))          # co-expression only
  expect_true(all(ed$arrows.to.enabled[ev %in% c(3, 4)]))
  expect_true(all(ed$arrows.to.type[ev == 3] == "arrow"))
  expect_true(all(ed$arrows.to.type[ev == 4] == "bar"))
  expect_equal(ed$width, 1 + 6 * e$weight[match(key(ed$from, ed$to), key(e$from, e$to))])
})

test_that("genes without a drawn edge are greyed out", {
  skip_if_no_widgets()
  obj <- nnet_fixture()
  # cell 1 has no edge above a 0.95 cutoff in the fixture
  w <- suppressWarnings(visualise.network(obj, 1, meta.network = FALSE,
                                          cutoff = 0.95, interactive = TRUE))
  expect_equal(nrow(w$x$edges), 0)
  expect_true(all(w$x$nodes$color.background == "#F2F2F2"))
})

test_that("the requested index is reported even when pathways reuse `i`", {
  skip_if_no_widgets()
  obj <- nnet_fixture()
  w <- suppressWarnings(visualise.network(obj, 2, meta.network = TRUE,
                                          cutoff = meta_cutoff(obj, 2),
                                          interactive = TRUE))
  expect_equal(w$x$main$text, "Meta-network 2")
})

test_that("widget data uses only nested vis.js colour and arrow options", {
  # A plain `color` / `arrows` column next to `color.border` /
  # `arrows.to.type` makes vis.js fail to render the widget.
  skip_if_no_widgets()
  obj <- nnet_fixture()
  w <- suppressWarnings(visualise.network(obj, 1, meta.network = TRUE,
                                          cutoff = meta_cutoff(obj),
                                          interactive = TRUE))
  expect_false("color" %in% names(w$x$nodes))
  expect_false("arrows" %in% names(w$x$edges))
})

test_that("the legend is added below the network", {
  skip_if_no_widgets()
  obj <- nnet_fixture()
  w <- suppressWarnings(visualise.network(obj, 1, meta.network = TRUE,
                                          cutoff = meta_cutoff(obj),
                                          interactive = TRUE))
  legend <- paste(vapply(w$append, as.character, character(1)), collapse = "")
  expect_match(legend, "TF cluster")
  expect_match(legend, "receptor")
  expect_match(legend, "prior knowledge: inhibition")
})

test_that("a missing visNetwork gives an informative error", {
  obj <- nnet_fixture()
  local_mocked_bindings(.has.package = function(pkg) pkg != "visNetwork")
  expect_error(
    suppressWarnings(  # igraph deprecation warnings from the shared view code
      visualise.network(obj, 1, meta.network = TRUE, cutoff = meta_cutoff(obj),
                        interactive = TRUE)),
    "install.packages\\(\"visNetwork\"\\)"
  )
})
