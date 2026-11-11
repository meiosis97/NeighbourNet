# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal: draw a network view as an interactive visNetwork widget
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Used by visualise.network(..., interactive = TRUE). Draws the same genes,
# edges, layout and clusters as the static plot from .network.view(), with:
#   - shapes by role (g2 genes: diamonds; g1 genes: circles; receptors:
#     triangles; intermediate signalling steps: small squares)
#   - colour = g2 cluster only (g2 genes filled, other genes neutral)
#   - edge width = weight; arrow / T-bar only for prior-supported edges
#   - grey genes: no edge above the cutoff
#   - tooltips, click-to-highlight, gene and role menus, an HTML legend

# Cluster colours: the static plot's palette (grey = no cluster, then Set1).
.nnet.cluster.cols <- function() c("#BDBDBD", RColorBrewer::brewer.pal(9, "Set1"))

.render.interactive <- function(view, height = "650px"){
  for(pkg in c("visNetwork", "htmlwidgets", "htmltools")){
    if(!requireNamespace(pkg, quietly = TRUE)){
      stop("interactive = TRUE needs the '", pkg, "' package. Install it with ",
           "install.packages(\"", pkg, "\").", call. = FALSE)
    }
  }
  tabs <- .view.tables(view)
  nd <- tabs$nodes
  ed <- tabs$edges

  cols <- .nnet.cluster.cols()
  col.of <- function(m) ifelse(is.na(m), cols[1], cols[pmin(m + 1L, length(cols))])
  role.g2 <- tabs$role.g2
  module.word <- paste(role.g2, "cluster")

  # Spread rings so that nodes and labels do not collide in large networks:
  # aim for ~70 px between neighbouring node centres on the most crowded
  # ring (node `size` is a radius, up to 20 px, plus room for the label).
  r <- round(sqrt(nd$x^2 + nd$y^2), 6)
  ring <- tapply(r, r, length)
  ring.r <- as.numeric(names(ring))
  px <- max(320, max(ring[ring.r > 0] * 70 / (2 * pi * ring.r[ring.r > 0])))

  # Tooltip text
  strip <- function(x) sub("^g[0-9]\\.", "", x)
  partners <- vapply(nd$name, function(n){
    e <- ed[ed$from == n | ed$to == n, , drop = FALSE]
    if(!nrow(e)) return("<i>none above the cutoff</i>")
    other <- ifelse(e$from == n, e$to, e$from)
    paste0(nd$label[match(other, nd$name)], " (", sprintf("%.2f", e$weight), ")",
           collapse = ", ")
  }, character(1))
  dup <- nd$label %in% nd$label[duplicated(nd$label) & nd$label != ""]
  isolated <- nd$degree == 0
  weight.txt <- ifelse(
    nd$role == "receptor",
    paste0("relative receptor activity: ", sprintf("%.2f", nd$weight)),
    ifelse(nd$role == "pathway step", "",
           paste0("strongest edge significance (any partner): ", sprintf("%.2f", nd$weight))))
  hidden.partner <- nd$role %in% c("target", "TF") & isolated & nd$weight > view$cutoff

  nodes <- data.frame(
    id    = nd$name,
    label = nd$label,
    group = nd$role,
    shape = unname(c(TF = "diamond", target = "dot", receptor = "triangle",
                     `pathway step` = "square")[nd$role]),
    size  = ifelse(nd$role == "pathway step" | isolated, 7, 9 + 11 * nd$weight),
    x     = nd$x * px,
    y     = -nd$y * px,                     # vis.js y axis points down
    title = paste0(
      "<b>", nd$label, "</b> &middot; ", nd$role,
      ifelse(nd$role == role.g2 & !is.na(nd$module),
             paste0(" &middot; ", module.word, " ", nd$module), ""),
      ifelse(dup, "<br><i>also drawn in its other layer</i>", ""),
      ifelse(weight.txt == "", "", paste0("<br>", weight.txt)),
      ifelse(hidden.partner,
             "<br><i>its strongest partner is not among the genes drawn</i>", ""),
      "<br>partners drawn: ", partners),
    stringsAsFactors = FALSE
  )
  # Only dotted (nested) vis.js columns: mixing a plain `color` column with
  # `color.border` stops the widget from rendering.
  nodes$color.background <- ifelse(isolated, "#F2F2F2",
                                   ifelse(nd$role == role.g2, col.of(nd$module),
                                          ifelse(nd$role == "pathway step", "#DDDDDD", "#FFFFFF")))
  nodes$color.border <- ifelse(isolated, "#BDBDBD", "#333333")
  nodes$borderWidth  <- ifelse(isolated, 1, 2)
  nodes$font.color   <- ifelse(isolated, "#9E9E9E",
                               ifelse(nd$role == "pathway step", "#555555", "#222222"))
  nodes$font.size    <- ifelse(nd$role == "pathway step", 13, 17)

  edges <- data.frame(
    from  = ed$from,
    to    = ed$to,
    width = 1 + 6 * ed$weight,
    color = col.of(ed$module),
    title = paste0(
      "<b>", nd$label[match(ed$from, nd$name)], " &ndash; ",
      nd$label[match(ed$to, nd$name)], "</b> (", ed$kind, ")",
      "<br>edge weight: ", sprintf("%.2f", ed$weight),
      ifelse(is.na(ed$significance), "",
             paste0("<br>significance: ", sprintf("%.2f", ed$significance))),
      "<br>prior knowledge: ",
      c(none = "none", activation = "activation", inhibition = "inhibition")[ed$prior]),
    stringsAsFactors = FALSE
  )
  edges$arrows.to.enabled <- ed$prior != "none"
  edges$arrows.to.type    <- ifelse(ed$prior == "inhibition", "bar", "arrow")

  title <- if(view$meta.network) paste("Meta-network", view$i) else paste("Cell", view$i)
  sub <- sprintf(paste0("Edges with significance > %.2f. Hover for details; ",
                        "click a gene to highlight its partners."), view$cutoff)

  mods <- sort(unique(stats::na.omit(nd$module[nd$role == role.g2])))
  legend <- .interactive.legend(mods, col.of, role.g2 = role.g2,
                                roles = unique(nd$role))

  w <- visNetwork::visNetwork(nodes, edges, main = title, submain = sub,
                              height = height, width = "100%")
  w <- visNetwork::visPhysics(w, enabled = FALSE)
  w <- visNetwork::visEdges(w, smooth = FALSE, selectionWidth = 2, hoverWidth = 1)
  w <- visNetwork::visNodes(w, font = list(face = "arial"))
  w <- visNetwork::visOptions(
    w,
    highlightNearest = list(enabled = TRUE, degree = 1, hover = FALSE,
                            algorithm = "hierarchical"),
    nodesIdSelection = list(enabled = TRUE, main = "Select a gene",
                            values = nd$name[order(nd$label)]),
    selectedBy = list(variable = "group", main = "Select a layer"))
  w <- visNetwork::visInteraction(w, hover = TRUE, tooltipDelay = 80,
                                  navigationButtons = TRUE, dragNodes = TRUE)
  # Keep the requested height instead of filling the window, so the legend
  # below the network stays visible.
  w$sizingPolicy$browser$fill <- FALSE
  w$sizingPolicy$viewer$fill  <- FALSE
  htmlwidgets::appendContent(w, legend)
}

# Nodes and edges of a network view as data frames (drawn elements only).
.view.tables <- function(view){
  graph <- view$graph
  l <- view$layout
  vn <- names(igraph::V(graph))
  type <- igraph::V(graph)$type
  role.g1 <- if(identical(view$as.g2, "responses")) "TF" else "target"
  role.g2 <- if(identical(view$as.g2, "responses")) "target" else "TF"
  role <- c(role.g1, role.g2, "pathway step", "receptor", "spacer")[type]

  pie <- view$pie.values[vn, , drop = FALSE]
  slices <- setdiff(colnames(pie), c("x", "y", "scale", "cluster", "N"))
  weight <- 1 - pie[, "N"]
  module <- ifelse(role == role.g2, igraph::V(graph)$cluster - 1L, NA_integer_)

  nodes <- data.frame(name = vn, label = view$vertex.label, role = role,
                      module = module, weight = weight,
                      x = l[vn, "x"], y = l[vn, "y"], stringsAsFactors = FALSE)
  nodes <- nodes[nodes$role != "spacer", ]

  e <- igraph::as_data_frame(graph, "edges")
  e <- e[e$evidence >= 2 & e$from %in% nodes$name & e$to %in% nodes$name, ]
  rl <- stats::setNames(nodes$role, nodes$name)
  coexp <- rl[e$from] %in% c("TF", "target") & rl[e$to] %in% c("TF", "target")
  edges <- data.frame(
    from = e$from, to = e$to,
    kind = ifelse(coexp, "co-expression", "signalling path"),
    weight = e$weight,
    significance = NA_real_,
    prior = c("none", "activation", "inhibition")[e$evidence - 1],
    module = ifelse(e$cluster > 1, e$cluster - 1L, NA_integer_),
    stringsAsFactors = FALSE)
  # g1-g2 significance (prob.network is g1 x g2 with prefixed names)
  pn <- view$prob.network
  g1.end <- ifelse(startsWith(edges$from, "g1."), edges$from, edges$to)
  g2.end <- ifelse(startsWith(edges$from, "g2."), edges$from, edges$to)
  ok <- coexp & g1.end %in% rownames(pn) & g2.end %in% colnames(pn)
  edges$significance[ok] <- pn[cbind(g1.end[ok], g2.end[ok])]

  nodes$degree <- vapply(nodes$name, function(n) sum(edges$from == n | edges$to == n),
                         integer(1))
  list(nodes = nodes, edges = edges, role.g2 = role.g2)
}

# Legend as plain HTML below the widget (fixed, readable text size).
.interactive.legend <- function(mods, col.of, role.g2, roles, font.px = 15){
  svg <- function(body, w = 20){
    sprintf(paste0('<svg width="%d" height="20" viewBox="0 0 %d 20" ',
                   'style="vertical-align:middle;">%s</svg>'), w, w, body)
  }
  node.svg <- function(shape, fill, stroke = "#333333"){
    body <- switch(shape,
      diamond  = '<polygon points="10,2 18,10 10,18 2,10"',
      dot      = '<circle cx="10" cy="10" r="7"',
      triangle = '<polygon points="10,2 18,17 2,17"',
      square   = '<rect x="5" y="5" width="10" height="10"')
    svg(sprintf('%s fill="%s" stroke="%s" stroke-width="1.5"/>', body, fill, stroke))
  }
  edge.svg <- function(end = c("none", "arrow", "bar")){
    end <- match.arg(end)
    line <- '<line x1="2" y1="10" x2="34" y2="10" stroke="#555555" stroke-width="3"/>'
    tip <- switch(end,
      none  = "",
      arrow = '<polygon points="38,10 29,5 29,15" fill="#555555"/>',
      bar   = '<line x1="37" y1="3" x2="37" y2="17" stroke="#555555" stroke-width="3"/>')
    svg(paste0(line, tip), w = 40)
  }
  item <- function(icon, label){
    sprintf(paste0('<span style="display:inline-block;margin:0 16px 6px 0;',
                   'white-space:nowrap;">%s<span style="vertical-align:middle;',
                   'margin-left:6px;">%s</span></span>'), icon, label)
  }
  row <- function(title, items){
    sprintf(paste0('<div style="display:block;width:100%%;">',
                   '<span style="display:inline-block;font-weight:bold;',
                   'margin:0 12px 6px 0;vertical-align:middle;">%s</span>%s</div>'),
            title, paste(items, collapse = " "))
  }
  shape.of <- c(TF = "diamond", target = "dot")
  role.g1 <- setdiff(c("TF", "target"), role.g2)
  genes <- c(
    item(node.svg(shape.of[[role.g2]], "#FFFFFF"), paste0(role.g2, " (colour = cluster)")),
    item(node.svg(shape.of[[role.g1]], "#FFFFFF"), role.g1),
    if("receptor" %in% roles) item(node.svg("triangle", "#FFFFFF"), "receptor"),
    if("pathway step" %in% roles) item(node.svg("square", "#DDDDDD"), "signalling step"),
    vapply(mods, function(m) item(node.svg(shape.of[[role.g2]], col.of(m)),
                                  paste(role.g2, "cluster", m)), character(1)),
    item(node.svg("dot", "#F2F2F2", "#BDBDBD"), "no edge above cutoff")
  )
  edges <- c(
    item(edge.svg("none"), "co-expression (width = weight)"),
    item(edge.svg("arrow"), "+ prior knowledge: activation"),
    item(edge.svg("bar"), "+ prior knowledge: inhibition")
  )
  htmltools::HTML(sprintf(
    paste0('<div style="font-family:Arial,Helvetica,sans-serif;font-size:%dpx;',
           'color:#222222;padding:8px 4px 4px;">%s%s</div>'),
    font.px, row("Genes", genes), row("Edges", edges)))
}
