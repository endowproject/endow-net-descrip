## ============================================================================
## 08_network_visualizations.R
## ----------------------------------------------------------------------------
## Network visualizations.
##
## Draws each community's network layers with ggraph, one PDF per network
## type (sites as subplots) and one PDF per site (network types as
## subplots), plus a few example aggregate-network illustrations.
##
## Outputs:
##   Figures/Networks/*_all_sites.pdf
##   Figures/Networks/*_all_nets.pdf
##   Figures/Networks/eg_agg_nets.pdf
##   Figures/Networks/eg_agg_nets_wealth.pdf
##
## Depends on: 00_setup.R
## ============================================================================

library(tidygraph)
library(ggraph)

source("net_descriptives/00_setup.R")

# ── Network Plots ────────────────────────────────────────────────
plot_net <- function(
  net,
  title,
  layout_style = "fr",
  edge_weight = FALSE,
  node_size_attr = NULL,
  node_colour = "royalblue"
) {
  if (ecount(net) == 0) edge_weight <- FALSE

  edge_width <- if (edge_weight) aes(width = weight) else aes()
  arr        <- if (is_directed(net)) arrow(length = unit(0.15, "cm")) else NULL
  graph      <- tidygraph::as_tbl_graph(net)

  if (!is.null(node_size_attr)) {
    vals <- igraph::vertex_attr(graph, node_size_attr)
    was_missing <- is.na(vals)
    if (any(was_missing)) {
      vals[was_missing] <- median(vals, na.rm = TRUE)
    }
    graph <- graph %>%
      tidygraph::activate(nodes) %>%
      tidygraph::mutate(
        !!node_size_attr := vals,
        was_missing       = was_missing
      )
  }

  p <- ggraph(graph, layout = layout_style) +
    geom_edge_loop(
      aes(span = 50, direction = 45, strength = 0.5),
      colour      = adjustcolor("black", alpha.f = 0.25),
      show.legend = FALSE
    ) +
    geom_edge_fan(
      edge_width,
      colour      = adjustcolor("black", alpha.f = 0.25),
      arrow       = arr,
      end_cap     = circle(1, "mm"),
      show.legend = FALSE
    ) +
    scale_edge_width_continuous(range = c(0.3, 2))

  p <- if (!is.null(node_size_attr)) {
    p +
      geom_node_point(aes(size = .data[[node_size_attr]], colour = was_missing),
                       show.legend = FALSE) +
      scale_size_continuous(range = c(0.5, 5), guide = "none") +
      scale_colour_manual(values = c(`FALSE` = node_colour, `TRUE` = "grey"))
  } else {
    p + geom_node_point(size = 1.5, colour = node_colour)
  }

  p +
    ggtitle(paste(strwrap(title, width = 20), collapse = "\n")) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title      = element_text(size = 9, hjust = 0.5)
    )
}

blank_plot <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste(strwrap(label, width = 20), collapse = "\n"),
             colour = "grey50", size = 3) +
    theme_void()
}

pad_to_full <- function(page, n) {
  while (length(page) < n) page <- c(page, list(ggplot() + theme_void()))
  page
}

write_net_pdf <- function(plot_list, path, plots_per_page, grid_rows, grid_cols, width, height) {
  pdf(path, width = width, height = height)
  for (page in split(plot_list, ceiling(seq_along(plot_list) / plots_per_page))) {
    print(wrap_plots(pad_to_full(page, plots_per_page), nrow = grid_rows, ncol = grid_cols))
  }
  dev.off()
}

# ---- One PDF per network type, sites as subplots ----
for (net_type in tie_meta$short[1:12]) {
  plot_list <- lapply(site_codes, function(site) {
    net <- su_nets_expanded[[site]][[net_type]]
    if (!is.null(net)) {
      if (net_type != "primary_kin") {
        plot_net(net, site, edge_weight = TRUE, node_colour = short_to_colour[net_type])
      } else {
        plot_net(net, site, edge_weight = FALSE, node_colour = short_to_colour[net_type])
      }
    } else {
      blank_plot(paste0(site, " (not asked)"))
    }
  })

  write_net_pdf(
    plot_list,
    path           = file.path(EndowDropbox, "Figures", "Networks",
                                paste0(net_type, "_all_sites.pdf")),
    plots_per_page = 16,
    grid_rows      = 4,
    grid_cols      = 4,
    width          = 10,
    height         = 10
  )
}

# ---- One PDF per site, network types as subplots ----
for (site in site_codes) {
  if(site != "MK") {
    plot_list <- lapply(c(1,3,5,7,9,2,4,6,8,10), function(i) {
      net_type <- tie_meta$short[i]
      net_label <- tie_meta$label[i]
      net <- su_nets_expanded[[site]][[net_type]]
      if (!is.null(net)) {
        plot_net(net, net_label, edge_weight = TRUE, node_colour = short_to_colour[net_type])
      } else {
        blank_plot(paste0(net_label, "\n(not asked)"))
      }
    })
  }

  if (site == "MK") {
    plot_list <-
    lapply(c(1,3,5,7,9,2,4,6,8,10), function(i) {
      net_type <- c("e1", "e2", "e3", "e4", "e56a", "e56b", "e78", "extra", "govt", "external")[i]
      net_label <- c(
        tie_meta$label[1:4],
        "Behavioral support (Agri. Team)",
        "Behavioral support (Market)",
        "Informational/affective support (Advice)",
        "Gender-specific informational/affective support",
        tie_meta$label[9:10]
      )[i]
      node_colour <- c(
        tie_colours[1:4],
        "#73C05B",
        "#73C05B",
        "orange",
        "orange",
        tie_colours[9:10]
      )[i]
      net <- su_nets$MK[[net_type]]
      if (!is.null(net)) {
        plot_net(net, net_label, edge_weight = TRUE, node_colour = node_colour)
      } else {
        blank_plot(paste0(net_label, "\n(not asked)"))
      }
    })
  }

  write_net_pdf(
    plot_list,
    path           = file.path(EndowDropbox, "Figures", "Networks",
                                paste0(site, "_all_nets.pdf")),
    plots_per_page = 10,
    grid_rows      = 2,
    grid_cols      = 5,
    width          = 10,
    height         = 5
  )
}

# ---- Individual agg nets as examples ----
eg_sites <- c("CR", "AZ", "PT")
if (!all(eg_sites %in% site_codes)) eg_sites <- head(site_codes, 3) # to work with faux data

eg_agg_nets <- lapply(eg_sites, function(site) {
    net <- su_nets_expanded[[site]]$sum
      plot_net(net, site, edge_weight = TRUE)
  })

write_net_pdf(
  eg_agg_nets,
  path           = file.path(EndowDropbox, "Figures", "Networks",
                              "eg_agg_nets.pdf"),
  plots_per_page = 3,
  grid_rows      = 1,
  grid_cols      = 3,
  width          = 10,
  height         = 5
)

eg_sites <- c("CP", "PT", "NI")
if (!all(eg_sites %in% site_codes)) eg_sites <- tail(site_codes, 3) # to work with faux data

eg_agg_nets_wealth <- lapply(eg_sites, function(site) {
    net <- su_nets_expanded[[site]]$sum
      plot_net(net, site, edge_weight = TRUE, node_size_attr = "su_wealth")
  })

write_net_pdf(
  eg_agg_nets_wealth,
  path           = file.path(EndowDropbox, "Figures", "Networks",
                              "eg_agg_nets_wealth.pdf"),
  plots_per_page = 3,
  grid_rows      = 1,
  grid_cols      = 3,
  width          = 10,
  height         = 5
)

