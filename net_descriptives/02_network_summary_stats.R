## ============================================================================
## 03_network_summary_stats.R
## ----------------------------------------------------------------------------
## Network-level summary statistics.
##
## Computes per-community, per-layer network summary statistics (density,
## degree, centralization, clustering, path length, etc.), saves them to
## Latex tables, and plots their cross-community distributions. Plus additional
## pieces looking at the use of edge weights, reciprocity across layers, etc.
##
## Outputs:
##   Tables/DescriptiveStats/network-summary-stats/net-summary-stat-*.tex
##   Tables/DescriptiveStats/network-summary-stats/global-net-summary-medians.tex
##   Figures/DescriptiveStats/*_distrib_violin.png
##   Figures/DescriptiveStats/density_vs_size.png
##   Figures/DescriptiveStats/selfloop_prop_distrib_violin.png
##   Figures/DescriptiveStats/compare_edge_weights.png
##   Figures/DescriptiveStats/obs_exp_apl.png
##   Figures/DescriptiveStats/recip_across_layers.png
##
## Depends on: 02_shared_definitions.R
## ============================================================================

source("net_descriptives/00_setup.R")

summarise_network <- function(g, site_name) {

  is_directed <- is_directed(g)

  # NOTE: edge weights reflect strength:
  # greater values imply more nominations of i by j.
  # For mean degree and transitivity, want to imply higher = stronger connection
  # Distance measures (av path length, diameter), assume greater values = more distant
  # So for those, doing 1/weight
  # for weighted transitivity, need undirected network without multiedges
  # collapsing both as mean value; could also justify "sum" for weight

  #if(is_directed) {
  g_und_mean <- as_undirected(
    g,
    mode = "collapse",
      edge.attr.comb = list(
      weight = "mean",
      weight_per_resp = "mean"
      )
    )

    # try also with "sum"
    g_und_sum <- as_undirected(
    g,
    mode = "collapse",
      edge.attr.comb = list(
      weight = "sum",
      weight_per_resp = "sum"
      )
    )
  #}

  data.frame(
    site              = site_name,
    n_nodes           = vcount(g),
    n_edges           = ecount(g),
    density           = edge_density(g, loops = FALSE),
    mean_degree       = mean(degree(g, mode = "in")), # had previously used "all" -- this is in + out, so will be double
    mean_degree_w     = if (is_directed) mean(strength(g, mode = "in", weights = igraph::edge_attr(g, "weight_per_resp")) ) else NA_real_,
    in_degree_range   = if (is_directed) paste( min(degree(g, mode = "in")), "-", max(degree(g, mode = "in")) ) else NA_character_,
    out_degree_range  = if (is_directed) paste( min(degree(g, mode = "out")), "-", max(degree(g, mode = "out")) ) else NA_character_,
    deg_centr         =  centr_degree(g, mode = "all")$centralization,
    in_deg_centr      = if (is_directed) centr_degree(g, mode = "in")$centralization else NA_real_,
    out_deg_centr     = if (is_directed) centr_degree(g, mode = "out")$centralization else NA_real_,
    deg_assort        = assortativity_degree(g, directed = is_directed),
    avg_path_length   = mean_distance(g, directed = is_directed, weights = NA),
    avg_path_length_w = mean_distance(g, directed = is_directed, weights = 1/igraph::edge_attr(g, "weight_per_resp")),
    diameter          = diameter(g, weights = NA), # without NA, uses "weight" attr by default
    diameter_w        = diameter(g, weights = 1/igraph::edge_attr(g, "weight_per_resp")),
    n_isolates        = sum(as.numeric(degree(g) == 0)),
    prop_isolates     = sum(as.numeric(degree(g) == 0)) / vcount(g),
    n_comp            = count_components(g, mode = "weak"),
    largest_comp      = max(components(g, mode = "weak")$csize) / vcount(g),
    global_clust      = transitivity(g, type = "global", weights = NA), # without NA, uses "weight" attr by default
    av_local_clust    = transitivity(g, type = "localaverage", weights = NA), # without NA, uses "weight" attr by default
    av_local_clust_w = if (is_directed) mean(transitivity(g_und_mean, type = "weighted", weights = igraph::edge_attr(g_und_mean, "weight_per_resp")), na.rm = TRUE) else NA_real_,
    av_local_clust_w_sum  = if (is_directed) mean(transitivity(g_und_sum, type = "weighted", weights = igraph::edge_attr(g_und_sum, "weight_per_resp")), na.rm = TRUE) else NA_real_,
    reciprocity       = if (is_directed) reciprocity(g) else NA_real_,
    ## number of *ordered* pairs with an edge -- the denominator of reciprocity,
    ## needed for the empirical-logit robustness specification
    recip_denom     = ecount(g),
    ## number of connected (directed) triples -- the denominator of global
    ## transitivity, needed likewise
    triple_denom    = sum(choose(degree(g_und_mean, mode = "all"), 2)),
    stringsAsFactors  = FALSE
  )
}

network_summaries <- lapply(tie_short[1:12], function(net_type) {

  rows <- lapply(site_codes, function(site) {
    g <- su_nets_expanded[[site]][[net_type]]
    if (is.null(g) || !inherits(g, "igraph")) {
      # Return an all-NA row if the network is missing for this site
      data.frame(
        site              = site,
        n_nodes           = NA_integer_,
        n_edges           = NA_integer_,
        density           = NA_real_,
        mean_degree       = NA_real_,
        mean_degree_w     = NA_real_,
        in_degree_range   = NA_character_,
        out_degree_range  = NA_character_,
        deg_centr         = NA_real_,
        in_deg_centr      = NA_real_,
        out_deg_centr      = NA_real_,
        deg_assort        = NA_real_,
        avg_path_length   = NA_real_,
        avg_path_length_w = NA_real_,
        diameter          = NA_real_,
        diameter_w        = NA_real_,
        n_isolates        = NA_real_,
        prop_isolates     = NA_real_,
        n_comp            = NA_real_,
        largest_comp      = NA_real_,
        global_clust      = NA_real_,
        av_local_clust    = NA_real_,
        av_local_clust_w  = NA_real_,
        av_local_clust_w_sum  = NA_real_,
        reciprocity       = NA_real_,
        recip_denom       = NA_real_,
        triple_denom      = NA_real_,
        stringsAsFactors  = FALSE
      )
    } else {
      summarise_network(g, site)
    }
  })
  do.call(rbind, rows)
})

names(network_summaries) <- tie_short[1:12]

network_summaries_long <- bind_rows(network_summaries, .id = "network") |>
  dplyr::mutate(
    network = factor(network,
      levels = tie_meta$short[1:12]
      ),
    network_label = factor(short_to_label[network],
      levels = rev(tie_meta$label[1:12])
    )
  ) |>
  dplyr::filter(!network %in% c("govt", "external", "primary_kin"))


## Save tables to Latex

net_summary_table_path <- file.path(EndowDropbox, "Tables", "DescriptiveStats", "network-summary-stats")
net_summary_csv_path <- file.path(csvpath, "network-summary-stats")

dir.create(net_summary_table_path, showWarnings = FALSE)

for (i in tie_short[1:12]) {
  base_name <- paste0("net-summary-stat-", i)
  filename <- paste0(net_summary_table_path, "/", base_name, ".tex")
  csv_filename <- paste0(net_summary_csv_path, "/", base_name, ".csv")
  message("Saving: ", base_name)

  if (i == "primary_kin") {

    # Undirected network: drop in/out-degree, in-degree centralization, and reciprocity
    tab <- network_summaries[[i]] |>
      dplyr::arrange(site) |>
      dplyr::select(
        site, n_nodes, n_edges, density, mean_degree, deg_centr, deg_assort,
        avg_path_length, diameter, n_isolates, prop_isolates,
        n_comp, largest_comp, global_clust, av_local_clust
      ) |>
      dplyr::mutate(dplyr::across(
        c(density, mean_degree, deg_centr, deg_assort,
          avg_path_length, diameter, prop_isolates,
          global_clust, av_local_clust, largest_comp),
        \(x) round(x, 2)
      )) |>
      dplyr::filter(!dplyr::if_all(where(is.numeric), is.na))

    csv_names <- c(
      "Site", "Nodes", "Edges", "Density", "Mean Deg.", "Deg. Centr.",
      "Deg. Assort.", "APL", "Diam.", "N Isol.", "Prop. Isol.", "N Comp.",
      "Largest Comp.", "Global Clust.", "Av. Loc. Clust."
    )

    tex_names <- c(
      "Site", "Nodes", "Edges", "Density",
      "\\shortstack{Mean\\\\ Deg.}",
      "\\shortstack{Deg.\\\\ Centr.}",
      "\\shortstack{Deg.\\\\ Assort.}",
      "APL", "Diam.", "N Isol.",
      "\\shortstack{Prop.\\\\ Isol.}",
      "N Comp.",
      "\\shortstack{Largest\\\\ Comp.}",
      "\\shortstack{Global\\\\ Clust.}",
      "\\shortstack{Av. Loc.\\\\ Clust.}"
    )

    caption <- "Network summary statistics by community for the close kin network. `Mean Deg.' shows the mean degree. `Deg. Centr.' shows the degree centralization. `Deg. Assort.' shows the degree assortativity. `APL' shows the average path length. `Diam.' shows the diameter. `N Isol.' shows the number of isolates, and `Prop. Isol.' their proportion of all nodes. `N Comp.' shows the number of weakly-connected components. `Largest Comp.' shows the proportion of nodes in the largest component. `Global Clust.' shows the global clustering coefficient. `Av. Local Clust.' shows the average local clustering coefficient."

  } else {

    tab <- network_summaries[[i]] |>
      dplyr::arrange(site) |>
      dplyr::select( # not including all of the weighted variants
        site, n_nodes, n_edges, density, mean_degree, mean_degree_w,
        in_degree_range, out_degree_range, deg_centr, in_deg_centr, deg_assort,
        avg_path_length, diameter, n_isolates, prop_isolates,
        n_comp, largest_comp, global_clust, av_local_clust, reciprocity
      ) |>
      dplyr::mutate(dplyr::across(
        c(density, mean_degree, mean_degree_w, deg_centr, in_deg_centr, deg_assort,
          avg_path_length, diameter, prop_isolates,
          global_clust, av_local_clust, largest_comp, reciprocity),
        \(x) round(x, 2)
      )) |>
      dplyr::filter(!dplyr::if_all(where(is.numeric), is.na))

    csv_names <- c(
      "Site", "Nodes", "Edges", "Density", "Mean Deg.", "Mean W. Deg.",
      "In-Deg.", "Out-Deg.", "Deg. Centr.", "In-Deg. Centr.", "Deg. Assort.",
      "APL", "Diam.", "N Isol.", "Prop. Isol.", "N Comp.", "Largest Comp.",
      "Global Clust.", "Av. Loc. Clust.", "Recip."
    )

    tex_names <- c(
      "Site", "Nodes", "Edges", "Density",
      "\\shortstack{Mean\\\\ Deg.}",
      "\\shortstack{Mean\\\\ W. Deg.}",
      "\\shortstack{In-\\\\ Deg.}",
      "\\shortstack{Out-\\\\ Deg.}",
      "\\shortstack{Deg.\\\\ Centr.}",
      "\\shortstack{In-Deg.\\\\ Centr.}",
      "\\shortstack{Deg.\\\\ Assort.}",
      "APL", "Diam.", "N Isol.",
      "\\shortstack{Prop.\\\\ Isol.}",
      "N Comp.",
      "\\shortstack{Largest\\\\ Comp.}",
      "\\shortstack{Global\\\\ Clust.}",
      "\\shortstack{Av. Loc.\\\\ Clust.}",
      "Recip."
    )

    caption <- paste0("Network summary statistics by community for the ", short_to_label[i], " layer. `Mean Deg.' shows the mean degree, while `Mean W. Deg.' uses edge weights. `In-Deg.' and `Out-Deg.' show the range. `(In-)Deg. Centr.' shows the (in-)degree centralization. `Deg. Assort.' shows the degree assortativity. `APL' shows the average path length. `Diam.' shows the diameter. `N Isol.' shows the number of isolates, and `Prop. Isol.' their proportion of all nodes. `N Comp.' shows the number of weakly-connected components. `Largest Comp.' shows the proportion of nodes in the largest component. `Global Clust.' shows the global clustering coefficient. `Av. Local Clust.' shows the average local clustering coefficient. `Recip.' shows reciprocity.")
  }

  # CSV
  tab |>
    setNames(csv_names) |>
    readr::write_csv(csv_filename)

  # LaTeX table
  tab |>
    dplyr::mutate(dplyr::across(where(is.numeric), \(x) ifelse(is.nan(x), "-", x))) |>
    kableExtra::kable(
      format    = "latex",
      booktabs  = TRUE,
      longtable = TRUE,
      escape    = FALSE,
      caption   = caption,
      label     = paste0("net-summary-stats-", gsub(" ", "-", i)),
      col.names = tex_names,
      linesep   = ""
    ) |>
    kableExtra::kable_styling(
      latex_options = "repeat_header",
      font_size     = 8
    ) |>
    kableExtra::save_kable(filename)
}

## Cross-site summary table showing median and then mean and SD for the main measures
stat_labels <- c(density = "Density",
                  mean_degree = "Mean Deg.", mean_degree_w = "Mean W. Deg.",
                  deg_centr = "Deg. Centr.", in_deg_centr = "In-Deg. Centr.",
                  deg_assort = "Deg. Assort.", avg_path_length = "APL",
                  diameter = "Diam.", n_isolates = "N Isol.", prop_isolates = "Prop. Isol.",
                  n_comp = "N Comp.", largest_comp = "Largest Comp.",
                  global_clust = "Global Clust.", av_local_clust = "Av. Local Clust.",
                  reciprocity = "Recip.")

xsite_summary_stats <- network_summaries_long |>
  group_by(network) |>
  dplyr::select(site, network, all_of(names(stat_labels))) |>
  dplyr::summarise(
    across(all_of(names(stat_labels)),
           list(
                median = \(x) median(x, na.rm = TRUE),
                mean   = \(x) mean(x, na.rm = TRUE),
                sd     = \(x) sd(x, na.rm = TRUE),
                cov    = \(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(network = factor(
    setNames(tie_meta$label_short[1:12], tie_meta$short[1:12])[network],
    levels = tie_meta$label_short[1:12]
  )) |>
  dplyr::arrange(network)

xsite_summary_long <- purrr::map_dfr(names(stat_labels), function(v) {
  tibble::tibble(
    network = xsite_summary_stats$network,
    stat    = stat_labels[[v]],
    Median      = as.character(round(xsite_summary_stats[[paste0(v, "_median")]], 2)),
    `Mean (SD)` = paste0(round(xsite_summary_stats[[paste0(v, "_mean")]], 2),
                         " (", round(xsite_summary_stats[[paste0(v, "_sd")]], 2), ")")
  )
})

xsite_summary_wide <- xsite_summary_long |>
  tidyr::pivot_longer(cols = c(Median, `Mean (SD)`), names_to = "type", values_to = "value") |>
  dplyr::mutate(
    type = factor(type, levels = c("Median", "Mean (SD)")),
    stat = factor(stat, levels = stat_labels)
  ) |>
  tidyr::pivot_wider(names_from = network, values_from = value) |>
  dplyr::arrange(stat, type) |>
  dplyr::select(-type)

xsite_summary_wide$stat <- as.character(xsite_summary_wide$stat)
xsite_summary_wide$stat[duplicated(xsite_summary_wide$stat)] <- ""

xsite_summary_wide |>
  kableExtra::kbl(
    format    = "latex",
    booktabs  = TRUE,
    longtable = TRUE,
    align     = c("l", rep("c", 9)),
    col.names = c("", tie_meta$label_short[c(1:8,12)]),
    caption   = "Median (above) and mean and standard deviation (below) of network summary statistics for each network layer. `Mean Deg.' shows the mean degree, while `Mean W. Deg.' uses edge weights. `(In-)Deg. Centr.' shows the (in-)degree centralization. `Deg. Assort.' shows the degree assortativity. `APL' shows the average path length. `Diam.' shows the diameter. `N Isol.' shows the number of isolates, and `Prop. Isol.' their proportion of all nodes. `N Comp.' shows the number of weakly-connected components. `Largest Comp.' shows the proportion of nodes in the largest component. `Global Clust.' shows the global clustering coefficient. `Av. Local Clust.' shows the average local clustering coefficient. `Recip.' shows reciprocity.",
    label     = "mean-median-net-summary-stats",
    linesep   = ""
  ) |>
  kableExtra::kable_styling(
    latex_options = c("hold_position", "repeat_header"),
    font_size     = 8
  ) |>
  writeLines(paste0(net_summary_table_path, "/global-net-summary-median-meansd.tex"))

# CSV version
xsite_summary_wide |>
  setNames(c("Statistic", tie_meta$label_short[c(1:8, 12)])) |>
  readr::write_csv(
    file.path(net_summary_csv_path, "global-net-summary-median-meansd.csv")
  )

## PLOTS

plot_vars <- list(
  density         = list(label = "density", bw = 0.02),
  mean_degree     = list(label = "mean degree", bw = 0.5),
  mean_degree_w   = list(label = "mean weighted degree", bw = 0.5),
  deg_centr       = list(label = "degree centralization", bw = 0.05),
  in_deg_centr    = list(label = "in-degree centralization", bw = 0.05),
  deg_assort      = list(label = "degree assortativity", bw = 0.05),
  avg_path_length = list(label = "average path length", bw = 0.2),
  avg_path_length_w = list(label = "average weighted path length", bw = 0.2),
  global_clust    = list(label = "global clustering coefficient", bw = 0.05),
  av_local_clust  = list(label = "average local clustering coefficient", bw = 0.05),
  av_local_clust_w  = list(label = "average weighted local clustering coefficient", bw = 0.05),
  reciprocity     = list(label = "reciprocity", bw = 0.05)
)

purrr::iwalk(plot_vars, \(spec, var) {

  p <- make_violin_plot(
    network_summaries_long, x_var = "network_label", y_var = var, group_var = "network_label",
    palette = setNames(tie_meta$colour, tie_meta$label),
    y_label = spec$label,
    hline = if (var == "deg_assort") 0 else NULL
  )

  ggsave(
    file.path(EndowDropbox, "Figures", "DescriptiveStats", paste0(var, "_distrib_violin.png")),
    plot = p, width = 9, height = 6
  )
  message("Saved: ", var)
})

combined_pairs <- list(
  clustering = c("global_clust", "av_local_clust"),
  av_local_clust = c("av_local_clust", "av_local_clust_w"),
  mean_degree = c("mean_degree", "mean_degree_w"),
  avg_path_length = c("avg_path_length", "avg_path_length_w")
)

purrr::iwalk(combined_pairs, \(vars, combo_name) {
  p_left  <- make_violin_plot(network_summaries_long, "network_label", vars[1], "network_label",
                               palette = setNames(tie_meta$colour, tie_meta$label),
                               y_label = plot_vars[[vars[1]]]$label)
  p_right <- make_violin_plot(network_summaries_long, "network_label", vars[2], "network_label",
                               palette = setNames(tie_meta$colour, tie_meta$label),
                               y_label = plot_vars[[vars[2]]]$label,
                               show_group_labels = FALSE)

  combined <- p_left + p_right + patchwork::plot_layout(nrow = 1)

  ggsave(
    file.path(EndowDropbox, "Figures", "DescriptiveStats", paste0(combo_name, "_combined_distrib_violin.png")),
    plot = combined, width = 14, height = 6
  )
  message("Saved combined: ", combo_name)
})

# purrr::iwalk(plot_vars, \(spec, var) {
#   p <- network_summaries_long |>
#     dplyr::mutate(
#       site_label = factor(site, levels = sort(unique(site)))
#     ) |>
#     ggplot(aes(x = .data[[var]], y = site_label, colour = network)) +
#     geom_segment(
#       data        = network_summaries_long |>
#                       dplyr::distinct(site) |>
#                       dplyr::mutate(site_label = factor(site, levels = sort(unique(site)))),
#       aes(x = -Inf, xend = Inf, y = site_label, yend = site_label),
#       colour    = "grey80",
#       linewidth = 0.3,
#       inherit.aes = FALSE
#     ) +
#     geom_point(
#       #position = ggplot2::position_jitter(height = 0.15, seed = 42),
#       size  = 5,
#       alpha = 0.9,
#       shape = 16
#     ) +
#     scale_colour_manual(
#       values = short_to_colour,
#       labels = setNames(tie_meta$label_short, tie_meta$short),
#       name   = "Network type"
#     ) +
#     labs(
#       x     = spec$label,
#       y     = NULL,
#       # title = paste0("Distribution of ", spec$label)
#     ) +
#     theme_classic(base_size = 13) +
#     theme(legend.position = "bottom")

#   ggsave(
#     file.path(EndowDropbox, "Figures", "DescriptiveStats", paste0(var, "_distrib_bysite.png")),
#     plot   = p,
#     width  = 12,
#     height = 10
#   )
#   message("Saved: ", var)
# })

density_vs_size <- network_summaries$sum |>
  ggplot(aes(x = n_nodes, y = density)) +
  # geom_smooth(method = "loess") +
  #geom_point() +
  #ggrepel::geom_text_repel(aes(label = site), size = 3, max.overlaps = 10) +
  geom_label(
    aes(label = site),
    family = "roboto",
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    y = "Edge Density in the Aggregate Network",
    x = "Number of Sharing Units"
  ) +
  theme_classic()

density_vs_size_log <- network_summaries$sum |>
  ggplot(aes(x = n_nodes, y = density)) +
  # geom_smooth(method = "loess") +
  #geom_point() +
  #ggrepel::geom_text_repel(aes(label = site), size = 3, max.overlaps = 10) +
  geom_point(
    alpha = 0.5
  ) +
  scale_x_log10(guide = "axis_logticks") +
  scale_y_log10(guide = "axis_logticks") +
  labs(
    y = " ",
    x = " "
  ) +
  theme_classic()

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "density_vs_size.png"),
  plot = density_vs_size + inset_element(density_vs_size_log, left = 0.4, bottom = 0.4, right = 1, top = 1),
  width = 5,
  height = 5
)

# Extract self-loop proportions
selfloop_df <- purrr::imap_dfr(su_nets_expanded_withloops, \(site_nets, site) {
  purrr::imap_dfr(site_nets, \(g, layer) {
    if (!igraph::is_igraph(g)) return(NULL)
    n_edges     <- igraph::ecount(g)
    n_selfloops <- sum(igraph::which_loop(g))
    dplyr::tibble(
      site     = site,
      network  = layer,
      prop_selfloop = if (n_edges > 0) n_selfloops / n_edges else NA_real_
    )
  })
}) |>
  dplyr::filter(network %in% tie_short, !network %in% "primary_kin") |>
  dplyr::mutate(network_label = factor(short_to_label[network], levels = rev(short_to_label[tie_short[-c(13,17)]]))) # -13 to drop "main_req", -17 to drop "any"

p <- make_violin_plot(
  selfloop_df, x_var = "network_label", y_var = "prop_selfloop", group_var = "network_label",
  palette = setNames(tie_meta$colour, tie_meta$label),
  y_label = "Proportion of ties that are self-loops"
)
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "selfloop_prop_distrib_violin.png"),
  plot   = p,
  width  = 9,
  height = 6
)


## comparing weighted and unweighted measures

# Map: unweighted variable name -> display label
var_labels <- tibble::tribble(
  ~var,              ~label,
  "avg_path_length", "Average Path Length",
  "diameter",        "Diameter",
  "mean_degree",     "Mean Degree",
  "av_local_clust",  "Average Local Clustering"
)

plots <- purrr::pmap(var_labels, function(var, label) {

  unweighted <- network_summaries$sum[[var]]
  weighted   <- network_summaries$sum[[paste0(var, "_w")]]

  r <- cor(unweighted, weighted, use = "complete.obs")

  network_summaries$sum |>
    mutate(unweighted = .data[[var]],
           weighted   = .data[[paste0(var, "_w")]]) |>
    ggplot(aes(x = unweighted, y = weighted)) +
    geom_abline(intercept = 0, slope = 1) +
    geom_label(
      aes(label = site),
      family = "roboto",
      size = 1.5,
      alpha = 0.5,
      linewidth = 0.25,
      label.padding = unit(0.3, "lines"),
      label.r = unit(0.6, "lines")
    ) +
    coord_cartesian(ylim = c(0, NA), xlim = c(0, NA)) +
    labs(
      title = sprintf("%s (r = %.2f)", label, r),
      x = "unweighted",
      y = "weighted"
    ) +
    theme_classic()
})

names(plots) <- var_labels$var

p <- wrap_plots(list(plots$mean_degree, plots$avg_path_length, plots$av_local_clust), ncol = 1)

showtext::showtext_opts(dpi = 300)

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "compare_edge_weights.png"),
  plot = p,
  width = 4,
  height = 8
)


# average path length: observed versus expected
p <- network_summaries$sum %>%
  mutate(
    ln_n           = log(n_nodes), # jackson's suggestion
    ln_n_over_ln_k = log(n_nodes) / log(mean_degree) # watts & strogatz
  ) %>%
  ggplot(aes(x = avg_path_length, y = ln_n_over_ln_k)) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  geom_label(
    aes(label = site),
    family = "roboto",
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  coord_cartesian(
    xlim = range(c(
      network_summaries$sum$avg_path_length,
      log(network_summaries$sum$n_nodes) / log(network_summaries$sum$mean_degree)
    )),
    ylim = range(c(
      network_summaries$sum$avg_path_length,
      log(network_summaries$sum$n_nodes) / log(network_summaries$sum$mean_degree)
    ))
  ) +
  labs(
    x     = "observed average path length",
    y     = "predicted average path length: ln(n) / ln(k)",
    title = paste(
      "Observed vs Predicted Average Path Length, Corr:",
      round(cor(
        network_summaries$sum$avg_path_length,
        log(network_summaries$sum$n_nodes) / log(network_summaries$sum$mean_degree)
      ), 2)
    )
  ) +
  theme_classic()

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "obs_exp_apl.png"),
  plot = p,
  width  = 5.5,
  height = 5.5
  )


## Reciprocity across layers
aggregate_support_layers <- c(
  "loan_req", "exchange_req", "f_behav", "m_behav", "f_info", "m_info"
)

dyad_list <- su_nets_expanded %>%
  map(~ {
    .x %>%
      keep(~ !is.null(.x)) %>%
      map(~ igraph::as_data_frame(.x)) %>%
      .[names(.) %in% tie_short[1:10]] %>%
      # Select from, to and add layer column with value 1
      imap(~ dplyr::select(.x, from, to) %>% dplyr::mutate(!!.y := 1)) %>%
      # Join all by from and to
      reduce(full_join, by = c("from", "to")) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        sup_get_sum = rowSums(dplyr::pick(any_of(aggregate_support_layers))), # include only main within-comm requests
        sup_giv_sum = rowSums(dplyr::pick(any_of(c("loan_giv", "exchange_giv")))) # include only double-sampled questions
        )
  })

# Now add alter nominations and reciprocity columns
dyad_list <- dyad_list %>%
  imap(
    function(df, site_name) {
      # Get layer columns
      exclude_cols <- c(
        "from",
        "to",
        "sup_get_sum",
        "sup_giv_sum"
      )
      layer_names <- setdiff(names(df), exclude_cols)

      # Create ego nomination columns (rename original layers)
      df_ego <- df %>%
        dplyr::rename_with(~ paste0(.x, "_egonom"), dplyr::all_of(layer_names))

      # Create alter nomination columns (reverse from/to, rename layers)
      df_alter <- df %>%
        dplyr::select(from, to, dplyr::all_of(layer_names)) %>%
        dplyr::rename(from = to, to = from) %>% # Swap from and to
        dplyr::rename_with(~ paste0(.x, "_alternom"), dplyr::all_of(layer_names))

      # Join back together
      result <- df_ego %>%
        dplyr::left_join(df_alter, by = c("from", "to")) %>%
        dplyr::mutate(across(ends_with("_alternom"), ~ replace_na(.x, 0)))

      # Add reciprocity columns for each layer
      for (layer in layer_names) {
        result <- result %>%
          dplyr::mutate(
            !!paste0(layer, "_recip") := case_when(
              .data[[paste0(layer, "_egonom")]] == 0 ~ NA_real_,
              .data[[paste0(layer, "_alternom")]] == 1 ~ 1,
              TRUE ~ 0
            )
          )
      }

      # Use the paper's constructed aggregate so site-specific prompt
      # combinations (especially MK's e56a/e78) are included consistently.
      aggregate_reverse <- igraph::as_data_frame(
        su_nets_expanded[[site_name]]$sum, what = "edges"
      ) %>%
        dplyr::select(from, to) %>%
        dplyr::rename(from = to, to = from) %>%
        dplyr::distinct(from, to) %>%
        dplyr::mutate(any_altersup = 1L)

      # any_alternom separately includes all available standard prompts.
      alternom_cols <- paste0(layer_names, "_alternom")
      result <- result %>%
        dplyr::left_join(aggregate_reverse, by = c("from", "to")) %>%
        dplyr::mutate(
          any_altersup = tidyr::replace_na(any_altersup, 0L),
          any_alternom = if_else(
            rowSums(dplyr::select(., all_of(alternom_cols))) > 0,
            1,
            0
          )
        )

      result
    }
  )

# now going long by layer, such that each observed edge gets its own row
edge_df <- bind_rows(dyad_list, .id = "site") %>%
  pivot_longer(
    cols = -c(
      site,
      from,
      to,
      sup_get_sum,
      sup_giv_sum,
      ends_with("_alternom"),
      ends_with("_recip"),
      any_altersup,
      any_alternom
    ),
    names_to = "layer",
    values_to = "edge_exists"
  ) %>%
  dplyr::mutate(layer = str_remove(layer, "_egonom")) %>%
  dplyr::filter(edge_exists == 1) %>%
  dplyr::relocate(layer, .after = site)


## ---- cross-layer reciprocity ----
layer_recip_summary <- edge_df %>%
  group_by(site, layer) %>%
  group_modify(
    ~ {
      lyr <- .y$layer
      sv <- .y$site

      # Helper: proportion of current group's dyads that exist in check_layer
      prop_egonom <- function(check_layer) {
        layer_pairs <- edge_df %>%
          dplyr::filter(site == sv, layer == check_layer) %>%
          dplyr::mutate(pair = paste(from, to)) %>%
          dplyr::pull(pair)
        mean(paste(.x$from, .x$to) %in% layer_pairs)
      }

      .x %>%
        dplyr::summarise(
          n = n(),
          mean_sup = mean(sup_get_sum),
          median_sup = median(sup_get_sum),
          sd_sup = sd(sup_get_sum),
          prop_recip = mean(.data[[paste0(lyr, "_recip")]], na.rm = TRUE),
          prop_egonom_loan_giv = prop_egonom("loan_giv"),
          prop_egonom_loan_req = prop_egonom("loan_req"),
          prop_egonom_exchange_giv = prop_egonom("exchange_giv"),
          prop_egonom_exchange_req = prop_egonom("exchange_req"),
          prop_egonom_f_behav = prop_egonom("f_behav"),
          prop_egonom_m_behav = prop_egonom("m_behav"),
          prop_egonom_f_info = prop_egonom("f_info"),
          prop_egonom_m_info = prop_egonom("m_info"),
          prop_alternom_loan_giv = mean(loan_giv_alternom, na.rm = TRUE),
          prop_alternom_loan_req = mean(loan_req_alternom, na.rm = TRUE),
          prop_alternom_exchange_giv = mean(
            exchange_giv_alternom,
            na.rm = TRUE
          ),
          prop_alternom_exchange_req = mean(
            exchange_req_alternom,
            na.rm = TRUE
          ),
          prop_alternom_f_behav = mean(f_behav_alternom, na.rm = TRUE),
          prop_alternom_m_behav = mean(m_behav_alternom, na.rm = TRUE),
          prop_alternom_f_info = mean(f_info_alternom, na.rm = TRUE),
          prop_alternom_m_info = mean(m_info_alternom, na.rm = TRUE),
          prop_any_altersup = mean(any_altersup, na.rm = TRUE)
        )
    }
  ) %>%
  ungroup()

alter_recip_df <- layer_recip_summary |>
  dplyr::rename(prop_alternom_any = prop_any_altersup) |>
  dplyr::select(site, layer, starts_with("prop_alternom_")) |>
  pivot_longer(
    starts_with("prop_alternom_"),
    names_to = "recip_layer",
    names_prefix = "prop_alternom_",
    values_to = "reciprocity"
  ) |>
  dplyr::mutate(
    network_label = factor(short_to_label_short[layer],
      levels = short_to_label_short[tie_short[-c(12, 13)]]),
    recip_label = factor(short_to_label[recip_layer],
      levels = rev(short_to_label[tie_short[-c(12, 13)]]))
    ) |>
  dplyr::filter(!layer %in% c("external", "govt"))

p <- ggplot(alter_recip_df, aes(x = recip_label, y = reciprocity, fill = recip_label, colour = recip_label)) +
  geom_violin(
      alpha = 0.4,
      linewidth = 0.5,
      scale = "width",
      draw_quantiles = 0.5
    ) +
  #geom_label(
  #  aes(label = site),
  #  family = "roboto",
  #  size = 0.5,
  #  alpha = 0.5,
  #  linewidth = 0.25,
  #  label.padding = unit(0.3, "lines"),
  #  label.r = unit(0.6, "lines")
  #) +
  geom_jitter(alpha = 0.25, width = 0.2) +
  #coord_cartesian(ylim = c(0,1)) +
  ggh4x::facet_wrap2(~ network_label, nrow = 2, dir = "v", axes = "x") +
  scale_fill_manual(values = tie_colours, guide = "none") +
  scale_colour_manual(values = tie_colours, guide = "none") +
    labs(
    x = "",
    y = "cross-layer reciprocity",
  ) +
  theme_classic() +
    theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    panel.spacing = unit(0.5, "lines")
  ) +
  coord_flip()

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "recip_across_layers.png"),
  plot = p,
  width  = 12,
  height = 8
  )

# Save network_summaries_long, so can be called in by respondent_coverage
save(network_summaries_long, file = file.path(EndowDerivedData, "network_summaries_long.rdata"))
