## ============================================================================
## 10_pca.R
## ----------------------------------------------------------------------------
## Principal component analysis of network layers.
##
## Runs PCA on the stacked (weighted and binarized) edge data across network
## layers, by community and pooled, and saves scree plots, correlation-circle
## biplots, and loadings tables.
##
## Outputs:
##   Tables/pca_loadings_weighted_pooled.tex
##   Tables/pca_loadings_weighted_sites.tex
##   Tables/pca_loadings_binary_sites.tex
##   Figures/Networks/pca_scree_weighted_sites.pdf
##   Figures/Networks/pca_biplot_weighted_sites.pdf
##   Figures/Networks/pca_biplot_weighted_23_sites.pdf
##   Figures/Networks/pca_scree_binary_sites.pdf
##   Figures/Networks/pca_biplot_binary_sites.pdf
##   Figures/Networks/pca_biplot_binary_23_sites.pdf
##   Figures/Networks/pca_scree_weighted_pooled.png
##   Figures/Networks/pca_biplot_weighted_pooled.png
##   Figures/Networks/pca_biplot_weighted_23_pooled.png
##   Figures/Networks/pca_scree_binary_pooled.png
##   Figures/Networks/pca_biplot_binary_pooled.png
##   Figures/Networks/pca_biplot_binary_23_pooled.png
##   Figures/Networks/pca_biplot_weighted_compare.png
##
## Depends on: 00_setup.R
## ============================================================================

source("net_descriptives/00_setup.R")

# ---- Plotting helpers ----
facets_per_page <- 20
n_pages <- ceiling(length(site_codes) / facets_per_page)

# ---- Principal Component Analysis --------------------------------------------

# turn adj matrices into dyad data frames, either weighted or binary
# (one row for each dyad, showing edges across the various network layers)
# compute PCA on this object

compute_layer_pca <- function(
  sites,
  net_types,
  nets_list,
  binary = FALSE,
  scale = TRUE,
  include_aggregate = FALSE,
  aggregate_site_label = "all"
) {
  run_pca <- function(site_vec) {
    stacked <- as.data.frame(setNames(
      lapply(net_types, function(nt) {
        vals <- unlist(lapply(site_vec, function(site) {
          all_nodes <- V(nets_list[[site]][[net_types[1]]])$name
          as.vector(get_adj_weighted(nets_list[[site]][[nt]], all_nodes))
        }))
        if (binary) (vals > 0) * 1L else vals
      }),
      net_types
    ))
    stacked <- stacked[stats::complete.cases(stacked), ]

    keep <- vapply(
      stacked,
      function(v) stats::var(v, na.rm = TRUE) > 0,
      logical(1)
    )
    dropped <- names(keep)[!keep]
    if (length(dropped) > 0) {
      message(
        "Dropping zero-variance layers: ",
        paste(dropped, collapse = ", ")
      )
    }
    stacked <- stacked[, keep, drop = FALSE]

    if (ncol(stacked) < 2) {
      message("Fewer than 2 variable layers remaining — skipping PCA")
      return(NULL)
    }

    stats::prcomp(stacked, scale. = scale, center = TRUE)
  }

  site_pcas <- setNames(
    lapply(sites, function(site) run_pca(site)),
    sites
  )

  if (!include_aggregate) {
    return(purrr::compact(site_pcas))
  }

  agg_pca <- run_pca(sites)

  purrr::compact(c(site_pcas, setNames(list(agg_pca), aggregate_site_label)))
}

layer_pca_weighted <- compute_layer_pca(
  sites = site_codes,
  net_types = tie_meta$short[1:8],
  nets_list = su_nets_expanded,
  binary = FALSE,
  include_aggregate = TRUE
)

layer_pca_binary <- compute_layer_pca(
  sites = site_codes,
  net_types = tie_meta$short[1:8],
  nets_list = su_nets_expanded,
  binary = TRUE,
  include_aggregate = TRUE
)

extract_pca_dfs <- function(layer_pca, pc_x = 1, pc_y = 2) {
  variance <- purrr::map_dfr(
    layer_pca,
    function(pca) {
      var_exp <- pca$sdev^2 / sum(pca$sdev^2)
      data.frame(
        pc = paste0("Principal Component ", seq_along(var_exp)),
        variance_explained = var_exp,
        cumulative = cumsum(var_exp)
      )
    },
    .id = "site"
  )

  loadings <- purrr::map_dfr(
    layer_pca,
    function(pca) {
      as.data.frame(pca$rotation) |>
        tibble::rownames_to_column("layer")
    },
    .id = "site"
  ) |>
    dplyr::mutate(
      layer_label = factor(
        short_to_label_short[layer],
        levels = short_to_label_short[net_type_order] # defined in 00_setup
      )
    )

  biplot <- purrr::map_dfr(
    layer_pca,
    function(pca) {
      rot <- as.data.frame(pca$rotation)
      var_exp <- pca$sdev^2 / sum(pca$sdev^2)
      pcx_col <- paste0("PC", pc_x) # index into prcomp object
      pcy_col <- paste0("PC", pc_y)
      pcx_label <- paste0("Principal Component ", pc_x) # display only
      pcy_label <- paste0("Principal Component ", pc_y)
      data.frame(
        layer = rownames(rot),
        PCx = rot[[pcx_col]],
        PCy = rot[[pcy_col]],
        pcx_var = round(var_exp[pc_x] * 100, 1),
        pcy_var = round(var_exp[pc_y] * 100, 1),
        pcx_label = pcx_label,
        pcy_label = pcy_label
      )
    },
    .id = "site"
  ) |>
    dplyr::mutate(layer_label = short_to_label_short[layer])

  list(variance = variance, loadings = loadings, biplot = biplot)
}


# ---- Extract data frames from both PCA objects ----
pca_weighted <- extract_pca_dfs(layer_pca_weighted)
pca_binary <- extract_pca_dfs(layer_pca_binary)
pca_weighted23 <- extract_pca_dfs(layer_pca_weighted, pc_x = 2, pc_y = 3)
pca_binary23 <- extract_pca_dfs(layer_pca_binary, pc_x = 2, pc_y = 3)

# ---- Split into site-wise and aggregate ----
pca_weighted_variance_sites <- dplyr::filter(
  pca_weighted$variance,
  site != "all"
)
pca_weighted_variance_agg <- dplyr::filter(pca_weighted$variance, site == "all")
pca_weighted_loadings_sites <- dplyr::filter(
  pca_weighted$loadings,
  site != "all"
)

pca_weighted_loadings_agg <- dplyr::filter(pca_weighted$loadings, site == "all")
pca_weighted_biplot_sites <- dplyr::filter(pca_weighted$biplot, site != "all")
pca_weighted_biplot_23_sites <- dplyr::filter(pca_weighted23$biplot, site != "all")
pca_weighted_biplot_agg <- dplyr::filter(pca_weighted$biplot, site == "all")
pca_weighted_biplot_23_agg <- dplyr::filter(pca_weighted23$biplot, site == "all")

pca_binary_variance_sites <- dplyr::filter(pca_binary$variance, site != "all")
pca_binary_variance_agg <- dplyr::filter(pca_binary$variance, site == "all")
pca_binary_loadings_sites <- dplyr::filter(pca_binary$loadings, site != "all")
pca_binary_loadings_agg <- dplyr::filter(pca_binary$loadings, site == "all")
pca_binary_biplot_sites <- dplyr::filter(pca_binary$biplot, site != "all")
pca_binary_biplot_23_sites <- dplyr::filter(pca_binary23$biplot, site != "all")
pca_binary_biplot_agg <- dplyr::filter(pca_binary$biplot, site == "all")
pca_binary_biplot_23_agg <- dplyr::filter(pca_weighted23$biplot, site == "all")

# sites for comparison:
compare_sites <- c("AZ", "TI", "NI", "PT")
if (!all(compare_sites %in% site_codes)) compare_sites <- site_codes # to work with faux data

pca_weighted_biplot_az <- dplyr::filter(pca_weighted$biplot, site == compare_sites[1])
pca_weighted_biplot_ti <- dplyr::filter(pca_weighted$biplot, site == compare_sites[2])
pca_weighted_biplot_ni <- dplyr::filter(pca_weighted$biplot, site == compare_sites[3])
pca_weighted_biplot_yi <- dplyr::filter(pca_weighted$biplot, site == compare_sites[4])

pca_weighted_biplot_23_az <- dplyr::filter(pca_weighted23$biplot, site == compare_sites[1])
pca_weighted_biplot_23_ti <- dplyr::filter(pca_weighted23$biplot, site == compare_sites[2])
pca_weighted_biplot_23_ni <- dplyr::filter(pca_weighted23$biplot, site == compare_sites[3])
pca_weighted_biplot_23_yi <- dplyr::filter(pca_weighted23$biplot, site == compare_sites[4])

# check on % of variance explained
summary(pca_weighted_variance_sites$variance_explained[
  pca_weighted_variance_sites$pc == "Principal Component 1"
])

summary(pca_binary_variance_sites$variance_explained[
  pca_binary_variance_sites$pc == "Principal Component 1"
])

# ---- Scree plot function ----
scree_plot <- function(df) {
  ggplot(df, aes(x = pc, y = variance_explained)) +
    geom_col(fill = "steelblue", width = 0.6) +
    geom_line(aes(group = 1), colour = "grey40", linewidth = 0.3) +
    geom_point(size = 0.8, colour = "grey40") +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    theme_classic(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      strip.text = element_text(size = 6, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    labs(x = NULL, y = "Variance explained")
}

make_scree_paginated <- function(page_data, page_sites, p, title_str) {
  scree_plot(page_data) +
    facet_wrap(~site, ncol = 4, nrow = 5) +
    labs(title = title_str)
}

# ---- Biplot function ----
biplot_plot <- function(df) {
  circle <- data.frame(
    x = cos(seq(0, 2 * pi, length.out = 200)),
    y = sin(seq(0, 2 * pi, length.out = 200))
  )
  axis_labels <- df |>
    dplyr::group_by(site) |>
    dplyr::slice(1) |>
    dplyr::select(site, pcx_var, pcy_var, pcx_label, pcy_label)

  layer_colours <- setNames(tie_meta$colour, tie_meta$label_short)

  ggplot(df, aes(x = PCx, y = PCy, colour = layer_label)) +
    geom_path(
      data = circle,
      aes(x = x, y = y),
      colour = "grey80",
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.2) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.2) +
    geom_segment(
      aes(x = 0, y = 0, xend = PCx, yend = PCy, colour = layer_label),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
      linewidth = 0.4,
      alpha = 0.6
    ) +
    ggrepel::geom_text_repel(
      aes(label = layer_label, colour = layer_label),
      size = 1.8,
      max.overlaps = Inf,
      segment.size = 0.3,
      segment.alpha = 0.6,
      box.padding = 0.3,
      point.padding = 0.1,
      force = 2,
      min.segment.length = 0
    ) +
    geom_text(
      data = axis_labels,
      aes(x = 0, y = 1.05, label = paste0(pcx_label, " (", pcx_var, "%)")),
      hjust = 0.5,
      vjust = 0,
      size = 1.8,
      colour = "grey40",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = axis_labels,
      aes(x = 1.05, y = 0, label = paste0(pcy_label, " (", pcy_var, "%)")),
      hjust = 0.5,
      vjust = 0,
      angle = -90,
      size = 1.8,
      colour = "grey40",
      inherit.aes = FALSE
    ) +
    scale_colour_manual(values = layer_colours) +
    scale_x_continuous(limits = c(-1.1, 1.1)) +
    scale_y_continuous(limits = c(-1.1, 1.1)) +
    coord_fixed(clip = "off") +
    theme_classic(base_size = 7) +
    theme(
      axis.text = element_text(size = 5),
      axis.title = element_text(size = 6),
      strip.text = element_text(size = 6, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    labs(x = axis_labels$pcx_label[1], y = axis_labels$pcy_label[1])
}

make_biplot_paginated <- function(page_data, page_sites, p, title_str) {
  biplot_plot(page_data) +
    facet_wrap(~site, ncol = 4, nrow = 5) +
    labs(title = title_str)
}

# ---- Site-wise PDFs: weighted ----
make_pdf(
  "pca_scree_weighted_sites.pdf",
  width = 14,
  height = ceiling(facets_per_page / 6) * 3.5,
  data_sites = pca_weighted_variance_sites,
  plot_fn = make_scree_paginated,
  title_str = NULL #"Variance explained by PC, weighted"
)

make_pdf(
  "pca_biplot_weighted_sites.pdf",
  width = 14,
  height = ceiling(facets_per_page / 6) * 3.5,
  data_sites = pca_weighted_biplot_sites,
  plot_fn = make_biplot_paginated,
  title_str = NULL #"PCA correlation circle, weighted"
)

make_pdf(
  "pca_biplot_weighted_23_sites.pdf",
  width = 14,
  height = ceiling(facets_per_page / 6) * 3.5,
  data_sites = pca_weighted_biplot_23_sites,
  plot_fn = make_biplot_paginated,
  title_str = NULL #"PCA correlation circle, weighted"
)

# ---- Site-wise PDFs: binary ----
make_pdf(
  "pca_scree_binary_sites.pdf",
  width = 14,
  height = ceiling(facets_per_page / 6) * 3.5,
  data_sites = pca_binary_variance_sites,
  plot_fn = make_scree_paginated,
  title_str = NULL #"Variance explained by PC, binary"
)

make_pdf(
  "pca_biplot_binary_sites.pdf",
  width = 14,
  height = ceiling(facets_per_page / 6) * 3.5,
  data_sites = pca_binary_biplot_sites,
  plot_fn = make_biplot_paginated,
  title_str = NULL #"PCA correlation circle, binary"
)

make_pdf(
  "pca_biplot_binary_23_sites.pdf",
  width = 14,
  height = ceiling(facets_per_page / 6) * 3.5,
  data_sites = pca_binary_biplot_23_sites,
  plot_fn = make_biplot_paginated,
  title_str = NULL #"PCA correlation circle, weighted"
)

# ---- Aggregate PNGs: weighted ----
make_png(
  "pca_scree_weighted_pooled.png",
  width = 5,
  height = 4,
  scree_plot(pca_weighted_variance_agg) +
    labs(title = NULL #"Variance explained by PC, weighted (pooled)"
  )
)

make_png(
  "pca_biplot_weighted_pooled.png",
  width = 4,
  height = 4,
  biplot_plot(pca_weighted_biplot_agg) +
    labs(title = NULL #"PCA correlation circle, weighted (pooled)"
  )
)

make_png(
  "pca_biplot_weighted_23_pooled.png",
  width = 4,
  height = 4,
  biplot_plot(pca_weighted_biplot_23_agg) +
    labs(
      title = NULL #"PCA correlation circle, weighted (pooled)"
    )
)

# ---- Aggregate PNGs: binary ----
make_png(
  "pca_scree_binary_pooled.png",
  width = 5,
  height = 4,
  scree_plot(pca_binary_variance_agg) +
    labs(title = NULL #"Variance explained by PC, binary (pooled)"
  )
)

make_png(
  "pca_biplot_binary_pooled.png",
  width = 4,
  height = 4,
  biplot_plot(pca_binary_biplot_agg) +
    labs(title = NULL #"PCA correlation circle, binary (pooled)"
  )
)

make_png(
  "pca_biplot_binary_23_pooled.png",
  width = 4,
  height = 4,
  biplot_plot(pca_binary_biplot_23_agg) +
    labs(
      title = NULL #"PCA correlation circle, binary (pooled)"
    )
)

## ---- Comparison of four selected sites in 2 x 2 ----
pca_weighted_biplot_compare <- dplyr::filter(pca_weighted$biplot, site %in% compare_sites)

p_compare <- biplot_plot(pca_weighted_biplot_compare) +
  facet_wrap(~ site, nrow = 2, ncol = 2)

make_png(
  "pca_biplot_weighted_compare.png",
  width = 5,
  height = 5,
  p_compare
)

# Loadings tables
pcs <- paste0("PC", 1:8)

pca_weighted$loadings |>
  dplyr::select(site, layer_label, dplyr::any_of(pcs)) |>
  dplyr::bind_rows(
    pca_weighted$variance |>
      dplyr::mutate(pc = sub("Principal Component ", "PC", pc)) |>
      dplyr::select(site, pc, variance_explained) |>
      tidyr::pivot_wider(names_from = pc, values_from = variance_explained) |>
      dplyr::mutate(layer_label = "Variance explained", .after = site)
  ) |>
  dplyr::mutate(
    site = as.character(site),
    dplyr::across(dplyr::any_of(pcs), \(x) round(x, 2))
  ) |>
  dplyr::arrange(site == "all", site, layer_label != "Variance explained") |>
  setNames(c("Community", "Network Layer or Variance Explained", pcs)) |>
  readr::write_csv(file.path(csvpath, "pca_loadings_weighted.csv"))

pca_binary$loadings |>
  dplyr::select(site, layer_label, dplyr::any_of(pcs)) |>
  dplyr::bind_rows(
    pca_binary$variance |>
      dplyr::mutate(pc = sub("Principal Component ", "PC", pc)) |>
      dplyr::select(site, pc, variance_explained) |>
      tidyr::pivot_wider(names_from = pc, values_from = variance_explained) |>
      dplyr::mutate(layer_label = "Variance explained", .after = site)
  ) |>
  dplyr::mutate(
    site = as.character(site),
    dplyr::across(dplyr::any_of(pcs), \(x) round(x, 2))
  ) |>
  dplyr::arrange(site == "all", site, layer_label != "Variance explained") |>
  setNames(c("Community", "Network Layer or Variance Explained", pcs)) |>
  readr::write_csv(file.path(csvpath, "pca_loadings_binary.csv"))


# ---- Helper to format loadings as kable ----
format_loadings_table <- function(df, caption = NULL, by_site = FALSE) {
  pcs <- paste0("PC", 1:8)

  sites_ordered <- sort(unique(df$site))

  selected <- purrr::map_dfr(sites_ordered, function(s) {
    df |>
      dplyr::filter(site == s) |>
      dplyr::select(layer_label, dplyr::any_of(pcs)) |>
      dplyr::mutate(dplyr::across(dplyr::any_of(pcs), ~ round(.x, 2)))
  })

  pc_cols <- setdiff(names(selected), "layer_label")

  default_caption <- paste0(
    caption,
    ". ",
    "Each entry is the correlation (loading) between a network layer and a principal component. ",
    "Loadings near $\\pm 1$ indicate the layer is strongly associated with that component; ",
    "loadings near 0 indicate little association. ",
    "Reading across a row shows how a layer's variance is distributed across components; ",
    "reading down a column reveals the substantive meaning of each component. ",
    "Components are ordered by descending variance explained."
  )

  tbl <- selected |>
    kableExtra::kbl(
      caption = default_caption,
      col.names = c("Layer", pc_cols),
      booktabs = TRUE,
      format = "latex",
      align = c("l", rep("r", length(pc_cols))),
      escape = FALSE,
      longtable = by_site # only need longtable when spanning many sites
    ) |>
    kableExtra::kable_styling(
      latex_options = c(
        #"striped",
        if (by_site) "repeat_header" else "hold_position"
      ),
      font_size = 8
    )

  if (!by_site) {
    return(tbl)
  }

  # Add site headings via pack_rows
  rows_per_site <- purrr::map_int(sites_ordered, function(s) {
    df |>
      dplyr::filter(site == s) |>
      dplyr::select(layer_label, dplyr::any_of(pcs)) |>
      nrow()
  })
  row_end <- cumsum(rows_per_site)
  row_start <- c(1, head(row_end, -1) + 1)

  purrr::reduce2(
    sites_ordered,
    seq_along(sites_ordered),
    function(tbl, s, i) {
      kableExtra::pack_rows(
        tbl,
        s,
        row_start[i],
        row_end[i],
        bold = TRUE,
        latex_gap_space = "0.5em"
      )
    },
    .init = tbl
  )
}

# ---- Aggregate ----
kableExtra::save_kable(
  format_loadings_table(pca_weighted_loadings_agg,
                        caption = "PCA loadings by layer, weighted (pooled)"),
  file = file.path(EndowDropbox, "Tables", "pca_loadings_weighted_pooled.tex")
)

# ---- Site-wise ----
kableExtra::save_kable(
  format_loadings_table(pca_weighted_loadings_sites,
                        caption = "PCA loadings by layer, weighted",
                        by_site = TRUE),
  file = file.path(EndowDropbox, "Tables", "pca_loadings_weighted_sites.tex")
)

kableExtra::save_kable(
  format_loadings_table(
    pca_binary_loadings_sites,
    caption = "PCA loadings by layer, binary",
    by_site = TRUE
  ),
  file = file.path(EndowDropbox, "Tables", "pca_loadings_binary_sites.tex")
)

## ALT Version to get MK in the edge overlap plots
## Could do this, but then would want to do throughout and just not worth it!
## Kept here and commented out as demo of possibility if want to revisit.

# su_nets_expanded_mk <- su_nets_expanded
# su_nets_expanded_mk$MK$f_behav <- su_nets$MK$e56a
# su_nets_expanded_mk$MK$m_behav <- su_nets$MK$e56a
# su_nets_expanded_mk$MK$f_info  <- su_nets$MK$e78
# su_nets_expanded_mk$MK$m_info  <- su_nets$MK$e78

# mk_full_overlap <- compute_layer_overlap(
#   sites = "MK",
#   net_types = net_type_order,
#   nets_list = su_nets_expanded_mk,
#   include_aggregate = FALSE
# )

# edge_overlap <- edge_overlap |>
#   dplyr::filter(site != "MK") |>
#   dplyr::bind_rows(mk_full_overlap)

# merge_groups <- list(
#   behav = c("f_behav", "m_behav"),
#   info  = c("f_info", "m_info")
# )

# site_merges <- list(MK = merge_groups)  # extend per-site as needed

# block_of <- function(tie, site) {
#   groups <- site_merges[[site]]
#   if (is.null(groups)) return(tie)
#   hit <- purrr::keep(names(groups), \(g) tie %in% groups[[g]])
#   if (length(hit) == 0) tie else hit[1]
# }

# block_positions <- function(site) {
#   tibble::tibble(
#     tie = net_type_order,
#     pos = seq_along(net_type_order)
#   ) |>
#     dplyr::mutate(block = purrr::map_chr(tie, block_of, site = site)) |>
#     dplyr::group_by(block) |>
#     dplyr::summarise(center = mean(pos), width = dplyr::n(), .groups = "drop")
# }

# collapse_blocks <- function(df, target_site) {
#   bp <- block_positions(target_site)

#   df |>
#     dplyr::filter(site == .env$target_site) |>
#     dplyr::mutate(
#       block1 = purrr::map_chr(layer1, block_of, site = .env$target_site),
#       block2 = purrr::map_chr(layer2, block_of, site = .env$target_site),
#       value  = dplyr::if_else(block1 == block2, NA_real_, value)
#     ) |>
#     dplyr::distinct(block1, block2, value, .keep_all = TRUE) |>
#     dplyr::left_join(bp, by = c("block2" = "block")) |>   # layer2/col -> x
#     dplyr::rename(x_center = center, x_width = width) |>
#     dplyr::left_join(bp, by = c("block1" = "block")) |>   # layer1/row -> y
#     dplyr::rename(y_center = center, y_width = width)
# }

# plot_df <- edge_overlap |>
#   dplyr::mutate(
#     row_idx = match(layer1, net_type_order),
#     col_idx = match(layer2, net_type_order),
#     value = dplyr::case_when(
#       row_idx == col_idx ~ NA_real_,
#       row_idx < col_idx ~ weighted,
#       row_idx > col_idx ~ binary
#     ),
#     layer1_label = factor(
#       short_to_label_short[layer1],
#       levels = short_to_label_short[net_type_order]
#     ),
#     layer2_label = factor(
#       short_to_label_short[layer2],
#       levels = short_to_label_short[net_type_order]
#     )
#   )

# plot_df_mk <- collapse_blocks(plot_df, "MK")

# heatmap_plot_blocks <- function(df, label_lookup) {
#   ggplot(df, aes(x = x_center, y = y_center, fill = value)) +
#     geom_tile(aes(width = x_width, height = y_width),
#               colour = "white", linewidth = 0.3) +
#     geom_text(
#       data = \(d) dplyr::filter(d, !is.na(value)),
#       aes(label = sprintf("%.2f", value)),
#       size = 1.8
#     ) +
#     scale_fill_gradient2(
#       low = "firebrick", mid = "white", high = "steelblue",
#       midpoint = 0, na.value = "grey92", name = "Corr.", limits = c(-1, 1)
#     ) +
#     scale_x_continuous(
#       breaks = label_lookup$center, labels = label_lookup$label,
#       expand = c(0, 0)
#     ) +
#     scale_y_continuous(
#       breaks = label_lookup$center, labels = label_lookup$label,
#       expand = c(0, 0)
#     ) +
#     theme_classic(base_size = 7) +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
#       axis.text.y = element_text(size = 5),
#       panel.grid = element_blank(),
#       legend.position = "bottom"
#     ) +
#     labs(x = NULL, y = NULL) +
#     coord_fixed()
# }

# make_heatmap_paginated <- function(
#   fill_scale  = NULL,
#   caption_str = NULL,
#   mk_df       = NULL,
#   mk_labels   = NULL,
#   ncol        = 4,
#   nrow        = 5
# ) {
#   function(page_data, page_sites, p, title_str) {

#     n_sites <- length(page_sites)

#     site_panel <- function(s, idx) {
#       row <- ceiling(idx / ncol)
#       col <- ((idx - 1) %% ncol) + 1

#       # bottom row = last actual row this page uses, OR last row of the grid
#       is_bottom <- (row == ceiling(n_sites / ncol)) || (idx + ncol > n_sites)
#       is_left   <- (col == 1)

#       if (!is.null(mk_df) && s == "MK") {
#         plt <- heatmap_plot_blocks(mk_df, mk_labels) + labs(title = "MK")
#       } else {
#         plt <- heatmap_plot(dplyr::filter(page_data, site == s), fill_scale = fill_scale) +
#           labs(title = s)
#       }

#       plt <- plt + theme(plot.title = element_text(size = 6, face = "bold"))

#       if (!is_left)   plt <- plt + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
#       if (!is_bottom) plt <- plt + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#       plt
#     }

#     panels <- purrr::imap(page_sites, site_panel)

#     patchwork::wrap_plots(panels, ncol = ncol, nrow = nrow) +
#       patchwork::plot_layout(guides = "collect") +
#       patchwork::plot_annotation(title = title_str, caption = caption_str)
#   }
# }

# # build once
# label_lookup_mk <- block_positions("MK") |>
#   dplyr::mutate(
#     label = dplyr::case_when(
#       block == "behav" ~ "Behav.",
#       block == "info"  ~ "Info.",
#       TRUE              ~ short_to_label_short[block]
#     )
#   )

# # ---- Site-wise and aggregate heatmaps ----
# old_theme <- theme_get()
# theme_set(theme_classic(base_size = 7) + theme(legend.position = "bottom"))
# make_pdf(
#   "edge_overlap_heatmap_alt.pdf",
#   width      = 12,
#   height     = ceiling(facets_per_page / 6) * 3.5,
#   data_sites = plot_df_sites,
#   plot_fn    = make_heatmap_paginated(
#     mk_df     = plot_df_mk,
#     mk_labels = label_lookup_mk
#     #caption_str = paste0("Upper triangle: correlation of edge weights  |  ",
#     #                     "Lower triangle: correlation of binary edges")
#   ),
#   title_str  = NULL#"Edge overlap across network layers"
# )
# theme_set(old_theme)
