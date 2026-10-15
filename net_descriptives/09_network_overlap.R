## ============================================================================
## 09_network_overlap.R
## ----------------------------------------------------------------------------
## Layer overlap (edge- and alter-level).
##
## Edge-level correlation (weighted and binary) between network layers, and
## alter-level (Jaccard) overlap of who gets nominated across layers -- both
## by community and pooled.
##
## Outputs:
##   Figures/Networks/edge_overlap_heatmap.pdf
##   Figures/Networks/edge_overlap_heatmap_pooled.png
##   Figures/Networks/alter_overlap_alterwise_heatmap.pdf
##   Figures/Networks/alter_overlap_alterwise_heatmap_pooled.png
##   Figures/Networks/alter_overlap_personwise_heatmap.pdf
##   Figures/Networks/alter_overlap_personwise_heatmap_pooled.png
##
## Depends on: 00_setup.R
## ============================================================================

source("net_descriptives/00_setup.R")

load(file.path(EndowDerivedData, "alters_summary.rdata"))
# this file is created by nomination_counts_composition; alternatively could source

# ---- Compute weighted and binary edge correlations across layers ----

compute_layer_overlap <- function(
  sites,
  net_types,
  nets_list,
  include_aggregate = FALSE,
  aggregate_site_label = "all"
) {
  site_results <- purrr::map_dfr(sites, function(site) {
    all_nodes <- V(nets_list[[site]][[net_types[1]]])$name
    stacked <- as.data.frame(setNames(
      lapply(net_types, function(nt) {
        as.vector(get_adj_weighted(nets_list[[site]][[nt]], all_nodes))
      }),
      net_types
    ))
    stacked_binary <- as.data.frame(lapply(stacked, function(v) (v > 0) * 1))

    pairs <- expand.grid(
      layer1 = net_types,
      layer2 = net_types,
      stringsAsFactors = FALSE
    )
    pairs$weighted <- mapply(
      function(l1, l2) {
        cor(stacked[[l1]], stacked[[l2]], use = "complete.obs")
      },
      pairs$layer1,
      pairs$layer2
    )
    pairs$binary <- mapply(
      function(l1, l2) {
        cor(stacked_binary[[l1]], stacked_binary[[l2]], use = "complete.obs")
      },
      pairs$layer1,
      pairs$layer2
    )
    pairs$site <- site
    pairs[, c("site", "layer1", "layer2", "weighted", "binary")]
  })

  if (!include_aggregate) {
    return(site_results)
  }

  agg_stacked <- as.data.frame(setNames(
    lapply(net_types, function(nt) {
      unlist(lapply(sites, function(site) {
        all_nodes <- V(nets_list[[site]][[net_types[1]]])$name
        as.vector(get_adj_weighted(nets_list[[site]][[nt]], all_nodes))
      }))
    }),
    net_types
  ))
  agg_binary <- as.data.frame(lapply(agg_stacked, function(v) (v > 0) * 1))

  pairs <- expand.grid(
    layer1 = net_types,
    layer2 = net_types,
    stringsAsFactors = FALSE
  )
  pairs$weighted <- mapply(
    function(l1, l2) {
      cor(agg_stacked[[l1]], agg_stacked[[l2]], use = "complete.obs")
    },
    pairs$layer1,
    pairs$layer2
  )
  pairs$binary <- mapply(
    function(l1, l2) {
      cor(agg_binary[[l1]], agg_binary[[l2]], use = "complete.obs")
    },
    pairs$layer1,
    pairs$layer2
  )
  pairs$site <- aggregate_site_label

  dplyr::bind_rows(
    site_results,
    pairs[, c("site", "layer1", "layer2", "weighted", "binary")]
  )
}

edge_overlap <- compute_layer_overlap(
  sites = site_codes,
  net_types = tie_meta$short[1:10],
  nets_list = su_nets_expanded,
  include_aggregate = TRUE
)

# MK exception, not doing anything with it now, as too awkward to incorporate
if ("MK" %in% names(su_nets)) {
  mk_edge_overlap <- compute_layer_overlap(
    sites = "MK",
    net_types = names(su_nets$MK)[1:7],
    nets_list = su_nets,
    include_aggregate = FALSE
  ) |>
    dplyr::mutate(
      dplyr::across(
        c(layer1, layer2),
        \(x) {
          dplyr::case_when(
            x == "e1" ~ "loan_req",
            x == "e2" ~ "loan_giv",
            x == "e3" ~ "exchange_req",
            x == "e4" ~ "exchange_giv",
            x == "e56a" ~ "agri",
            x == "e56b" ~ "market",
            x == "e78" ~ "advice",
            .default = x
          )
        }
      )
    )
}

# ---- Prepare plot data ----
# Lower triangle (row_idx < col_idx): same direction
# Upper triangle (row_idx > col_idx): reverse direction

plot_df <- edge_overlap |>
  dplyr::mutate(
    row_idx = match(layer1, net_type_order), # defined in 00_setup
    col_idx = match(layer2, net_type_order),
    value = dplyr::case_when(
      row_idx == col_idx ~ NA_real_,
      row_idx < col_idx ~ weighted,
      row_idx > col_idx ~ binary
    ),
    layer1_label = factor(
      short_to_label_short[layer1],
      levels = short_to_label_short[net_type_order]
    ),
    layer2_label = factor(
      short_to_label_short[layer2],
      levels = short_to_label_short[net_type_order]
    )
  )

plot_df_sites <- dplyr::filter(plot_df, site != "all")
plot_df_agg <- dplyr::filter(plot_df, site == "all")

# CSV:
plot_df |>
  dplyr::filter(layer1 != layer2) |>
  dplyr::select(site, layer1_label, layer2_label, weighted, binary) |>
  dplyr::mutate(dplyr::across(c(weighted, binary), \(x) round(x, 3))) |>
  dplyr::arrange(site != "all", site) |>
  setNames(c("Community", "Layer 1", "Layer 2", "Pearsons r Weighted", "Pearsons r Binary")) |>
  readr::write_csv(file.path(csvpath,"edge_overlap.csv"))

# ---- Heatmap ----
heatmap_plot <- function(df, fill_scale = NULL) {
  if (is.null(fill_scale)) {
    fill_scale <- scale_fill_gradient2(
      low = "firebrick",
      mid = "white",
      high = "steelblue",
      midpoint = 0,
      na.value = "grey92",
      name = "Corr.",
      limits = c(-1, 1)
    )
  }
  ggplot(df, aes(x = layer2_label, y = layer1_label, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(
      data = \(d) dplyr::filter(d, !is.na(value)),
      aes(label = sprintf("%.2f", value)),
      size = 1.8
    ) +
    fill_scale +
    theme_classic(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      strip.text = element_text(size = 6, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.caption = element_text(size = 6, colour = "grey40")
    ) +
    labs(x = NULL, y = NULL) +
    coord_fixed()
}

# ---- Plotting helpers ----
facets_per_page <- 20
n_pages <- ceiling(length(site_codes) / facets_per_page)

# ---- Plot functions ----
make_heatmap_paginated <- function(fill_scale = NULL, caption_str = NULL) {
  function(page_data, page_sites, p, title_str) {
    heatmap_plot(page_data, fill_scale = fill_scale) +
      facet_wrap(~site, ncol = 4, nrow = 5) +
      labs(title = title_str, caption = caption_str)
  }
}

# ---- Site-wise and aggregate heatmaps ----
make_pdf(
  "edge_overlap_heatmap.pdf",
  width      = 12,
  height     = ceiling(facets_per_page / 6) * 3.5,
  data_sites = plot_df_sites,
  plot_fn    = make_heatmap_paginated(
    #caption_str = paste0("Lower triangle: correlation of edge weights  |  ",
    #                     "Upper triangle: correlation of binary edges")
  ),
  title_str  = NULL#"Edge overlap across network layers"
)

make_png("edge_overlap_heatmap_pooled.png", width = 4, height = 4,
  heatmap_plot(plot_df_agg) +
    labs(
      title   = NULL, #"Edge overlap across network layers (pooled)",
      #caption = paste0("Lower triangle: correlation of edge weights  |  ",
      #                 "Upper triangle: correlation of binary edges")
    )
)


## Similar approach, but now for all nominations, not SU-SU nets

## Two ways of doing it:
# - for the *typical person*, how much do their two layers overlap? (person wise)
# - of all alters named in either layer, what share were named in both? (alter wise)
## The former weights each respondent equally, the latter effectively
## gives more weight to respondents who make more nominations.

# ---- Jaccard overlap for a pair of id vectors ----
jaccard <- function(a, b) {
  a <- unlist(a)
  b <- unlist(b)
  a <- a[!is.na(a)]
  b <- b[!is.na(b)]
  if (length(a) == 0 && length(b) == 0) {
    return(NA_real_)
  }
  length(intersect(a, b)) / length(union(a, b))
}

jaccard_scale <- scale_fill_gradient(
  low = "white",
  high = "steelblue",
  na.value = "grey92",
  name = "Jaccard",
  limits = c(0, 1)
)

# First approach:

# ---- Compute person-level Jaccard for every layer pair ----
compute_alter_overlap <- function(site_data, id_col) {
  site_data <- site_data |>
    dplyr::mutate(tie_label = tie_meta$label_short[match(tie, tie_meta$id)])

  wide <- site_data |>
    dplyr::select(personid, tie_label, {{ id_col }}) |>
    tidyr::nest(data = {{ id_col }}) |>
    dplyr::mutate(ids = purrr::map(data, ~ unlist(.x))) |>
    dplyr::select(personid, tie_label, ids) |>
    tidyr::pivot_wider(names_from = tie_label, values_from = ids)

  net_types <- setdiff(names(wide), "personid")

  pairs <- expand.grid(
    layer1 = net_types,
    layer2 = net_types,
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(layer1 != layer2)

  pairs |>
    dplyr::rowwise() |>
    dplyr::mutate(
      overlap = mean(
        purrr::map2_dbl(
          wide[[layer1]],
          wide[[layer2]],
          ~ if (!is.null(.x) && !is.null(.y)) jaccard(.x, .y) else NA_real_
        ),
        na.rm = TRUE
      )
    ) |>
    dplyr::ungroup()
}

# ---- Run across all sites ----
site_alter_overlap_pw <- alters_summary |>
  dplyr::filter(tie != 12) |>
  dplyr::group_by(site) |>
  tidyr::nest() |>
  dplyr::mutate(
    raw_overlap = purrr::map(data, compute_alter_overlap, id_col = alter_ids),
    su_overlap = purrr::map(data, compute_alter_overlap, id_col = alter_su_ids)
  ) |>
  dplyr::select(site, raw_overlap, su_overlap) |>
  tidyr::unnest(c(raw_overlap, su_overlap), names_sep = "_") |>
  dplyr::rename(
    layer1 = raw_overlap_layer1,
    layer2 = raw_overlap_layer2,
    overlap_raw = raw_overlap_overlap,
    overlap_su = su_overlap_overlap
  ) |>
  dplyr::select(site, layer1, layer2, overlap_raw, overlap_su)

# ---- Pooled aggregate: stack all persons across sites ----
pooled_alter_overlap_pw <- alters_summary |>
  dplyr::filter(tie != 12) |>
  dplyr::group_by(site) |>
  tidyr::nest() |>
  dplyr::pull(data) |>
  dplyr::bind_rows() |>
  (\(d) {
    dplyr::mutate(
      compute_alter_overlap(d, id_col = alter_ids) |>
        dplyr::rename(overlap_raw = overlap),
      overlap_su = compute_alter_overlap(d, id_col = alter_su_ids)$overlap
    )
  })() |>
  dplyr::mutate(site = "all") |>
  dplyr::select(site, layer1, layer2, overlap_raw, overlap_su)

alter_overlap_pw <- dplyr::bind_rows(site_alter_overlap_pw, pooled_alter_overlap_pw)

# CSV:
alter_overlap_pw |>
  dplyr::filter(layer1 != layer2) |>
  dplyr::select(site, layer1, layer2, overlap_raw, overlap_su) |>
  dplyr::mutate(dplyr::across(c(overlap_raw, overlap_su), \(x) round(x, 3))) |>
  dplyr::arrange(site != "all", site) |>
  setNames(c("Community", "Layer 1", "Layer 2", "Jaccard Indivs", "Jaccard SUs")) |>
  readr::write_csv(file.path(csvpath,"alter_overlap_personwise.csv"))

# ---- Plot data ----
label_order <- tie_meta$label_short[1:10]

complete_grid <- tidyr::expand_grid(
  site = c(site_codes, "all"),
  layer1 = label_order,
  layer2 = label_order
)

plot_df_alter_pw <- complete_grid |>
  dplyr::left_join(alter_overlap_pw, by = c("site", "layer1", "layer2")) |>
  dplyr::mutate(
    row_idx = match(layer1, label_order),
    col_idx = match(layer2, label_order),
    value = dplyr::case_when(
      row_idx == col_idx ~ NA_real_,
      row_idx < col_idx ~ overlap_raw,
      row_idx > col_idx ~ overlap_su
    ),
    layer1_label = factor(layer1, levels = label_order),
    layer2_label = factor(layer2, levels = label_order)
  )

plot_df_alter_sites_pw <- dplyr::filter(plot_df_alter_pw, site != "all")
plot_df_alter_agg_pw <- dplyr::filter(plot_df_alter_pw, site == "all")

# Second approach:

# ---- Compute summed intersection and union counts for each layer pair ----
compute_alter_overlap_counts <- function(site_data, id_col) {
  site_data <- site_data |>
    dplyr::mutate(tie_label = tie_meta$label_short[match(tie, tie_meta$id)])

  wide <- site_data |>
    dplyr::select(personid, tie_label, {{ id_col }}) |>
    tidyr::nest(data = {{ id_col }}) |>
    dplyr::mutate(ids = purrr::map(data, ~ unlist(.x))) |>
    dplyr::select(personid, tie_label, ids) |>
    tidyr::pivot_wider(names_from = tie_label, values_from = ids)

  net_types <- setdiff(names(wide), "personid")

  expand.grid(layer1 = net_types, layer2 = net_types,
              stringsAsFactors = FALSE) |>
    dplyr::filter(layer1 != layer2) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      keep  = list(!purrr::map_lgl(wide[[layer1]], is.null) &
                   !purrr::map_lgl(wide[[layer2]], is.null)),
      n_int = sum(purrr::map2_int(wide[[layer1]][keep], wide[[layer2]][keep],
                                  ~ length(intersect(.x, .y)))),
      n_uni = sum(purrr::map2_int(wide[[layer1]][keep], wide[[layer2]][keep],
                                  ~ length(union(.x, .y))))
    ) |>
    dplyr::ungroup() |>
    dplyr::select(layer1, layer2, n_int, n_uni)
}

# ---- Run across all sites ----
site_counts_aw <- alters_summary |>
  dplyr::filter(tie != 12) |>
  dplyr::group_by(site) |>
  tidyr::nest() |>
  dplyr::mutate(
    raw = purrr::map(data, compute_alter_overlap_counts, id_col = alter_ids),
    su  = purrr::map(data, compute_alter_overlap_counts, id_col = alter_su_ids)
  ) |>
  dplyr::select(site, raw, su) |>
  tidyr::unnest(c(raw, su), names_sep = "_") |>
  dplyr::transmute(
    site,
    layer1 = raw_layer1, layer2 = raw_layer2,
    int_raw = raw_n_int, uni_raw = raw_n_uni,
    int_su  = su_n_int,  uni_su  = su_n_uni
  )

site_alter_overlap_aw <- site_counts_aw |>
  dplyr::mutate(overlap_raw = int_raw / uni_raw,
                overlap_su  = int_su  / uni_su) |>
  dplyr::select(site, layer1, layer2, overlap_raw, overlap_su)

# ---- Pooled: sum the counts across sites ----
pooled_alter_overlap_aw <- site_counts_aw |>
  dplyr::group_by(layer1, layer2) |>
  dplyr::summarise(
    overlap_raw = sum(int_raw) / sum(uni_raw),
    overlap_su  = sum(int_su)  / sum(uni_su),
    .groups = "drop"
  ) |>
  dplyr::mutate(site = "all") |>
  dplyr::select(site, layer1, layer2, overlap_raw, overlap_su)

alter_overlap_aw <- dplyr::bind_rows(site_alter_overlap_aw, pooled_alter_overlap_aw)

# CSV:
alter_overlap_aw |>
  dplyr::filter(layer1 != layer2) |>
  dplyr::select(site, layer1, layer2, overlap_raw, overlap_su) |>
  dplyr::mutate(dplyr::across(c(overlap_raw, overlap_su), \(x) round(x, 3))) |>
  dplyr::arrange(site != "all", site) |>
  setNames(c("Community", "Layer 1", "Layer 2", "Jaccard Indivs", "Jaccard SUs")) |>
  readr::write_csv(file.path(csvpath,"alter_overlap_alterwise.csv"))

# ---- Plot data ----
plot_df_alter_aw <- complete_grid |>
  dplyr::left_join(alter_overlap_aw, by = c("site", "layer1", "layer2")) |>
  dplyr::mutate(
    row_idx = match(layer1, label_order),
    col_idx = match(layer2, label_order),
    value = dplyr::case_when(
      row_idx == col_idx ~ NA_real_,
      row_idx < col_idx ~ overlap_raw,
      row_idx > col_idx ~ overlap_su
    ),
    layer1_label = factor(layer1, levels = label_order),
    layer2_label = factor(layer2, levels = label_order)
  )

plot_df_alter_sites_aw <- dplyr::filter(plot_df_alter_aw, site != "all")
plot_df_alter_agg_aw <- dplyr::filter(plot_df_alter_aw, site == "all")

# The two approaches are very close:
# In text, will present the alter-wise version.

cmp <- dplyr::inner_join(
  alter_overlap_pw, alter_overlap_aw,
  by = c("site", "layer1", "layer2"),
  suffix = c("_pw", "_aw")
) |>
  dplyr::filter(layer1 < layer2) |>   # drop mirrored pairs (Jaccard is symmetric)
  tidyr::pivot_longer(
    dplyr::starts_with("overlap_"),
    names_to = c("id_type", ".value"),
    names_pattern = "overlap_(raw|su)_(pw|aw)"
  ) |>
  dplyr::mutate(diff = aw - pw)   # positive = alterwise higher

cmp |>
  dplyr::group_by(id_type, pooled = site == "all") |>
  dplyr::summarise(
    cor_pearson  = cor(pw, aw, use = "complete.obs"),
    cor_spearman = cor(pw, aw, use = "complete.obs", method = "spearman"),
    mean_diff    = mean(diff, na.rm = TRUE),
    mean_abs     = mean(abs(diff), na.rm = TRUE),
    max_abs      = max(abs(diff), na.rm = TRUE),
    .groups = "drop"
  )


# ---- Alter overlap plots ----
make_pdf(
  "alter_overlap_personwise_heatmap.pdf",
  width = 14,
  height = ceiling(facets_per_page / 6) * 3.5,
  data_sites = plot_df_alter_sites_pw,
  plot_fn = make_heatmap_paginated(
    fill_scale = jaccard_scale#,
    #caption_str = paste0(
    #  "Lower triangle: overlap of individual IDs  |  ",
    #  "Upper triangle: overlap of sharing unit IDs"
    #)
  ),
  title_str = NULL #"Alter overlap across network layers"
)

make_png(
  "alter_overlap_personwise_heatmap_pooled.png",
  width = 4,
  height = 4,
  heatmap_plot(plot_df_alter_agg_pw, fill_scale = jaccard_scale) +
    labs(
      title = NULL
    )
)

make_pdf(
  "alter_overlap_alterwise_heatmap.pdf",
  width = 14,
  height = ceiling(facets_per_page / 6) * 3.5,
  data_sites = plot_df_alter_sites_aw,
  plot_fn = make_heatmap_paginated(
    fill_scale = jaccard_scale#,
    #caption_str = paste0(
    #  "Lower triangle: overlap of individual IDs  |  ",
    #  "Upper triangle: overlap of sharing unit IDs"
    #)
  ),
  title_str = NULL #"Alter overlap across network layers"
)

make_png(
  "alter_overlap_alterwise_heatmap_pooled.png",
  width = 4,
  height = 4,
  heatmap_plot(plot_df_alter_agg_aw, fill_scale = jaccard_scale) +
    labs(
      title = NULL
    )
)
