## ============================================================================
## 11_kinship.R
## ----------------------------------------------------------------------------
## Kinship plots.
##
## Relates support-network ties to kinship: what proportion of support
## partners are close kin (and vice versa), how kin counts relate to
## reachability/network degree, and how this varies with local exogamy.
##
## Outputs:
##   Figures/DescriptiveStats/kin_support_violin.png
##   Figures/DescriptiveStats/prop_kin_vs_cores_kin.png
##   Figures/DescriptiveStats/supp_vs_cores_kin.png
##   Figures/DescriptiveStats/prop_extended_kin.png
##   Figures/DescriptiveStats/prop_kin_by_endogamy.png
##
## Depends on: 00_setup.R
## ============================================================================

source("net_descriptives/00_setup.R")

require(data.table)

# ── Kinship Plots ────────────────────────────────────────────────

## Note: this is created in the kinship.R script
kin_supp <- read.csv(file.path(
  EndowDerivedData,
  "su_kin_support.csv"
))

# recode 'collapse_sum' to 'sum' to align
kin_supp <- kin_supp |>
  dplyr::mutate(
    network = dplyr::recode(
      network,
      collapse_sum = "sum"
    )
  )

kin_supp_summary <- kin_supp |>
  filter(!network %in% c("behav", "info", "govt", "external")) |>
  group_by(site, network) |>
  dplyr::summarise(
    mean_prop_kin = mean(prop_connected_kin, na.rm = TRUE),
    mean_prob_kin = mean(prob_kin, na.rm = TRUE),
    .groups = "drop_last"
  ) |>
  dplyr::mutate(
    network_label = factor(
      short_to_label[network],
      levels = rev(short_to_label[tie_short[-c(13,17)]])
    ),
    network_label_short = factor(
      short_to_label_short[network],
      levels = short_to_label_short[tie_short[-c(13,17)]]
    )
  )

# CSV:
kin_supp_summary |>
  dplyr::select(site, network_label_short, mean_prop_kin, mean_prob_kin) |>
  dplyr::mutate(dplyr::across(c(mean_prop_kin, mean_prob_kin), \(x) round(x, 3))) |>
  dplyr::arrange(site, network_label_short) |>
  setNames(c(
    "Community", "Network Layer",
    "Average Proportion of Support Partners Who Are Kin",
    "Average Proportion of Kin Who Are Support Partners")) |>
  readr::write_csv(file.path(csvpath,"kin_support_summary.csv"))

kin_supp_summary_long <- kin_supp_summary |>
  dplyr::select(-network_label_short) |>
  tidyr::pivot_longer(
    cols = c(mean_prop_kin, mean_prob_kin),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    metric_label = factor(
      metric,
      levels = c("mean_prop_kin", "mean_prob_kin"),
      labels = c(
        "Support Partners Who Are Kin",
        "Kin Who Are Support Partners"
      )
    )
  )

p <- make_violin_plot(
  kin_supp_summary_long, x_var = "network_label", y_var = "value", group_var = "network_label",
  palette = setNames(tie_meta$colour, tie_meta$label),
  facet_var = "metric_label",
  y_label = "Average Proportion"
)
ggsave(
  file.path(
    EndowDropbox,
    "Figures",
    "DescriptiveStats",
    "kin_support_violin.png"
  ),
  plot = p,
  width = 12,
  height = 6
)

kin_count_summary <- kin_supp |>
  filter(network == "sum") |>
  group_by(site) |>
  dplyr::summarise(
    median_kin = median(kins, na.rm = TRUE),
    mean_kin = mean(kins, na.rm = TRUE),
    sd_kins = sd(kins, na.rm = TRUE),
    .groups = "drop")

summary(kin_count_summary)

p <- kin_supp |>
  filter(network == "sum") |>
  ggplot(aes(x = kins, y = prop_connected_kin)) +
  geom_point(
    alpha = 0.2,
    colour = "grey"
  ) +
  geom_smooth(
    method = lm,
    formula = y ~ splines::bs(x, 3),
    se = TRUE,
    alpha = 0.2,
    fill = "blue"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(.~ site, scales = "free_x") +
  scale_x_continuous(breaks = integer_breaks()) +
  labs(
    title    = NULL,
    x        = "Count of Close Kin Sharing Units",
    y        = "Proportion of Support Partners Who Are Kin",
  ) +
  theme_classic(base_size = 11)
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "prop_kin_vs_cores_kin.png"),
  plot = p,
  width  = 9,
  height = 10
  )

p <- kin_supp |>
  filter(network == "sum") |>
  ggplot(aes(x = kins, y = connection_network_degree)) +
  geom_point(
    alpha = 0.2,
    colour = "grey"
  ) +
  geom_smooth(
    method = lm,
    formula = y ~ splines::bs(x, 3),
    se = TRUE,
    alpha = 0.2,
    fill = "blue"
  ) +
  facet_wrap(.~ site, scales = "free") +
  scale_x_continuous(breaks = integer_breaks()) +
  scale_y_continuous(breaks = integer_breaks()) +
  labs(
    title    = NULL,
    x        = "Count of Close Kin Sharing Units",
    y        = "Count of Support Partners",
  ) +
  theme_classic(base_size = 11)
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "supp_vs_cores_kin.png"),
  plot = p,
  width  = 9,
  height = 10
  )


prop_sus_reachable <- function(g, unit_attr = "name") {
  #g_und <- igraph::as.undirected(g, mode = "collapse")

  sus <- igraph::vertex_attr(g, unit_attr)
  total_sus <- unique(sus)

  comp <- igraph::components(g, mode = "weak")
  membership <- comp$membership

  purrr::map_dfr(seq_len(igraph::vcount(g)), function(v) {
    own_su <- sus[v]
    own_comp <- membership[v]

    reachable_nodes <- which(membership == own_comp)
    reachable_sus <- unique(sus[reachable_nodes])
    reachable_sus <- setdiff(reachable_sus, own_su)

    n_other_sus_total <- length(setdiff(total_sus, own_su))

    tibble::tibble(
      su_id = own_su,
      n_reachable_sus = length(reachable_sus),
      n_other_sus_total = n_other_sus_total,
      prop_sus_reachable = ifelse(
        n_other_sus_total == 0, NA_real_,
        length(reachable_sus) / n_other_sus_total
      )
    )
  })
}

reach_by_network <- purrr::map_dfr(
  names(su_nets_expanded),
  function(site_name) {
    purrr::map_dfr(
      tie_short,
      function(net_name) {
        g <- su_nets_expanded[[site_name]][[net_name]]

        if (is.null(g)) return(NULL)

        prop_sus_reachable(g, unit_attr = "name") |>
          dplyr::mutate(site = site_name, network = net_name)
      }
    )
  }
)

p <- reach_by_network |>
  filter(network == "primary_kin") |>
  ggplot(aes(prop_sus_reachable)) +
  geom_histogram(bins = 30) +
#  coord_cartesian(xlim = c(0, 1)) +
 facet_wrap(.~ site, scales = "free_y") +
  labs(
    title    = NULL,
    x        = "Proportion of Sharing Units that are Extended Kin",
    y        = "Count",
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "prop_extended_kin.png"),
  plot = p,
  width  = 9,
  height = 10
  )


dfl <- read.csv(file.path(EndowDerivedData, "DataForLASSO.csv"))

p <- kin_supp_summary_long |>
  dplyr::left_join(
    dfl |> dplyr::select(site_name, LocalExogamy),
    by = c("site" = "site_name")
  ) |>
  filter(metric == "mean_prop_kin") |>
  filter(network == "sum") |>
  ggplot(aes(value, LocalExogamy, fill = LocalExogamy)) +
  geom_violin(
    alpha = 0.4,
    linewidth = 0.5,
    scale = "width",
    draw_quantiles = 0.5
  ) +
  geom_label(
    aes(label = site),
    family = "roboto",
    position = ggplot2::position_jitter(height = 0.2),
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    #    colour = "white",
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
  scale_fill_viridis_d(guide = "none") +
  labs(
    title = NULL,
    x = "Average Proportion of Support Partners Who Are Kin",
    y = NULL,
  ) +
  theme_classic(base_size = 11)
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "prop_kin_by_endogamy.png"),
  plot = p,
  width  = 9,
  height = 6
  )

p <- kin_count_summary |>
  dplyr::left_join(
    dfl |> dplyr::select(site_name, LocalExogamy),
    by = c("site" = "site_name")
  ) |>
  ggplot(aes(mean_kin, LocalExogamy, fill = LocalExogamy)) +
  geom_violin(
    alpha = 0.4,
    linewidth = 0.5,
    scale = "width",
    draw_quantiles = 0.5
  ) +
  geom_label(
    aes(label = site),
    position = ggplot2::position_jitter(height = 0.2),
    size = 2,
    alpha = 0.6,
    linewidth = 0.3,
#    colour = "white",
    label.padding = unit(0.15, "lines"),
    label.r = unit(0.5, "lines")
  ) +
  scale_fill_viridis_d(guide = "none") +
  labs(
    title    = NULL,
    x        = "Average Count of Co-Resident Kin",
    y        = NULL,
  ) +
  theme_classic(base_size = 11)


# create dyad data frame
nets_to_use <- tie_meta$short[c(11, 1:8, 12)]

dyad_list <- lapply(site_codes, function(site) {
  site_nets <- su_nets_expanded[[site]]

  g0 <- site_nets[[nets_to_use[1]]]
  A0 <- as_adjacency_matrix(g0)
  ids <- rownames(A0)

  df <- data.frame(site = site, i = ids[row(A0)], j = ids[col(A0)])

  for (net in nets_to_use) {
    if (is.null(site_nets[[net]])) {
      next
    } # skip missing networks
    df[[net]] <- as.vector(as_adjacency_matrix(site_nets[[net]]))
  }

  # drop self-loops
  df <- df[df$i != df$j, ]

  # as upstream in 00_setup drop all unsampled SUs, this isn't necessary
  # mask layers to NA where i is not su_sampled
  #unsampled <- su_meta[[site]]$su_id[su_meta[[site]]$su_sampled == FALSE]

  #layer_cols <- intersect(nets_to_use, names(df))
  #df[df$i %in% unsampled, layer_cols] <- NA

  df
})

dyads_df <- dplyr::bind_rows(dyad_list)

dyads_dt <- data.table::as.data.table(dyads_df)

kin_cov <- dyads_dt[,
  .(
    n_su = data.table::uniqueN(i),
    n_pairs = sum(!is.na(primary_kin)),
    n_kin_pairs = sum(primary_kin, na.rm = TRUE),
    prop_dyads_kin = mean(primary_kin, na.rm = TRUE),
    n_su_with_any_kin = data.table::uniqueN(i[primary_kin %in% 1])
  ),
  by = site
]
kin_cov[, prop_su_with_any_kin := n_su_with_any_kin / n_su]
kin_cov[,
  kin_data_usable := n_kin_pairs >= 20 & prop_su_with_any_kin >= 0.10
]
# only one deemed 'unusable' by this count (FF), so ignoring; it's close enough.

## reliance quantities per sharing unit, for the aggregate support network
rel <- dyads_dt[
  !is.na(sum),
  .(
    n_alters_possible = .N,
    n_kin_available = sum(primary_kin),
    n_kin_partners = sum(sum * primary_kin),
    prop_available_kin = mean(primary_kin),
    out_degree = sum(sum),
    n_nonkin_available = sum(1 - primary_kin),
    n_nonkin_partners = sum(sum * (1 - primary_kin))
  ),
  by = .(site, i)
]
setnames(rel, "i", "su_id")
rel[, `:=`(
  ## (1) share of unit's support partners who are primary kin
  prop_partners_kin = ifelse(
    out_degree > 0,
    n_kin_partners / out_degree,
    NA_real_
  ),
  ## (2) share of unit's available primary kin units that are support partners
  prop_avail_kin_used = ifelse(
    n_kin_available > 0,
    n_kin_partners / n_kin_available,
    NA_real_
  ),
  ## comparison quantity: same for non-kin
  prop_avail_nonkin_used = ifelse(
    n_nonkin_available > 0,
    n_nonkin_partners / n_nonkin_available,
    NA_real_
  )
)]

site_layer_summary <- rbindlist(lapply(nets_to_use[-1], function(network) {
  d <- dyads_dt[!is.na(get(network))]
  if (nrow(d) == 0) {
    return(NULL)
  }
  out <- d[,
    {
      n_pairs <- .N
      n_kin <- sum(primary_kin)
      n_edge <- sum(get(network))
      n_kin_edge <- sum(get(network) * primary_kin)
      p_kin <- if (n_kin > 0) n_kin_edge / n_kin else NA_real_
      p_nonkin <- if (n_pairs - n_kin > 0) {
        (n_edge - n_kin_edge) / (n_pairs - n_kin)
      } else {
        NA_real_
      }
      .(
        n_pairs = n_pairs,
        n_kin_pairs = n_kin,
        n_edges = n_edge,
        n_kin_edges = n_kin_edge,
        avail_kin_share = n_kin / n_pairs,
        partners_kin_share = if (n_edge > 0) n_kin_edge / n_edge else NA_real_,
        p_edge_kin = p_kin,
        p_edge_nonkin = p_nonkin,
        enrichment = if (n_edge > 0 && n_kin > 0) {
          (n_kin_edge / n_edge) / (n_kin / n_pairs)
        } else {
          NA_real_
        },
        risk_ratio = if (!is.na(p_kin) && !is.na(p_nonkin) && p_nonkin > 0) {
          p_kin / p_nonkin
        } else {
          NA_real_
        },
        odds_ratio = if (
          !is.na(p_kin) &&
            !is.na(p_nonkin) &&
            p_kin > 0 &&
            p_kin < 1 &&
            p_nonkin > 0 &&
            p_nonkin < 1
        ) {
          (p_kin / (1 - p_kin)) / (p_nonkin / (1 - p_nonkin))
        } else {
          NA_real_
        }
      )
    },
    by = site
  ]
  out[, layer := network]
  out
}))

site_layer_summary[,
  domain := short_to_pair[layer]
]

tabC1 <- site_layer_summary[,
  {
    ## community-level (equal weight per community)
    .(
      n_communities = .N,
      avail_kin_share_mean = mean(avail_kin_share),
      avail_kin_share_med = median(avail_kin_share),
      partners_kin_share_mean = mean(partners_kin_share, na.rm = TRUE),
      partners_kin_share_med = median(partners_kin_share, na.rm = TRUE),
      p_edge_kin_mean = mean(p_edge_kin, na.rm = TRUE),
      p_edge_kin_med = median(p_edge_kin, na.rm = TRUE),
      p_edge_nonkin_mean = mean(p_edge_nonkin, na.rm = TRUE),
      p_edge_nonkin_med = median(p_edge_nonkin, na.rm = TRUE),
      enrichment_med = median(enrichment, na.rm = TRUE),
      enrichment_q25 = quantile(enrichment, .25, na.rm = TRUE),
      enrichment_q75 = quantile(enrichment, .75, na.rm = TRUE),
      risk_ratio_mean = mean(risk_ratio, na.rm = TRUE),
      risk_ratio_med = median(risk_ratio, na.rm = TRUE),
      odds_ratio_mean = mean(odds_ratio, na.rm = TRUE),
      odds_ratio_med = median(odds_ratio, na.rm = TRUE),
      ## pooled over dyads (large communities dominate)
      pooled_avail_kin_share = sum(n_kin_pairs) / sum(n_pairs),
      pooled_partners_kin_share = sum(n_kin_edges) / sum(n_edges),
      pooled_p_edge_kin = sum(n_kin_edges) / sum(n_kin_pairs),
      pooled_p_edge_nonkin = (sum(n_edges) - sum(n_kin_edges)) /
        (sum(n_pairs) - sum(n_kin_pairs))
    )
  },
  by = layer
]
tabC1[, pooled_enrichment := pooled_partners_kin_share / pooled_avail_kin_share]
tabC1[, pooled_risk_ratio := pooled_p_edge_kin / pooled_p_edge_nonkin]
tabC1[, layer_label := short_to_label[layer]]
setcolorder(tabC1, c("layer", "layer_label"))


## Cross-layer summary table showing median and then mean and SD for the kin measures
kin_stat_labels <- c(
  #avail_kin_share = "Co-Res. Kin Share",
  partners_kin_share = "Prop Partners Kin",
  p_edge_kin = "P(Tie $\\mid$ Kin)",
  p_edge_nonkin = "P(Tie $\\mid$ Non-kin)",
  #enrichment = "Enrichment",
  risk_ratio = "Risk Ratio",
  odds_ratio = "Odds Ratio"
)

xlayer_summary_stats <- site_layer_summary |>
  group_by(layer) |>
  dplyr::select(layer, all_of(names(kin_stat_labels))) |>
  dplyr::summarise(
    across(
      all_of(names(kin_stat_labels)),
      list(
        median = \(x) median(x, na.rm = TRUE),
        mean = \(x) mean(x, na.rm = TRUE),
        sd = \(x) sd(x, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    layer = factor(
      short_to_label[layer],
      levels = short_to_label[tie_meta$short[c(11, 1:8, 12)]]
    )
  ) |>
  dplyr::arrange(layer)

xlayer_summary_long <- purrr::map_dfr(names(kin_stat_labels), function(v) {
  tibble::tibble(
    layer = xlayer_summary_stats$layer,
    stat = kin_stat_labels[[v]],
    Median = as.character(round(
      xlayer_summary_stats[[paste0(v, "_median")]],
      2
    )),
    `Mean (SD)` = paste0(
      round(xlayer_summary_stats[[paste0(v, "_mean")]], 2),
      " (",
      round(xlayer_summary_stats[[paste0(v, "_sd")]], 2),
      ")"
    )
  )
})

xlayer_summary_wide <- xlayer_summary_long |>
  tidyr::pivot_longer(
    cols = c(Median, `Mean (SD)`),
    names_to = "type",
    values_to = "value"
  ) |>
  dplyr::mutate(
    type = factor(type, levels = c("Median", "Mean (SD)")),
    stat = factor(stat, levels = kin_stat_labels)
  ) |>
  tidyr::pivot_wider(names_from = layer, values_from = value) |>
  dplyr::arrange(stat, type) |>
  dplyr::select(-type)

xlayer_summary_wide$stat <- as.character(xlayer_summary_wide$stat)
xlayer_summary_wide$stat[duplicated(xlayer_summary_wide$stat)] <- ""

xlayer_summary_wide |>
  kableExtra::kbl(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    escape = FALSE,
    align = c("l", rep("c", 9)),
    col.names = c("", tie_meta$label_short[c(1:8, 12)]),
    caption = "Median (above) and mean and standard deviation (below) of kin availability and tie probability measures for each network layer. `Prop Partners Kin' is the proportion of connected dyads (i.e., partners) that are kin. `P(Tie $\\mid$ Kin)' is the probability that a kin dyad is connected in the layer. `P(Tie $\\mid$ Non-kin)' is the probability that a non-kin dyad is connected in the layer. `Risk Ratio' is the ratio of P(Tie $\\mid$ Kin) to P(Tie $\\mid$ Non-kin). `Odds Ratio' shows the odds of a tie existing between kin, relative to the odds of an tie existing between non-kin.",
    label = "global-kin-support",
    linesep = ""
  ) |>
  kableExtra::kable_styling(
    latex_options = c("hold_position", "repeat_header"),
    font_size = 8
  ) |>
  writeLines(file.path(
    tablepath,
    "global-kin-support.tex"
  ))


## ---------------------------------------------------------------------------
##  Figure: Availability and composition are not the same thing
## ---------------------------------------------------------------------------
agg <- site_layer_summary |> dplyr::filter(layer == "sum")

p1 <- ggplot(agg, aes(x = avail_kin_share, y = partners_kin_share)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_label(
    aes(label = site),
    family = "roboto",
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  labs(
    x = "Availability: share of dyads that are close kin",
    y = "Composition: share of support ties that are close kin ties"#,
    #title = "A. Composition far exceeds availability in every community",
    #subtitle = "Dashed line = random mixing with respect to kinship"
  ) +
  theme_classic()

p2 <- ggplot(agg, aes(x = p_edge_nonkin, y = p_edge_kin)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_label(
    aes(label = site),
    family = "roboto",
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  labs(
    x = "P(support tie | non-kin dyad)",
    y = "P(support tie | close kin dyad)"#,
    #title = "B. Reliance: kin dyads carry support far more often",
    #subtitle = "Dashed line = no kin/non-kin difference"
  ) +
  theme_classic()

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "kin_availability_vs_reliance.png"),
  plot = p1 + p2,
  width = 11,
  height = 5
)


## ---------------------------------------------------------------------------
##  Figure: Distribution of kin availability across communities and
##  sharing units
## ---------------------------------------------------------------------------
fold <- max(agg$avail_kin_share) / min(agg$avail_kin_share)
p3 <- ggplot(agg, aes(x = avail_kin_share)) +
  geom_histogram(bins = 18, fill = "grey", colour = "white") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Share of dyads that are close kin",
    y = "Communities",
    #   title = sprintf(
    #     paste(
    #       "A. Kin availability spans a %.0f-fold range across communities",
    #       "(%.1f%% to %.1f%% of possible pairs)"
    #     ),
    #     fold,
    #     100 * min(agg$avail_kin_share),
    #     100 * max(agg$avail_kin_share)
    #   )
  ) +
  theme_classic()

p4 <- rel %>%
  dplyr::mutate(site = forcats::fct_reorder(site, n_kin_available, median)) %>%
  ggplot(aes(x = site, y = n_kin_available)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, fill = "grey92") +
  labs(
    x = "Community",
    y = "Count of co-resident close kin sharing units"#,
    #title = "B. Within-community spread in the number of locally available close kin units"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))

ggsave(
  file.path(
    EndowDropbox,
    "Figures",
    "DescriptiveStats",
    "kin_availability.png"
  ),
  plot = p3 / p4,
  width = 10,
  height = 7
)

## ---------------------------------------------------------------------------
##  Table: sharing-unit level availability and use
## ---------------------------------------------------------------------------
tabC3 <- rel[, .(
  n_sharing_units = .N,
  mean_kin_available = mean(n_kin_available, na.rm = TRUE),
  median_kin_available = median(n_kin_available, na.rm = TRUE),
  pct_with_no_close_kin = 100 * mean(n_kin_available == 0, na.rm = TRUE),
  mean_prop_partners_kin = mean(prop_partners_kin, na.rm = TRUE),
  mean_prop_avail_kin_used = mean(prop_avail_kin_used, na.rm = TRUE),
  mean_prop_avail_nonkin_used = mean(prop_avail_nonkin_used, na.rm = TRUE)
)]
tabC3_bysite <- rel[,
  .(
    n_sharing_units = .N,
    mean_kin_available = mean(n_kin_available, na.rm = TRUE),
    median_kin_available = median(n_kin_available, na.rm = TRUE),
    pct_with_no_close_kin = 100 * mean(n_kin_available == 0, na.rm = TRUE),
    mean_prop_partners_kin = mean(prop_partners_kin, na.rm = TRUE),
    mean_prop_avail_kin_used = mean(prop_avail_kin_used, na.rm = TRUE),
    mean_prop_avail_nonkin_used = mean(prop_avail_nonkin_used, na.rm = TRUE)
  ),
  by = site
][order(site)]

tabC3[, site := "Pooled"]

tabC3_combined <- data.table::rbindlist(
  list(tabC3, tabC3_bysite),
  fill = TRUE
)

setcolorder(tabC3_combined, c("site", setdiff(names(tabC3_combined), "site")))


tabC3_combined |>
  dplyr::mutate(
    kin_avail = paste0(sprintf("%.2f", round(mean_kin_available, 2)), " (", median_kin_available, ")")
  ) |>
  dplyr::relocate(kin_avail, .after = n_sharing_units) |>
  dplyr::select(-mean_kin_available, -median_kin_available) |>
  kableExtra::kbl(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    escape = FALSE,
    align = c("l", rep("c", 6)),
    digits = 2,
    col.names = c("Community", "N", "Co-Res. Kin", "No Kin", "Partners Kin", "Named Kin", "Named Non-Kin"),
    caption = "Availability of close kin and ties to them in the aggregate network, across all communities (`Pooled') and separately by community. `N' is the number of sharing units. `Co-Res. Kin' is the mean and median (in brackets) of number of co-residing sharing units with close kin that each sharing unit has. `No Kin' is the percent of sharing units with no co-residing close kin. `Partners Kin' is the mean proportion of partners in the aggregate network who are close kin. `Named Kin' is the mean proportion of co-resident kin sharing units that are partners in the aggregate network. `Named Non-Kin' is equivalently the proportion of non-kin sharing units that are named as partners in the aggregate network.",
    label = "sitewise-kin-support",
    linesep = ""
  ) |>
  kableExtra::kable_styling(
    latex_options = c("hold_position", "repeat_header"),
    font_size = 8
  ) |>
  writeLines(file.path(
    tablepath,
    "sitewise-kin-support.tex"
  ))

# CSV
tabC3_combined |>
  dplyr::relocate(mean_kin_available, median_kin_available, .after = n_sharing_units) |>
  dplyr::mutate(dplyr::across(is.numeric, \(x) round(x, 3))) |>
  setNames(c(
    "Community", "N", "Co-Res. Kin Mean", "Co-Res. Kin Median",
    "No Kin", "Partners Kin", "Named Kin", "Named Non-Kin"
  )) |>
  readr::write_csv(file.path(csvpath, "sitewise_kin_support.csv"))
