## ============================================================================
## 07_degree_distributions.R
## ----------------------------------------------------------------------------
## Degree distributions.
##
## CCDF degree distributions by network layer (all communities overlaid), and
## for the aggregate network alone.
##
## Outputs:
##   Figures/DescriptiveStats/degree_distributions.png
##   Figures/DescriptiveStats/degree_distribution_aggregate.png
##
## Depends on: 00_setup.R
## ============================================================================

source("net_descriptives/00_setup.R")

# ── Degree Distribution ────────────────────────────────────────────────

degree_distributions <- purrr::imap_dfr(su_nets_expanded, \(nets, site) {
  purrr::imap_dfr(nets[tie_short], \(g, net_type) {
    if (is.null(g) || !inherits(g, "igraph") || igraph::vcount(g) == 0L) return(NULL)
    data.frame(
      site     = site,
      network  = net_type,
      degree   = igraph::degree(g, mode = "all"),
      stringsAsFactors = FALSE
    )
  })
})

p <- degree_distributions |>
  dplyr::filter(!network %in% c("govt", "external", "primary_kin")) |>
  dplyr::mutate(
    network_label = factor(short_to_label[network], levels = short_to_label[tie_short[-c(13,17)]])
  ) |>
  dplyr::group_by(network_label, site) |>
  dplyr::arrange(degree, .by_group = TRUE) |>
  dplyr::mutate(ccdf = 1 - (rank(degree, ties.method = "max") / dplyr::n())) |>
  dplyr::ungroup() |>
  dplyr::filter(degree > 0, ccdf > 0) |>
  ggplot(aes(x = degree, y = ccdf, colour = site, linetype = site, group = site)) +
  geom_line(linewidth = 0.4, alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks(sides = "bl") +
  scale_colour_manual(values = site_colour_map) +
  scale_linetype_manual(values = site_linetype_map) +
  facet_wrap(~ network_label, scales = "free", ncol = 3) +
  labs(
    x     = "Degree",
    y     = "P(X ≥ x)",
    #title = "Degree distributions (CCDF)"
  ) +
  guides(colour = guide_legend(ncol = 20)) +
  theme_classic(base_size = 11) +
  theme(strip.text = element_text(size = 8), legend.position = "bottom")

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "degree_distributions.png"),
  plot = p,
  width  = 12,
  height = 8
)

p <- degree_distributions |>
  dplyr::filter(network == "sum") |>
  dplyr::group_by(site) |>
  dplyr::arrange(degree, .by_group = TRUE) |>
  dplyr::mutate(ccdf = 1 - (rank(degree, ties.method = "max") / dplyr::n())) |>
  dplyr::ungroup() |>
  dplyr::filter(degree > 0, ccdf > 0) |>
  ggplot(aes(x = degree, y = ccdf, colour = site, linetype = site, group = site)) +
  geom_line(linewidth = 0.4, alpha = 0.7) +
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 20, 50, 100),
    labels = scales::label_number(),
    limits = c(1,100)
  ) +
  scale_y_log10(labels = scales::label_number()) +
  annotation_logticks(sides = "bl") +
  scale_colour_manual(values = site_colour_map) +
  scale_linetype_manual(values = site_linetype_map) +
  labs(
    x     = "Degree",
    y     = "P(X ≥ x)",
    #title = "Degree distributions (CCDF)"
    colour = NULL,
    linetype = NULL
  ) +
  guides(colour = guide_legend(ncol = 9)) +
  theme_classic(base_size = 11) +
  theme(strip.text = element_text(size = 8), legend.position = "bottom")

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "degree_distribution_aggregate.png"),
  plot = p,
  width  = 6,
  height = 6
)