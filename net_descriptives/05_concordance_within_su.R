## ============================================================================
## 05_concordance_within_su.R
## ----------------------------------------------------------------------------
## Within-sharing-unit concordance.
##
## For sharing units with more than one surveyed member, measures how much
## agreement there is between co-resident respondents' nominations (both at
## the individual and sharing-unit level), overall and by gender.
##
## Outputs:
##   Figures/DescriptiveStats/su_concordance_layerwise_personid.png
##   Figures/DescriptiveStats/su_concordance_layerwise_suid.png
##   Figures/DescriptiveStats/su_concordance_aggregate_both.png
##   Figures/DescriptiveStats/su_concordance_aggregate_gender_both.png
##
## Depends on: 00_setup.R
## ============================================================================

source("net_descriptives/00_setup.R")

load(file.path(EndowDerivedData, "alters_summary.rdata"))
# this file is created by nomination_counts_composition; alternatively could source

# ── Concordance across SU members ────────────────────────────────────────────

# derive multi_member_sus from the complete alters_summary
multi_member_sus <- alters_summary |>
  dplyr::filter(!site %in% c("MA", "KO", "MY")) |>
  dplyr::filter(!is.na(su_id)) |>
  dplyr::group_by(site, tie, su_id) |>
  dplyr::summarise(n_su_members_asked = dplyr::n_distinct(personid), .groups = "drop") |>
  dplyr::filter(n_su_members_asked > 1)

# restrict alters_summary to multi-member SUs
alters_summary_multi <- alters_summary |>
  dplyr::semi_join(multi_member_sus, by = c("site", "tie", "su_id"))

concordance <- alters_summary_multi |>
  dplyr::group_by(site, tie, su_id) |>
  dplyr::mutate(
    # --- alter personid concordance ---
    others_alter_ids    = purrr::map(personid, ~ unique(unlist(alter_ids[personid != .x]))),
    n_shared_ids        = purrr::map2_int(alter_ids, others_alter_ids, ~ length(intersect(.x, .y))),
    prop_shared_ids     = dplyr::if_else(alter_count > 0, n_shared_ids / alter_count, NA_real_),

    # --- alter su_id concordance ---
    others_alter_su_ids = purrr::map(personid, ~ unique(unlist(alter_su_ids[personid != .x]))),
    n_shared_su_ids     = purrr::map2_int(alter_su_ids, others_alter_su_ids, ~ length(intersect(.x, .y))),
    prop_shared_su_ids  = dplyr::if_else(alter_su_count > 0, n_shared_su_ids / alter_su_count, NA_real_)
  ) |>
  dplyr::ungroup() |>
  dplyr::left_join(multi_member_sus, by = c("site", "tie", "su_id")) |>
  dplyr::select(site, tie, personid, gender, su_id, n_su_members_asked,
                alter_count, n_shared_ids, prop_shared_ids,
                alter_su_count, n_shared_su_ids, prop_shared_su_ids)


concordance_cols <- c("n_su_members_asked", "alter_count", "n_shared_ids",
                      "prop_shared_ids", "alter_su_count", "n_shared_su_ids",
                      "prop_shared_su_ids")

# make sure we're dropping cases that aren't valid
# most already are
concordance <- concordance |>
  dplyr::mutate(
    excluded = (site %in% two_nine_ten & tie %in% c(2, 9, 10)) |
               (site %in% nine_ten & tie %in% c(9, 10)) |
               (site %in% nine & tie %in% c(9)) |
               (site %in% c("MA", "KO") & tie %in% 1:10) |
               (site == "MK" & tie %in% c(5, 6, 7, 8)) |
               (site == "BM")
  ) |>
  dplyr::filter(!excluded) |>
  dplyr::select(-excluded)


concordance_plot_data <- concordance |>
  dplyr::mutate(tie = factor(tie, levels = c(1:10, 12), labels = tie_label_short[c(1:10, 12)]))


concordance_tables <- concordance |>
  dplyr::group_by(tie, site) |>
  dplyr::summarise(
    n_multi_su          = dplyr::n_distinct(su_id),
    mean_n_respondents  = mean(n_su_members_asked),
    median_prop_ids     = median(prop_shared_ids, na.rm = TRUE),
    mean_prop_ids       = mean(prop_shared_ids, na.rm = TRUE),
    median_prop_su_ids  = median(prop_shared_su_ids, na.rm = TRUE),
    mean_prop_su_ids    = mean(prop_shared_su_ids, na.rm = TRUE),
    .groups = "drop"
  )

site_order <- concordance_tables |>
  dplyr::filter(tie == 12) |>
  dplyr::group_by(site) |>
  dplyr::arrange(median_prop_su_ids) |>
  dplyr::pull(site)

p1 <- concordance_plot_data |>
  dplyr::filter(!is.na(prop_shared_ids)) |>
  dplyr::filter(tie != "Aggregate") |>
  dplyr::mutate(site = factor(site, levels = site_order)) |>
  dplyr::semi_join(
    concordance_tables |>
      dplyr::filter(n_multi_su >= 10) |>
      dplyr::mutate(tie = factor(tie, levels = c(1:10, 12), labels = tie_label_short[c(1:10, 12)])),
    by = c("site", "tie")
  ) |>
  dplyr::mutate(
    tie = factor(tie, levels = tie_label_short[c(1:10, 12)])
  ) |>
  ggplot(aes(x = prop_shared_ids, y = site, fill = tie, colour = tie)) +
  geom_boxplot(
    alpha         = 0.4,
    linewidth     = 0.3,
    outlier.size  = 0.4,
    outlier.alpha = 0.3
  ) +
  facet_wrap(~ tie, nrow = 2, dir = "v") +
  scale_fill_manual(values = tie_colours_short[c(1:10, 12)], guide = "none") +
  scale_colour_manual(values = tie_colours_short[c(1:10, 12)], guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(
    x = "Proportion concordant individuals",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0),
    panel.spacing    = unit(0.5, "lines"),
    axis.text.y      = element_text(size = 6)
  )

p2 <- concordance_plot_data |>
  dplyr::filter(!is.na(prop_shared_ids)) |>
  dplyr::filter(tie == "Aggregate") |>
  dplyr::mutate(site = factor(site, levels = site_order)) |>
  dplyr::semi_join(
    concordance_tables |>
      dplyr::filter(n_multi_su >= 10) |>
      dplyr::mutate(tie = factor(tie, levels = c(1:10, 12), labels = tie_label_short[c(1:10, 12)])),
    by = c("site", "tie")
  ) |>
  dplyr::mutate(
    tie = factor(tie, levels = tie_label_short[c(1:10, 12)])
  ) |>
  ggplot(aes(x = prop_shared_ids, y = site, fill = tie, colour = tie)) +
  geom_boxplot(
    alpha         = 0.4,
    linewidth     = 0.3,
    outlier.size  = 0.4,
    outlier.alpha = 0.3
  ) +
  facet_wrap(~ tie, nrow = 2, dir = "v") +
  scale_fill_manual(values = tie_colours_short[c(1:10, 12)], guide = "none") +
  scale_colour_manual(values = tie_colours_short[c(1:10, 12)], guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(
    x = "Proportion concordant individuals",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0),
    panel.spacing    = unit(0.5, "lines"),
    axis.text.y      = element_text(size = 6)
  )

p <- p1 + p2 + patchwork::plot_layout(widths = c(5, 1), axes = "collect")

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "su_concordance_layerwise_personid.png"),
  plot = p,
  width  = 12,
  height = 8
)

p1 <- concordance_plot_data |>
  dplyr::filter(!is.na(prop_shared_su_ids)) |>
  dplyr::filter(tie != "Aggregate") |>
  dplyr::mutate(site = factor(site, levels = site_order)) |>
  dplyr::semi_join(
    concordance_tables |>
      dplyr::filter(n_multi_su >= 10) |>
      dplyr::mutate(tie = factor(tie, levels = c(1:10, 12), labels = tie_label_short[c(1:10, 12)])),
    by = c("site", "tie")
  ) |>
  dplyr::mutate(
    tie = factor(tie, levels = tie_label_short[c(1:10, 12)])
  ) |>
  ggplot(aes(x = prop_shared_su_ids, y = site, fill = tie, colour = tie)) +
  geom_boxplot(
    alpha         = 0.4,
    linewidth     = 0.3,
    outlier.size  = 0.4,
    outlier.alpha = 0.3
  ) +
  facet_wrap(~ tie, nrow = 2, dir = "v") +
  scale_fill_manual(values = tie_colours_short[c(1:10, 12)], guide = "none") +
  scale_colour_manual(values = tie_colours_short[c(1:10, 12)], guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(
    x = "Proportion concordant sharing units",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0),
    panel.spacing    = unit(0.5, "lines"),
    axis.text.y      = element_text(size = 6)
  )

p2 <- concordance_plot_data |>
  dplyr::filter(!is.na(prop_shared_su_ids)) |>
  dplyr::filter(tie == "Aggregate") |>
  dplyr::mutate(site = factor(site, levels = site_order)) |>
  dplyr::semi_join(
    concordance_tables |>
      dplyr::filter(n_multi_su >= 10) |>
      dplyr::mutate(tie = factor(tie, levels = c(1:10, 12), labels = tie_label_short[c(1:10, 12)])),
    by = c("site", "tie")
  ) |>
  dplyr::mutate(
    tie = factor(tie, levels = tie_label_short[c(1:10, 12)])
  ) |>
  ggplot(aes(x = prop_shared_su_ids, y = site, fill = tie, colour = tie)) +
  geom_boxplot(
    alpha         = 0.4,
    linewidth     = 0.3,
    outlier.size  = 0.4,
    outlier.alpha = 0.3
  ) +
  facet_wrap(~ tie, nrow = 2, dir = "v") +
  scale_fill_manual(values = tie_colours_short[c(1:10, 12)], guide = "none") +
  scale_colour_manual(values = tie_colours_short[c(1:10, 12)], guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(
    x = "Proportion concordant sharing units",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0),
    panel.spacing    = unit(0.5, "lines"),
    axis.text.y      = element_text(size = 6)
  )

p <- p1 + p2 + patchwork::plot_layout(widths = c(5, 1), axes = "collect")

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "su_concordance_layerwise_suid.png"),
  plot = p,
  width  = 12,
  height = 8
)


p_shared_ids_aggregate <- concordance_plot_data |>
  dplyr::filter(!is.na(prop_shared_ids)) |>
  dplyr::filter(tie == "Aggregate") |>
  dplyr::mutate(site = factor(site, levels = site_order)) |>
  dplyr::semi_join(
    concordance_tables |>
      dplyr::filter(n_multi_su >= 10) |>
      dplyr::mutate(tie = factor(tie, levels = c(1:10, 12), labels = tie_labels[c(1:10, 12)])),
    by = c("site", "tie")
  ) |>
  ggplot(aes(x = site, y = prop_shared_ids)) +
  geom_boxplot(outlier.size = 0.5, fill = "grey") +
  #geom_jitter(alpha = 0.1, width = 0.3) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(
    title    = NULL,
    subtitle = "Individuals",
    x        = NULL,
    y        = "Proportion concordant"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))

p_shared_su_ids_aggregate <- concordance_plot_data |>
  dplyr::filter(!is.na(prop_shared_su_ids)) |>
  dplyr::filter(tie == "Aggregate") |>
  dplyr::mutate(site = factor(site, levels = site_order)) |>
  dplyr::semi_join(
    concordance_tables |>
      dplyr::filter(n_multi_su >= 10) |>
      dplyr::mutate(tie = factor(tie, levels = c(1:10, 12), labels = tie_labels[c(1:10, 12)])),
    by = c("site", "tie")
  ) |>
  ggplot(aes(x = site, y = prop_shared_su_ids)) +
  geom_boxplot(outlier.size = 0.5, fill = "grey") +
  #geom_jitter(alpha = 0.1, width = 0.3) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(
    title    = NULL,
    subtitle = "Sharing units",
    x        = NULL,
    y        = "Proportion concordant"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))

p <- patchwork::wrap_plots(p_shared_ids_aggregate | p_shared_su_ids_aggregate) #+
#  patchwork::plot_annotation(
#    title = "Concordance in alter nominations within sharing units"
#  )
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "su_concordance_aggregate_both.png"),
  plot = p,
  width  = 9,
  height = 6
  )


# CSV to match general info of boxplots
# Same filtering as the plots: only site x tie combinations with >= 10 multi-respondent SUs
keep_combos <- concordance_tables |>
  dplyr::filter(n_multi_su >= 10) |>
  dplyr::mutate(
    tie = factor(tie, levels = c(1:10, 12), labels = tie_label_short[c(1:10, 12)])
  )

su_concordance_summary <- concordance_plot_data |>
  dplyr::semi_join(keep_combos, by = c("site", "tie")) |>
  dplyr::mutate(
    tie  = factor(tie, levels = tie_label_short[c(1:10, 12)])
  ) |>
  dplyr::group_by(tie, site) |>
  dplyr::summarise(
    #n_su_total        = dplyr::n(),
    #n_ids             = sum(!is.na(prop_shared_ids)),
    mean_shared_ids   = mean(prop_shared_ids, na.rm = TRUE),
    median_shared_ids = median(prop_shared_ids, na.rm = TRUE),
    sd_shared_ids     = sd(prop_shared_ids, na.rm = TRUE),
    #n_sus             = sum(!is.na(prop_shared_su_ids)),
    mean_shared_sus   = mean(prop_shared_su_ids, na.rm = TRUE),
    median_shared_sus = median(prop_shared_su_ids, na.rm = TRUE),
    sd_shared_sus     = sd(prop_shared_su_ids, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 3))) |>
  # groups with no non-missing values give NaN; make them proper NAs
  dplyr::mutate(dplyr::across(where(is.numeric), \(x) dplyr::if_else(is.nan(x), NA_real_, x))) |>
  dplyr::arrange(tie, site) |>
  dplyr::rename(
    `Network Prompt` = tie, Site = site,
    #`N SU` = n_su_total,
    #`N SU (Indivs)` = n_ids,
    `Mean Shared Indivs` = mean_shared_ids,
    `Median Shared Indivs` = median_shared_ids,
    `SD Shared Indivs` = sd_shared_ids,
    #`N SU (SUs)` = n_sus,
    `Mean Shared SUs` = mean_shared_sus,
    `Median Shared SUs` = median_shared_sus,
    `SD Shared SUs` = sd_shared_sus
  )

# CSV (numeric, full precision)
readr::write_csv(
  su_concordance_summary,
  file.path(csvpath, "within_su_concordance.csv")
)

## broken out by gender
p_shared_ids_aggregate_gender <- concordance_plot_data |>
  dplyr::filter(!is.na(prop_shared_ids)) |>
  dplyr::filter(tie == "Aggregate") |>
  dplyr::mutate(site = factor(site, levels = site_order)) |>
  dplyr::semi_join(
    concordance_tables |>
      dplyr::filter(n_multi_su >= 10) |>
      dplyr::mutate(tie = factor(tie, levels = c(1:10, 12), labels = tie_labels[c(1:10, 12)])),
    by = c("site", "tie")
  ) |>
  ggplot(aes(x = site, y = prop_shared_ids, fill = gender)) +
  geom_boxplot(outlier.size = 0.5) +
  #geom_jitter(alpha = 0.1, width = 0.1) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(3,"PRGn")[3], RColorBrewer::brewer.pal(3,"PRGn")[1]))  +
  labs(
    title    = NULL,
    subtitle = "Individuals",
    x        = NULL,
    y        = "Proportion concordant"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 6, angle = 45, hjust = 1)) +
  guides(colour = guide_legend(title = "Gender"))


p_shared_su_ids_aggregate_gender <- concordance_plot_data |>
  dplyr::filter(!is.na(prop_shared_su_ids)) |>
  dplyr::filter(tie == "Aggregate") |>
  dplyr::mutate(site = factor(site, levels = site_order)) |>
  dplyr::semi_join(
    concordance_tables |>
      dplyr::filter(n_multi_su >= 10) |>
      dplyr::mutate(tie = factor(tie, levels = c(1:10, 12), labels = tie_labels[c(1:10, 12)])),
    by = c("site", "tie")
  ) |>
  ggplot(aes(x = site, y = prop_shared_su_ids, fill = gender)) +
  geom_boxplot(outlier.size = 0.5) +
  #geom_jitter(alpha = 0.1, width = 0.1) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(3,"PRGn")[3], RColorBrewer::brewer.pal(3,"PRGn")[1]))  +
  labs(
    title    = NULL,
    subtitle = "Sharing units",
    x        = NULL,
    y        = "Proportion concordant"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 6, angle = 45, hjust = 1)) +
  guides(colour = guide_legend(title = "Gender"))

p <- (p_shared_ids_aggregate_gender | p_shared_su_ids_aggregate_gender) / patchwork::guide_area() +
  patchwork::plot_layout(guides = "collect", heights = c(20, 1)) #+
#  patchwork::plot_annotation(
#    title = "Concordance in alter nominations within sharing units"
#  )
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "su_concordance_aggregate_gender_both.png"),
  plot = p,
  width  = 9,
  height = 6
  )
