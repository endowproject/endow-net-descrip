## ============================================================================
## 03_respondent_coverage.R
## ----------------------------------------------------------------------------
## Respondent coverage, and its consequences.
##
## Site-level respondent/sharing-unit coverage, a handful of site-specific
## data-cleaning fixes to `people_observations` and `networks` (NP residency
## filter, MY joint-interview handling, MK question recoding), and the
## respondent-coverage Latex table.
## With that in hand, plots of how sampling influences key unweighted and
## weighted measures.
##
## Outputs:
##   Tables/DescriptiveStats/respondents_by_site.tex
##   Figures/DescriptiveStats/out_deg_by_num_surv.png
##   Figures/DescriptiveStats/weighted_out_deg_by_num_surv.png
##   Figures/DescriptiveStats/degree_vs_sampling.png
##   Figures/DescriptiveStats/weighted_degree_vs_sampling.png
##   Figures/DescriptiveStats/local_clust_by_num_surv.png
##   Figures/DescriptiveStats/weighted_local_clust_by_num_surv.png
##   Figures/DescriptiveStats/av_clust_vs_sampling.png
##
## Depends on: 00_setup.R
## ============================================================================

source("net_descriptives/00_setup.R")

su_df <- readRDS(file.path(EndowDerivedData, "su_master.rds"))
# for sampling coverage

load(file.path(EndowDerivedData, "network_summaries_long.rdata"))
# this file is created by network_summary_stats; alternatively could source

# ── Respondent counts ────────────────────────────────────────

su_summary <- su_meta %>%
  # Some sites have multiple family records per SU (e.g. CR). Count each
  # sharing unit once, including when averaging the sampling controls.
  purrr::map(~ dplyr::distinct(.x, su_id, .keep_all = TRUE)) %>%
  imap_dfr(
    ~ tibble(
      site = .y,
      su_count = nrow(.x),
      su_sampled = sum(.x$su_sampled, na.rm = TRUE),
      mean_su_size = mean(.x$su_size, na.rm = TRUE),
      mean_num_surveyed = mean(.x$num_surv, na.rm = TRUE),
      mean_num_surveyed_fq = mean(.x$num_surv_female_q, na.rm = TRUE),
      mean_num_surveyed_mq = mean(.x$num_surv_male_q, na.rm = TRUE),
      prop_sampled = su_sampled / su_count
    )
  )

# Save, so can be called in by variance decomp file
save(su_summary, file = file.path(EndowDerivedData, "su_sampling.rdata"))

respondent_counts <- imap_dfr(people_observations, ~ {
  if (.y %in% c("TN")) {
    tibble(site = .y)
  } else {
    .x |>
      dplyr::filter(surveyed == TRUE) |>
      dplyr::count(gender) |>
      dplyr::mutate(total = sum(n)) |>
      tidyr::pivot_wider(names_from = gender, values_from = n, values_fill = 0) |>
      dplyr::mutate(site = .y)
  }
}) |>
  relocate(site) |>
  dplyr::mutate(across(where(is.numeric), ~ if_else(site %in% c("MY","TN"), NA_integer_, replace_na(., 0L))))

respondent_counts <- respondent_counts |>
  left_join(su_summary, by = "site") |>
  dplyr::mutate(perc_sampled = prop_sampled * 100) |>
  dplyr::select(site, su_count, mean_su_size, perc_sampled, mean_num_surveyed, mean_num_surveyed_fq, mean_num_surveyed_mq, total, female, male) |>
arrange(site)

# intended respondents, based on collection-entry-metadata files
respondent_counts <- respondent_counts |>
  dplyr::mutate(intended_resp = case_when(
    site %in% c("AV", "AZ", "BY", "DJ", "IH", "NI", "NP", "PC", "TE", "TI", "TZ", "UP", "WH", "WL", "PQ") ~ "all adult residents",
    site %in% c("AR", "BD", "CP", "CR", "EC", "EK", "FF", "FJ", "KA", "KM", "KS", "MG", "MN", "RA", "TS", "YN") ~ "female head",
    site %in% c("AP", "BM", "MC", "MK", "MY", "PE", "PS", "SN", "TP", "TM", "ZS") ~ "male and female heads",
    site %in% c("AH", "HE", "HI", "LB", "SH", "VT", "YI") ~ "male head",
    site %in% c("KT", "PT") ~ "male or female head",
    site %in% c("KO", "MA") ~ "male and female heads (different prompts)",
    site %in% c("AN", "BF", "BH", "BZ", "DA", "OT", "SP", "SY") ~ "unclear",
    .default = NA_character_
  ))

respondent_counts <- respondent_counts[!respondent_counts$site %in% c(SITES_TO_DROP, "TN"),]

digits_vector <- ifelse(grepl("mean", colnames(respondent_counts)), 2, 0)

tab <- respondent_counts |>
  dplyr::select(-mean_num_surveyed_fq, -mean_num_surveyed_mq)

for (i in 2:8) {
  tab[[i]] <- round(tab[[i]], digits_vector[i])
}

csv_names <- c(
  "Community", "# SU", "Mean # Members", "% Sampled", "Mean Surveyed per SU",
  "# Surveyed", "# Women", "# Men", "Intended Respondents"
)

# CSV
tab |>
  setNames(csv_names) |>
  readr::write_csv(file.path(csvpath, "respondents_by_site.csv"))

# LaTeX table
tab |>
  kableExtra::kable(
  format = "latex",
  booktabs = TRUE,
  longtable = TRUE,
  caption = "Summary of the sharing unit coverage and respondents for the social support network survey in each community.",
  label = "respondents_by_site",
  col.names = c(
    "Community", "\\# SU", "Mean \\# Members", "\\% Sampled", "Mean Surveyed per SU", "\\# Surveyed",  "\\# Women", "\\# Men", "Intended Respondents"
  ),
  digits = digits_vector,
  linesep = "",
  escape = FALSE
) |>
  kableExtra::kable_styling(
      latex_options = "repeat_header",
      font_size     = 8
    ) |>
  kableExtra::column_spec(9, width = "5cm") |>
  kableExtra::save_kable(file.path(EndowDropbox, "Tables", "DescriptiveStats", "respondents_by_site.tex"))

# Note: code at the top to filter for CR and code here to filter for NP are a result of noting discrepancies between counts gathered without these.


# ── Sampling consequences ────────────────────────────────────────

p <- su_df |>
  dplyr::filter(num_surveyed != 0) |>
  dplyr::filter(network == "sum") |>
  ggplot(aes(x = num_surveyed, y = out_degree, colour = num_surveyed, group = num_surveyed)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(
    width = 0.2,
    alpha = 0.1,
  ) +
  facet_wrap(.~ site, scales = "free") +
  scale_x_continuous(breaks = integer_breaks()) +
  scale_y_continuous(breaks = integer_breaks()) +
  viridis::scale_color_viridis() +
  viridis::scale_fill_viridis() +
  labs(
    title    = NULL,
    x        = "Count of Respondents",
    y        = "Out Degree in Aggregate Network",
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "out_deg_by_num_surv.png"),
  plot = p,
  width  = 9,
  height = 10
  )

p <- su_df |>
  dplyr::filter(num_surveyed != 0) |>
  dplyr::filter(network == "sum") |>
  ggplot(aes(x = num_surveyed, y = out_degree_weighted, colour = num_surveyed, group = num_surveyed)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(
    width = 0.2,
    alpha = 0.1,
  ) +
  facet_wrap(.~ site, scales = "free") +
  scale_x_continuous(breaks = integer_breaks()) +
  scale_y_continuous(breaks = integer_breaks()) +
  viridis::scale_color_viridis() +
  viridis::scale_fill_viridis() +
  labs(
    title    = NULL,
    x        = "Count of Respondents",
    y        = "Weighted Out Degree in Aggregate Network",
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "weighted_out_deg_by_num_surv.png"),
  plot = p,
  width  = 9,
  height = 10
  )


respondent_counts <- respondent_counts |>
  mutate(
    intended_resp = factor(
      intended_resp,
      levels = c(
        "male or female head",
        "female head",
        "male head",
        "male and female heads (different prompts)",
        "male and female heads",
        "all adult residents"
      )
    )
  ) |>
  left_join(
    network_summaries_long |>
      dplyr::filter(network == "sum") |>
      dplyr::select(
        site,
        mean_degree,
        mean_degree_w,
        density,
        av_local_clust,
        global_clust,
        reciprocity,
        avg_path_length
      ),
    by = "site"
  )

out_deg_by_num_surv <- respondent_counts |>
  ggplot(aes(x = intended_resp, y = mean_degree)) +
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
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
  labs(
    x = "Intended Respondents",
    y = "Mean Degree in Aggregate Network"
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))

mean_deg_by_sampling <- respondent_counts |>
  ggplot(aes(x = mean_num_surveyed, y = mean_degree)) +
  geom_smooth(
    alpha = 0.2,
    fill = "blue",
    method = "lm"
  ) +
  geom_label(
    aes(label = site),
    family = "roboto",
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
  coord_cartesian(
    ylim = c(
      min(respondent_counts$mean_degree),
      max(respondent_counts$mean_degree)
    )
  ) +
  labs(
    x = "Mean Respondents Sampled Per Sharing Unit",
    y = "Mean Degree in Aggregate Network"
  ) +
  theme_classic()

p <- (out_deg_by_num_surv | mean_deg_by_sampling) + plot_layout(axis_titles = "collect")
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "degree_vs_sampling.png"),
  plot = p,
  width  = 12,
  height = 6
  )

w_out_deg_by_num_surv <- respondent_counts |>
  ggplot(aes(x = intended_resp, y = mean_degree_w)) +
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
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
  labs(
    x = "Intended Respondents",
    y = "Mean Weighted Degree in Aggregate Network"
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))

mean_w_deg_by_sampling <- respondent_counts |>
  ggplot(aes(x = mean_num_surveyed, y = mean_degree_w)) +
  geom_smooth(
    alpha = 0.2,
    fill = "blue",
    method = "lm"
  ) +
  geom_label(
    aes(label = site),
    family = "roboto",
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
  coord_cartesian(
    ylim = c(
      min(respondent_counts$mean_degree_w),
      max(respondent_counts$mean_degree_w)
    )
  ) +
  labs(
    x = "Mean Respondents Sampled Per Sharing Unit",
    y = "Mean Weighted Degree in Aggregate Network"
  ) +
  theme_classic()

p <- (w_out_deg_by_num_surv | mean_w_deg_by_sampling) + plot_layout(axis_titles = "collect")
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "weighted_degree_vs_sampling.png"),
  plot = p,
  width  = 12,
  height = 6
  )


## other standard measures
## (not currently saved)
p <- respondent_counts |>
  ggplot(aes(x = mean_num_surveyed, y = density)) +
  geom_smooth(
    alpha = 0.2,
    fill = "blue",
    method = "lm"
  ) +
  geom_label(
    aes(label = site),
    family = "roboto",
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
 labs(
   x = "Mean Respondents Sampled Per Sharing Unit",
   y = "Density in Aggregate Network"
 ) +
  theme_classic()

p <- respondent_counts |>
  ggplot(aes(x = mean_num_surveyed, y = reciprocity)) +
  geom_smooth(
    alpha = 0.2,
    fill = "blue",
    method = "lm"
  ) +
  geom_label(
    aes(label = site),
    family = "roboto",
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
   labs(
     x = "Mean Respondents Sampled Per Sharing Unit",
     y = "Reciprocity in Aggregate Network"
   ) +
  theme_classic()

p <- respondent_counts |>
  ggplot(aes(x = mean_num_surveyed, y = global_clust)) +
  geom_smooth(
    alpha = 0.2,
    fill = "blue",
    method = "lm"
  ) +
  geom_label(
    aes(label = site),
    family = "roboto",
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
   labs(
     x = "Mean Respondents Sampled Per Sharing Unit",
     y = "Global Clustering Coefficient in Aggregate Network"
   ) +
  theme_classic()


clust_list <- vector("list", length(su_nets_expanded))
names(clust_list) <- names(su_nets_expanded)

for (site in names(su_nets_expanded)) {

  g <- su_nets_expanded[[site]]$sum

  g_und <- as_undirected(
    g,
    mode = "collapse",
    edge.attr.comb = list(
      weight = "mean",          # default is "sum" -- avoid double-counting
      weight_per_resp = "mean"  # same
    )
  )

  local_clust   <- transitivity(g_und, type = "local", weights = NA)
  local_clust_w <- transitivity(g_und, type = "weighted")

  su_id <- V(g_und)$name
  if (is.null(su_id)) su_id <- as.character(seq_len(vcount(g_und)))

  clust_list[[site]] <- tibble(
    site         = site,
    su_id        = su_id,
    local_clust  = local_clust,
    local_clust_w = local_clust_w
  )
}

clust_df <- bind_rows(clust_list)

clust_df <- su_df |>
  dplyr::filter(num_surveyed != 0) |>
  dplyr::filter(network == "sum") |>
  dplyr::left_join(clust_df, by = c("site", "su_id"))

p <- clust_df |>
  ggplot(aes(x = num_surveyed, y = local_clust, colour = num_surveyed, group = num_surveyed)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(
    width = 0.2,
    alpha = 0.1,
  ) +
  facet_wrap(.~ site, scales = "free") +
  scale_x_continuous(breaks = integer_breaks()) +
  scale_y_continuous(limits = c(0,1)) +
  viridis::scale_color_viridis() +
  viridis::scale_fill_viridis() +
  labs(
    title    = NULL,
    x        = "Count of Respondents",
    y        = "Local Clustering in Aggregate Network",
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "local_clust_by_num_surv.png"),
  plot = p,
  width  = 9,
  height = 10
  )

p <- clust_df |>
  ggplot(aes(x = num_surveyed, y = local_clust_w, colour = num_surveyed, group = num_surveyed)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(
    width = 0.2,
    alpha = 0.1,
  ) +
  facet_wrap(.~ site, scales = "free") +
  scale_x_continuous(breaks = integer_breaks()) +
  scale_y_continuous(limits = c(0,1)) +
  viridis::scale_color_viridis() +
  viridis::scale_fill_viridis() +
  labs(
    title    = NULL,
    x        = "Count of Respondents",
    y        = "Weighted Local Clustering in Aggregate Network",
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "weighted_local_clust_by_num_surv.png"),
  plot = p,
  width  = 9,
  height = 10
  )

av_clust_by_num_surv <- respondent_counts |>
  ggplot(aes(x = intended_resp, y = av_local_clust)) +
  geom_violin(
    alpha = 0.4,
    linewidth = 0.5,
    scale = "width",
    draw_quantiles = 0.5
  ) +
  geom_label(
    aes(label = site),
    family = "roboto",
    position = ggplot2::position_jitter(height = 0.02),
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
    coord_cartesian(
    ylim = c(0,1)
  ) +
  labs(
    x = "Intended Respondents",
    y = "Average Local Clustering in Aggregate Network"
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))

av_clust_by_sampling <- respondent_counts |>
  ggplot(aes(x = mean_num_surveyed, y = av_local_clust)) +
  geom_smooth(
    alpha = 0.2,
    fill = "blue",
    method = "lm"
  ) +
  geom_label(
    aes(label = site),
    family = "roboto",
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
  coord_cartesian(
    ylim = c(0,1)
  ) +
  labs(
    x = "Mean Respondents Sampled Per Sharing Unit",
    y = "Average Local Clustering in Aggregate Network"
  ) +
  theme_classic()

p <- (av_clust_by_num_surv | av_clust_by_sampling) + plot_layout(axis_titles = "collect")
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "av_clust_vs_sampling.png"),
  plot = p,
  width  = 12,
  height = 6
  )
