## ============================================================================
## 04_nomination_counts_composition.R
## ----------------------------------------------------------------------------
## Nomination counts and composition.
##
## NOTE: the `alters` object below was computed immediately after the name-
## generator-prompt tables in the original script (right before the "Network
## summary statistics" section began), but it is only used here, for the
## prop_external_violin plot -- so it has been moved down to this script
## rather than left stranded in 01_name_generator_prompts.R.
##
## Builds per-respondent nomination counts (individuals and sharing units),
## internal/external and within-SU composition, and gender/age patterns.
##
## Outputs:
##   Tables/DescriptiveStats/nomination-summary-stats/global-nominations-summary.tex
##   Tables/DescriptiveStats/nomination-summary-stats/nominations-summary-*.tex
##   Figures/DescriptiveStats/prop_external_violin.png
##   Figures/DescriptiveStats/nomination_counts_boxplot.png
##   Figures/DescriptiveStats/nomination_counts_density.png
##   Figures/DescriptiveStats/nomination_counts_violin.png
##   Figures/DescriptiveStats/internal_external_violin.png
##   Figures/DescriptiveStats/count_alters_both.png
##   Figures/DescriptiveStats/count_alters_gender_both.png
##   Figures/DescriptiveStats/count_alters_gender_age.png
##   Figures/DescriptiveStats/nomination_gender_diffs.png
##
## Depends on: 00_setup.R
## ============================================================================

source("net_descriptives/00_setup.R")

alters <- purrr::imap_dfr(su_alters, function(x, site_name) {
  x |>
    ## drop clearly "extra" name generators and the external kinship ones
    dplyr::filter(!tie %in%
      c(11:25, ".", "familyrelation", "farmers", "internal_support",
      "lineage", "related", "relation", "relations", "relations-affine",
      "relations-close", "relative", "relbin", "unrelated",
      "watch_house_f", "watch_house_m")) |>
    dplyr::group_by(tie) |>
    dplyr::summarise(
      tot_alt = mean(alter_count),
      tot_ext = {
        ext <- externals_count
        ext[is.na(ext) & alter_count > 0] <- 0
        mean(ext)
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      site     = site_name,
      prop_ext = tot_ext / tot_alt
    )
}) |>
  dplyr::select(site, tie, tot_alt, tot_ext, prop_ext)

# Save alters, so can be called in by variance decomp file
save(alters, file = file.path(EndowDerivedData, "external_alters.rdata"))


# ---- Proportion of nominations of externals, based on SU networks ----

alters1 <- alters |>
  dplyr::filter(tie %in% c(1:10, "main_req"), !is.na(prop_ext)) |>
  dplyr::mutate(tie_name = factor(tie_labels[tie], levels = rev(tie_labels[-c(12, 17)]))) # -12 and 17 to cut the extra "Aggregate" entries

p <- make_violin_plot(
  alters1, x_var = "prop_ext", y_var = "tie_name", group_var = "tie_name",
  palette = tie_colours,
  flip = FALSE,
  jitter_width = 0.05,
  x_limits = c(0, 1),
  x_breaks = c(0, 0.25, 0.5, 0.75, 1),
  x_label = "Proportion of external alters"
)
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "prop_external_violin.png"),
  plot = p,
  width = 9,
  height = 6
)

# ── Respondent-level nomination counts ────────────────────────────────────────

# NOTE: the su_alters object (created in process_data) also has information on alters,
# aggregated to the SU. this includes things like internal/external status.
# should add external/internal status, to show the # or % of nominations being dropped for the external name generators, in particular
# do with something like:
# residents <- lapply(people_observations, function(df) {
#   df[df$is_resident == 1,]
# })
# and then mark alters as res or not

## To get nomination counts from which we can construct some summary measures
## when working off of the edge lists, we need to account for the fact that
## people who nominated no one would not have any entry, so we have to account
## for those real zeroes; and we also have to account for which questions people
## were presented with -- hence the exceptions below for those sites where
## gender-specific questions were only asked of gender-aligned respondents,
## and sites where some prompts weren't asked.
## Those cases are marked off in 00_setup as are also called by next script.


# lookup. Note: presumes each person has ONE sharing unit
# there are currently two exceptions to this in KS! (ks0003 and ks0143)
su_lookup <- purrr::imap_dfr(
  people_observations,
  function(.x, .y) {
    if (.y == "PE") {
      .x |>
        dplyr::mutate(su_id = su_id_1) |>
        dplyr::select(personid, su_id) |>
        dplyr::mutate(site = .y)
    } else {
      .x |>
        dplyr::select(personid, su_id) |>
        dplyr::mutate(site = .y)
    }
  }
)

alters_summary <- purrr::imap_dfr(
  networks[!names(networks) %in% c("TN")],
  function(df, site_name) {

    persons <- people_observations[[site_name]] |>
      dplyr::filter(surveyed == TRUE) |>
      dplyr::select(personid, gender, age, dplyr::any_of("marital_status"))

    # determine valid ties for the site
    excluded_ties <- if (site_name %in% two_nine_ten) {
      c(2, 9, 10)
    } else if (site_name %in% nine_ten) {
      c(9, 10)
    } else if (site_name %in% nine) {
      c(9)
    } else if (site_name == "MK") {
      c(6, 8) # treating 56a recoded to 5 as sole behav and 78 recoded to 7 as sole info
    }
    else {
      integer(0)
    }

    valid_ties <- setdiff(1:10, excluded_ties)

    if (site_name %in% gender_specific_sites) {
      valid_combinations <- dplyr::bind_rows(
        persons |>
          dplyr::filter(gender == "male") |>
          dplyr::select(personid) |>
          tidyr::crossing(tie = setdiff(valid_ties, c(5, 7))),
        persons |>
          dplyr::filter(gender == "female") |>
          dplyr::select(personid) |>
          tidyr::crossing(tie = setdiff(valid_ties, c(6, 8)))
      )
    } else if (site_name %in% c("KO", "MA")) {
      valid_combinations <- dplyr::bind_rows(
        persons |>
          dplyr::filter(gender == "male") |>
          dplyr::select(personid) |>
          tidyr::crossing(tie = ko_ma_men),
        persons |>
          dplyr::filter(gender == "female") |>
          dplyr::select(personid) |>
          tidyr::crossing(tie = ko_ma_women)
      )
      # empirically add cross-gender responses where they exist (KO)
      exception_ties <- df |>
        dplyr::mutate(tie = as.integer(tie)) |>
        dplyr::inner_join(
          persons |> dplyr::select(personid, gender),
          by = "personid"
        ) |>
        dplyr::filter(
          (gender == "male" & tie %in% ko_ma_women) |
            (gender == "female" & tie %in% ko_ma_men)
        ) |>
        dplyr::distinct(personid, tie)

      valid_combinations <- dplyr::bind_rows(
        valid_combinations,
        exception_ties
      ) |>
        dplyr::distinct(personid, tie)
    } else if (site_name == "MG") {
      valid_combinations <- dplyr::bind_rows(
        persons |>
          dplyr::filter(gender == "male") |>
          dplyr::select(personid) |>
          tidyr::crossing(tie = mg_men),
        persons |>
          dplyr::filter(gender == "female", marital_status %in% "married") |>
          dplyr::select(personid) |>
          tidyr::crossing(tie = mg_married_women),
        persons |>
          dplyr::filter(gender == "female", !(marital_status %in% "married")) |>
          dplyr::select(personid) |>
          tidyr::crossing(tie = mg_other_women)
      )
    } else {
      valid_combinations <- tidyr::crossing(
        persons |> dplyr::select(personid),
        tie = valid_ties
      )
    }

    actual_alters <- df |>
      dplyr::mutate(tie = as.integer(tie)) |>
      dplyr::filter(tie %in% valid_ties) |>
      dplyr::left_join(
        su_lookup |>
          dplyr::filter(site == site_name) |>
          dplyr::select(-site) |>
          dplyr::rename(alterid = personid, alter_su_id = su_id),
        by = "alterid"
      ) |>
      dplyr::group_by(personid, tie) |>
      dplyr::summarise(
        alter_ids      = list(unique(alterid)),
        alter_ids_all  = list(alterid), # for use with KM below
        alter_count    = length(unique(alterid)),
        alter_count_all    = length(alterid), # for when alter named multiple times, for use with alter_su_ids_all
        alter_su_ids_all   = list(alter_su_id),
        alter_su_ids   = list(unique(alter_su_id[!is.na(alter_su_id)])),
        alter_su_count = length(unique(alter_su_id[!is.na(alter_su_id)])),
        .groups        = "drop"
      )

    ## KM alter entries are already SUIDs, though for externals some are IDs
    ## so need to shift those and retain only actual SUIDs; given sampling here
    ## using those only in people_observations will exclude unsampled SUs
    ## so instead, keying off the fact that external alters have "n" in their name,
    if (site_name == "KM") {
      #valid_su_ids <- setdiff(unique(people_observations$KM$su_id), NA)
      actual_alters <- actual_alters |>
        dplyr::mutate(
          #alter_su_ids = purrr::map(alter_ids, ~ .x[.x %in% valid_su_ids]),
          alter_su_ids = purrr::map(
            alter_ids,
            ~ grep("n", .x, value = TRUE, invert = TRUE)
          ),
          alter_su_count = purrr::map_int(alter_su_ids, length),
          alter_count = NA_integer_
        )
    }

    valid_combinations |>
      dplyr::left_join(actual_alters, by = c("personid", "tie")) |>
      dplyr::left_join(
        persons |> dplyr::select(personid, gender, age),
        by = "personid"
      ) |>
      dplyr::mutate(
        alter_ids = purrr::map(
          alter_ids,
          ~ if (is.null(.x)) character(0) else .x
        ),
        alter_su_ids = purrr::map(
          alter_su_ids,
          ~ if (is.null(.x)) character(0) else .x
        ),
        alter_count = if (site_name == "KM") {
          NA_integer_
        } else {
          tidyr::replace_na(alter_count, 0L)
        },
        alter_count_all = tidyr::replace_na(alter_count_all, 0L),
        alter_su_count = tidyr::replace_na(alter_su_count, 0L),
        site = site_name
      )
  }
) |>
  dplyr::left_join(su_lookup, by = c("site", "personid")) |>
  dplyr::mutate(
    alter_same_su_count = purrr::map2_int(
      alter_su_ids_all, su_id,
      ~ sum(.x == .y, na.rm = TRUE)
    ),
    prop_same_su = dplyr::if_else(
      alter_count > 0, alter_same_su_count / alter_count, NA_real_
    ),
    external_count = dplyr::if_else(
      site == "KM",
      purrr::map_int(alter_ids_all, ~ length(grep("n", .x))),
      purrr::map_int(alter_su_ids_all, ~ sum(is.na(.x)))
    ),
    prop_external = dplyr::if_else(
      alter_count_all > 0,
      external_count / purrr::map_int(alter_ids_all, ~ length(.x)),
      NA_real_
    )
  )

alters_summary_aggregate <- alters_summary |>
  dplyr::filter(tie %in% c(1, 3, 5, 6, 7, 8)) |>
  dplyr::group_by(site, personid, gender, age, su_id) |>
  dplyr::summarise(
    alter_ids = list(unique(unlist(alter_ids))),
    alter_count = length(unique(unlist(alter_ids))),
    alter_su_ids = list(unique(unlist(alter_su_ids))),
    alter_su_count = length(unique(unlist(alter_su_ids))),
    alter_same_su_count = sum(unlist(alter_su_ids) == su_id[1], na.rm = TRUE),
    alter_ids_all = list(unlist(alter_ids_all)),
    alter_su_ids_all = list(unlist(alter_su_ids_all)),
    alter_count_all = length(unlist(alter_ids_all)),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    alter_count = dplyr::if_else(site == "KM", NA_integer_, alter_count)
  ) |>
  dplyr::mutate(
    prop_same_su = dplyr::if_else(
      alter_count > 0, alter_same_su_count / alter_count, NA_real_
    ),
    external_count = dplyr::if_else(
      site == "KM",
      purrr::map_int(alter_ids_all, ~ length(grep("n", .x))),
      purrr::map_int(alter_su_ids_all, ~ sum(is.na(.x)))
    ),
    prop_external = dplyr::if_else(
      alter_count_all > 0,
      external_count / purrr::map_int(alter_ids_all, ~ length(.x)),
      NA_real_
    )
  ) |>
  dplyr::mutate(tie = 12L)

alters_summary <- dplyr::bind_rows(alters_summary, alters_summary_aggregate)

## Now that have aggregate measures, dropping MK ties 5 and 7
alters_summary <- alters_summary |>
  dplyr::filter(!(site == "MK" & tie %in% c(5:8)))

## overall proportion of nominations to external alters:
## excluding KM, as most 'alters' are actually sharing units
sum(alters_summary$external_count[
  alters_summary$tie != 12 & alters_summary$site != "KM"
]) /
  sum(alters_summary$alter_count[
    alters_summary$tie != 12 & alters_summary$site != "KM"
  ])

## specific for IH, for line in data paper
sum(alters_summary$external_count[
  alters_summary$tie != 12 & alters_summary$site == "IH"
]) /
  sum(alters_summary$alter_count[
    alters_summary$tie != 12 & alters_summary$site == "IH"
  ])

## overall proportion of nominations to own SU
 sum(alters_summary$alter_same_su_count[
  alters_summary$tie != 12 & alters_summary$site != "KM"
]) /
  sum(alters_summary$alter_count[
    alters_summary$tie != 12 & alters_summary$site != "KM"
  ])

## Summary tables
nom_csv_path <- file.path(csvpath, "nomination-summary-stats")

## Overall
global_nom_table <- alters_summary |>
  dplyr::mutate(
    tie_name = factor(tie_label_short[as.character(tie)], levels = tie_label_short[-c(13,17)]) # -13 to drop "main_req", -17 to drop "any"
  ) |>
  dplyr::group_by(tie_name) |>
  dplyr::summarise(
    prop_zero = mean(alter_count_all == 0, na.rm = TRUE), # using _all so can include KM; given def doesn't make a difference otherwise
    dplyr::across(
      c(alter_count, alter_su_count, prop_external, prop_same_su),
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd   = \(x) sd(x,   na.rm = TRUE)
      ),
      .names = "{.col}__{.fn}"
    ),
    .groups = "drop"
  ) |>
  dplyr::transmute(
    `Name Generator` = tie_name,
    `Prop. Zero`     = sprintf("%.2f", prop_zero),
    `N Indiv.`       = sprintf("%.2f (%.2f)", alter_count__mean,    alter_count__sd),
    `N SU`           = sprintf("%.2f (%.2f)", alter_su_count__mean, alter_su_count__sd),
    `Prop. Same SU`  = sprintf("%.2f (%.2f)", prop_same_su__mean,   prop_same_su__sd),
    `Prop. Ext.`     = sprintf("%.2f (%.2f)", prop_external__mean,  prop_external__sd)
  )

# CSV
readr::write_csv(
  global_nom_table,
  file.path(nom_csv_path, "global-nominations-summary.csv")
)

# LaTeX table
global_nom_table |>
  kableExtra::kbl(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    caption = paste0(
      "Mean (and standard deviation) values for nominations made across all respondents, for each name generator. ",
      "`Prop. Zero' is the proportion of respondents nominating no alters. ",
      "`N Indiv' is the number of individuals nominated. ",
      "`N SU' is the number of sharing units. ",
      "`Prop. Same SU' is the proportion of individuals nominated who are members of the same sharing unit as the respondent. ",
      "`Prop. Ext' is the proportion of individuals nominated who are external to the community."
    ),
    label = "nominations-summary-global"
  ) |>
  kableExtra::kable_styling(latex_options = c("hold_position")) |>
  kableExtra::save_kable(file.path(
    EndowDropbox, "Tables", "DescriptiveStats",
    "nomination-summary-stats", "global-nominations-summary.tex"
  ))

## Site-wise
alters_summary_table <- alters_summary |>
  #dplyr::filter(tie != 12) |>
  dplyr::mutate(
    tie_name = factor(tie_short[as.character(tie)], levels = tie_short[-c(13,17)]),  # -13 to drop "main_req", -17 to drop "any"
  ) |>
  dplyr::group_by(tie_name, site) |>
  dplyr::summarise(
    n_respondents = dplyr::n(),
    prop_zero = mean(alter_count_all == 0, na.rm = TRUE), # using _all so can include KM; given def doesn't make a difference otherwise
    dplyr::across(
      c(alter_count, alter_su_count, prop_external, prop_same_su),
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd   = \(x) sd(x,   na.rm = TRUE)
      ),
      .names = "{.col}__{.fn}"
    ),
    .groups = "drop"
  ) |>
  dplyr::transmute(
    `Name Generator` = tie_name,
    `Site`           = site,
    `N Resp.`        = n_respondents,
    `Prop. Zero`     = sprintf("%.2f", prop_zero),
    `N Indiv.`       = sprintf("%.2f (%.2f)", alter_count__mean,    alter_count__sd),
    `N SU`           = sprintf("%.2f (%.2f)", alter_su_count__mean, alter_su_count__sd),
    `Prop. Same SU`  = sprintf("%.2f (%.2f)", prop_same_su__mean,   prop_same_su__sd),
    `Prop. Ext.`     = sprintf("%.2f (%.2f)", prop_external__mean,  prop_external__sd)
  ) |>
  dplyr::mutate(dplyr::across(where(is.character), \(x) dplyr::if_else(grepl("^NaN|\\(NaN\\)|\\(NA\\)", x), "---", x))) |>
  dplyr::arrange(`Name Generator`, `Site`)

for (i in tie_meta$short[c(1:10, 12)]) {
  base_name <- paste0("nominations-summary-", i)
  filename <- file.path(EndowDropbox, "Tables", "DescriptiveStats",
                            "nomination-summary-stats", paste0(base_name, ".tex"))
  csv_filename <- file.path(nom_csv_path, paste0(base_name, ".csv"))
  message("Saving: ", base_name)

  tab <- alters_summary_table |>
      dplyr::filter(`Name Generator` == i) |>
      dplyr::select(-`Name Generator`) |>
      dplyr::arrange(Site)

    # CSV: turn the "---" placeholders back into proper missing values
    tab |>
      dplyr::mutate(dplyr::across(where(is.character), \(x) dplyr::na_if(x, "---"))) |>
      readr::write_csv(csv_filename)

    # LaTeX table
    tab |>
      kableExtra::kable(
      format    = "latex",
      booktabs  = TRUE,
      longtable = TRUE,
      escape    = FALSE,
      caption   = paste0(
        "Mean (and standard deviation) values for nominations made across all respondents, for the ",
        short_to_label[i], " name generator. ",
        "`N Resp.' is the number of respondents. ",
        "`Prop. Zero' is the proportion of respondents nominating no alters. ",
        "`N Indiv' is the number of individuals nominated. ",
        "`N SU' is the number of sharing units. ",
        "`Prop. Same SU' is the proportion of individuals nominated who are members of the same sharing unit as the respondent. ",
        "`Prop. Ext' is the proportion of individuals nominated who are external to the community."
      ),
      label     = paste0("nomination-summary-", i),
      linesep   = ""
    ) |>
    kableExtra::kable_styling(
      latex_options = "repeat_header",
      font_size     = 8
    ) |>
    kableExtra::save_kable(filename)
}


## PLOT: Boxplot of distributions of number of alters per question per site
p1 <- alters_summary |>
  dplyr::filter(tie != 12) |>
  dplyr::mutate(
    tie_name = factor(tie_label_short[as.character(tie)], levels = tie_label_short[-c(13,17)])  # -13 to drop "main_req", -17 to drop "any"
  ) |>
  ggplot(aes(x = alter_count, y = site, fill = tie_name, colour = tie_name)) +
  geom_boxplot(
    alpha         = 0.4,
    linewidth     = 0.3,
    outlier.size  = 0.4,
    outlier.alpha = 0.3
  ) +
  facet_wrap(~ tie_name, nrow = 2, dir = "v") +
  scale_fill_manual(values = tie_colours_short, guide = "none") +
  scale_colour_manual(values = tie_colours_short, guide = "none") +
  coord_cartesian(xlim = c(0, NA)) +
  geom_vline(xintercept = c(5, 10), linewidth = 0.4, linetype = rep(c("dashed", "dotted"), 10)) +
  labs(x = "Nominations per respondent", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0),
    panel.spacing    = unit(0.5, "lines"),
    axis.text.y      = element_text(size = 6)
  )

p2 <- alters_summary |>
  dplyr::filter(tie == 12) |>
  dplyr::mutate(
    tie_name = factor(tie_label_short[as.character(tie)], levels = tie_label_short[-c(13,17)])  # -13 to drop "main_req", -17 to drop "any"
  ) |>
  ggplot(aes(x = alter_count, y = site, fill = tie_name, colour = tie_name)) +
  geom_boxplot(
    alpha         = 0.4,
    linewidth     = 0.3,
    outlier.size  = 0.4,
    outlier.alpha = 0.3
  ) +
  facet_wrap(~ tie_name, nrow = 2, dir = "v") +
  scale_fill_manual(values = tie_colours_short, guide = "none") +
  scale_colour_manual(values = tie_colours_short, guide = "none") +
  coord_cartesian(xlim = c(0, NA)) +
  geom_vline(xintercept = c(5, 10), linewidth = 0.4, linetype = c("dashed", "dotted")) +
  labs(x = "Nominations per respondent", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(hjust = 0),
    panel.spacing    = unit(0.5, "lines"),
    axis.text.y      = element_text(size = 6)
  )

p <- p1 + p2 + patchwork::plot_layout(widths = c(5, 1), axes = "collect")

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "nomination_counts_boxplot.png"),
  plot = p,
  width = 12,
  height = 8
)

safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

nomination_summary <- alters_summary |>
  dplyr::group_by(site, tie) |>
  dplyr::summarise(
    n_respondents          = dplyr::n_distinct(personid),
    n_respondents_male     = dplyr::n_distinct(personid[gender == "male"]),
    n_respondents_female   = dplyr::n_distinct(personid[gender == "female"]),

    # alter_count summaries
    min_nominations        = safe_min(alter_count),
    max_nominations        = safe_max(alter_count),
    mean_nominations       = mean(alter_count,   na.rm = TRUE),
    median_nominations     = median(alter_count, na.rm = TRUE),
    sd_nominations         = sd(alter_count,     na.rm = TRUE),
    ttest_ids = tryCatch(
      list(broom::tidy(t.test(
        alter_count[gender == "male"],
        alter_count[gender == "female"]
      ))),
      error = function(e) list(tibble::tibble(
        estimate1 = NA_real_, estimate2 = NA_real_, statistic = NA_real_,
        p.value = NA_real_, parameter = NA_real_, conf.low = NA_real_,
        conf.high = NA_real_, method = NA_character_, alternative = NA_character_
      ))
    ),

    # alter_count ~ age + age^2
    lm_ids = list(tryCatch(
      broom::tidy(lm(alter_count ~ age + I(age^2))) |>
        dplyr::mutate(term = dplyr::recode(term,
          "(Intercept)" = "intercept",
          "age"         = "age",
          "I(age^2)"    = "age2"
        )) |>
        tidyr::pivot_wider(
          names_from  = term,
          values_from = c(estimate, std.error, statistic, p.value)
        ),
      error = function(e) tibble::tibble(
        estimate_intercept = NA_real_, std.error_intercept = NA_real_,
        statistic_intercept = NA_real_, p.value_intercept = NA_real_,
        estimate_age = NA_real_, std.error_age = NA_real_,
        statistic_age = NA_real_, p.value_age = NA_real_,
        estimate_age2 = NA_real_, std.error_age2 = NA_real_,
        statistic_age2 = NA_real_, p.value_age2 = NA_real_
      )
    )),

    # alter_su_count summaries
    min_nominations_su     = safe_min(alter_su_count),
    max_nominations_su     = safe_max(alter_su_count),
    mean_nominations_su    = mean(alter_su_count,   na.rm = TRUE),
    median_nominations_su  = median(alter_su_count, na.rm = TRUE),
    sd_nominations_su      = sd(alter_su_count,     na.rm = TRUE),
    ttest_su_ids = tryCatch(
      list(broom::tidy(t.test(
        alter_su_count[gender == "male"],
        alter_su_count[gender == "female"]
      ))),
      error = function(e) list(tibble::tibble(
        estimate1 = NA_real_, estimate2 = NA_real_, statistic = NA_real_,
        p.value = NA_real_, parameter = NA_real_, conf.low = NA_real_,
        conf.high = NA_real_, method = NA_character_, alternative = NA_character_
      ))
    ),

    # alter_su_count ~ age + age^2
    lm_su_ids = list(tryCatch(
      broom::tidy(lm(alter_su_count ~ age + I(age^2))) |>
        dplyr::mutate(term = dplyr::recode(term,
          "(Intercept)" = "intercept",
          "age"         = "age",
          "I(age^2)"    = "age2"
        )) |>
        tidyr::pivot_wider(
          names_from  = term,
          values_from = c(estimate, std.error, statistic, p.value)
        ),
      error = function(e) tibble::tibble(
        estimate_intercept = NA_real_, std.error_intercept = NA_real_,
        statistic_intercept = NA_real_, p.value_intercept = NA_real_,
        estimate_age = NA_real_, std.error_age = NA_real_,
        statistic_age = NA_real_, p.value_age = NA_real_,
        estimate_age2 = NA_real_, std.error_age2 = NA_real_,
        statistic_age2 = NA_real_, p.value_age2 = NA_real_
      )
    )),

    # prop_same_su summaries
    mean_prop_same_su      = mean(prop_same_su,   na.rm = TRUE),
    prop_any_same_su       = mean(prop_same_su > 0, na.rm = TRUE),
    mean_prop_external      = mean(prop_external,   na.rm = TRUE),
    prop_any_external       = mean(prop_external > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  tidyr::unnest(ttest_ids,    names_sep = "_") |>
  tidyr::unnest(lm_ids,       names_sep = "_") |>
  tidyr::unnest(ttest_su_ids, names_sep = "_") |>
  tidyr::unnest(lm_su_ids,    names_sep = "_") |>
  dplyr::rename(
    mean_nominations_male          = ttest_ids_estimate1,
    mean_nominations_female        = ttest_ids_estimate2,
    ttest_statistic                = ttest_ids_statistic,
    ttest_p_value                  = ttest_ids_p.value,

    age_effect_intercept            = lm_ids_estimate_intercept,
    age_effect_age                  = lm_ids_estimate_age,
    age_effect_age2                 = lm_ids_estimate_age2,
    age_effect_age_se               = lm_ids_std.error_age,
    age_effect_age2_se              = lm_ids_std.error_age2,
    age_effect_age_stat             = lm_ids_statistic_age,
    age_effect_age2_stat            = lm_ids_statistic_age2,
    age_effect_age_pvalue           = lm_ids_p.value_age,
    age_effect_age2_pvalue          = lm_ids_p.value_age2,

    mean_nominations_su_male       = ttest_su_ids_estimate1,
    mean_nominations_su_female     = ttest_su_ids_estimate2,
    ttest_su_statistic              = ttest_su_ids_statistic,
    ttest_su_p_value                = ttest_su_ids_p.value,

    age_effect_su_intercept         = lm_su_ids_estimate_intercept,
    age_effect_su_age               = lm_su_ids_estimate_age,
    age_effect_su_age2              = lm_su_ids_estimate_age2,
    age_effect_su_age_se            = lm_su_ids_std.error_age,
    age_effect_su_age2_se           = lm_su_ids_std.error_age2,
    age_effect_su_age_stat          = lm_su_ids_statistic_age,
    age_effect_su_age2_stat         = lm_su_ids_statistic_age2,
    age_effect_su_age_pvalue        = lm_su_ids_p.value_age,
    age_effect_su_age2_pvalue       = lm_su_ids_p.value_age2
  ) |>
  dplyr::select(
    site, tie,
    n_respondents, n_respondents_male, n_respondents_female,
    min_nominations, max_nominations, mean_nominations, median_nominations, sd_nominations,
    mean_nominations_male, mean_nominations_female,
    ttest_statistic, ttest_p_value,
    age_effect_age, age_effect_age2, age_effect_age_se, age_effect_age2_se,
    age_effect_age_stat, age_effect_age2_stat, age_effect_age_pvalue, age_effect_age2_pvalue,
    min_nominations_su, max_nominations_su, mean_nominations_su, median_nominations_su, sd_nominations_su,
    mean_nominations_su_male, mean_nominations_su_female,
    ttest_su_statistic, ttest_su_p_value,
    age_effect_su_age, age_effect_su_age2, age_effect_su_age_se, age_effect_su_age2_se,
    age_effect_su_age_stat, age_effect_su_age2_stat, age_effect_su_age_pvalue, age_effect_su_age2_pvalue,
    mean_prop_same_su, prop_any_same_su,
    mean_prop_external, prop_any_external
  )

## manually correcting MY men, interviewed jointly with wives
## NOTE: this presumes that if recorded as male_hh (or female_hh), were present!
nomination_summary$n_respondents_male[
  nomination_summary$site == "MY"
] <- length(setdiff(unique(networks$MY$male_hh), NA))

## PLOT: density plot of the average number of people nominated for each
## question by each respondent, for each site
## Combines pairs of questions together (double, sampled, gender, external)
p <- nomination_summary |>
  dplyr::filter(tie != 12) |>
  dplyr::mutate(
    tie_name = factor(tie_labels[as.character(tie)], levels = tie_meta$label[-c(13,17)]),  # -13 to drop "main_req", -17 to drop "any"
    pair     = factor(tie_pairs[tie_name], levels = pair_levels)
  ) |>
  ggplot(aes(x = mean_nominations, fill = tie_name, colour = tie_name)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  scale_fill_manual(values = tie_colours, name = "Name generator") +
  scale_colour_manual(values = tie_colours, name = "Name generator") +
  scale_x_continuous(breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 1)) +
  facet_wrap(~ pair, ncol = 2, scales = "free_y") +
  labs(
    x     = "Mean nominations per respondent",
    y     = "Density",
    # title = "Distribution of mean nominations across communities by name generator"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "nomination_counts_density.png"),
  plot = p,
  width = 9,
  height = 6
)


## PLOT: Violin plot version of the above, with both indiv and su noms
df_nom1 <- nomination_summary |>
  dplyr::filter(tie != 12) |>
  dplyr::mutate(tie_name = factor(tie_labels[as.character(tie)], levels = rev(tie_labels[-c(13,17)]))) |>  # -13 to drop "main_req", -17 to drop "any"
  tidyr::pivot_longer(cols = c(mean_nominations, mean_nominations_su),
                       names_to = "measure", values_to = "value") |>
  dplyr::mutate(measure = dplyr::recode(measure,
    mean_nominations = "Individuals", mean_nominations_su = "Sharing units"))

p1 <- make_violin_plot(
  df_nom1, x_var = "value", y_var = "tie_name", group_var = "tie_name",
  palette = tie_colours,
  flip = FALSE,
  jitter_width = 0.1,
  x_breaks = seq(0, 10, by = 1),
  coord_xlim = c(0, 10),
  facet_var = "measure",
  x_label = "Mean nominations per respondent"
)

df_nom2 <- nomination_summary |>
  dplyr::filter(tie == 12) |>
  dplyr::mutate(tie_name = factor(tie_labels[as.character(tie)])) |>
  tidyr::pivot_longer(cols = c(mean_nominations, mean_nominations_su),
                       names_to = "measure", values_to = "value") |>
  dplyr::mutate(measure = dplyr::recode(measure,
    mean_nominations = "Individuals", mean_nominations_su = "Sharing units"))

p2 <- make_violin_plot(
  df_nom2, x_var = "value", y_var = "tie_name", group_var = "tie_name",
  palette = tie_colours,
  flip = FALSE,
  jitter_width = 0.1,
  x_breaks = seq(0, 20, by = 2),
  coord_xlim = c(0, 20),
  facet_var = "measure",
  x_label = "Mean nominations per respondent"
)

p2 <- p2 + ggplot2::theme(strip.text = ggplot2::element_blank())

p <- p1 / p2 + patchwork::plot_layout(heights = c(8, 1), axes = "collect")

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "nomination_counts_violin.png"),
  plot = p,
  width = 12,
  height = 6
)

df_nom_intext <- nomination_summary |>
  #dplyr::filter(tie != 12) |>
  dplyr::mutate(tie_name = factor(tie_labels[as.character(tie)], levels = rev(tie_labels[-c(13,17)]))) |> # -13 to drop "main_req", -17 to drop "any"
  tidyr::pivot_longer(cols = c(mean_prop_same_su, mean_prop_external), # alternative: use prop_any_external
                       names_to = "measure", values_to = "value") |>
  dplyr::mutate(measure = dplyr::recode(measure,
    mean_prop_same_su = "Same Sharing Unit", mean_prop_external = "External to Community"))

p <- make_violin_plot(
  df_nom_intext, x_var = "value", y_var = "tie_name", group_var = "tie_name",
  palette = tie_colours,
  flip = FALSE,
  jitter_width = 0.01,
  coord_xlim = c(0, 1),
  facet_var = "measure",
  x_label = "Proportion"
)
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "internal_external_violin.png"),
  plot = p,
  width = 12,
  height = 6
)

# Summary table, by median
nomination_summary |>
  dplyr::group_by(tie) |>
  dplyr::summarise(across(where(is.numeric), \(x) median(x, na.rm = TRUE))) |>
  dplyr::mutate(network = factor(
    setNames(tie_meta$label_short, tie_meta$id)[tie],
    levels = tie_meta$label_short[-c(13,17)] # -13 to drop "main_req", -17 to drop "any"
  )) |>
  dplyr::relocate(network)


## PLOT: number of nominations
site_order <- alters_summary |>
  dplyr::filter(tie == 12) |>
  dplyr::group_by(site) |>
  dplyr::summarise(median_alter_su_count = median(alter_su_count, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(median_alter_su_count) |>
  dplyr::pull(site)

count_ids_aggregate <- alters_summary |>
  dplyr::filter(!site %in% c("KM", "KO", "MA")) |>
  dplyr::filter(tie == 12) |>
  dplyr::mutate(site = factor(site, levels = site_order)) |>
  ggplot(aes(x = site, y = alter_count)) +
  geom_boxplot(outlier.size = 0.5, fill = "grey") +
  labs(
    title    = NULL,
    subtitle = "Individuals",
    x        = NULL,
    y        = "Count",
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))


count_suids_aggregate <- alters_summary |>
  dplyr::filter(!site %in% c("KO", "MA")) |>
  dplyr::filter(tie == 12) |>
  dplyr::mutate(site = factor(site, levels = site_order)) |>
  ggplot(aes(x = site, y = alter_su_count)) +
  geom_boxplot(outlier.size = 0.5, fill = "grey") +
  labs(
    title = NULL,
    subtitle = "Sharing units",
    x = NULL,
    y = "Count",
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))
p <- (count_ids_aggregate | count_suids_aggregate) / patchwork::guide_area() +
  patchwork::plot_layout(guides = "collect", heights = c(20, 1))
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "count_alters_both.png"),
  plot = p,
  width  = 12,
  height = 6
  )

## PLOT: number of nominations by gender
count_ids_aggregate_gender <- alters_summary |>
  dplyr::filter(!site %in% c("BM", "KM", "KO", "MA", "MY")) |>
  dplyr::filter(tie == 12) |>
  dplyr::mutate(site = factor(site, levels = site_order)) |>
  ggplot(aes(x = site, y = alter_count, fill = gender)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(3,"PRGn")[3], RColorBrewer::brewer.pal(3,"PRGn")[1])) +
  labs(
    title    = NULL,
    subtitle = "Individuals",
    x        = NULL,
    y        = "Count",
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 6, angle = 45, hjust = 1)) +
  guides(colour = guide_legend(title = "Gender"))

count_suids_aggregate_gender <- alters_summary |>
  dplyr::filter(!site %in% c("BM", "KO", "MA", "MY")) |>
  dplyr::filter(tie == 12) |>
  dplyr::mutate(site = factor(site, levels = site_order)) |>
  ggplot(aes(x = site, y = alter_su_count, fill = gender)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(3,"PRGn")[3], RColorBrewer::brewer.pal(3,"PRGn")[1])) +
  labs(
    title    = NULL,
    subtitle = "Sharing units",
    x        = NULL,
    y        = "Count",
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 6, angle = 45, hjust = 1)) +
  guides(colour = guide_legend(title = "Gender"))
p <- (count_ids_aggregate_gender | count_suids_aggregate_gender) / patchwork::guide_area() +
  patchwork::plot_layout(guides = "collect", heights = c(20, 1))
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "count_alters_gender_both.png"),
  plot = p,
  width  = 12,
  height = 6
  )


## PLOT: alter count by gender *and age*
alter_base_data <- alters_summary |>
  # remove MA and KO because asked different questions;
  # currently DOESN'T deal with other gender-specific nom cases;
  # have BM and MY show entries only for wife's age because pairs asked jointly;
  # they will be coloured grey/purple
  dplyr::filter(!site %in% c("MA", "KO")) |>
  dplyr::filter(tie == 12) |>
  dplyr::filter(!site %in% c("BM", "MY") | gender == "female")

# only plot curve if have at least 10 responses for a gender; still plot all points
alter_smooth_data <- alter_base_data |>
  dplyr::group_by(site, gender) |>
  dplyr::filter(n() >= 10) |>
  ungroup()

p <- alter_base_data |>
  ggplot(aes(age, alter_count, colour = gender)) +
  geom_point(
    data = dplyr::filter(alter_base_data, !site %in% c("BM", "MY")),
    alpha = 0.2
  ) +
  geom_point(
    data = dplyr::filter(alter_base_data, site %in% c("BM", "MY")),
    alpha = 0.2,
    colour = "grey"
  ) +
  geom_smooth(
    data = dplyr::filter(alter_smooth_data, !site %in% c("BM", "MY")),
    method = lm, formula = y ~ splines::bs(x, 3), se = FALSE
  ) +
  geom_smooth(
    data = dplyr::filter(alter_smooth_data, site %in% c("BM", "MY")),
    method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, colour = "grey"
  ) +
  scale_colour_manual(values = c(RColorBrewer::brewer.pal(3,"PRGn")[3], RColorBrewer::brewer.pal(3,"PRGn")[1]))  +
  facet_wrap(. ~ site, scales = "free") +
  labs(
    title    = NULL,
    x        = "Age",
    y        = "Count",
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(title = "Gender"))
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "count_alters_gender_age.png"),
  plot = p,
  width  = 9,
  height = 10
  )

## PLOT: differences in nominations by gender

## manually prep NAs for nine_ten, two_nine_ten, nine cases so they show as grey!
na_cases <- bind_rows(
  expand_grid(site = two_nine_ten, tie = c(2, 9, 10)),
  expand_grid(site = nine_ten, tie = c(9, 10)),
  expand_grid(site = nine, tie = c(9)),
  expand_grid(site = "MK", tie = 5:8)
)

gender_noms <- bind_rows(nomination_summary, na_cases)

gender_noms <- gender_noms |>
  # remove BM and MY because pairs asked jointly;
  # remove MA and KO because asked different questions;
  # remove BD and KS because only one gender responded
  # remove MG, because men and women asked different sets of questions,
  # with married women asked about both self and husb inflating responses
  dplyr::filter(!site %in% c("BM", "MA", "KO", "BD", "KS", "MY", "MG")) |>
  dplyr::mutate(
    tie_name = factor(tie_labels[as.character(tie)], levels = tie_meta$label[-c(13,17)]), # -13 to drop "main_req", -17 to drop "any"
    diff_indiv = mean_nominations_male - mean_nominations_female,
    diff_su = mean_nominations_su_male - mean_nominations_su_female,
    sig_indiv = dplyr::case_when(
      ttest_p_value < 0.001 ~ "***",
      ttest_p_value < 0.01 ~ "**",
      ttest_p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    sig_su = dplyr::case_when(
      ttest_su_p_value < 0.001 ~ "***",
      ttest_su_p_value < 0.01 ~ "**",
      ttest_su_p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# gender
table(gender_noms$ttest_p_value < 0.01, gender_noms$tie)
table(abs(gender_noms$diff_indiv) > 2, gender_noms$tie)

table(gender_noms$ttest_su_p_value < 0.01, gender_noms$tie)
table(abs(gender_noms$diff_su) > 2, gender_noms$tie)

# age
table(nomination_summary$age_effect_age_pvalue < 0.01, nomination_summary$tie)
table(nomination_summary$age_effect_age2_pvalue < 0.01, nomination_summary$tie)

table(nomination_summary$age_effect_su_age_pvalue < 0.01, nomination_summary$tie)
table(nomination_summary$age_effect_su_age2_pvalue < 0.01, nomination_summary$tie)


# PLOT Heatmap: site x tie, fill = male - female difference
shared_limits <- range(
  gender_noms$diff_indiv,#[gender_noms$tie != 12],
  gender_noms$diff_su,#[gender_noms$tie != 12],
  na.rm = TRUE)
shared_limits <- c(-max(abs(shared_limits)), max(abs(shared_limits)))

gender_noms_indiv <- gender_noms |>
  #dplyr::filter(tie != 12) |>
  ggplot(aes(y = tie_name, x = site, fill = diff_indiv)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_indiv), size = 3, vjust = 0.75) +
  scale_fill_distiller(
    palette = "PRGn",
    na.value = "#cccccc",
    limits = shared_limits,
    name = "Male − Female\nmean nominations"
  ) +
  labs(
    title = "Individuals",
    #subtitle = "Asterisks indicate t-test p < 0.05 / 0.01 / 0.001",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1)
  )

gender_noms_su <- gender_noms |>
  #dplyr::filter(tie != 12) |>
  ggplot(aes(y = tie_name, x = site, fill = diff_su)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_su), size = 3, vjust = 0.75) +
  scale_fill_distiller(
    palette = "PRGn",
    na.value = "#cccccc",
    limits = shared_limits,
    name = "Male − Female\nmean nominations"
  ) +
  labs(
    title = "Sharing units",
    #subtitle = "Asterisks indicate t-test p < 0.05 / 0.01 / 0.001",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1)
  )
p <- (gender_noms_indiv / gender_noms_su) / patchwork::guide_area() +
  patchwork::plot_layout(guides = "collect", heights = c(10, 10, 1))
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "nomination_gender_diffs.png"),
  plot = p,
  width = 10,
  height = 10
)

# Save alters_summary, so can be called in by concordance_within_su
save(alters_summary, file = file.path(EndowDerivedData, "alters_summary.rdata"))