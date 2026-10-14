## ============================================================================
## 00_setup.R
## ----------------------------------------------------------------------------
## Shared setup for the ENDOW network descriptive-statistics pipeline.
##
## Loads key files and constructs shared lookups and an oft-used plot function
##
## Run this first (or simply source it from any downstream script).
## ============================================================================

## Network Descriptive Statistics
##
## Generates summary statistics of social support nominations and networks
##
## Outputs:
##
## -- Name generator prompts --
##   Tables/DescriptiveStats/name-generator-prompts/*.tex      - prompt wording used for each name generator, by community
##
## -- Network-level summary statistics --
##   Tables/DescriptiveStats/network-summary-stats/net-summary-stat-*.tex      - per-layer tables of basic network summary statistics (density, degree, centralization, clustering, etc.), by community
##   Tables/DescriptiveStats/network-summary-stats/global-net-summary-medians.tex - median summary statistics across communities, one column per network layer
##   Figures/DescriptiveStats/*_distrib_violin.png             - violin plots showing the cross-community distribution of each network summary statistic, by layer
##   Figures/DescriptiveStats/density_vs_size.png               - aggregate-network density vs. number of sharing units, by community
##   Figures/DescriptiveStats/selfloop_prop_distrib_violin.png  - violin plot of the proportion of self-loop ties, by layer
##
## -- Respondent coverage --
##   Tables/DescriptiveStats/respondents_by_site.tex            - sharing-unit and respondent coverage by community
##
## -- Nomination counts and composition --
##   Tables/DescriptiveStats/nomination-summary-stats/global-nominations-summary.tex   - mean/SD nomination counts pooled across all communities, by name generator
##   Tables/DescriptiveStats/nomination-summary-stats/nominations-summary-*.tex          - mean/SD nomination counts by community, one table per name generator
##   Figures/DescriptiveStats/prop_external_violin.png          - distribution of the proportion of nominations to external alters, by name generator
##   Figures/DescriptiveStats/nomination_counts_boxplot.png      - per-respondent nomination counts, faceted by name generator and community
##   Figures/DescriptiveStats/nomination_counts_density.png      - density of mean nominations per respondent, faceted by paired name generator
##   Figures/DescriptiveStats/nomination_counts_violin.png       - mean nominations per respondent (individuals vs. sharing units), by name generator
##   Figures/DescriptiveStats/internal_external_violin.png       - proportion of nominations that are within-sharing-unit vs. external to the community
##   Figures/DescriptiveStats/own_su_prop_violin.png              - proportion of respondents naming their own sharing unit, by name generator
##   Figures/DescriptiveStats/count_alters_both.png                - aggregate nomination counts (individuals and sharing units) by community
##   Figures/DescriptiveStats/count_alters_gender_both.png          - aggregate nomination counts by community, split by respondent gender
##   Figures/DescriptiveStats/count_alters_gender_age.png            - aggregate nomination counts vs. respondent age, by gender and community
##   Figures/DescriptiveStats/nomination_gender_diffs.png            - heatmap of male-minus-female mean nominations, by community and name generator
##
## -- Within-sharing-unit concordance --
##   Figures/DescriptiveStats/su_concordance_layerwise_personid.png  - concordance in individual-level nominations among co-resident sharing-unit members, by layer
##   Figures/DescriptiveStats/su_concordance_layerwise_suid.png       - concordance in sharing-unit-level nominations among co-resident members, by layer
##   Figures/DescriptiveStats/su_concordance_aggregate_both.png        - aggregate-network concordance (individuals and sharing units), by community
##   Figures/DescriptiveStats/su_concordance_aggregate_gender_both.png - aggregate-network concordance by community, split by gender
##
## -- Request/give double-sampling concordance --
##   Tables/DescriptiveStats/double_sampled_concordance.tex     - proportion of concordant ties between paired request and give networks (loan, exchange), by community
##   Figures/DescriptiveStats/double_sampled_concordance.png    - distribution of concordance proportions, by paired name generator
##
## -- Degree distributions --
##   Figures/DescriptiveStats/degree_distributions.png          - CCDF degree distributions by layer, all communities overlaid
##   Figures/DescriptiveStats/degree_distribution_aggregate.png  - CCDF degree distribution for the aggregate network, all communities overlaid
##
## -- Network visualizations --
##   Figures/Networks/*_all_sites.pdf       - one network-layer visualization per community, for each network type
##   Figures/Networks/*_all_nets.pdf         - all network-layer visualizations for a single community, one PDF per community
##   Figures/Networks/eg_agg_nets.pdf         - example aggregate-network visualizations for selected communities
##
## -- Layer overlap (edge- and alter-level) --
##   Figures/Networks/edge_overlap_heatmap.pdf          - edge-weight and binary-edge correlation across layers, by community
##   Figures/Networks/edge_overlap_heatmap_pooled.png    - edge-weight and binary-edge correlation across layers, pooled across communities
##   Figures/Networks/alter_overlap_heatmap.pdf           - Jaccard similarity of nominated individuals/sharing units across layers, by community
##   Figures/Networks/alter_overlap_heatmap_pooled.png     - Jaccard similarity of nominated individuals/sharing units across layers, pooled
##
## -- Principal component analysis of network layers --
##   Tables/pca_loadings_weighted_pooled.tex            - PCA loadings on edge weights, pooled across communities
##   Tables/pca_loadings_weighted_sites.tex              - PCA loadings on edge weights, by community
##   Tables/pca_loadings_binary_sites.tex                 - PCA loadings on binarized edges, by community
##   Figures/Networks/pca_scree_weighted_sites.pdf        - variance explained by PC (weighted), by community
##   Figures/Networks/pca_biplot_weighted_sites.pdf        - PC1/PC2 correlation circle (weighted), by community
##   Figures/Networks/pca_biplot_weighted_23_sites.pdf      - PC2/PC3 correlation circle (weighted), by community
##   Figures/Networks/pca_scree_binary_sites.pdf             - variance explained by PC (binary), by community
##   Figures/Networks/pca_biplot_binary_sites.pdf             - PC1/PC2 correlation circle (binary), by community
##   Figures/Networks/pca_biplot_binary_23_sites.pdf            - PC2/PC3 correlation circle (binary), by community
##   Figures/Networks/pca_scree_weighted_pooled.png              - variance explained by PC (weighted), pooled
##   Figures/Networks/pca_biplot_weighted_pooled.png               - PC1/PC2 correlation circle (weighted), pooled
##   Figures/Networks/pca_biplot_weighted_23_pooled.png             - PC2/PC3 correlation circle (weighted), pooled
##   Figures/Networks/pca_scree_binary_pooled.png                    - variance explained by PC (binary), pooled
##   Figures/Networks/pca_biplot_binary_pooled.png                    - PC1/PC2 correlation circle (binary), pooled
##   Figures/Networks/pca_biplot_binary_23_pooled.png                  - PC2/PC3 correlation circle (binary), pooled
##   Figures/Networks/pca_biplot_weighted_az.png, _ti.png, _ni.png       - PC1/PC2 correlation circles for selected individual communities

# to prevent function conflicts!
if ("package:reshape2" %in% search()) {
  detach("package:reshape2", unload = TRUE)
}

if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
}

library(dplyr)
library(ggplot2)
library(plm)
library(tidyverse)
library(gridExtra)
library(igraph)
library(patchwork)
library(ggforce)

## specifying here sites to drop:
SITES_TO_DROP <- c("BZ", "PS")
#  "BZ" # not fully processed through, and doesn't want to be included
#  "PS" # poor coverage of sharing units (53%)

# NOTE: If this sources versions of the underlying data files generated for
# the main ENDOW paper *you will be missing some sites that we can include
# here* -- sites that are dropped there because of poor _wealth_ coverage.
# If so, you need to go back to 00_run_all.R and regenerate the main files.

# These three variables must be defined in your .Rprofile (see README.md).
# EndowGitHub  — path to this repository (ENDOW-wealth-corr)
# EndowDatabase — path to the endow-database repository
# EndowDropbox — path to your local ENDOW Dropbox folder

if (!exists("EndowGitHub") || !is.character(EndowGitHub)) {
  stop(
    "EndowGitHub is not defined. ",
    "Follow the R Setup section of README.md to configure your .Rprofile."
  )
}

if (!exists("EndowDatabase") || !is.character(EndowDatabase)) {
  stop(
    "EndowDatabase is not defined. ",
    "Follow the R Setup section of README.md to configure your .Rprofile."
  )
}

if (!exists("EndowDropbox") || !is.character(EndowDropbox)) {
  stop(
    "EndowDropbox is not defined. ",
    "Follow the R Setup section of README.md to configure your .Rprofile."
  )
}

# Derived path used throughout — all intermediate data files live here
# EndowDerivedData is defined in .Rprofile on Windows; on Mac we fall back to
# the DerivedData subfolder of this repository.
if (!exists("EndowDerivedData")) {
  EndowDerivedData <- file.path(EndowGitHub, "DerivedData")
}

figpath <- file.path(EndowDropbox, "Figures", "DescriptiveStats")
tablepath <- file.path(EndowDropbox, "Tables", "DescriptiveStats")

# Load key data files
load(file.path(EndowDerivedData, "su_alters.rdata"))
load(file.path(EndowDerivedData,"su_nets.rdata")) # mostly for MK
load(file.path(EndowDerivedData,"su_nets_expanded.rdata"))
load(file.path(EndowDerivedData,"su_nets_expanded_withloops.rdata")) # for self-loop analysis
load(file.path(EndowDerivedData, "su_meta.rdata"))
load(file.path(EndowDerivedData, "people_observations.rdata"))

networks <- Map(read.csv, Sys.glob(file.path(EndowDatabase, "networks", "*_networks.csv"))) |>
  setNames(regmatches(
    Sys.glob(file.path(EndowDatabase, "networks", "*_networks.csv")),
    regexpr("[A-Z]{2}(?=.*\\.csv)", Sys.glob(file.path(EndowDatabase, "networks", "*_networks.csv")), perl = TRUE)
  ))

site_codes <- sort(names(su_alters))

## Minor site-specific cleaning

## specifically for NP, for the counts here, need to drop people who are not residents
# some of them did the survey but the ethnog says they shouldn't be included
people_observations$NP <- people_observations$NP[people_observations$NP$is_resident == 1,]

## For MY, M and F household heads surveyed jointly, vars are female_hh, male_hh
# NOTE: this presumes that if recorded as male_hh (or female_hh), were present!
people_observations$MY$surveyed <- people_observations$MY$personid %in%
  c(networks$MY$female_hh, networks$MY$male_hh)
# treating female_hh as single reporter for our purposes to get personid variable
networks$MY$personid <- networks$MY$female_hh

## For MK, questions 5/6 asked without gender specification of everyone; two
## questions posed (so entries in networks of 5/6a and 5/6b; in agg network, use
## 5/6a). Also question 7/8 was singular and not gender specific
## DOING JANKY RECODING HERE, SO WORKS WITH STANDARD CODE! Aligning here, so
## that for 'aggregate' nominations, end up with same set used in main analyses
## but then for things using gender-specific questions, removing those questions
## (see below following delineation of alters_summary object)
networks$MK <- networks$MK |>
  dplyr::mutate(tie_orig = tie) |>
  # Expand "7/8" rows into two rows before recoding
  #(\(d) {
    #dplyr::bind_rows(
    #  d,
    #  d |> dplyr::filter(tie == "7/8") |> dplyr::mutate(tie = "8")
    #)
  #})() |>
  dplyr::mutate(
    tie = dplyr::case_when(
      tie == "5/6a" ~ "5",
      #tie == "5/6b" ~ "6",
      tie == "7/8" ~ "7",
      TRUE ~ tie
    ),
    tie = as.numeric(tie)
  )

## split AZ and TE from TN
networks$AZ <- networks$TN |>
  dplyr::filter(personid %in% (
    people_observations$TN |>
      dplyr::filter(location == "alakapuram") |>
      dplyr::pull(personid)
  ))

networks$TE <- networks$TN |>
  dplyr::filter(personid %in% (
    people_observations$TN |>
      dplyr::filter(location == "tenpatti") |>
      dplyr::pull(personid)
  ))

## split MC into MC and AP
networks$AP <- networks$MC |>
  dplyr::filter(personid %in% (
    people_observations$MC |>
      dplyr::filter(location == "2") |>
      dplyr::pull(personid)
  ))

networks$MC <- networks$MC |>
  dplyr::filter(personid %in% (
    people_observations$MC |>
      dplyr::filter(location == "1") |>
      dplyr::pull(personid)
  ))

## Defining some sets of sites based on who was asked which questions, used in
## nomination_counts_composition and concordance_within_su.
## In MA and KO, *all* questions seem to have been asked *either* of men or women
## specifically: men asked 1, 2, 6, 8, 9, 10; women asked 3, 4, 5, 7
## no exceptions to this in MA, a few exceptions to this in KO --  WRONG to be
## listing both genders in as having been presented all questions, so below we
## turn these two sites to NAs (meaning we ignore the handful of exceptions in KO)
## In MK, behav and info questions were asked of all respondents (56a, 56b, 78), but
## questions weren't gender-specific, so wouldn't fit in with those name gens.
## They did have questions which are used in the aggregate, but we can't label
## them within the men or women-specific questions. So used to construct agg, then
## removed from outputs using gender-specific name gen/layer
## In MG, men were asked Qs now coded as: 1, 2, 3, 4, 6, 8, 9, 10 (no female Qs)
## In MG, *married* women were asked Qs now coded as: 1:10 (for 1:4, 9, 10 this
## collapses questions posed about them and about their husbands)
## In MG, other women were asked 1, 2, 3, 4, 5, 7, 9, 10

gender_specific_sites <- c("IH", "WH", "PC", "WL", "AV", "DJ", "PE", "TZ", "AP", "MC", "PQ")
two_nine_ten <- c("PC", "IH", "WH", "WL")
nine_ten <- c("AV", "DJ", "MK")
nine <- c("BY")
ko_ma_men <- c(1, 2, 6, 8, 9, 10)
ko_ma_women <- c(3, 4, 5, 7)
mg_men <- c(1, 2, 3, 4, 6, 8, 9, 10)
mg_married_women <- 1:10
mg_other_women <- c(1, 2, 3, 4, 5, 7, 9, 10)
# for MG, also have mg_resp object created in process_data.R

sysfonts::font_add_google("Roboto Mono", "roboto", regular.wt = 400)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 300)

# ── Common look ups for plots and tables ────────────────────────────────────────────────

tie_meta <- data.frame(
  id        = c(1:12, "main_req", "req", "main_all", "all", "any"), # 11 = primary kin; 12 = "aggregate"; "main_req" = 1,3,5,6,7,8; "req" = 1,3,5,6,7,8,9,10; "main_all" = 1:10, "all" = truly all in edge list -- these last ones made in process_data and including exceptions like 5/6a, 7/8a and KM external IDs
  short     = c("loan_req", "loan_giv", "exchange_req", "exchange_giv",
                "f_behav", "m_behav", "f_info", "m_info", "govt", "external",
                "primary_kin", "sum",
                "main_req", "req", "main_all", "all", "any"), # these last are for su_alters and su_externals, generated in process_data
  label     = c("Borrow money", "Lend money", "Borrow goods", "Lend goods",
                "Behavioral support (women)", "Behavioral support (men)",
                "Informational/affective support (women)",
                "Informational/affective support (men)",
                "Government and NGO ties", "Support outside the community",
                "Primary kin", "Aggregate",
                "Aggregate", "Aggregate (with external NGs)", "Aggregate (with double-sampled and external NGs)", "Aggregate (all)", "Aggregate"),
  label_short = c("Borrow money", "Lend money", "Borrow goods", "Lend goods",
                "Women behav.", "Men behav.",
                "Women info.", "Men info.",
                "Gov/NGO", "External support",
                "Primary kin", "Aggregate",
                "Aggregate", "Aggregate (w/ ext)", "Aggregate (w/ DS and ext)", "Aggregate (all)", "Aggregate"),
  pair      = c("Monetary loan", "Monetary loan",
                "Household exchange", "Household exchange",
                "Behavioral support", "Behavioral support",
                "Informational/affective support", "Informational/affective support",
                "External support", "External support",
                "Primary kin", "Aggregate",
                "Aggregate", "Aggregate", "Aggregate", "Aggregate", "Aggregate"),
  colour    = c("#1f78b4", "#a6cee3", "#e31a1c", "#fb9a99",
                "#33a02c", "#b2df8a", "#ff7f00", "#fdbf6f",
                "#6a3d9a", "#cab2d6", "#8c510a", "#969696",
                "#969696", "#969696", "#969696", "#969696", "#969696"),
  stringsAsFactors = FALSE
)

pair_levels <- c("Monetary loan", "Household exchange",
                 "Behavioral support", "Informational/affective support",
                 "External support", "Primary kin", "Aggregate")

# integer id  -> label / short / pair / colour
tie_labels <- setNames(tie_meta$label, as.character(tie_meta$id))
tie_label_short  <- setNames(tie_meta$label_short, as.character(tie_meta$id))
tie_short <- setNames(tie_meta$short, as.character(tie_meta$id))

# label       -> pair / colour
tie_pairs <- setNames(tie_meta$pair, tie_meta$label)
tie_colours <- setNames(tie_meta$colour, tie_meta$label)
tie_colours_short <- setNames(tie_meta$colour, tie_meta$label_short)

# short name  -> label / colour (extra mappings, same palette)
short_to_label  <- setNames(tie_meta$label, tie_meta$short)
short_to_label_short <- setNames(tie_meta$label_short, tie_meta$short)
short_to_colour <- setNames(tie_meta$colour, tie_meta$short)
short_to_pair <- setNames(tie_meta$pair, tie_meta$short)

net_type_order  <- tie_meta$short[1:10]

# site colors/lines
site_colours <- c(
  "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
  "#a65628", "#f781bf", "#999999", "#1b9e77", "#d95f02"
)

site_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")

# Lookup for all sites present
n_sites   <- length(site_codes)

site_style <- data.frame(
  site     = site_codes,
  colour   = site_colours[  ((seq_len(n_sites) - 1) %% length(site_colours))   + 1],
  linetype = site_linetypes[((seq_len(n_sites) - 1) %/% length(site_colours))  + 1],
  stringsAsFactors = FALSE
)

site_colour_map   <- setNames(site_style$colour,   site_style$site)
site_linetype_map <- setNames(site_style$linetype, site_style$site)

## create new weighted edges that for each layer align with what we do in the aggregate "sum" network
layer_divisor_map <- list(
  loan_req     = "num_surv",
  loan_giv     = "num_surv",
  exchange_req = "num_surv",
  exchange_giv = "num_surv",
  f_behav      = "num_surv_female_q",
  f_info       = "num_surv_female_q",
  m_behav      = "num_surv_male_q",
  m_info       = "num_surv_male_q",
  govt         = "num_surv",
  external     = "num_surv"
)

for (site_name in names(su_nets_expanded)) {
  for (layer_name in names(layer_divisor_map)) {

    g <- su_nets_expanded[[site_name]][[layer_name]]
    if (is.null(g)) next

    divisor_attr <- layer_divisor_map[[layer_name]]

    divisor_vals <- igraph::vertex_attr(g, divisor_attr)

    if(site_name == "MG") {
      if(layer_name %in% c("loan_req", "loan_giv", "govt", "external")) {
        # for these questions, married women were additionally
        # asked second name generators *about their husbands*, giving them two chances to
        # name partners, so while there is only 1 respondent, then number of opportunities
        # to answer for those questions were 2, not one.
        # doing equivalent here as in process_data.R and analyse_data_and_create_networks.R
        mg_resp <- people_observations$MG |>
          dplyr::filter(personid %in% networks$MG$personid) |>
          dplyr::select(su_id, personid, gender, marital_status) |>
          dplyr::mutate(
            surveyed_loan_ext = dplyr::case_when(
              gender == "male" ~ 1,
              gender == "female" & !marital_status %in% "married" ~ 1,
              gender == "female" & marital_status %in% "married" ~ 2,
              TRUE ~ 0
            )
          )
        divisor_vals <- mg_resp$surveyed_loan_ext[match(
          igraph::vertex_attr(g, "name"),
          mg_resp$su_id
          )]
      }
    }

      divisor_vals[is.na(divisor_vals)] <- 0

    # Get the source (tail) node index for each edge
    edge_ends <- igraph::ends(g, igraph::E(g), names = FALSE)
    tail_idx <- edge_ends[, 1]

    divisor_per_edge <- divisor_vals[tail_idx]

    new_weight <- igraph::edge_attr(g, "weight") / divisor_per_edge

    igraph::edge_attr(g, "weight_per_resp") <- new_weight

    su_nets_expanded[[site_name]][[layer_name]] <- g
  }
}

# for ease, adding similarly named edge attr for sum network
for (site_name in names(su_nets_expanded)) {
    g <- su_nets_expanded[[site_name]]$sum
    if (is.null(g)) next
    igraph::edge_attr(g, "weight_per_resp") <- igraph::edge_attr(g, "weight")
    su_nets_expanded[[site_name]]$sum <- g
}

# DROP SHARING UNITS NOT SAMPLED FOR NETWORK QUESTIONS
for (site_name in names(su_nets_expanded)) {
  ref_g <- su_nets_expanded[[site_name]][[1]]
  sampled <- igraph::vertex_attr(ref_g, "sampled")

  drop_idx <- which(sampled == FALSE)  # NA-safe: NAs are excluded, not dropped
  if (length(drop_idx) == 0) next

  for (layer_name in names(su_nets_expanded[[site_name]])) {
    g <- su_nets_expanded[[site_name]][[layer_name]]
    if (is.null(g)) next
    su_nets_expanded[[site_name]][[layer_name]] <- igraph::delete_vertices(g, drop_idx)
  }
  # and equiv for one with self-loops
  for (layer_name in names(su_nets_expanded_withloops[[site_name]])) {
    g <- su_nets_expanded_withloops[[site_name]][[layer_name]]
    if (is.null(g)) next
    su_nets_expanded_withloops[[site_name]][[layer_name]] <- igraph::delete_vertices(g, drop_idx)
  }
}

## Function for common violin plot
make_violin_plot <- function(
  data,
  x_var,
  y_var,
  group_var,                    # categorical var driving fill/colour (always ends up on the visual y-axis)
  palette,                      # named vector, names matching levels of group_var
  label_var      = "site",
  jitter_width   = 0.2,
  flip           = TRUE,        # TRUE: x_var is categorical, coord_flip() applied. FALSE: y_var is already categorical
  x_label        = NULL,
  y_label        = NULL,
  x_limits       = NULL,        # passed to scale_x_continuous (drops data outside range)
  x_breaks       = NULL,        # vector or function, passed to scale_x_continuous
  coord_xlim     = NULL,        # passed to coord_cartesian(xlim = ...) (clips view, keeps data)
  facet_var      = NULL,
  hline          = NULL,        # geom_hline(yintercept = ...), pre-flip continuous axis
  vline          = NULL,        # geom_vline(xintercept = ...), for non-flipped continuous axis
  show_group_labels = TRUE,     # set FALSE to strip axis text/title on the categorical (visual y) axis
  base_size      = 13,
  legend         = FALSE
) {

  p <- data |>
    ggplot(aes(
      x = .data[[x_var]],
      y = .data[[y_var]],
      fill = .data[[group_var]],
      colour = .data[[group_var]]
    )) +
    geom_violin(
      alpha = 0.4,
      linewidth = 0.5,
      scale = "width",
      draw_quantiles = 0.5
    ) +
    geom_label(
      aes(label = .data[[label_var]]),
      family = "roboto",
      position = ggplot2::position_jitter(width = jitter_width),
      size = 1.5,
      alpha = 0.5,
      linewidth = 0.25,
      colour = "white",
      label.padding = unit(0.3, "lines"),
      label.r = unit(0.6, "lines")
    ) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_colour_manual(values = palette, guide = "none") +
    labs(x = x_label, y = y_label) +
    theme_classic(base_size = base_size) +
    theme(legend.position = if (legend) "right" else "none")

  if (flip) p <- p + coord_flip()

  if (!is.null(x_limits) || !is.null(x_breaks)) {
  p <- p + scale_x_continuous(limits = x_limits, breaks = x_breaks, oob = scales::oob_squish)
  }

  if (!is.null(coord_xlim)) {
    p <- p + coord_cartesian(xlim = coord_xlim)
  }

  if (!is.null(facet_var)) {
    p <- p + facet_wrap(stats::as.formula(paste0("~", facet_var)))
  }

  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline)
  if (!is.null(vline)) p <- p + geom_vline(xintercept = vline)

  if (!show_group_labels) {
    p <- p + theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
  }

  p
}

integer_breaks <- function(n = 5) {
  function(x) {
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2])) return(numeric(0))
    if (rng[1] == rng[2]) return(rng[1])

    raw_breaks <- pretty(rng, n)
    step <- diff(raw_breaks)[1]

    if (step < 1) {
      breaks <- seq(ceiling(rng[1]), floor(rng[2]), by = 1)
      if (length(breaks) > n + 2) {
        by <- ceiling(length(breaks) / n)
        breaks <- breaks[seq(1, length(breaks), by = by)]
      }
    } else {
      breaks <- unique(round(raw_breaks))
    }

    breaks
  }
}

## functions used for network_overlap and pca:

get_adj_weighted <- function(net, all_nodes) {
  if (is.null(net)) {
    return(matrix(0, length(all_nodes), length(all_nodes),
                  dimnames = list(all_nodes, all_nodes)))
  }
  attr  <- if ("weight" %in% edge_attr_names(net)) "weight" else NULL
  m     <- as.matrix(igraph::as_adjacency_matrix(net, attr = attr, sparse = FALSE))
  diag(m) <- 0
  m
}

make_png <- function(filename, width, height, expr) {
  png(
    file.path(EndowDropbox, "Figures", "Networks", filename),
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  print(expr)
  dev.off()
}

make_pdf <- function(filename, width, height, data_sites, plot_fn, title_str) {
  sites_in_data <- sort(unique(data_sites$site))
  n_pages_local <- ceiling(length(sites_in_data) / facets_per_page)

  pdf(
    file.path(EndowDropbox, "Figures", "Networks", filename),
    width = width,
    height = height
  )
  for (p in seq_len(n_pages_local)) {
    page_sites <- sites_in_data[
      ((p - 1) * facets_per_page + 1):min(
        p * facets_per_page,
        length(sites_in_data)
      )
    ]
    page_data <- dplyr::filter(data_sites, site %in% page_sites)
    print(plot_fn(page_data, page_sites, p, title_str))
  }
  dev.off()
}