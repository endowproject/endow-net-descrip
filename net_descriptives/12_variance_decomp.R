## ============================================================================
## 12_variance_decomp.R
## ----------------------------------------------------------------------------
##  ANALYSIS: descriptive variance decomposition of cross-community
##  heterogeneity in network structure.
##
##  ESTIMAND
##  --------
##  For a network statistic y measured on the binary network of community s in
##  support layer l, we decompose the observed variation in y into (i) a part
##  associated with the support layer, (ii) a part shared by all layers within a
##  community, and (iii) residual community x layer variation. The target
##  quantities are the between-community standard deviation of y net of
##  systematic layer differences, the share of non-layer variation attributable
##  to the community, and how much of that between-community variance survives
##  adjustment for observable measurement and sampling differences. These are
##  descriptive variance shares for the 50 communities observed, not causal
##  attributions and not estimates for a hypothetical super-population of
##  cultures.
##
##  MODELS (per statistic, on the eight constituent prompt layers)
##  --------------------------------------------------------------
##   M1  baseline      y ~ layer + (1 | site)
##   M2  adjusted      y ~ layer + measurement controls + (1 | site)
##
##  Layer enters as a saturated fixed effect (8 levels). The community enters as
##  a random intercept because communities are observed repeatedly across
##  layers, which is what lets community-level and residual variation be
##  separated. Each of the eight prompts retains its own layer fixed effect.
##
##  THE AGGREGATE NETWORK is analysed separately. Its
##  edge set is by construction the union of six of the constituent prompts,
##  so its edge count and every density-like statistic are mechanically
##  larger. Entering it as a ninth "layer" would let a purely mechanical
##  construction effect load onto the layer fixed effects and inflate the
##  residual. We therefore report its between-community dispersion as a
##  single-level quantity and compare it with the multilevel results.
##
##  Outputs: extensions/output/tables/A_tab*, extensions/output/figures/A_fig*
##
## Depends on: 00_setup.R
## ============================================================================

library(scales)
library(lme4)

set.seed(20260915)
NSIM_BOOT <- 1000 # parametric bootstrap replicates for variance-component CIs

output_dir <- file.path(EndowDropbox, "Output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data

source("net_descriptives/00_setup.R")

load(file.path(EndowDerivedData, "network_summaries_long.rdata"))
# this file is created by network_summary_stats; alternatively could source

load(file.path(EndowDerivedData, "su_sampling.rdata"))
# this file is created by respondent_coverage; alternatively could source
# called su_summary once loaded!

load(file.path(EndowDerivedData, "external_alters.rdata"))
# this file is created by nomination_counts_composition; alternatively could source
# called alters once loaded!

# convenience functions
emp_logit <- function(p, n) {
  k <- p * n
  log((k + 0.5) / (n - k + 0.5))
}

write_table <- function(df, name, digits = 3, caption = NULL) {
  write.csv(df, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  if (!is.null(caption)) {
    message(sprintf("  [%s] %s", name, caption))
  }
  invisible(df)
}

save_fig <- function(plot, name, width, height) {
  ggsave(
    file.path(output_dir, paste0(name, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
  invisible(plot)
}

panel_layers <- levels(network_summaries_long$network)[c(1:8, 12)]

layer_resp_col_map <- c(
  loan_req = "mean_num_surveyed",
  loan_giv = "mean_num_surveyed",
  exchange_req = "mean_num_surveyed",
  exchange_giv = "mean_num_surveyed",
  f_behav = "mean_num_surveyed_fq",
  f_info = "mean_num_surveyed_fq",
  m_behav = "mean_num_surveyed_mq",
  m_info = "mean_num_surveyed_mq",
  sum = "mean_num_surveyed"
)

resp_long <- su_summary |>
  dplyr::select(
    site,
    mean_num_surveyed,
    mean_num_surveyed_fq,
    mean_num_surveyed_mq
  ) |>
  tidyr::pivot_longer(
    -site,
    names_to = "resp_col",
    values_to = "mean_resp_layer"
  )

alters <- alters |>
  dplyr::mutate(network = tie_short[alters$tie]) |>
  dplyr::mutate(network = recode(network, "main_req" = "sum"))


panel <- network_summaries_long |>
  dplyr::filter(network %in% panel_layers) |>
  dplyr::mutate(
    layer = as.character(network),
    is_aggregate = as.character(network) == "sum",
    keep = "all"
  ) |>
  dplyr::mutate(resp_col = layer_resp_col_map[layer]) |>
  dplyr::left_join(resp_long, by = c("site", "resp_col")) |>
  dplyr::select(-resp_col) |>
  dplyr::left_join(
    su_summary |> dplyr::select(site, n_su_sampled = su_sampled, prop_sampled),
    by = "site"
  ) |>
  dplyr::left_join(alters, by = c("site", "layer" = "network"))

dat <- panel %>% filter(!is_aggregate)
agg <- panel %>% filter(is_aggregate)

dat <- dat %>%
  dplyr::mutate(
    layer = factor(layer, levels = unique(dat$layer)),
    site = factor(site),
    ## transformations, see OUTCOMES below
    log_mean_degree = log(mean_degree),
    log_mean_degree_w = log(mean_degree_w),
    logit_recip = emp_logit(reciprocity, recip_denom),
    logit_globalclust = emp_logit(global_clust, triple_denom),
    logit_localclust = emp_logit(av_local_clust, triple_denom),
    logit_indeg = emp_logit(in_deg_centr, n_nodes),
    ## centred measurement controls (centring makes the intercept and the
    ## layer contrasts refer to a community of average measurement profile)
    c_log_n = scale(log(n_su_sampled), scale = FALSE)[, 1],
    c_prop_samp = scale(prop_sampled, scale = FALSE)[, 1],
    c_log_resp = scale(log(mean_resp_layer), scale = FALSE)[, 1],
    ## prop_ext (was prop_ext_req) -- see header note: now varies by
    ## layer as well as by site.
    c_prop_ext = scale(prop_ext, scale = FALSE)[, 1]
  )

## ---------------------------------------------------------------------------
##  OUTCOMES
##  -------
##  Primary (binary networks):
##    log_mean_degree  mean out-degree on the log scale. Degree is a positive,
##                     right-skewed, scale-varying count average; on the log
##                     scale layer fixed effects are multiplicative and the
##                     community random intercept is a proportional community
##                     effect, which is both better behaved and easier to read.
##    reciprocity      share of tied ordered pairs that are reciprocated.
##    gloabl_clust     global clustering coefficient
##    av_local_clust   mean node-level local clustering
##    indeg_centr      Freeman in-degree centralization
##    degree_assort    Degree assortativity (ranges [-1,1])
##
##  Robustness for the four [0,1] bounded statistics: the same models on the
##  empirical-logit scale, which respects the (0,1) support and is finite for
##  the exact zeros.
##
##  Robustness for degree and clustering: the weighted analogues
##  (mean out-strength; Barrat weighted local clustering), where weights are the
##  layers' own nomination counts.
##
##  NOT used as primary outcomes:
##    density        = mean_degree / (n - 1), mechanically size-dependent, so it
##                     would largely re-express the size control.
##    av_path_length   the prompt layers are often disconnected, and igraph
##                     averages over reachable pairs only,
##                     so the statistic is not comparable across communities.
##                     Reported for the aggregate network only, where the median
##                     reachable share is 0.93.
## ---------------------------------------------------------------------------

OUTCOMES <- tibble::tribble(
  ~key               , ~var                , ~label                                  , ~family , ~role        ,
  "mean_degree"      , "log_mean_degree"   , "Mean degree"                           , "gauss" , "primary"    ,
  "reciprocity"      , "reciprocity"       , "Reciprocity"                           , "gauss" , "primary"    ,
  "global_clust"     , "global_clust"      , "Global clustering"                     , "gauss" , "primary"    ,
  "av_local_clust"   , "av_local_clust"    , "Average local clustering"              , "gauss" , "primary"    ,
  "in_deg_centr"     , "in_deg_centr"      , "In-degree centralization"              , "gauss" , "primary"    ,
  "deg_assort"       , "deg_assort"        , "Degree assortativity"                  , "gauss" , "primary"    ,
  "reciprocity_lgt"  , "logit_recip"       , "Reciprocity (emp. logit)"              , "gauss" , "scale_rob"  ,
  "globalclust_lgt"  , "logit_globalclust" , "Global clustering (emp. logit)"        , "gauss" , "scale_rob"  ,
  "localclust_lgt"   , "logit_localclust"  , "Average local clustering (emp. logit)" , "gauss" , "scale_rob"  ,
  "indeg_lgt"        , "logit_indeg"       , "In-degree centralization (emp. logit)" , "gauss" , "scale_rob"  ,
  "mean_degree_w"    , "log_mean_degree_w" , "Mean degree (weighted)"                , "gauss" , "weight_rob" ,
  "av_local_clust_w" , "av_local_clust_w" , "Average local clustering (weighted)"    , "gauss" , "weight_rob"
)

CONTROLS <- "c_log_n + c_log_resp + c_prop_ext" # CUTTING PROP_SAMP BECAUSE NET IS ONLY OF SAMPLED SUs

## n_su_sampled: larger community size can decrease density (obviously), but also reciprocity, transitivity, etc
## prop_sampled: NOT used as everything here is on the network of SAMPLED sharing units (and our coverage is very high)
## mean_resp_layer: in communities where more people were surveyed, more opportunities to name alters and increase density, etc
## prop_ext: in communities where more nominations went to those outside the community, will have greater sparsity within community

## ---------------------------------------------------------------------------
##  Fitting helpers
## ---------------------------------------------------------------------------
fit_lmm <- function(formula, d) {
  suppressMessages(suppressWarnings(
    lmer(
      formula,
      data = d,
      REML = TRUE,
      control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
    )
  ))
}

# Variance components from a fitted lmer, as a named vector of variances.
vc_vec <- function(m) {
  v <- as.data.frame(VarCorr(m))
  setNames(v$vcov, ifelse(v$grp == "Residual", "residual", as.character(v$grp)))
}

# Between-community variance share (ICC): community variance divided by total
# non-layer (i.e. random + residual) variance.
icc_site <- function(m) {
  v <- vc_vec(m)
  unname(v[["site"]] / sum(v))
}

# Parametric-bootstrap CI for an arbitrary scalar function of the fit.
# Parallelised across cores; lmer refits are cheap individually but there are
# thousands of them.
NCPUS <- max(1L, parallel::detectCores() - 2L)

## Make our own cluster so we can export what FUN needs. Without this, every
## bootstrap replicate silently errors ("could not find function 'vc_vec'")
## and boot_ci() ends up taking quantile() of zero finite values -- the [NA,
## NA] interval above. The point estimate is unaffected since icc1/sd_site
## are computed directly on m1, not through the bootstrap.
cl <- parallel::makeCluster(NCPUS)
parallel::clusterEvalQ(cl, library(lme4))
parallel::clusterExport(cl, varlist = c("vc_vec", "icc_site"))

boot_ci <- function(m, fn, nsim = NSIM_BOOT, probs = c(0.025, 0.975)) {
  bb <- suppressMessages(suppressWarnings(
    bootMer(
      m,
      FUN = fn,
      nsim = nsim,
      type = "parametric",
      use.u = FALSE,
      parallel = "snow",
      cl = cl
    )
  ))
  v <- bb$t[, 1]
  v <- v[is.finite(v)]
  c(
    lo = unname(quantile(v, probs[1])),
    hi = unname(quantile(v, probs[2])),
    n_ok = length(v)
  )
}

# Convergence / boundary diagnostics for one fit.
diag_fit <- function(m) {
  msgs <- m@optinfo$conv$lme4$messages
  data.frame(
    singular = isSingular(m, tol = 1e-4),
    n_warn = length(msgs),
    warn = if (length(msgs)) paste(msgs, collapse = "; ") else "",
    n_obs = nobs(m),
    n_site = length(unique(m@frame$site)),
    stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------------
##  Main loop
## ---------------------------------------------------------------------------
results <- list()
diagnostics <- list()
site_effects <- list()

# May take some time (~20 min) if on computer with few cores
for (i in seq_len(nrow(OUTCOMES))) {
  o <- OUTCOMES[i, ]
  d <- dat %>% filter(is.finite(.data[[o$var]]))
  n_drop <- nrow(dat) - nrow(d)

  f_base <- as.formula(paste(o$var, "~ layer + (1 | site)"))
  f_adj <- as.formula(paste(o$var, "~ layer +", CONTROLS, "+ (1 | site)"))
  m1 <- fit_lmm(f_base, d)
  m2 <- fit_lmm(f_adj, d)

  v1 <- vc_vec(m1)
  v2 <- vc_vec(m2)

  icc1 <- icc_site(m1)
  icc2 <- icc_site(m2)
  ci_sd1 <- boot_ci(m1, function(x) sqrt(vc_vec(x)[["site"]]))
  ci_icc1 <- boot_ci(m1, icc_site)
  ci_sd2 <- boot_ci(m2, function(x) sqrt(vc_vec(x)[["site"]]))
  ci_icc2 <- boot_ci(m2, icc_site)

  results[[o$key]] <- data.frame(
    key = o$key,
    statistic = o$label,
    role = o$role,
    n_obs = nrow(d),
    n_dropped = n_drop,
    n_sites = length(unique(d$site)),
    ## M1 baseline
    sd_site = sqrt(v1[["site"]]),
    sd_site_lo = ci_sd1[["lo"]],
    sd_site_hi = ci_sd1[["hi"]],
    sd_resid = sqrt(v1[["residual"]]),
    icc = icc1,
    icc_lo = ci_icc1[["lo"]],
    icc_hi = ci_icc1[["hi"]],
    ## M2 measurement-adjusted
    sd_site_adj = sqrt(v2[["site"]]),
    sd_site_adj_lo = ci_sd2[["lo"]],
    sd_site_adj_hi = ci_sd2[["hi"]],
    sd_resid_adj = sqrt(v2[["residual"]]),
    icc_adj = icc2,
    icc_adj_lo = ci_icc2[["lo"]],
    icc_adj_hi = ci_icc2[["hi"]],
    ## attenuation of the between-community VARIANCE after adjustment
    pct_var_attenuation = 100 * (1 - v2[["site"]] / v1[["site"]]),
    stringsAsFactors = FALSE
  )

  diagnostics[[o$key]] <- bind_rows(
    cbind(
      key = o$key,
      statistic = o$label,
      model = "M1 baseline",
      diag_fit(m1)
    ),
    cbind(
      key = o$key,
      statistic = o$label,
      model = "M2 adjusted",
      diag_fit(m2)
    )
  )

  ## community random intercepts (BLUPs) for the influence diagnostic / figure
  re <- ranef(m1)$site
  site_effects[[o$key]] <- data.frame(
    key = o$key,
    statistic = o$label,
    site = rownames(re),
    blup = re[, 1],
    stringsAsFactors = FALSE
  )

  message(sprintf(
    "Model: %-38s sd_site=%.3f  ICC=%.2f [%.2f,%.2f]  ICC_adj=%.2f  atten=%.0f%%",
    o$label,
    sqrt(v1[["site"]]),
    icc1,
    ci_icc1[["lo"]],
    ci_icc1[["hi"]],
    icc2,
    100 * (1 - v2[["site"]] / v1[["site"]])
  ))
}

res <- bind_rows(results)
dg <- bind_rows(diagnostics)
se <- bind_rows(site_effects)

write.csv(
  res,
  file.path(output_dir, "A_variance_components_all.csv"),
  row.names = FALSE
)
write.csv(
  dg,
  file.path(output_dir, "A_model_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  se,
  file.path(output_dir, "A_site_random_intercepts.csv"),
  row.names = FALSE
)

if (exists("cl")) {
  parallel::stopCluster(cl)
}

role_labels <- c(
  primary = "Primary (binary networks)",
  scale_rob  = "Robustness: empirical logit",
  weight_rob = "Robustness: weighted edges")

res_long <- dplyr::bind_rows(
  res |>
    dplyr::transmute(
      statistic, role,
      order = dplyr::row_number(),
      model = "Base",
      sd_summary = sprintf("%.2f [%.2f, %.2f]", sd_site, sd_site_lo, sd_site_hi),
      resid = sd_resid,
      icc_summary = sprintf("%.2f [%.2f, %.2f]", icc, icc_lo, icc_hi)
    ),
  res |>
    dplyr::transmute(
      statistic, role,
      order = dplyr::row_number(),
      model = "Adj.",
      sd_summary = sprintf("%.2f [%.2f, %.2f]", sd_site_adj, sd_site_adj_lo, sd_site_adj_hi),
      resid = sd_resid_adj,
      icc_summary = sprintf("%.2f [%.2f, %.2f]", icc_adj, icc_adj_lo, icc_adj_hi)
    )
) |>
  dplyr::arrange(order, factor(model, levels = c("Base", "Adj."))) |>
  dplyr::mutate(statistic = ifelse(duplicated(statistic), "", statistic))

tbl <- res_long |>
  dplyr::select(statistic, model, sd_summary, resid, icc_summary) |>
  kableExtra::kbl(
    format = "latex",
    booktabs = TRUE,
    longtable = FALSE,
    escape = FALSE,
    align = c("l", "l", "c", "c", "c"),
    digits = 2,
    col.names = c("", "", "Between-community SD [95\\% CI]", "Residual SD", "ICC [95\\% CI]"),
    caption = paste(
      "Variance decomposition of network statistics across the ENDOW",
      "communities and the eight main network layers. Specifically,",
      "eight primary (binary) statistics and two robustness sets covering",
      "the same statistics where applicable -- on the bounded-scale",
      "empirical logit, and using their weighted-edge analogues.",
      "For each statistic, the `Base' row reports a model with layer fixed",
      "effects and a community random intercept; the `Adj.' row reports the",
      "same model additionally controlling for (centred) log community size,",
      "log mean respondents per prompt,",
      "and the share of nominations directed outside the community.",
      "`Between-community SD' is the standard deviation of the community random",
      "intercept, net of layer differences; `Residual SD' is the standard",
      "deviation of the community-by-layer residual; `ICC' is the",
      "intraclass correlation coefficient, capturing the",
      "between-community share of their combined variance. Bracketed values",
      "are 95\\% parametric-bootstrap intervals."
    ),
    label = "variance_decomp",
    linesep = ""
  ) |>
  kableExtra::kable_styling(
    latex_options = c("hold_position", "repeat_header"),
    font_size = 8
  )

## adding headers for the different model types
role_group_sizes <- res |>
  dplyr::count(role, .drop = FALSE) |>
  dplyr::mutate(n_rows = n * 2)

row_end   <- cumsum(role_group_sizes$n_rows)
row_start <- c(1, head(row_end, -1) + 1)

for (i in seq_along(role_labels)) {
  tbl <- kableExtra::pack_rows(tbl, role_labels[i], row_start[i], row_end[i],
                                bold = TRUE, latex_gap_space = "0.5em")
}

writeLines(tbl, file.path(tablepath, "variance_decomp.tex"))


## ---------------------------------------------------------------------------
##  Aggregate network, analysed separately
## ---------------------------------------------------------------------------
agg_stats <- data.frame(
  statistic = c(
    "Mean degree",
    "Reciprocity",
    "Global clustering",
    "Average local clustering",
    "In-degree centralization",
    "Degree assortativity",
    "Average path length"
  ),
  mean = c(
    mean(log(agg$mean_degree)),
    mean(agg$reciprocity),
    mean(agg$global_clust),
    mean(agg$av_local_clust),
    mean(agg$in_deg_centr),
    mean(agg$deg_assort),
    mean(agg$avg_path_length)
  ),
  sd_between_community = c(
    sd(log(agg$mean_degree)),
    sd(agg$reciprocity),
    sd(agg$global_clust),
    sd(agg$av_local_clust),
    sd(agg$in_deg_centr),
    sd(agg$deg_assort),
    sd(agg$avg_path_length)
  ),
  p10 = c(
    quantile(log(agg$mean_degree), .1),
    quantile(agg$reciprocity, .1),
    quantile(agg$global_clust, .1),
    quantile(agg$av_local_clust, .1),
    quantile(agg$in_deg_centr, .1),
    quantile(agg$deg_assort, .1),
    quantile(agg$avg_path_length, .1)
  ),
  p90 = c(
    quantile(log(agg$mean_degree), .9),
    quantile(agg$reciprocity, .9),
    quantile(agg$global_clust, .9),
    quantile(agg$av_local_clust, .9),
    quantile(agg$in_deg_centr, .9),
    quantile(agg$deg_assort, .9),
    quantile(agg$avg_path_length, .9)
  ),
  stringsAsFactors = FALSE
)

## Measurement-adjusted residual SD for the aggregate: OLS on the same
## controls, one observation per community.
agg_adj <- agg %>% mutate(log_mean_degree = log(mean_degree))
agg_resid_sd <- sapply(
  c(
    "log_mean_degree",
    "reciprocity",
    "global_clust",
    "av_local_clust",
    "in_deg_centr",
    "deg_assort",
    "avg_path_length"
  ),
  function(v) {
    f <- as.formula(paste(
      v,
      "~ log(n_su_sampled) + log(mean_resp_layer) + prop_ext" # CUTTING PROP_SAMPLED, SINCE MEASURES CALCULATED ON NET OF SAMPLED SUs
    ))
    m <- lm(f, data = agg_adj)
    summary(m)$sigma
  }
)
agg_stats$sd_after_controls <- c(agg_resid_sd)
agg_stats$pct_var_attenuation <- 100 *
  (1 -
    (agg_stats$sd_after_controls /
      agg_stats$sd_between_community)^2)

message("Aggregate network, between-community dispersion")
print(agg_stats, row.names = FALSE, digits = 3)

agg_stats |>
  dplyr::mutate(
    range_summary = sprintf("%.2f [%.2f, %.2f]", mean, p10, p90)
  ) |>
  dplyr::select(
    statistic,
    range_summary,
    sd_between_community,
    sd_after_controls
  ) |>
  kableExtra::kbl(
    format = "latex",
    booktabs = TRUE,
    longtable = FALSE,
    escape = FALSE,
    align = c("l", "c", "c", "c"),
    digits = c(NA, NA, 2, 2),
    col.names = c("", "Mean [P10, P90]", "Between-community SD",
                  "SD after controls"),
    caption = paste(
      "Between-community dispersion of the aggregate network,",
      "one observation per community.",
      "`Mean [P10, P90]' summarises the empirical distribution of the",
      "statistic across communities (showing the 10th to 90th percentile",
      "in brackets). `Between-community SD' is its standard",
      "deviation. `SD After Controls' is the residual standard deviation",
      "from an OLS regression of the statistic on (log) community size,",
      "(log) mean respondents per prompt, and the share of nominations",
      "directed outside the community."
    ),
    label = "agg_variance_decomp",
    linesep = ""
  ) |>
  kableExtra::kable_styling(
    latex_options = c("hold_position", "repeat_header"),
    font_size = 8
  ) |>
  writeLines(file.path(tablepath, "agg_variance_decomp.tex"))


## ---------------------------------------------------------------------------
##  Influence diagnostic: leave-one-community-out refit of the ICC
## ---------------------------------------------------------------------------
loo_path <- file.path(output_dir, "A_loo_icc.csv")
loo <- map_dfr(OUTCOMES$key[OUTCOMES$role == "primary"], function(k) {
  o <- OUTCOMES[OUTCOMES$key == k, ]
  d <- dat %>% filter(is.finite(.data[[o$var]]))
  f <- as.formula(paste(o$var, "~ layer + (1 | site)"))
  full <- icc_site(fit_lmm(f, d))
  map_dfr(unique(as.character(d$site)), function(s) {
    m <- fit_lmm(f, d %>% filter(site != s))
    data.frame(
      key = k,
      statistic = o$label,
      dropped_site = s,
      icc = icc_site(m),
      icc_full = full,
      stringsAsFactors = FALSE
    )
  })
})

write.csv(loo, loo_path, row.names = FALSE)
loo_summary <- loo %>%
  group_by(key, statistic) %>%
  summarise(
    icc_full = first(icc_full),
    icc_loo_min = min(icc),
    icc_loo_max = max(icc),
    most_influential_site = dropped_site[which.max(abs(icc - icc_full))],
    max_abs_change = max(abs(icc - icc_full)),
    .groups = "drop"
  )

message("Leave-one-community-out influence on the ICC")
print(as.data.frame(loo_summary), row.names = FALSE, digits = 3)

## ---------------------------------------------------------------------------
##  Sparsity gradient: sensitivity of the baseline ICC to sparse layers
## ---------------------------------------------------------------------------
##  Sampling error in statistics computed on sparse layers can inflate the
##  residual and attenuate the estimated community share if that error is not
##  shared across layers within a community. Refit M1 on progressively denser
##  community-layer cells to examine sensitivity to this potential source of
##  noise.
##
##  This is a diagnostic, not a preferred estimate: restricting on tie count
##  selects on both the community and the prompt, so the higher thresholds
##  describe a different (and smaller) set of communities. Report the number
##  of retained communities and whether the fit is singular.
## ---------------------------------------------------------------------------
sparsity_grad <- map_dfr(OUTCOMES$key[OUTCOMES$role == "primary"], function(k) {
  o <- OUTCOMES[OUTCOMES$key == k, ]
  map_dfr(c(0, 50, 100, 150), function(thr) {
    d <- dat %>% filter(is.finite(.data[[o$var]]), n_edges >= thr)
    if (length(unique(as.character(d$site))) < 20) {
      return(NULL)
    }
    m <- fit_lmm(as.formula(paste(o$var, "~ layer + (1 | site)")), d)
    data.frame(
      statistic = o$label,
      min_ties = thr,
      n_cells = nrow(d),
      n_communities = length(unique(as.character(d$site))),
      singular = isSingular(m, tol = 1e-4),
      icc_two_level = icc_site(m),
      stringsAsFactors = FALSE
    )
  })
})

message("Sparsity gradient in the baseline ICC")
print(sparsity_grad, row.names = FALSE, digits = 2)
write_table(
  sparsity_grad,
  "A_tab5_sparsity_gradient_decomposition",
  digits = 3,
  caption = "Baseline ICC refitted on progressively denser community-layer cells"
)

## Companion figure: how the baseline ICC changes as sparse cells are dropped.
pGrad <- sparsity_grad %>%
  ggplot(aes(x = min_ties, y = icc_two_level)) +
  geom_line(linewidth = 0.5, colour = "#1f4e79") +
  geom_point(size = 2.1, colour = "#1f4e79") +
  facet_wrap(~statistic, nrow = 1) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = "Minimum ties per community-layer cell",
    y = "Between-community variance share (ICC)",
    title = "Sensitivity of the baseline ICC to sparse network layers",
    subtitle = paste(
      "M1 refitted at each threshold; the ICC divides community variance by",
      "community plus residual variance.\nHigher thresholds change the mix of",
      "communities and prompts retained."
    )
  ) +
  theme_classic()
save_fig(pGrad, "A_fig4_sparsity_gradient", width = 12, height = 4.6)

## ---------------------------------------------------------------------------
##  Tables
## ---------------------------------------------------------------------------

## Superceded with table to main folder above
tab_main <- res %>%
  filter(role == "primary") %>%
  transmute(
    Statistic = statistic,
    `Sites` = n_sites,
    `Site-layer obs.` = n_obs,
    `Between-community SD` = sd_site,
    `Between-comm. SD 95% CI low` = sd_site_lo,
    `Between-comm. SD 95% CI high` = sd_site_hi,
    `Residual SD` = sd_resid,
    `ICC (community share)` = icc,
    `ICC 95% CI low` = icc_lo,
    `ICC 95% CI high` = icc_hi,
    `Adjusted between-comm. SD` = sd_site_adj,
    `Adjusted ICC` = icc_adj,
    `Variance attenuation (%)` = pct_var_attenuation
  )
write_table(
  tab_main,
  "A_tab1_variance_components_primary",
  digits = 3,
  caption = "Variance decomposition of binary network statistics, 8 prompt layers"
)

tab_rob <- res %>%
  filter(role != "primary") %>%
  transmute(
    Statistic = statistic,
    Type = ifelse(
      role == "scale_rob",
      "Scale (empirical logit)",
      "Weighted edges"
    ),
    `Between-community SD` = sd_site,
    `Residual SD` = sd_resid,
    ICC = icc,
    `ICC 95% CI low` = icc_lo,
    `ICC 95% CI high` = icc_hi,
    `Adjusted ICC` = icc_adj,
    `Variance attenuation (%)` = pct_var_attenuation
  )
write_table(
  tab_rob,
  "A_tab3_robustness_scale_and_weights",
  digits = 3,
  caption = "Robustness: bounded-scale transformation and weighted edge definitions"
)

write_table(
  agg_stats,
  "A_tab4_aggregate_network_separately",
  digits = 3,
  caption = "Aggregate (composite) support network: between-community dispersion"
)
write_table(
  loo_summary,
  "A_tabS2_leave_one_community_out_icc",
  digits = 3,
  caption = "Leave-one-community-out influence on the ICC"
)
write_table(
  dg,
  "A_tabS3_model_diagnostics",
  digits = 0,
  caption = "Convergence and boundary diagnostics"
)

## ---------------------------------------------------------------------------
##  Figure A1 (manuscript): how much heterogeneity is between communities?
## ---------------------------------------------------------------------------
lev <- res %>% filter(role == "primary") %>% arrange(icc) %>% pull(statistic)

p1 <- res %>%
  filter(role == "primary") %>%
  mutate(statistic = factor(statistic, levels = lev)) %>%
  select(
    statistic,
    Baseline = icc,
    baseline_lo = icc_lo,
    baseline_hi = icc_hi,
    Adjusted = icc_adj,
    adj_lo = icc_adj_lo,
    adj_hi = icc_adj_hi
  ) %>%
  pivot_longer(c(Baseline, Adjusted), names_to = "spec", values_to = "icc") %>%
  mutate(
    lo = ifelse(spec == "Baseline", baseline_lo, adj_lo),
    hi = ifelse(spec == "Baseline", baseline_hi, adj_hi),
    spec = factor(spec, levels = c("Baseline", "Adjusted"))
  ) %>%
  ggplot(aes(x = icc, y = statistic, colour = spec)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_linerange(
    aes(xmin = lo, xmax = hi),
    position = position_dodge(width = 0.55),
    linewidth = 0.6
  ) +
  geom_point(position = position_dodge(width = 0.55), size = 2.4) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(
    values = c(Baseline = "#1f4e79", Adjusted = "#c0504d"),
    name = NULL
  ) +
  labs(
    x = "Between-community variance share (ICC)",
    y = NULL,
    #title = "Between-community share of variation",
    #subtitle = "Two-level model; baseline vs measurement-adjusted"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "between_site_variance.png"),
  plot = p1,
  width  = 7,
  height = 4
  )

save_fig(
  p1,
  "A_fig1_between_community_heterogeneity",
  width = 7,
  height = 4
)

## ---------------------------------------------------------------------------
##  Figure A2 (robustness): scale and weighted-edge definitions
## ---------------------------------------------------------------------------
pair_map <- tibble::tribble(
  ~base             , ~alt              , ~pair_label                , ~kind,
  "mean_degree"     , "mean_degree_w"   , "Mean degree"              , "Binary vs weighted",
  "in_deg_centr"    , "indeg_lgt"       , "In-degree centralization" , "Raw vs empirical logit",
  "av_local_clust"  , "av_local_clust_w", "Average local clustering" , "Binary vs weighted",
  "av_local_clust"  , "localclust_lgt"  , "Average local clustering" , "Raw vs empirical logit",
  "global_clust"    , "globalclust_lgt" , "Global clustering"        , "Raw vs empirical logit",
  "reciprocity"     , "reciprocity_lgt" , "Reciprocity"              , "Raw vs empirical logit"
)

rob_long <- pair_map %>%
  mutate(
    icc_base = res$icc[match(base, res$key)],
    lo_base = res$icc_lo[match(base, res$key)],
    hi_base = res$icc_hi[match(base, res$key)],
    icc_alt = res$icc[match(alt, res$key)],
    lo_alt = res$icc_lo[match(alt, res$key)],
    hi_alt = res$icc_hi[match(alt, res$key)]
  ) %>%
  select(pair_label, kind, icc_base, lo_base, hi_base, icc_alt, lo_alt, hi_alt) %>%
  pivot_longer(
    -c(pair_label, kind),
    names_to = c(".value", "spec"),
    names_pattern = "(.*)_(base|alt)"
  ) %>%
  mutate(
    spec = ifelse(
      spec == "base",
      "Baseline",
      "Alternative"
    ),
    spec = factor(
      spec,
      levels = c("Baseline", "Alternative")
    ),
    pair_label = factor(
      pair_label,
      levels = lev
    )
  )

pRob <- ggplot(rob_long, aes(x = icc, y = pair_label, colour = spec)) +
  geom_linerange(
    aes(xmin = lo, xmax = hi),
    position = position_dodge(width = 0.5),
    linewidth = 0.6
  ) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(position = position_dodge(width = 0.5), size = 2.2) +
  facet_grid(kind ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = c("#1f4e79", "#e8a33d"), name = NULL) +
  labs(
    x = "Between-community variance share (ICC)",
    y = NULL,
    #title = "Robustness of the between-community share to measure definition"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "between_site_variance_alt.png"),
  plot = pRob,
  width  = 7,
  height = 5
  )

save_fig(pRob, "A_fig2_robustness_weighted_and_scale", width = 8, height = 6)

## ---------------------------------------------------------------------------
##  Figure A3 (SI): community random intercepts (caterpillar), key statistics
## ---------------------------------------------------------------------------
pCat <- se %>%
  filter(key %in% c("mean_degree", "reciprocity", "av_local_clust")) %>%
  group_by(key) %>%
  mutate(rank = rank(blup)) %>%
  ungroup() %>%
  ggplot(aes(x = rank, y = blup)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(size = 1.4, colour = "#1f4e79") +
  facet_wrap(~statistic, scales = "free", nrow = 1) +
  labs(
    x = "Community (ranked)",
    y = "Shrunken community effect",
    title = "Community-level deviations after removing systematic layer differences"
  ) +
  theme_classic()
save_fig(pCat, "A_fig3_community_effects_caterpillar", width = 10, height = 3.4)

message("Done.")
