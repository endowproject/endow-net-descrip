## ============================================================================
## 06_concordance_req_give.R
## ----------------------------------------------------------------------------
## Request/give double-sampling concordance.
##
## Compares the paired "request" and (transposed) "give" networks for loans
## and household exchange, to see how often a tie reported by the requester
## is corroborated by the giver.
##
## Outputs:
##   Tables/DescriptiveStats/double_sampled_concordance.tex
##   Figures/DescriptiveStats/double_sampled_concordance.png
##
## Depends on: 00_setup.R
## ============================================================================

source("net_descriptives/00_setup.R")

# ── Concordance between req and t(giv) networks ────────────────────────────────────────────────

# incoming = *_req, outgoing = t(*_giv)
# A tie i->j in req means i says they requested from j.
# A tie i->j in giv means i says they gave to j, i.e. j->i in the req direction.
# So we compare req[i,j] against t(giv)[i,j] = giv[j,i].

compute_concordance <- function(req_net, giv_net, site, pair) {
  adj_req   <- as.matrix(igraph::as_adjacency_matrix(req_net))
  adj_giv_t <- t(as.matrix(igraph::as_adjacency_matrix(giv_net)))

  # Align vertex order to req_net (it surely does, but still)
  vnames    <- rownames(adj_req)
  adj_giv_t <- adj_giv_t[vnames, vnames]

  n       <- nrow(adj_req)
  summats <- adj_req + adj_giv_t

  present <- sum(summats == 2)                    # tie in both
  null    <- sum(summats == 0) - n                # absent in both, minus diagonal
  net1    <- sum(adj_req)    - present            # only in req
  net2    <- sum(adj_giv_t)  - present            # only in t(giv)
  total   <- present + net1 + net2

  tibble::tibble(
    site              = site,
    pair              = pair,
    present           = present,
    null              = null,
    net1              = net1,          # edges only in *_req
    net2              = net2,          # edges only in t(*_giv)
    concord_prop      = present / total,
    net1_concord_prop = present / (present + net1),
    net2_concord_prop = present / (present + net2),
    net1_prop         = net1 / total,
    net2_prop         = net2 / total
  )
}
# note: concord_prop + net1_prop + net2_prop = 1

concordance_pairs <- list(
  loan     = c(req = "loan_req",     giv = "loan_giv"),
  exchange = c(req = "exchange_req", giv = "exchange_giv")
)

concordance_results <- purrr::imap_dfr(su_nets_expanded, function(site_nets, site) {
  purrr::imap_dfr(concordance_pairs, function(pair, pair_name) {
    req_net <- site_nets[[pair["req"]]]
    giv_net <- site_nets[[pair["giv"]]]

    if (is.null(req_net) || is.null(giv_net)) {
      message(glue::glue("Skipping {site} / {pair_name}: missing network(s)"))
      return(NULL)
    }

    compute_concordance(req_net, giv_net, site = site, pair = pair_name)
  })
})

ds_concordance_table <- concordance_results |>
  dplyr::select(site, pair, concord_prop, net1_concord_prop, net2_concord_prop) |>
  tidyr::pivot_wider(
    names_from  = pair,
    values_from = c(concord_prop, net1_concord_prop, net2_concord_prop)
  ) |>
  dplyr::select(site, dplyr::ends_with("_loan"), dplyr::ends_with("_exchange"))

withr::with_options(
  list(knitr.kable.NA = "---"),
  ds_concordance_table |>
    dplyr::arrange(site) |>
    kableExtra::kbl(
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      digits = 3,
      col.names = c("Community", rep(c("Overall", "Request", "Give"), 2)),
      caption = "Proportion of concordant ties in the sharing-unit networks for each double-sampled relationship, in each community. ``Request'' refers to the network based on the name generator in which respondents reporting receiving a type of support, while ``Give'' refers to respondents reporting giving support to others. The column ``Overall'' reports the proportion of concordant ties in the union of the request and (transposed) give networks.",
      label = "ds_concordance",
    ) |>
    kableExtra::add_header_above(c(" " = 1, "Loan" = 3, "Exchange" = 3)) |>
    kableExtra::kable_styling(latex_options = c("hold_position")) |>
    kableExtra::save_kable(file.path(
      EndowDropbox,
      "Tables",
      "DescriptiveStats",
      "double_sampled_concordance.tex"
    ))
)

# CSV:
csv_names <- c(
  "Community",
  paste(rep(c("Loan", "Exchange"), each = 3), rep(c("Overall", "Request", "Give"), 2))
)

ds_concordance_table |>
  dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 3))) |>
  dplyr::arrange(site) |>
  setNames(csv_names) |>
  readr::write_csv(
    file.path(csvpath, "double_sampled_concordance.csv")
  )

p <- concordance_results |>
  dplyr::mutate(pair = forcats::fct_relevel(pair, "loan", "exchange")) |>
  ggplot(aes(y = concord_prop, x = pair, colour = pair, fill = pair)) +
  geom_violin(
    alpha = 0.4,
    linewidth = 0.5,
    scale = "width",
    draw_quantiles = 0.5
  ) +
  geom_label(
    aes(label = site),
    family = "roboto",
    position = ggplot2::position_jitter(width = 0.2),
    size = 1.5,
    alpha = 0.5,
    linewidth = 0.25,
    colour = "white",
    label.padding = unit(0.3, "lines"),
    label.r = unit(0.6, "lines")
  ) +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(
    y = "Proportion Concordant",
    x = "Name Generator"
  ) +
  theme_classic() +
  theme(legend.position = "none")
ggsave(
  file.path(EndowDropbox, "Figures", "DescriptiveStats", "double_sampled_concordance.png"),
  plot = p,
  width  = 4,
  height = 4
  )
