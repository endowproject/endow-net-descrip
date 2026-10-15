## ============================================================================
## 01_name_generator_prompts.R
## ----------------------------------------------------------------------------
## Name generator prompts.
##
## Pulls the wording used for each name generator, by community, from each
## site's private nodelist-metadata file, and saves one Latex table per name
## generator.
##
## NOTE: this requires access to each site's private GitHub repository! Those
## without access will not be able to execute this script.
##
## Outputs:
##   Tables/DescriptiveStats/name-generator-prompts/*.tex
##
## Depends on: 00_setup.R
## ============================================================================

source("net_descriptives/00_setup.R")

# ── Name Generator Prompt Tables ────────────────────────────────────────────────

# NOTE: This requires access to each site's private GitHub repository!
# Those without access will not be able to execute.

site_codes_for_nl <- setdiff(site_codes, c("AZ", "TE", "AP"))
site_codes_for_nl <- c(site_codes_for_nl, "TN")

NODELIST_PATTERN <- regex("nodelist[-_]?metadata", ignore_case = TRUE)

# ── Helper: find and load nodelist-metadata for one site ──────────────────────

find_nodelist_metadata <- function(site) {

  meta_dir <- paste0("../endow-", site, "/primary-sources/wave-1/metadata")

  if (!dir.exists(meta_dir))
    stop("Metadata folder not found for site '", site, "': ", meta_dir)

  csv_files <- list.files(meta_dir, pattern = "\\.csv$",
                          full.names = TRUE, ignore.case = TRUE)
  hits <- csv_files[stringr::str_detect(basename(csv_files), NODELIST_PATTERN)]

  if (length(hits) == 0)
    stop("No nodelist-metadata file found for site '", site, "' in: ", meta_dir)

  if (length(hits) > 1) {
    chosen <- hits[which.min(nchar(basename(hits)))]
    warning("Multiple matches for '", site, "'; using: ", chosen)
    hits <- chosen
  }

  first_line <- readLines(hits, n = 1, encoding = "latin1", warn = FALSE)
  delim <- if (lengths(regmatches(first_line, gregexpr(";", first_line))) >
               lengths(regmatches(first_line, gregexpr(",", first_line)))) ";" else ","

  readr::read_delim(hits, delim = delim, show_col_types = FALSE,
                    locale = readr::locale(encoding = "latin1"))
}

recode_map <- c(
  "individual-level money-lending question (ingoing)" = "standardized money-lending question",
  "Requesting loans" = "standardized money-lending question",
  "individual-level money-lending question (outgoing)" = "double sampled money-lending question",
  "Giving loans" = "double sampled money-lending question",
  "individual-level resource exchange question (ingoing)" = "household-level exchange question",
  "Borrowing" = "household-level exchange question",
  "individual-level resource exchange question (outgoing)" = "double sampled household-level exchange question",
  "Lending" = "double sampled household-level exchange question",
  "Female-oriented K-network (informational)" = "female-oriented K network",
  "Male-oriented K-network (informational)" = "male-oriented K network",
  "Female-oriented Q-network (behavioral)" = "female-oriented Q network",
  "Male-oriented Q-network (behavioral)" = "male-oriented Q network",
  "Other important external ties" = "other important external ties",
  "External organization ties (government and NGO)" = "external organizational ties"
)

# Messy list of symbols to try to remove/replace
symbol_recode_map <- c(
  "\u2018" = "'",   # left single quotation mark
  "\u0092" = "'",
  "\u0060" = "'",
  "\u2019" = "'",   # right single quotation mark
  "\u201C" = "",   # left double quotation mark
  "\u201D" = "",   # right double quotation mark
  "\u0094" = "",   # ”
  "�" = "",         # replacement character
  "\u2013" = "-",   # en dash
  "\u2014" = "-",   # em dash
  "\u00e2" = "",    # â (adjust replacement as needed)
  "\u00ca" = "",     # Ê
  "\u00d5" = "",     # Õ
  "\u00ed" = "",     # í
  "\u00ee" = "",     # î
  "\u00e2" = "",     # â
  "\u00e3" = "",     # ã
  "\u00f2" = "",     # ò
  "\u00f3" = "",     # ó
  "\u00f4" = "",     # ô
  "\u00f5" = "",     # õ
  "\u00f6" = "",     # ö
  "\u00ea" = "",     # ê
  "â" = "",
  "â" = "",
  "Ò" = "",
  "Ó" = "",
  "Õ" = "",
  "Ê" = "",
  "î" = "",
  "í" = "'",
  "’" = "'",
  "1.	" = "",
  "2.	" = "",
  "3.	" = "",
  "4.	" = "",
  "4a.	" = "",
  "5.	" = "",
  "6.	" = "",
  "7.	" = "",
  "8.	" = "",
  "9.	" = "",
  "10.	" = ""
)

# Numbered list prefixes handled separately with fixed matching
list_prefixes <- paste0(c("1\\. ", "2\\. ", "3\\. ", "4a?\\. ",
                           paste0(5:10, "\\. ")), collapse = "|")

nodelist_metadata <- purrr::map(
  purrr::set_names(site_codes_for_nl),
  find_nodelist_metadata
)

nodelist_metadata <- purrr::map(nodelist_metadata, \(su_df)
  dplyr::mutate(su_df, standard_variable_description = dplyr::recode(
    standard_variable_description, !!!recode_map
  ))
)

# Duplicate TN into TE/AZ, and AP from MC
nodelist_metadata[["TE"]] <- nodelist_metadata[["TN"]]
nodelist_metadata[["AZ"]] <- nodelist_metadata[["TN"]]
nodelist_metadata[["TN"]] <- NULL
if("MC" %in% site_codes) {
  if("AP" %in% site_codes) {
    nodelist_metadata[["AP"]] <- nodelist_metadata[["MC"]]
  }
}
# align one description
nodelist_metadata[["HE"]]$standard_variable_description[9] <- "external organizational ties"

# clarify MG prompts
nodelist_metadata[["MG"]] <- nodelist_metadata[["MG"]] |>
  dplyr::mutate(
    prompt_phrasing_english = stringr::str_replace_all(
      prompt_phrasing_english,
      stringr::fixed("your HUSBAND"),
      "you (if man) / your husband (if married woman, otherwise not asked)"
    )
  ) |>
  dplyr::mutate(
    prompt_phrasing_english = stringr::str_replace_all(
      prompt_phrasing_english,
      stringr::fixed("does you"),
      "do you"
    )
  ) |>
  dplyr::mutate(
    prompt_phrasing_english = dplyr::if_else(
      as.character(standard_variable_name) %in%
        c("5", "7") &
        !is.na(prompt_phrasing_english),
      paste(
        prompt_phrasing_english,
        "[Asked only of women respondents; men respondents not asked about their wives.]"
      ),
      prompt_phrasing_english
    )
  )

# ── Helper: resolve prompt column ──────────────────────────────────────────────
get_prompt_col <- function(df, site) {
  cols <- names(df)
  if ("english_variable_phrasing" %in% cols) return("english_variable_phrasing")
  if ("prompt_phrasing_english"   %in% cols) return("prompt_phrasing_english")
  stop("No prompt column found for site '", site,
       "': expected 'english_variable_phrasing' or 'prompt_phrasing_english'")
}

main_name_gens <- tibble::tribble(
  ~ng_name, ~caption,
  "standardized money-lending question", "Name generator prompts for borrowing money",
  "double sampled money-lending question", "Name generator prompts for lending money",
  "household-level exchange question", "Name generator prompts for borrowing goods",
  "double sampled household-level exchange question", "Name generator prompts for lending goods",
  "female-oriented Q network", "Name generator prompts for female-oriented labor/behavioral support",
  "male-oriented Q network", "Name generator prompts for male-oriented labor/behavioral support",
  "female-oriented K network", "Name generator prompts for female-oriented informational/affective support",
  "male-oriented K network", "Name generator prompts for male-oriented informational/affective support",
  "other important external ties", "Name generator prompts for ties to government or NGO employees",
  "external organizational ties", "Name generator prompts for support outside the community"
)

name_gens <- purrr::pmap(main_name_gens, \(ng_name, caption) {
  purrr::imap(nodelist_metadata, \(df, site) {
    prompt_col <- get_prompt_col(df, site)
    row <- dplyr::filter(df, standard_variable_description == ng_name)
      if (nrow(row) == 0)
        return(tibble::tibble(site = site, prompt = ""))
      if (nrow(row) > 1)
        warning("Multiple rows with standard_variable_description == ", ng_name,
                " for site '", site, "'; using first")
      tibble::tibble(site = site, prompt = row[[prompt_col]][1])
    }) |>
      purrr::list_rbind() |>
        dplyr::mutate(prompt = stringr::str_replace_all(prompt, symbol_recode_map)) |>
        dplyr::mutate(prompt = stringr::str_remove_all(prompt, list_prefixes)) |>
        dplyr::mutate(prompt = stringr::str_remove_all(prompt, '["\u201C\u201D]')) |>
        dplyr::mutate(prompt = paste0(toupper(stringr::str_sub(prompt, 1, 1)),
                stringr::str_sub(prompt, 2))) |>
        dplyr::arrange(site)
  }
) |>
  purrr::set_names(main_name_gens$ng_name)

## Deal with main site name generator exceptions
#  "IH", "PC", "WH", "WL" sites are missing 2, 9, and 10
#  "MK" is missing 9 and 10, and has collapsed 5/6 and 7/8 questions
# "AV", "DJ" didn't ask 9 or 10
# "BY" prompt for 9 (govt/NGO) not actually on govt/NGO

name_gens <- purrr::imap(name_gens, \(df, i) {
  df <- dplyr::mutate(df, prompt = dplyr::recode(prompt, `>> not found <<` = "(not asked)"))
  df <- dplyr::mutate(df, prompt = dplyr::recode(prompt, `NANA` = ""))
# also recode NA to ""
  df <- dplyr::mutate(df, prompt = dplyr::if_else(is.na(prompt), "", prompt))

  if (i == "double sampled money-lending question") {
    dplyr::mutate(df, prompt = dplyr::if_else(
      site %in% c("IH", "PC", "WH", "WL"),
      "(not asked)",
      prompt
    ))
  } else if (i == "external organizational ties") {
    dplyr::mutate(df, prompt = dplyr::if_else(
      site %in% c("AV", "DJ", "IH", "PC", "WH", "WL", "BY", "MK"),
      "(not asked)",
      prompt
    ))
  } else if (i == "other important external ties") {
    dplyr::mutate(df, prompt = dplyr::if_else(
      site %in% c("AV", "DJ", "IH", "PC", "WH", "WL", "MK"),
      "(not asked)",
      prompt
    ))
  } else if (i %in% c("female-oriented K network", "male-oriented K network")) {
    dplyr::mutate(df, prompt = dplyr::if_else(
      site == "MK",
      "(One joint question asked) If you need advice, who would you go to?",
      prompt
    ))
  } else if (i %in% c("female-oriented Q network", "male-oriented Q network")) {
    dplyr::mutate(df, prompt = dplyr::if_else(
      site == "MK",
      "(One joint question asked) If there was an agricultural wage labor opportunity that required a team of workers, who would you invite to join you?",
      prompt
    ))
  } else {
    df
  }
})

# ── Save results ───────────────────────────────────────────────────────────────
ng_table_path <- file.path(EndowDropbox, "Tables", "DescriptiveStats", "name-generator-prompts")
ng_csv_path <- file.path(csvpath, "name-generator-prompts")

# ensure folders exist
dir.create(ng_table_path, showWarnings = FALSE)

for (i in main_name_gens$ng_name) {
  base_name <- gsub(" ", "-", i)
  filename <- paste0(ng_table_path, "/", base_name, ".tex")
  csv_filename <- paste0(ng_csv_path, "/", base_name, ".csv")
  message("Saving: ", base_name)
  kableExtra::kable(
    name_gens[[i]],
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    caption = main_name_gens$caption[main_name_gens$ng_name == i],
    label = paste0("prompts-", gsub(" ", "-", i)),
    col.names = c("Community", "Prompt"),
    linesep = ""
    ) |>
    kableExtra::kable_styling(
      latex_options = "repeat_header",
      font_size = 8
    ) |>
    kableExtra::column_spec(2, width = "15cm") |>
    kableExtra::save_kable(filename)

  # CSV version
  csv_data <- name_gens[[i]]
  colnames(csv_data) <- c("Community", "Prompt")
  readr::write_csv(csv_data, csv_filename)
}
