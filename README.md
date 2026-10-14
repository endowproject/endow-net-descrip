# ENDOW Network Descriptives — Replication Kit

This repository contains the analysis code for:

> **"A cross-cultural dataset on the structure of social networks: The ENDOW project"**
> The ENDOW Team

The paper presents descriptive statistics on nomination patterns across a series of 10 social support network prompts and network summary statistics for sharing unit-level networks derived from those network prompts.

This repo provides code to generate all of the figures and tables, using faux data. With the use of faux data, the resulting tables/figures will not reproduce the actual figures/tables in the paper but instead fully show the code.

The repo also includes the summary tables included in the paper as a series of CSVs in the ActualDataOutputs folder.

Installation of the relevant R packages should be straightforward and take a matter of minutes.

With the faux data, the full script should take XXX MINUTES to run.

---

## Table of Contents

1. [Repository Layout](#repository-layout)
2. [Prerequisites](#prerequisites)
3. [R Setup](#r-setup)
4. [Running the Replication Kit](#running-the-replication-kit)
5. [Pipeline Overview](#pipeline-overview)
6. [Key Data Files](#key-data-files)
7. [Output Locations](#output-locations)
8. [R Package Dependencies](#r-package-dependencies)

---

## Repository Layout

```
ENDOW-net-descrip/
├── Scripts/
│   ├── 00_setup.R
│   ├── 01_name_generator_prompts.R
│   ├── 02_network_summary_stats.R
│   ├── 03_respondent_coverage.R
│   ├── 04_nomination_counts_composition.R
│   ├── 05_concordance_within_su.R
│   ├── 06_concordance_req_give.R
│   ├── 07_degree_distributions.R
│   ├── 08_network_visualizations.R
│   ├── 09_network_overlap.R
│   ├── 10_pca.R
│   ├── 11_kinship.R
│   └── 12_variance_decomp.R
│
├── DerivedData/
│   ├── ... interim files created by scripts
│   ├──
│   ├──
│   └──
│
├── FauxData/
│   ├── ... things coming from rep repo that are called into 00_setup
│   ├──
│   ├──
│   └──
│
├── ActualDataOutputs/
│   ├── ... actual data CSVs
│   ├──
│   ├──
│   └──

```

---

## Prerequisites

You need three things set up before running any code:

ELLY: UPDATE THIS TO PULL FROM FauxData FOLDER where I can dump the derived data files from main rep-code repo

**1. GitHub repositories (local clones)**
- `ENDOW-wealth-corr` — this repository

**2. ENDOW database**
- `endow-database` — the faux ENDOW database files; not actually needed here!

**3. A Dropbox**
- A local folder — all figures and tables are written here

**4. R 4.x** — see [R Package Dependencies](#r-package-dependencies) below

---

## R Setup

The scripts find files using path variables that you define in your `.Rprofile`. See [this guide](https://stackoverflow.com/questions/46819684/how-to-access-and-edit-rprofile) if you are not familiar with `.Rprofile`.

Add the following to your `~/.Rprofile`, adjusting paths to match your machine, _or_ set them at the start of any session before executing 00_run_all.R:

**Mac:**
```r
EndowGitHub <- file.path("~", "Documents", "GitHub", "endow-net-descrip")
EndowDatabase <- file.path("~", "Documents", "GitHub", "endow-net-descrip", "FauxData")
EndowDropbox <- file.path("~", "Dropbox", "ENDOW")
```

**Windows:**
```r
EndowGitHub <- file.path("C:/Users/tomru", "Dropbox (Personal)", "GitHub", "endow-net-descrip")
EndowDatabase <- file.path("C:/Users/tomru", "OneDrive - Stanford", "Documents", "GitHub", "endow-net-descrip", "FauxData")
EndowDropbox <- file.path("C:/Users/tomru", "Dropbox (Personal)", "ENDOW")
```

After editing `.Rprofile`, restart R to pick up the changes. The replication kit will stop with a clear error message if any required variable is missing.

---

## Running the Replication Kit

Open a **fresh R session**, set your working directory to the repo, and run any file directly — its `source()` calls will pull in what it needs:

```r
source("Scripts/11_pca.R")
```

Or run the whole thing:

```r
for (f in file.path("Scripts", sprintf("%02d_*.R", 1:12))) { for (m in Sys.glob(f)) source(m) }
```

---

## Pipeline Overview

The code to generate the faux input data used here is on another public repository:
  - https://github.com/endowproject/endow-wealth-corr-rep-code/blob/main/simulate_endow_database.R generates simulated data for four communities
  - https://github.com/endowproject/endow-wealth-corr-rep-code/blob/main/process_data.R cleans the input data and generates basic networks
  - https://github.com/endowproject/endow-wealth-corr-rep-code/blob/main/analyse_data_and_create_networks.R combines layers and computes descriptives
These generate the objects included in the FauxData folder.

### 00_setup.R


---

## Output Locations

| Output type | Location |
|-------------|----------|
| All figures | `EndowDropbox/Figures/` (including exploratory subfolders) |
| All tables | `EndowDropbox/Tables/` |
| Intermediate datasets | `EndowGitHub/DerivedData/` |

---

## R Package Dependencies

Install all required packages before running:

```r
# From CRAN
install.packages(c(
  "data.table", "dplyr", "tidyr", "ggplot2", "ggpubr", "ggrepel",
  "gridExtra", "ggstance", "viridis", "reshape2", "purrr", "tibble",
  "igraph", "DescTools", "ineq",
  "stargazer", "plm", "glmnet", "lme4",
  "stringr", "janitor", "readr", "egg",
  "wbwdi", "pipr", "kableExtra",
  "mice", "splitstackshape", "tidyverse",
  "fastDummies", "readxl", "latex2exp",
  "knitr", "ggridges", "magick", "ggforce",
  "sf", "rnaturalearth", "rnaturalearthdata",
  "cNORM", "acid", "randomForest", "reshape",
  "colorRamps", "ggcorrplot", "ppcor", "tidytext",
  "showtext", "sysfonts", "patchwork", "metafor"
))

```

The analysis has been tested on **R 4.6.1** on macOS.
