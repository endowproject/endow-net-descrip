# ENDOW Network Descriptives — Replication Kit & Data Outputs

This repository contains the analysis code and resulting data outputs for:

> **"A cross-cultural dataset on the structure of social networks: The ENDOW project"**
> The ENDOW Team

The paper presents descriptive statistics on nomination patterns across a series of 10 social support network prompts and network summary statistics for sharing unit-level networks derived from those prompts.

This repo provides code to generate all of the figures and tables, using **faux data** for four simulated communities (EG, EX, SI, SM). Because the data are simulated, the resulting tables and figures will not reproduce those in the paper; they demonstrate that the code runs and show what each output looks like.

In addition to the code, the summary tables underlying the paper (i.e., the key data files we are able to share) are included as CSVs in the `ActualDataOutputs/` folder. Running the code on the faux data writes to `FauxDataOutputs/` and does **not** overwrite them.

Installation of the relevant R packages should be straightforward and take a matter of minutes. With the faux data, the full pipeline should take just a few minutes to run.

---

## Table of Contents

1. [Repository Layout](#repository-layout)
2. [Actual Data Outputs](#actual-data)
2. [Prerequisites](#prerequisites)
3. [R Setup](#r-setup)
4. [Running the Replication Kit](#running-the-replication-kit)
5. [Pipeline Overview](#pipeline-overview)
6. [Key Data Files](#key-data-files)
7. [Output Locations](#output-locations)
8. [Known Limitations](#known-limitations)
9. [R Package Dependencies](#r-package-dependencies)

---

## Repository Layout

```
endow-net-descrip/
├── .Rprofile                       # defines the path variables used by all scripts
├── ActualDataOutputs/              # summary-table CSVs from the actual data, as presented in the paper
│   ├── name-generator-prompts/       # 10 files, one per name generator
│   ├── network-summary-stats/        # 13 files, one per network layer + a cross-community summary
│   ├── nomination-summary-stats/     # 12 files, one per name generator + a global summary
│   └── *.csv                         # coverage, concordance, overlap, PCA, and kinship tables
│
├── FauxData/                       # simulated inputs (see Key Data Files)
│   ├── networks/                   #   <SITE>_networks.csv (edgelists)
│   └── *.rdata, *.rds, *.csv       #   objects built from the simulated inputs; intermediate
│                                   #   files created by the scripts are also saved here
│
├── FauxDataOutputs/                # created when you run the scripts
│   ├── Figures/
│   ├── Tables/
│   └── ...                         # CSV summaries and other outputs
│
└── net_descriptives/               # scripts to generate all data outputs
    ├── 00_setup.R
    ├── 01_name_generator_prompts.R # included for completeness, but cannot be run (requires access to private repos)
    ├── 02_network_summary_stats.R
    ├── 03_respondent_coverage.R
    ├── 04_nomination_counts_composition.R
    ├── 05_concordance_within_su.R
    ├── 06_concordance_req_give.R
    ├── 07_degree_distributions.R
    ├── 08_network_visualizations.R
    ├── 09_network_overlap.R
    ├── 10_pca.R
    ├── 11_kinship.R
    └── 12_variance_decomp.R        # Included for completeness, but will throw some errors

```

---

## Actual Data Outputs

`ActualDataOutputs/` contains the summary tables underlying the paper, generated from the full (private) ENDOW data. Each CSV is the machine-readable version of a table or figure in the paper; the script that produces its faux-data counterpart (written to `FauxDataOutputs/` under the same name) is given below.

### Network layers

Many files are organized by network prompt / network layer. The layers are the ten name generators plus two derived layers:

| Layer code | Layer (as labelled in tables) |
|---|---|
| `loan_req` / `loan_giv` | Borrow money / Lend money |
| `exchange_req` / `exchange_giv` | Borrow goods / Lend goods |
| `f_behav` / `m_behav` | Behavioral support (women) / (men) |
| `f_info` / `m_info` | Informational/affective support (women) / (men) |
| `govt` | Government and NGO ties |
| `external` | Support outside the community |
| `primary_kin` | Primary (close) kin |
| `sum` | Aggregate of the name-generator layers |

### Top-level files

| File | Produced by | Description |
|---|---|---|
| `respondents_by_site.csv` | `03_respondent_coverage.R` | Sharing-unit coverage and respondents for the social support network survey, one row per community. Columns: number of sharing units (`# SU`), mean number of members per SU, percent of SUs sampled, mean number surveyed per SU, total surveyed, number of women and men surveyed, and the intended respondents (e.g., female head, male and female heads). |
| `within_su_concordance.csv` | `05_concordance_within_su.R` | Agreement in nominations among co-resident respondents, for sharing units with more than one surveyed member. For each respondent, concordance is the proportion of the individuals they nominated (`Indivs`), or of the sharing units those individuals belong to (`SUs`), who were also nominated by another respondent in the same SU. Columns give the mean, median and SD of that proportion by network prompt (ten name generators plus the aggregate) and community. Only community-by-prompt combinations with at least 10 multi-respondent SUs are included. |
| `double_sampled_concordance.csv` | `06_concordance_req_give.R` | Concordance between "request" and "give" networks for the two double-sampled relationships (monetary loans and household exchange), one row per community. Ties reported by the requester are compared with the transposed give network. `Overall` is the proportion of concordant ties in the union of the two networks; `Request` is the proportion of request ties also reported by the giver; `Give` is the proportion of give ties also reported by the requester. |
| `edge_overlap.csv` | `09_network_overlap.R` | Edge-level overlap between each pair of the ten name-generator layers, by community and pooled across communities (`all`). Gives the Pearson correlation between the layers' dyad-level edge values, using edge weights (`Pearsons r Weighted`) and binarized ties (`Pearsons r Binary`). |
| `alter_overlap_personwise.csv` | `09_network_overlap.R` | Alter-level overlap between each pair of layers, computed **per respondent** and then averaged across respondents who answered both prompts. `Jaccard Indivs` is the Jaccard index of the sets of individuals nominated; `Jaccard SUs` is the same for the sharing units nominated. Community `all` pools respondents across communities. |
| `alter_overlap_alterwise.csv` | `09_network_overlap.R` | Same measures as the personwise file, but computed **across alters**: intersection and union counts are summed over respondents before taking the ratio. Community `all` sums the counts across communities. |
| `pca_loadings_weighted.csv` | `10_pca.R` | Principal component loadings from a PCA of the eight main network layers (borrow/lend money, borrow/lend goods, behavioral support for women/men, informational/affective support for women/men), using **edge weights**. The data are stacked dyads and variables are standardized. One block per community and a pooled block (`all`); each block lists each layer's loading on PC1–PC8 and a final `Variance explained` row. Layers with no variance in a community are dropped. |
| `pca_loadings_binary.csv` | `10_pca.R` | As above, but using **binarized** ties (presence/absence of a tie). |
| `kin_support_summary.csv` | `11_kinship.R` | For each community and network layer, the average proportion of a sharing unit's support partners who are close kin, and the average proportion of its close kin who are support partners. Averages are across sharing units. |
| `sitewise_kin_support.csv` | `11_kinship.R` | Availability of close kin and ties to them in the aggregate network, pooled across communities (`Pooled`, first row) and by community. Columns: `N` sharing units; mean and median number of co-residing close-kin sharing units (`Co-Res. Kin`); percent of SUs with no co-residing close kin (`No Kin`); mean proportion of aggregate-network partners who are close kin (`Partners Kin`); mean proportion of co-resident kin SUs named as partners (`Named Kin`); and the equivalent proportion for non-kin SUs (`Named Non-Kin`). |

### `name-generator-prompts/`

*Produced by `01_name_generator_prompts.R` (not runnable on the faux data).* Each file has two columns, `Community` and `Prompt`, giving the English wording of the name-generator prompt used in each community. An empty prompt means the question was not asked there.

| File | Name generator |
|---|---|
| `standardized-money-lending-question.csv` | Borrowing money |
| `double-sampled-money-lending-question.csv` | Lending money |
| `household-level-exchange-question.csv` | Borrowing goods |
| `double-sampled-household-level-exchange-question.csv` | Lending goods |
| `female-oriented-Q-network.csv` / `male-oriented-Q-network.csv` | Female- / male-oriented labor (behavioral) support |
| `female-oriented-K-network.csv` / `male-oriented-K-network.csv` | Female- / male-oriented informational/affective support |
| `other-important-external-ties.csv` | Ties to government or NGO employees |
| `external-organizational-ties.csv` | Support outside the community |

### `network-summary-stats/`

*Produced by `02_network_summary_stats.R`.*

**`net-summary-stat-<layer>.csv`** (10 files, one per layer for the main network layers (excluding the two focused on external ties), plus the primary kin network and the aggregate network) contain network summary statistics for the sharing-unit network of that layer, one row per community in which the layer was collected. Columns: `Nodes`, `Edges`, `Density`, `Mean Deg.` (mean degree), `Mean W. Deg.` (mean weighted degree), `In-Deg.` and `Out-Deg.` (ranges), `Deg. Centr.` and `In-Deg. Centr.` (degree and in-degree centralization), `Deg. Assort.` (degree assortativity), `APL` (average path length), `Diam.` (diameter), `N Isol.` and `Prop. Isol.` (number and proportion of isolates), `N Comp.` (number of weakly-connected components), `Largest Comp.` (proportion of nodes in the largest component), `Global Clust.` (global clustering coefficient), `Av. Loc. Clust.` (average local clustering coefficient) and `Recip.` (reciprocity). The `primary_kin` network is undirected, so its file omits the weighted, in-/out-degree, in-degree-centralization and reciprocity columns.

**`global-net-summary-median-meansd.csv`** summarizes those statistics across communities for each layer: for every statistic, one row gives the median across communities and the next gives the mean and standard deviation. Columns are the eight main network layers and the aggregate network.

### `nomination-summary-stats/`

*Produced by `04_nomination_counts_composition.R`.* Entries of the form `mean (SD)` summarize nominations made by respondents.

**`global-nominations-summary.csv`**: one row per name generator, across all respondents in all communities.

**`nominations-summary-<layer>.csv`** (11 files: the ten name generators and `sum`, the aggregate): the same measures, one row per community. These include the additional column `N Resp.`, the number of respondents.

Columns: `Prop. Zero` (proportion of respondents nominating no alters), `N Indiv.` (number of individuals nominated), `N SU` (number of sharing units nominated), `Prop. Same SU` (proportion of nominated individuals who are in the respondent's own sharing unit) and `Prop. Ext.` (proportion of nominated individuals external to the community).

---

## Prerequisites for Executing Scripts

1. **R 4.x** (tested on R 4.6.1 on macOS) and the packages listed under [R Package Dependencies](#r-package-dependencies).
2. **A clone of this repository.** Everything the scripts need is in the repo (with the exception of `01_name_generator_prompts.R`, which cannot run).
3. **An internet connection on first run**, because `00_setup.R` downloads the *Roboto Mono* font from Google Fonts via `sysfonts::font_add_google()`.

---

## R Setup

The scripts find files using path variables (`EndowGitHub`, `EndowDatabase`, `EndowDerivedData`, `EndowDropbox`). This repository includes a project-level `.Rprofile` in its root that defines them relative to the repo, so no manual configuration is needed:

```r
EndowGitHub      <- normalizePath(getwd())                    # repo root
EndowDatabase    <- file.path(EndowGitHub, "FauxData")        # networks/*.csv
EndowDerivedData <- file.path(EndowGitHub, "FauxData")        # .rdata / .rds objects
EndowDropbox     <- file.path(EndowGitHub, "FauxDataOutputs") # Figures/ and Tables/ land here
```

R runs this file automatically at startup, provided R is **launched from the repository root**. The easiest way is to open the repo as an RStudio project. If you start R another way, set the working directory to the repo root and restart R before running `00_setup.R`. See [this guide](https://stackoverflow.com/questions/46819684/how-to-access-and-edit-rprofile) if you are not familiar with `.Rprofile`.

A few things to be aware of:

- A project-level `.Rprofile` is used **instead of** your personal `~/.Rprofile` when R starts in this folder. If you rely on settings in your personal file, add `source("~/.Rprofile")` at the top of the project `.Rprofile`.
- Scripts are run from the repo root and call `source("net_descriptives/00_setup.R")`, so the scripts folder must keep the name `net_descriptives`.
- If you would rather define the paths yourself, edit the `.Rprofile` in the repo root, or set the four variables at the start of a session before running `00_setup.R`.

The scripts will stop with a clear error message if any required variable is missing.

---

## Running the Replication Kit

Open a **fresh R session** from the repo root and run any script directly; its `source()` call will pull in `00_setup.R`:

```r
source("net_descriptives/10_pca.R")
```

Several scripts depend on intermediate files written by earlier ones (see [Pipeline Overview](#pipeline-overview)), so if you run scripts individually, run them in numerical order the first time.

To run the whole pipeline, skipping the two scripts that are not expected to run on the faux data (see [Known Limitations](#known-limitations)):

```r
scripts <- Sys.glob(file.path("net_descriptives", "[0-1][0-9]_*.R"))
scripts <- scripts[!grepl("00_setup|01_name_generator|12_variance", scripts)]
for (f in scripts) source(f)
```

---

## Pipeline Overview

The code that generated the faux input data is in another public repository:

- [`simulate_endow_database.R`](https://github.com/endowproject/endow-wealth-corr-rep-code/blob/main/simulate_endow_database.R) generates simulated data for four communities (including what is contained in networks/)
- [`process_data.R`](https://github.com/endowproject/endow-wealth-corr-rep-code/blob/main/process_data.R) cleans the input data and generates basic networks
- [`analyse_data_and_create_networks.R`](https://github.com/endowproject/endow-wealth-corr-rep-code/blob/main/analyse_data_and_create_networks.R) combines layers and computes descriptives
- [`analyse_data_and_create_networks.R`](https://github.com/endowproject/endow-wealth-corr-rep-code/blob/main/analyse_data_and_create_networks.R) combines layers and computes descriptives
- [`kinship.R`](https://github.com/endowproject/endow-wealth-corr-rep-code/blob/main/kinship.R) creates kinship-support variables

These generate the objects included in `FauxData/`. (Two figures and one table in this paper are also generated in code in the [`descriptivestats.R`](https://github.com/endowproject/endow-wealth-corr-rep-code/blob/main/descriptivestats.R) script in that repository: the overall world map with community locations (`sites_map.png`), the table showing primary and secondary subsistence type (`subsistence.png`), and the table giving details on each ENDOW community including participants, data lead(s), etc. (`su_net_summary.tex`).)

The scripts in `net_descriptives/` then produce the descriptive statistics, tables and figures of this paper:

| Script | What it does | Reads from earlier scripts | Writes for later scripts |
|---|------|---|---|
| `00_setup.R` | Loads data, applies site-specific cleaning, builds shared lookups and plotting functions. Sourced by every other script. | | |
| `01_name_generator_prompts.R` | Tables of the prompt wording used for each name generator, by community. *Not runnable here; see Known Limitations.* | | |
| `02_network_summary_stats.R` | Network-level summary statistics (density, degree, centralization, clustering, etc.) by layer and community; cross-community distributions. | | `network_summaries_long.rdata` |
| `03_respondent_coverage.R` | Sharing-unit and respondent coverage by community; degree and clustering vs. sampling. | `network_summaries_long.rdata` | `su_sampling.rdata` |
| `04_nomination_counts_composition.R` | Nomination counts and their composition (internal vs. external, gender, age) by name generator and community. | | `external_alters.rdata`, `alters_summary.rdata` |
| `05_concordance_within_su.R` | Concordance in nominations among members of the same sharing unit. | `alters_summary.rdata` | |
| `06_concordance_req_give.R` | Concordance between paired request and give (double-sampled) networks. | | |
| `07_degree_distributions.R` | Degree distributions (CCDFs) by layer. | | |
| `08_network_visualizations.R` | Network visualizations by layer and community. | | |
| `09_network_overlap.R` | Edge- and alter-level overlap across layers. | `alters_summary.rdata` | |
| `10_pca.R` | Principal component analysis of network layers, by community and pooled. | | |
| `11_kinship.R` | Kin composition of support networks. | | |
| `12_variance_decomp.R` | Between-community variance decomposition (multilevel models). *Will run only with errors; see Known Limitations.* | `network_summaries_long.rdata`, `su_sampling.rdata`, `external_alters.rdata` | |

---

## Key Data Files

All inputs are in `FauxData/`:

| File(s) | Used by |
|---|---|
| `networks/<SITE>_networks.csv` | `00_setup.R`|
| `people_observations.rdata`, `su_alters.rdata`, `su_meta.rdata` | `00_setup.R` |
| `su_nets.rdata`, `su_nets_expanded.rdata`, `su_nets_expanded_withloops.rdata` | `00_setup.R` and downstream scripts |
| `su_master.rds` | `03_respondent_coverage.R` |
| `su_kin_support.csv`, `DataForLASSO.csv` | `11_kinship.R` |

Two-letter codes in the `<SITE>` filenames identify the four simulated communities (EG, EX, SI, SM).

---

## Output Locations

| Output type | Location |
|-------------|----------|
| Figures | `FauxDataOutputs/Figures/` (`DescriptiveStats/` and `Networks/`) |
| Tables (LaTeX) | `FauxDataOutputs/Tables/` |
| CSV summaries of tables | `FauxDataOutputs/` |
| Intermediate datasets | `FauxData/` |
| Tables from the actual data (paper) | `ActualDataOutputs/` |

The output folders are created automatically by `00_setup.R`.

---

## Known Limitations

Because this repository uses simulated data for four communities rather than the full private dataset:

- **`01_name_generator_prompts.R` will not run.** It reads each community's name-generator metadata from private per-site repositories. It is included for completeness, to document how the prompt tables in the paper were produced.
- **`12_variance_decomp.R` will throw errors.** It fits multilevel models with a random intercept for community, which is not meaningful with four communities. It is included for completeness.
- Some scripts contain site-specific handling for communities not in the faux data (for example, special cases for particular sites in `00_setup.R` and `04_nomination_counts_composition.R`). These are skipped or have no effect with the faux data.

---

## R Package Dependencies

Install all required packages before running:

```r
# From CRAN
install.packages(c(
  # core
  "tidyverse", "igraph", "data.table", "patchwork", "scales",
  "kableExtra", "ggrepel", "ggh4x", "RColorBrewer", "viridis",
  "broom", "glue", "withr",
  # network plotting (08_network_visualizations.R)
  "tidygraph", "ggraph",
  # fonts (00_setup.R)
  "showtext", "sysfonts",
  # multilevel models (12_variance_decomp.R)
  "lme4"
))
```

The analysis has been tested on **R 4.6.1** on macOS.