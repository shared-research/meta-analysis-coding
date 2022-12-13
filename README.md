
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The Cognitive Effects of Computational Thinking: A Systematic Review and Meta-Analytic Study

<!-- badges: start -->

[<img alt="alt_text" src="https://img.shields.io/badge/OSF-DOI%2010.17605%2FOSF.IO%2FUVBCD-blue" />](https://osf.io/uvbcd/)
<!-- badges: end -->

This repository contains all the scripts, functions and data to fully
reproduce the statistical analysis of the *The cognitive effects of
computational thinking: A systematic review and meta-analytic study*
paper by Montuori C., Gambarota F., Altoè G. and Arfè B.

For a detailed description of the analysis method see the
**supplementary materials**
([HTML](docs/supplementary/supplementary.html),
[PDF](docs/supplementary/supplementary.pdf))

## Setup

In order to reproduce the analysis is necessary to clone or download
this repository and open the `metaCoding.Rproj`. Then the `renv` package
will install all required packages with the correct version. After
installing all packages, the `main_script.R` or individual scripts
(`scripts/*`) can be used to run each analysis step.

## Folders organization

- `data/`: contains raw and cleaned data in `.rds` and `.xlsx`/`.csv`
  format
- `docs/`: contains scripts to create the supplementary materials
- `main_script.R`: is a script for easily managing all analysis steps
- `objects/`: contains R objects in `.rds` format created from the
  analysis scripts
- `R/`: contains all custom functions for the project. Functions are
  automatically loaded when the project is loaded. For reloading use
  `devtools::load_all()`
- `scripts/`: main scripts for pre-processing, statistical analysis and
  creating tables/figures
  - `1_analysis.R`: Scripts for calculating the effect size and
    computing the meta-analysis models
  - `2_tables_figures.R`: Script to create tables and figures for the
    paper
  - `3_supplementary.R`: Script to create supplementary materials
    objects
- `tables/`: contains tables created with the `2_tables_figures.R` in
  `.docx` format
- `tests/`: is for testing the effect size computation functions

## Coding style

The analysis project is organized as an R package where functions within
the `R/` folder are automatically available into the global environment
when the project is activated (`metaCoding.Rproj`). Is possible to
manually load the functions using `devtools::load_all()`. The coding
style is based on the [*tidyverse*](https://www.tidyverse.org/) using
also [*metaprogramming*](https://adv-r.hadley.nz/meta-big-picture.html).
For managing multiple meta-analysis models we used *nested tibbles* as
data structure (see [here](https://r4ds.had.co.nz/many-models.html)). In
particular the `objects/dat_meta.rds` contains all processing steps and
models as different columns.

# Session Info

    #>  setting  value
    #>  version  R version 4.2.2 (2022-10-31 ucrt)
    #>  os       Windows 10 x64 (build 19044)
    #>  system   x86_64, mingw32
    #>  ui       RTerm
    #>  language (EN)
    #>  collate  English_United States.utf8
    #>  ctype    English_United States.utf8
    #>  tz       Europe/Berlin
    #>  date     2022-12-13
    #>  pandoc   2.19.2 @ C:/Program Files/RStudio/bin/quarto/bin/tools/ (via rmarkdown)

# R Packages

- `devtools`
- `testthat`
- `rmarkdown`
- `bookdown`
- `knitr`
- `flextable`
- `here`
- `tidyverse`
- `broom.mixed`
- `metafor`
- `purrr`
- `rlang`
- `cowplot`
- `ggh4x`
- `grid`
- `latex2exp`
- `ftExtra`
- `cli`
- `renv`
- `sessioninfo`
- `officer`
- `readxl`
- `dplyr`
- `magrittr`
