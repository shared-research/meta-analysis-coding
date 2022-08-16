---
title: "Meta analysis coding title"
subtitle: "Supplementary materials"
author:
    - Chiara Montuori^[Department of Developmental and Social Psychology, chiara.montuori@phd.unipd.it]
    - Filippo Gambarota^[Department of Developmental and Social Psychology, filippo.gambarota@phd.unipd.it]
    - Gianmarco Altoè^[Department of Developmental and Social Psychology, gianmarco.altoe@unipd.it]
    - Barbara Arfè^[Department of Developmental and Social Psychology, barbara.arfe@unipd.it]
bibliography: ../files/references.bib
csl: ../files/apa7.csl
header-includes:
    - \usepackage{fullpage}
    - \usepackage{pdflscape}
    - \usepackage{afterpage}
output:
    bookdown::pdf_document2:
        keep_tex: true
        latex_engine: lualatex
    bookdown::html_document2:
        keep_md: yes
        toc: true
        toc_float: true
---













\pagebreak

# Introduction

This document explain in details the statistical approach of the main paper. All the analysis are conducted with R [@r-lang] and all scripts are available on the [Open Science Framework repository](https://osf.io/uvbcd/).

- Section \@ref(dataset) present the full dataset description with all extracted variables. Despite only a subset of these variables are used to compute the meta-analysis model, we reported all extracted study-level information.
- Section \@ref(effectsize) describe the computation of the effect size measure and the sampling variance.
- Section \@ref(models) describe the main meta-analysis model presented in the paper together with modeling alternatives
- Section \@ref(correlations) describe the correlations imputation that is necessary to compute the meta-analysis model
- Section \@ref(sens) present the sensitivity analysis for different correlations values and different meta-analytic models 

## Meta-analysis

The meta-analysis is the statistical procedure to combine information from already conducted studies. The core aspect of a meta-analysis is to weight each included effect according to the amount of information (i.e., precision) that the study provide. In order to compute a meta-analysis model we need:

- The effect size measure
- The effect size variance (i.e., the inverse of the precision)
- The meta-analysis model

## Papers

Almost all published papers reported enough information to calculate the effect size. Di Lieto and colleagues [-@Di_Lieto2020-tc] provided us raw data in order to calculate all relevant parameters. They use a design with 3 time points (T0 = pre, T1 = post and T2 = follow-up) for both the experimental and control group. Furthermore they measured several different outcomes measures. In order to include the paper we performed these pre-processing steps:

- We selected only the T0 and T1 in order to have a single pre-post measure
- We included all participants that have both time points (T0 and T1) in at least one outcome measure. For example, if a child has T0 and T1 for the outcome $x$ but only T0 for the outcome $y$ we keep the child in the analysis. This create a situation where different outcomes have different sample sizes but maximize the amount of available information.

# Dataset description {#dataset}

The main dataset can be found on the online OSF repository under `data/raw/meta_table_cleaned.csv`. The following list describes the meaning of each column: 



- **paper**: Unique id for every paper 
 - **paper_id**: Unique id for every paper with authors and year 
 - **author**: Paper's authors 
 - **year**: Publication year 
 - **Participants**: Participants' characteristics: ST (typical development), SA (atypical development) 
 - **n_EX**: Sample size for the experimental group 
 - **n_CT**: Sample size for the control group 
 - **male_EX**: Number of males for the experimental group 
 - **female_EX**: Number of females for the experimental group 
 - **male_CT**: Number of males for the control group 
 - **female_CT**: Number of females for the control group 
 - **Training_length**: Lenght of the training (in weeks) 
 - **Minute_Training_length**: Lenght of a single training session (in minutes) 
 - **Training_Mode**: Type of training activity: CV (Virtual Coding), ER (Educational Robotics) 
 - **Training_CV_Type**: Types of virtual coding activities: structured and unstructured 
 - **M_SES_EX**: Mean socio-economic status (SES) for the experimental group 
 - **SD_SES_EX**: Standard deviation of socio-economic status (SES) for the experimental group 
 - **M_SES_CT**: Mean socio-economic status (SES) for the control group 
 - **SD_SES_CT**: Standard deviation of socio-economic status (SES) for the control group 
 - **M_age_CT**: Mean age for the control group 
 - **M_age_EX**: Mean age for the experimental group 
 - **outcome**: The outcome measure/test used 
 - **M_pre_EX**: Mean pre-training score for the experimental group 
 - **M_post_EX**: Mean post-training score for the experimental group 
 - **SD_pre_EX**: Standard deviation pre-training score for the experimental group 
 - **SD_post_EX**: Standard deviation post-training score for the experimental group 
 - **M_pre_CT**: Mean pre-training score for the control group 
 - **M_post_CT**: Mean post-training score for the control group 
 - **SD_pre_CT**: Standard deviation pre-training score for the control group 
 - **SD_post_CT**: Standard deviation post-training score for the control group 
 - **outcome2**: The recoded outcome variable 
 - **Flip_NoFlip**: Whether the effect size need to be flipped in order to have the same direction (i.e., positive values, the training improve performance) 
 - **eff_size**: The computed $dpcc_2$ 
 - **eff_size_var**: The computed $dpcc_2$ variance 
 - **eff_size_se**: The computed $dpcc_2$ standard error 

## Pre-processing steps

The raw dataset is available on the online [Open Science Framework repository](https://osf.io/uvbcd/). We performed few minimal pre-processing steps:

- renaming relevant columns
- re-coding the `outcome` variable into the `outcome2` variable
- separating the Arfè et al. (2019) paper into two disticint papers (see below)
- converting outcomes names in English

### Arfè et al. (2019)

Arfè et al. (2019) is the only paper that contain multiple sub-studies creating a multilevel structure. In order to reduce the dataset complexity and given that the two sub-studies refers to a different pool of subjects, we decided to consider this studies as two separate papers, creating Arfè et al. (2019a) and Arfè et al. (2019b).

## Multiple effects for the same outcome{#aggregation}

We decided to recode the `outcome` into `outcome2` in order to have a small set of outcomes according to the underlying psychological construct. For example, if the test $x$ and the test $y$ are both considered a Working Memory measure we decided to recode both outcomes as "Working Memory".
This reduce the dataset complexity but create a situation where a single paper has multiple effects belonging to the same outcome. We decided to aggregate multiple effects of the same `outcome2` using the approach suggested by Borenstein et al. [-@Borenstein2009-mo, pp.225-233] implemented using the `metafor::aggregate.escalc()` function (https://wviechtb.github.io/metafor/reference/aggregate.escalc.html). 
To note, we use a slightly different approach compared to Borenstein et al.[-@Borenstein2009-mo, pp.225-233] computing an *inverse-variance weighted average* instead of an *un-weighted average*. Essentially, the weighted average combine multiple effect sizes taking into account the precision (i.e., the inverse of the variance).

We used the `metafor::aggregate.escalc()` as follows:


```r
aggregate_effects <- function(data, rho, split_by = paper, weighted = TRUE){
    split_by <- rlang::enexpr(split_by)
    data <- metafor::escalc(yi = eff_size, vi = eff_size_var, 
                   var.names = c("eff_size", "eff_size_var"),
                   data = data)
    dat_split <- split(data, pull(data, !!split_by))
    dat_split <- purrr::map_dfr(dat_split, function(x) metafor::aggregate.escalc(x, cluster = outcome2, rho = rho, weighted = weighted))
    bind_rows(dat_split) %>% tibble()
}
```

This is a wrapper of the `metafor::aggregate.escalc()` function that split the dataset according to the `paper` and then aggregate within each `outcome2`.

# Effect size computation {#effectsize}

All included studies used a Pretest-Posttest-Control Group (PPC) design where the treated group is compared with an age-matched control group before and after a certain treatment. Morris [-@Morris2008-pe] described different effect size indexes for the PPC design. We decided to use the $d_{ppc2}$ index because it provide an unbiased estimation of the true effect size and the sampling variance can be analytically computed. The $d_{ppc2}$ is calculated as indicated in Equation \@ref(eq:dpcc).

\begin{equation} 
  d_{pcc_2} = c_p \frac{(M_{T,post} - M_{T,pre}) - (M_{C,post} - M_{C,pre})}{SD_{pooled,pre}}
  (\#eq:dpcc)
\end{equation}

The numerator describe the actual mean difference between treated and control group of the respective *pre-post* scores. The difference is standardized using only a pooled standard deviation of pre-test scores. The variance is calculated in Equation \@ref(eq:dpcc-sd):

\begin{equation} 
  \sigma^2(d_{pcc_2}) = c^2_p(1-\rho)(\frac{n_t+n_c}{n_t n_c})(\frac{n_t+n_c-2}{n_t+n_c-4})(\frac{1 + \Delta^2}{2(1-\rho)(\frac{n_t+n_c}{n_t n_c})})
  (\#eq:dpcc-sd)
\end{equation}

Both the $d_{ppc2}$ and the variance are adjusted using a small sample correction factor $c_p$ calculated in Equation \@ref(eq:cp)

\begin{equation} 
  c_p = 1 - \frac{3}{4(2n_t + 2n_c - 4) - 1}
  (\#eq:cp)
\end{equation}

Positive values of the effect size means that the treatment increased the performance of the experimental group. However for some measures, higher values corresponds to worse performance (e.g., number of errors). In this case we computed the effect size as described in the previous section and then we flipped the sign obtaining the same interpretation regardless the original measure.

For the actual computation we wrote two functions following the documentation of the `metafor` package: (http://www.metafor-project.org/doku.php/analyses:morris2008), `get_dpcc2()`:


```r
get_dpcc2 <- function(mt_pre, mt_post, mc_pre, mc_post,
                          st_pre, st_post, sc_pre, sc_post,
                          nt, nc){
    
    mdiff <- (mt_post - mt_pre) - (mc_post - mc_pre)
    poolvar <- sqrt(((((nt - 1) * st_pre^2) + ((nc - 1) * sc_pre^2)) / (nt + nc - 2)))
    cp <- 1 - (3/(4*(nt + nc - 2) - 1))
    
    dpcc2 <- (mdiff / poolvar) * cp
    
    return(dpcc2)
    
}
```

and `get_dpcc2_var()`:


```r
get_dpcc2_var <- function(dpcc2, nt, nc, rho){
    
    dpcc2_var <- 2*(1-rho) * (1/nt + 1/nc) + dpcc2^2 / (2*(nt + nc))
    
    return(dpcc2_var)
    
}
```

# Meta-analysis models {#models}

## Univariate vs multivariate model {#uni-vs-multi}

When multiple variables are collected from the same pool of participants (e.g., multiple outcomes studies) and/or there are multiple experiments within the same paper, effects size cannot be considered independent [@Cheung2019-po; @Cheung2014-fg]. Instead of pooling multiple effects into a single measure or keeping only a single effect per study, we computed a multivariate meta-analysis model taking into account the dependency (i.e., the correlation) between multiple effects size calculated on the same pool of participants. In order to compute the multivariate model we need to include or impute the variance-covariance matrix that describe how multiple effects correlate within each study. The section \@ref(correlations) explain our imputation approach.

## Random-effect vs fixed-effect model

Another important choice for the meta-analysis is between the fixed-effect and the random-effect model. Under the fixed-effect model we make the assumption of a single *true* effect size that is estimated by a sample of studies. For this reason, the variability observed in the empirical meta-analysis is only caused by the sampling error.

On the other side, the *random-effect* model assume a *distribution* of *true* effects. The model estimate the average effect (the *mean* of the distribution) and the between-study heterogeneity $\tau^2$ (i.e., the *variance* of the distribution). The critical assumption is that the observed heterogeneity is composed by sampling variability (as the *fixed-effect* model) and *true* between-study variability. The latter can be caused by different experimental designs, participants features or other study-level variables and can be explained using moderators (i.e., meta-regression).

The goal of the *random-effect* model is to generalize the meta-analytic findings at the population level while the *fixed-effect* model is focused on a specific pool of studies. 

We decided to use a *fixed-effect* model for several reasons. Under the *random-effect* framework, the multivariate model estimate the average effect and the variability ($\tau^2$) for each included outcome. Firstly, one of the main reason to estimate $\tau^2$ (i.e., the focus of the *random-effect* model) is to explain the observed heterogeneity with study-level variables. Given our limited pool of studies, we did not included a meta-regression analysis.
Crucially, with a limited number of studies the $\tau^2$ estimation can be strongly biased [@Veroniki2016-nw] influencing also the pooled effect estimation [@Borenstein2009-mo, pp.73-75]. Finally, the *fixed-effect* is less complex in terms of model parameters beacuse we need to estimate only the average effect of each included outcome, taking into account the statistical dependency (see \@ref(uni-vs-multi)).

## Modelling functions

We tested each model parameter (i.e., average effect for a specific outcome) using the Wald *z-test* with $\alpha = 0.05$.  
For the *univariate fixed-effect* model we use the `metafor::rma()` function:


```r
fit_uni_fixed <- function(data){
    rma(yi = eff_size, 
        vi = eff_size_var,
        method = "FE",
        data = data)
}
```

For the *univariate random-effect* model we use the `metafor::rma()` function and the REML estimator for $\tau^2$:


```r
fit_uni_random <- function(data){
    rma(yi = eff_size, 
        vi = eff_size_var,
        method = "REML",
        data = data)
}
```

For the *multivariate fixed-effect* model we use the `metafor::rma.mv()` function:


```r
fit_multi_fixed <- function(data, cov_matrix){
    rma.mv(
        yi = eff_size,
        V = cov_matrix,
        mods = ~ 0 + outcome2,
        data = data)
}
```

In this case we used `outcome2` as moderator with a *cell-mean parametrization* [@Schad2020-ht] (i.e., removing the intercept `0 + outcome2`). In this way model parameters and statistical tests correspond directly to the average effect for each outcome.

For the *multivariate random-effect* model we use the `metafor::rma.mv()` function:


```r
fit_multi_random <- function(data, cov_matrix){
    rma.mv(
        yi = eff_size,
        V = cov_matrix,
        mods = ~ 0 + outcome2,
        random = ~ outcome2|paper_id,
        struct = "UN",
        method = "ML",
        data = data)
}
```

As explained in the previous section, this model estimate also a $\tau^2$ for each outcome (`random = ~ outcome2|paper_id`).

### Variance-covariance matrix

For the *multivariate* model is necessary to include a block variance-covariance matrix. Essentially, each study has $n$ different outcomes and we need a $n \times n$ variance-covariance matrix where the diagonal is the $d_{pcc2}$ variance for a specific outcome and off-diagonal elements are the covariance between pairs of outcomes. The covariance between two outcomes is calculated as the product of the correlation and the variances:

\begin{equation} 
  Cov_{O_1,O_2} = \rho_{(O_1,O_2)}\sigma_{O_1}\sigma_{O_2}
  (\#eq:vcov)
\end{equation}

Combining all study-level matrices we obtain a full block-variance-covariance matrix to use within the `rma.mv()` function. The `metafor::vcalc()` function allow to automatically create the full matrix specifying the assumed correlation and a clustering variable:


```r
get_block_cov_matrix <- function(data, rho){
    vcalc(eff_size_var, cluster = paper_id, obs = effect_id, data = data, rho = rho)
}
```

# Correlations inputation{#correlations}

As explained in previous sections, we decided to use a *fixed-effect multivariate* model. In order to compute the model we need to include 3 correlations measures:

- For the effect size computation, we need the correlation between pre-post (*pre-post correlation*) scores for both groups
- For the aggregation of multiple effects within the same paper (*aggregation correlation*) (see Section \@ref(aggregation)) we need the correlation between the effects
- For the actual multivariate model we need the full variance-covariance matrix of different outcomes. To create the matrix we need to include the correlation between outcomes within each study (*multivariate correlation*).

Given that correlations are rarely reported in published papers, we decided to impute these values and assessing the impact with a multiverse-like approach [@Steegen2016-lz]. In particular we used:

- A $\rho_{pre-post}$ of 0.5, 0.7 and 0.9
- A $\rho_{agg}$ of 0.3, 0.5, 0.7
- A $\rho_{multi}$ of 0.3, 0.5 and 0.7

# Meta-analysis table

Table \@ref(tab:metatab) describe all included papers with pre-post scores for the experimental and control group along with the computed effect size. Furthermore, the table is organized grouping the affects according to the considered outcome in order to highlight the multivariate data structure.

```{=html}
<template id="cf04a4cd-1588-426b-b67d-7a80fd89a272"><style>
.tabwid table{
  border-spacing:0px !important;
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 1.275em;
  margin-bottom: 1.275em;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-a3c4621e{}.cl-a3b5e32e{font-family:'Arial';font-size:9pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-a3b5e32f{font-family:'Arial';font-size:5.4pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;position: relative;top:2.7pt;}.cl-a3b5e330{font-family:'Arial';font-size:5.4pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;position: relative;bottom:2.7pt;}.cl-a3b5e331{font-family:'Arial';font-size:7pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-a3b60a5c{margin:0;text-align:center;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-a3b6a67e{width:21.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6a67f{width:21.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 1pt solid rgba(102, 102, 102, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6a680{width:36pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6a681{width:21.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 1pt solid rgba(102, 102, 102, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6a682{width:36pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 1pt solid rgba(102, 102, 102, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6a683{width:72pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 1pt solid rgba(102, 102, 102, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6a684{width:72pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6a685{width:21.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 1pt solid rgba(102, 102, 102, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6a686{width:21.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 1pt solid rgba(102, 102, 102, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6a687{width:21.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6a688{width:72pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 1pt solid rgba(102, 102, 102, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6cd98{width:36pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 1pt solid rgba(102, 102, 102, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6cd99{width:72pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6cd9a{width:36pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6cd9b{width:21.6pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6cd9c{width:21.6pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 1pt solid rgba(102, 102, 102, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6cd9d{width:72pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 1pt solid rgba(102, 102, 102, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6cd9e{width:36pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 1pt solid rgba(102, 102, 102, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6cd9f{width:36pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6cda0{width:72pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6cda1{width:21.6pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 1pt solid rgba(102, 102, 102, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6cda2{width:21.6pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6f49e{width:36pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6f49f{width:36pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 1pt solid rgba(102, 102, 102, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6f4a0{width:21.6pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 1pt solid rgba(102, 102, 102, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6f4a1{width:72pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 1pt solid rgba(102, 102, 102, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6f4a2{width:72pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-a3b6f4a3{width:21.6pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 1pt solid rgba(102, 102, 102, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-a3c4621e'>
```
<caption class="Table Caption">

(\#tab:metatab)Pre-post scores (M = mean, SD = standard deviation, N = sample size) for the control and experimental group. The last two columns represent the computed effect size and the variance. The effects are organized according to the outcome value and the corresponding paper.

</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td  colspan="2"class="cl-a3b6f4a2"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e"></span></p></td><td  colspan="5"class="cl-a3b6f4a0"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">Experimental Group</span></p></td><td  colspan="5"class="cl-a3b6f4a0"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">Control Group</span></p></td><td  colspan="2"class="cl-a3b6f49f"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">Effect size</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6f4a2"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">Outcome</span></p></td><td class="cl-a3b6f4a1"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">Paper</span></p></td><td class="cl-a3b6f4a0"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">N</span></p></td><td class="cl-a3b6cda2"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">M</span><span class="cl-a3b5e32f">pre</span></p></td><td class="cl-a3b6cda2"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">M</span><span class="cl-a3b5e32f">post</span></p></td><td class="cl-a3b6cda2"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">SD</span><span class="cl-a3b5e32f">pre</span></p></td><td class="cl-a3b6f4a3"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">SD</span><span class="cl-a3b5e32f">post</span></p></td><td class="cl-a3b6f4a0"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">N</span></p></td><td class="cl-a3b6cda2"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">M</span><span class="cl-a3b5e32f">pre</span></p></td><td class="cl-a3b6cda2"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">M</span><span class="cl-a3b5e32f">post</span></p></td><td class="cl-a3b6cda2"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">SD</span><span class="cl-a3b5e32f">pre</span></p></td><td class="cl-a3b6f4a3"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">SD</span><span class="cl-a3b5e32f">post</span></p></td><td class="cl-a3b6f49f"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">d</span><span class="cl-a3b5e32f">ppc2</span></p></td><td class="cl-a3b6f49e"><p class="cl-a3b60a5c"><span class="cl-a3b5e32e">σ</span><span class="cl-a3b5e330">2</span><span class="cl-a3b5e32e"> d</span><span class="cl-a3b5e32f">ppc2</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td  rowspan="2"class="cl-a3b6a684"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Cognitive Flexibility Acc.</span></p></td><td class="cl-a3b6a683"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Di Lieto et al. (2020a)</span></p></td><td class="cl-a3b6a681"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">96</span></p></td><td class="cl-a3b6a67e"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">7.22</span></p></td><td class="cl-a3b6a67e"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">8.39</span></p></td><td class="cl-a3b6a67e"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.06</span></p></td><td class="cl-a3b6a67f"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.69</span></p></td><td class="cl-a3b6a681"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">91</span></p></td><td class="cl-a3b6a67e"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">7.19</span></p></td><td class="cl-a3b6a67e"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">8.27</span></p></td><td class="cl-a3b6a67e"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.10</span></p></td><td class="cl-a3b6a67f"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.66</span></p></td><td class="cl-a3b6a682"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.04</span></p></td><td class="cl-a3b6a680"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.01</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Di Lieto et al. (2020b)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">18</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">5.72</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">6.61</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.37</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.12</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">18</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">7.05</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">7.05</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.15</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.41</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.38</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.07</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  rowspan="5"class="cl-a3b6cd99"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Inhibition Acc.</span></p></td><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Arfè et al. (2019a)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">44</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.63</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.79</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.45</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.01</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">36</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">5.76</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.16</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.36</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.72</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.21</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.02</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Arfè et al. (2019b)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">19</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.03</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.58</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.80</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.81</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">19</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.74</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.82</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.79</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.27</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.46</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.05</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Arfè et al. (2020)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">88</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">5.71</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.97</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">6.01</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.23</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">91</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">5.66</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.12</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.26</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.50</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.44</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.01</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Di Lieto et al. (2020a)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">96</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">5.72</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.39</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.62</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.42</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">91</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.61</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.02</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.53</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.13</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.14</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.01</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Di Lieto et al. (2020b)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">18</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">8.16</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">6.50</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">7.64</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">6.14</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">18</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">6.58</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">5.38</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.39</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">6.66</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.06</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.05</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  rowspan="3"class="cl-a3b6cd99"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Planning Acc.</span></p></td><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Arfè et al. (2019a)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">44</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">7.13</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">9.84</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.46</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.45</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">36</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">5.75</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">7.57</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.91</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.80</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.30</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.02</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Arfè et al. (2019b)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">19</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">8.39</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">11.42</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.51</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.60</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">19</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">8.92</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">8.55</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.23</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.29</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.99</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.06</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Arfè et al. (2020)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">88</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">6.02</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">9.59</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.88</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.39</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">91</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.95</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">7.05</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.34</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.18</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.61</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.01</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  rowspan="9"class="cl-a3b6cd99"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Problem Solving</span></p></td><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Akcaoglu &amp; Koehler (2014)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">20</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">-1.11</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">-0.07</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.23</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.11</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">24</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">-1.17</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">-1.29</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.91</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.15</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.05</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.05</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Arfè et al. (2019a)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">44</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.31</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">6.12</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.46</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.06</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">36</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.09</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.68</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.60</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.92</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.79</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.03</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Arfè et al. (2019b)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">19</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">6.05</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">7.16</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.08</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.96</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">19</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">5.58</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">5.21</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.17</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.08</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.29</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.08</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Arfè et al. (2020)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">88</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.90</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">5.86</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.45</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.87</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">91</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.57</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.74</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.08</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.08</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.75</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.02</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Erol &amp; Cirak (2021)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">16</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">93.25</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">99.12</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">16.32</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">12.39</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">18</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">89.11</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">86.72</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">10.79</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">10.88</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.59</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.08</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">La Paglia et al. (2017)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">30</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">54.19</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">58.09</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">7.52</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">6.89</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">30</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">56.29</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">56.86</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">11.58</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">12.22</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.34</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.04</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Nam et al. (2010)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">30</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">9.55</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">13.36</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.50</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.55</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">30</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">10.09</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">11.50</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.16</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.84</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.71</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.04</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Pardamean et al. (2011)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">43</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">9.98</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">12.05</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.08</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.90</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">42</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">9.35</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">10.76</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.02</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.23</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.26</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.02</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">WonNam et al. (2019)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">25</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">25.72</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">37.16</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">6.29</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.88</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">28</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">26.79</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">29.75</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">5.07</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">6.33</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.47</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.07</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  rowspan="2"class="cl-a3b6cd99"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Working Memory Acc.</span></p></td><td class="cl-a3b6a688"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Di Lieto et al. (2020a)</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">96</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.34</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.93</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.92</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.31</span></p></td><td class="cl-a3b6a685"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">91</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.15</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">4.23</span></p></td><td class="cl-a3b6a687"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.76</span></p></td><td class="cl-a3b6a686"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.90</span></p></td><td class="cl-a3b6cd98"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.26</span></p></td><td class="cl-a3b6cd9a"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.01</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-a3b6cd9d"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">Di Lieto et al. (2020b)</span></p></td><td class="cl-a3b6cda1"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">18</span></p></td><td class="cl-a3b6cd9b"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.36</span></p></td><td class="cl-a3b6cd9b"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.34</span></p></td><td class="cl-a3b6cd9b"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.25</span></p></td><td class="cl-a3b6cd9c"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.94</span></p></td><td class="cl-a3b6cda1"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">18</span></p></td><td class="cl-a3b6cd9b"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.11</span></p></td><td class="cl-a3b6cd9b"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">3.42</span></p></td><td class="cl-a3b6cd9b"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">1.14</span></p></td><td class="cl-a3b6cd9c"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">2.01</span></p></td><td class="cl-a3b6cd9e"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">-0.21</span></p></td><td class="cl-a3b6cd9f"><p class="cl-a3b60a5c"><span class="cl-a3b5e331">0.04</span></p></td></tr></tbody></table></div></template>
<div class="flextable-shadow-host" id="7b140cfe-1373-4d82-93fa-657a488347b8"></div>
<script>
var dest = document.getElementById("7b140cfe-1373-4d82-93fa-657a488347b8");
var template = document.getElementById("cf04a4cd-1588-426b-b67d-7a80fd89a272");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;text-align:center;";
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

```


# Multiverse analysis{#sens}

We assessed the impact of our modelling assumptions using a sensitivity analysis approach. In particular we considered:

- 4 meta-analysis models:
    - *univariate fixed-effect*
    - *univariate random-effect*
    - *multivariate fixed-effect*
    - *multivariate random-effect*
- All correlations combinations ($3 \rho_{pre-post} \times 3 \rho_{aggregation} \times 3 \rho_{multivariate}$)

In this way we computed 72 meta-analysis directly assessing the impact on the parameters estimation and p-values^[For *multivariate* models (*fixed* and *random*) we have 3x3x3 correlation combinations while for *univariate* models (*fixed* and *random*) only 3x3 ($\rho_{pre-post}$ and $\rho_{agg}$)]. All plots in the following pages depicted the sensitivity analysis for each outcome. Triangles indicates parameters with $p < 0.05$ while points indicates parameters with $p > 0.05$. Red shapes indicates the correlations/model combination included in the paper. Clearly, *random effect models* (both *univariate* and *multivariate*) is associated with wider confidence intervals (i.e., less precise estimation) compared to *fixed effect models*. To note, given that the $\rho_{multi}$ is not relevant for *univariate* models, the lower part of the figure is repeated.





<img src="supplementary_files/figure-html/cognitive-flexibility-acc-plot-1.svg" width="864" style="display: block; margin: auto;" />

<img src="supplementary_files/figure-html/inhibition-accuracy-plot-1.svg" width="864" style="display: block; margin: auto;" />

<img src="supplementary_files/figure-html/planning-accuracy-plot-1.svg" width="864" style="display: block; margin: auto;" />

<img src="supplementary_files/figure-html/problem-solving-plot-1.svg" width="864" style="display: block; margin: auto;" />

<img src="supplementary_files/figure-html/working-memory-accuracy-plot-1.svg" width="864" style="display: block; margin: auto;" />



# References
