# Packages ----------------------------------------------------------------

library(tidyverse)
library(latex2exp)
library(ggh4x)
library(here)

# Loading Data ------------------------------------------------------------

dat_meta <- read_rds(here("objects", "dat_meta.rds"))
dat_meta_post <- read_rds(here("objects", "dat_meta_post.rds"))

# Dataframe metadata ------------------------------------------------------

# Loading a table with column description

df_metadata <- readxl::read_xlsx(here("data", "metadata", "df_metadata.xlsx"))

df_metadata %>% 
    drop_na() %>% 
    filter(description != "no") %>% 
    mutate(knitr = sprintf("- **%s**: %s \n", colname, description)) %>% 
    saveRDS(., here("objects", "df_metadata.rds"))

# Sensitivity Analysis Plot -----------------------------------------------

# Prepare data

dat_sens <- dat_meta_post %>% 
    select(ends_with("cor"), tidy_meta_uni_fixed, tidy_meta_multi_random, tidy_meta_multi_fixed, tidy_meta_uni_random) %>% 
    pivot_longer(starts_with("tidy"), names_to = "model", values_to = "data") %>% 
    unnest(data) %>% 
    filter(type == "summary") %>% 
    mutate(outcome = ifelse(is.na(outcome), term, outcome),
           meta_type = ifelse(str_detect(model, "fixed"), "Fixed", "Random"),
           model_type = ifelse(str_detect(model, "multi"), "Multivariate", "Univariate"),
           morris_cor_parse = factor(morris_cor, labels = create_cor_labels(., morris_cor, "pre-post")),
           boren_cor_parse = factor(boren_cor, labels = create_cor_labels(., boren_cor, "agg")),
           multi_cor_parse = factor(multi_cor, labels = create_cor_labels(., multi_cor, "multi"))) %>% 
    clean_outcome_names() %>% 
    prep_data_sens()

# Plotting

dat_sens <- dat_sens %>% 
    group_nest(outcome, keep = TRUE) %>% 
    mutate(plot_sens_estimate = map(data, sensitivity_plot_estimate),
           plot_sens_pvalue = map(data, sensitivity_plot_pvalue))

# Save Plots

saveRDS(dat_sens, "objects/dat_sens.rds")