# Environment -------------------------------------------------------------

library(tidyverse)
library(metafor)
library(broom.mixed)
library(here)

devtools::load_all()

# Prepare data for models

prep_data <- function(data){
    data %>% 
        # This flip effects size for having positive effect = treatment works
        mutate(eff_size_before_flip = eff_size,
               eff_size = ifelse(Flip_NoFlip == "Flip",
                                 eff_size * -1,
                                 eff_size),
               # Converting to factors
               year = factor(year),
               paper = factor(paper),
               paper_id = factor(paper_id))
}

# Importing Data ----------------------------------------------------------

dat <- read_rds(here("data", "raw", "meta_table_fit.rds"))

# Computing Effect Size ---------------------------------------------------

# We compute the effect size variance for a different range of pre-post
# correlation

morris_cor <- c(0.5, 0.7, 0.9) # range of correlations for morris
boren_cor <- c(0.3, 0.5, 0.7) # range of correlations for borenstein
multi_cor <- c(0.3, 0.5, 0.7) # range of correlations for multivariate model

dat_grid <- expand_grid(
    morris_cor,
    boren_cor,
    multi_cor
)

# Creating a list of datasets with the effect size variance computed with different
# correlations values

dat_meta <- dat_grid %>% 
    mutate(data = map(morris_cor, function(rho) compute_effect_size(dat, rho)),
           data = map(data, prep_data),
           data_agg = map2(data, boren_cor, aggregate_effects),
           data_agg = map(data_agg, clean_data),
           data_uni = map(data_agg, make_dat_univariate),
           cov_mat = map2(data_agg, multi_cor, get_block_cov_matrix), # the full cov mat
           cov_mat_list = map2(cov_mat, data_agg, get_splitted_cov_mat)) # splitted cov mat

# Models ------------------------------------------------------------------

# fitting all models
dat_meta <- dat_meta %>% 
    mutate(meta_uni_fixed = map_depth(data_uni, fit_uni_fixed, .depth = 2),
           meta_uni_random = map_depth(data_uni, fit_uni_random, .depth = 2),
           meta_multi_random = map2(data_agg, cov_mat, fit_multi_random),
           meta_multi_fixed = map2(data_agg, cov_mat, fit_multi_fixed))

#  Post-processing --------------------------------------------------------

# extracting all relevant parameters
dat_meta_post <- dat_meta %>% 
    mutate(tidy_meta_uni_fixed = map_depth(meta_uni_fixed, 
                                           tidy_meta, 
                                           type = "fixed", .depth = 2),
           tidy_meta_uni_fixed = map(tidy_meta_uni_fixed, bind_rows, .id = "outcome"),
           tidy_meta_uni_random = map_depth(meta_uni_random, 
                                           tidy_meta, 
                                           type = "random", .depth = 2),
           tidy_meta_uni_random = map(tidy_meta_uni_random, bind_rows, .id = "outcome"),
           tidy_meta_multi_random = map(meta_multi_random, tidy_meta, type = "multivar"),
           tidy_meta_multi_fixed = map(meta_multi_fixed, tidy_meta, type = "multivar")) %>% 
    select(ends_with("cor"), data_agg, starts_with("tidy"), starts_with("meta_"))

# Plots -------------------------------------------------------------------

dat_meta_post <- dat_meta_post %>% 
    mutate(plot_multi_fixed = map2(meta_multi_fixed, tidy_meta_multi_fixed, prep_data_forest_multi),
           plot_multi_fixed = map(plot_multi_fixed, forest_plot_multi),
           plot_multi_random = map2(meta_multi_random, tidy_meta_multi_random, prep_data_forest_multi),
           plot_multi_random = map(plot_multi_random, forest_plot_multi))

# Saving ------------------------------------------------------------------

saveRDS(dat_meta, here("objects", "dat_meta.rds"))
saveRDS(dat_meta_post, here("objects", "dat_meta_post.rds"))
