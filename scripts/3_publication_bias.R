# Packages ----------------------------------------------------------------

library(tidyverse)
library(here)
library(metafor)

# Loading Data ------------------------------------------------------------

dat_meta <- read_rds("objects/dat_meta.rds")

# problem solving
prob_solving <- get_model(dat_meta, meta_uni_fixed)$`Problem Solving`

# egger regression with the univariate fixed-effect model
res <- regtest(prob_solving)

# Saving ------------------------------------------------------------------

saveRDS(res, "objects/publication_bias.rds")