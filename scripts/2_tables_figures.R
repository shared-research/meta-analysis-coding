# Packages ----------------------------------------------------------------

library(tidyverse)
library(flextable)
library(officer)
library(here)

devtools::load_all()

# Loading Data ------------------------------------------------------------

dat_meta <- read_rds(here("objects", "dat_meta.rds"))
dat_meta_post <- read_rds(here("objects/dat_meta_post.rds"))

# Tables General Settings -------------------------------------------------

sect_properties <- prop_section(
    page_size = page_size(orient = "landscape"),
    type = "continuous",
    page_margins = page_mar(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
)

# Meta-analysis Table -----------------------------------------------------

## Column names ad paragraph using flextable::as_paragraph()

par_list <- list(
    # General
    outcome2 = as_paragraph("Outcome"),
    paper_id = as_paragraph("Paper"),
    # Experimental group
    #M_age_EX = as_paragraph("Age"),
    n_EX = as_paragraph("N"),
    M_pre_EX = as_paragraph("M", as_sub("pre")),
    M_post_EX = as_paragraph("M", as_sub("post")),
    SD_pre_EX = as_paragraph("SD", as_sub("pre")),
    SD_post_EX = as_paragraph("SD", as_sub("post")),
    # Control Group
    #M_age_CT = as_paragraph("Age"),
    n_CT = as_paragraph("N"),
    M_pre_CT = as_paragraph("M", as_sub("pre")),
    M_post_CT = as_paragraph("M", as_sub("post")),
    SD_pre_CT = as_paragraph("SD", as_sub("pre")),
    SD_post_CT = as_paragraph("SD", as_sub("post")),
    # Effect size
    #eff_size = as_paragraph("d", as_sub("ppc2")),
    #eff_size_var = as_paragraph("\u03C3", as_sup(2L), " d", as_sub("ppc2"))
    dpcc2 = as_paragraph("d", as_sub("ppc2"), "(\u03C3)")
)

## Characters for creating the upper header row

header_list <- c(
    rep("", 2),
    rep("Experimental Group", 5),
    rep("Control Group", 5),
    rep("", 1)
)

## Creating table

dat_table_meta <- dat_meta %>% 
    get_data(data_agg) %>% 
    mutate(n_EX = as.integer(n_EX),
           n_CT = as.integer(n_CT)) %>% 
    select(outcome2, paper_id,
           n_EX, M_pre_EX, M_post_EX, SD_pre_EX, SD_post_EX,
           n_CT, M_pre_CT, M_post_CT, SD_pre_CT, SD_post_CT,
           eff_size, eff_size_var) %>% 
    arrange(outcome2) %>% 
    mutate(dpcc2 = sprintf("%.2f (%.2f)", eff_size, eff_size_var)) %>% 
    select(-eff_size, -eff_size_var)

meta_table <- dat_table_meta %>% 
    flextable() %>% 
    colformat_double(digits = 2) %>% 
    add_par_header(par_list = par_list) %>% 
    autofit() %>% 
    theme_vanilla() %>% 
    add_header_row(values = header_list) %>% 
    merge_h(part = "header") %>% 
    merge_v(part = "body", j = 1) %>% 
    align(part = "header", align = "center") %>% 
    align(part = "body", align = "center") %>% 
    vline(j = c(2, 7, 12), part = "all") %>% 
    bold(i = 1:2, part = "header") %>% 
    fontsize(size = 7, part = "body") %>% 
    fontsize(size = 9, part = "header") %>% 
    width(j = 1:2, width = 1, unit = "in") %>% 
    width(j = 3:12, width = 0.3, unit = "in") %>% 
    width(j = 13, width = 1, unit = "in") %>% 
    fit_to_width(max_width = sect_properties$page_size$width)

## Saving as .docx file

save_as_docx(meta_table, 
             path = here("tables", "meta_table.docx"), 
             pr_section = sect_properties)

# Model Table -------------------------------------------------------------

model_table <- dat_meta_post %>% 
    get_data(tidy_meta_multi_fixed) %>% 
    multivariate_table(fit = get_model(dat_meta, meta_multi_fixed))

save_as_docx(model_table, 
             path = here("tables", "model_table.docx"), 
             pr_section = sect_properties)

# Forest Plot -------------------------------------------------------------

model_to_plot <- get_model(dat_meta, meta_multi_fixed)
data_to_plot <- get_data(dat_meta_post, tidy_meta_multi_fixed)

dat_forest_multi <- prep_data_forest_multi(model_to_plot, data_to_plot)

## Creating the polygons coordinates for the avg effect

get_diamonds_coord <- function(data){
    data %>% 
        mutate(x = c(.upper[1], eff_size[1], .lower[1], eff_size[1]),
               y = c(1, 1.15, 1, 0.85)) # hardcoding y coordinates
}

## Adding the coordinates to each effect

get_avg_diamonds <- function(data){
    data %>% 
        filter(paper_id == "Average") %>% 
        group_by(outcome2) %>% 
        expand_grid(id = 1:4) %>% 
        split(., .$outcome2) %>% 
        map(., get_diamonds_coord) %>% 
        bind_rows()
}

forest_multi <- dat_forest_multi %>% 
    forest_plot_multi() +
    geom_polygon(data = get_avg_diamonds(dat_forest_multi),
                 aes(x = x, y = y),
                 fill = "firebrick3")

# Paper Objects -----------------------------------------------------------

paper_objects <- list(
    fit = get_model(dat_meta, meta_multi_fixed), 
    data = get_data(dat_meta, data),
    data_agg = get_data(dat_meta, data_agg),
    dat_table_meta = dat_table_meta,
    meta_table = meta_table,
    model_table = model_table,
    forest_multi = forest_multi
)

# Saving ------------------------------------------------------------------

saveRDS(paper_objects, file = here("objects", "paper_objects.rds"))
ggsave("figures/forest_multi.pdf", paper_objects$forest_multi, width = 10, height = 6)
ggsave("figures/forest_multi.png", paper_objects$forest_multi, width = 10, height = 6)
