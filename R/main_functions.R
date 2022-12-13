# --- EFFECT SIZE UTILS --- #

# get_dppc2 -----------------------------------------------------------

# This function compute the dppc2

get_dppc2 <- function(mt_pre, mt_post, mc_pre, mc_post,
                      st_pre, st_post, sc_pre, sc_post,
                      nt, nc){
    
    mdiff <- (mt_post - mt_pre) - (mc_post - mc_pre)
    poolsd <- sqrt(((((nt - 1) * st_pre^2) + ((nc - 1) * sc_pre^2)) / (nt + nc - 2)))
    #cp <- 1 - (3/(4*(nt + nc - 2) - 1))
    cp <- metafor:::.cmicalc(nt + nc - 2)
    
    dppc2 <- cp * (mdiff / poolsd)
    
    return(dppc2)
    
}

# get_dppc2_var -----------------------------------------------------------

# This function compute the dppc2 variance. The rho argument is the pre-post
# correlation

get_dppc2_var <- function(dppc2, nt, nc, rho){
    
    dppc2_var <- 2*(1-rho) * (1/nt + 1/nc) + dppc2^2 / (2*(nt + nc))
    
    return(dppc2_var)
    
}

# compute_effect_size -----------------------------------------------------

compute_effect_size <- function(data, rho){
    data %>% 
        mutate(eff_size = get_dppc2(mt_pre = M_pre_EX, mt_post = M_post_EX, 
                                    mc_pre = M_pre_CT, mc_post = M_post_CT,
                                    st_pre = SD_pre_EX, st_post = SD_post_EX, 
                                    sc_pre = SD_pre_CT, sc_post = SD_post_CT,
                                    nt = n_EX, nc = n_CT),
               eff_size_var = get_dppc2_var(eff_size, 
                                            nt = n_EX, nc = n_CT, 
                                            rho = rho))
}

# --- TIDY META UTILS --- #

# tidy_meta_study ---------------------------------------------------------

# This function return the estimation of all effects from a fitted rma* model
# Is mainly used for plotting the forest plot

tidy_meta_study <- function(fit, conf_level = 0.95){
    broom.mixed::tidy(fit,
                      conf.int = TRUE,
                      conf.level = conf_level,
                      include_studies = TRUE)
}

# tidy_meta_params --------------------------------------------------------

tidy_meta_params <- function(fit, conf_level = 0.95){
    broom.mixed::tidy(fit,
                      conf.int = TRUE,
                      conf.level = conf_level)
}

# tidy_meta_tau -----------------------------------------------------------

# This function return the tau2 estimation from fitted models in a tidy way.
# Is ready to be binded to the parameters summary in order to have a dataframe
# with all estimations

# random = univariate random effect
# multilev = multilevel random effect model (aka 3 level model)
# mutlivar = multivariate random model

# The argument need to be specified because the tau estimations are logged in
# different ways depending on the model

tidy_meta_tau <- function(fit, type = c("fixed", "random", "multilev", "multivar")){
    
    type <- match.arg(type)
    
    if(type == "multivar"){
        meta_data <- list(
            term = rownames(fit$b),
            tau_value = fit$tau2
        )
    }else if(type == "multilev"){
        meta_data <- list(
            term = fit$s.names,
            tau_value = fit$sigma2
        )
    }else if(type == "random"){
        meta_data <- list(
            term = "tau",
            tau_value = fit$tau2
        )
    }
    
    tibble(
        term = meta_data$term,
        type = "tau",
        estimate = meta_data$tau_value,
        std.error = NA,
        statistic = NA,
        p.value = NA,
        conf.low = NA,
        conf.high = NA
    )
    
}

# tidy_meta ---------------------------------------------------------------

# This is a simple wrapper for tidy_meta_tau and tidy_meta_params in order to
# create the final dataframe

tidy_meta <- function(fit, conf_level = 0.95, type = c("fixed", "random", "multilev", "multivar")){
    
    type <- match.arg(type)
    
    out <- tidy_meta_params(fit, conf_level = conf_level)
    
    if(type != "fixed"){
        tau <- tidy_meta_tau(fit, type) 
        out <- bind_rows(out, tau)
    }
    
    out <- out %>% clean_param_names()
    
    return(out)
}

# --- AGGREGATION AND COV MATRIX UTILS --- #

# aggregate_effects -------------------------------------------------------

# This function aggregate effects sizes using the Borenstein (2009) approach
# the rho argument is the correlation between different measures. The function
# split data according to the study_i in order to combine the same outcome into
# a single effect size

aggregate_effects <- function(data, rho, split_by = paper, weighted = TRUE){
    split_by <- rlang::enexpr(split_by)
    data <- metafor::escalc(yi = eff_size, vi = eff_size_var, 
                   var.names = c("eff_size", "eff_size_var"),
                   data = data)
    dat_split <- split(data, pull(data, !!split_by))
    dat_split <- purrr::map_dfr(dat_split, function(x) metafor::aggregate.escalc(x, 
                                                                                 cluster = outcome2, 
                                                                                 rho = rho, 
                                                                                 weighted = weighted))
    bind_rows(dat_split) %>% 
        tibble() %>% 
        mutate(eff_size_se = sqrt(eff_size_var))
}

# get_block_cov_matrix ----------------------------------------------------

# This function is a simple wrapper for the metafor::vcalc() that create
# block covariance matrix to be used into the rma.mv() function to fit
# a multivariate model

get_block_cov_matrix <- function(data, rho){
    vcalc(eff_size_var, cluster = paper_id, obs = outcome2, data = data, rho = rho)
}

# --- CLEANING UTILS --- #

# put_names_cov_mat -------------------------------------------------------

put_names_cov_mat <- function(v, names){
    rownames(v) <- names
    colnames(v) <- paste0("O", 1:ncol(v)) # generic outcome name
    v <- data.frame(v)
    return(v)
}

# get_splitted_cov_mat ----------------------------------------------------

get_splitted_cov_mat <- function(V, data, split_by = paper_id, col_outcome = outcome2){
    
    # lazy eval
    split_by <- rlang::enexpr(split_by)
    col_outcome <- rlang::enexpr(col_outcome)
    
    # List of matrices
    V_list <- V_list <- blsplit(V, pull(data, !!split_by))
    
    # List of rownames
    dat_list <- split(data, pull(data, !!split_by))
    names_list <- map(dat_list, ~pull(.x, !!col_outcome))
    
    # Put row and columns names
    map2(V_list, names_list, put_names_cov_mat)
}

# make_dat_univariate -----------------------------------------------------

# This function split the main dataset according to the outcome in order
# to compute the univariate meta-analysis model

make_dat_univariate <- function(data){
    split(data, data$outcome2)
}

# clean_data --------------------------------------------------------------

# This function apply some data cleaning steps before aggregating data and
# fitting the models

clean_data <- function(data){
    data %>% 
        group_by(paper) %>% 
        mutate(effect_id = 1:n()) %>%
        ungroup() %>% 
        mutate(id = 1:nrow(.)) %>% 
        select(paper, paper_id, effect_id, id, everything())
}

# clean_param_names -------------------------------------------------------

# This function clean the tibble names after the tidy_meta function for better
# plotting and printing

clean_param_names <- function(data){
    data %>% 
        mutate(term = str_remove_all(term, "outcome2"))
}


# clean_column ------------------------------------------------------------

# Simple utils for removing patterns from specific columns. Use lazy-evaluation
# with rlang

clean_column <- function(data, column){
    
    column <- rlang::enexpr(column)
    
    data %>% 
        mutate(!!column := str_replace_all(!!column, "_", " "))
}


# clean_outcome_names ------------------------------------------------------

clean_outcome_names <- function(data, col = outcome){
    col <- rlang::enexpr(col)
    data %>% 
        mutate(!!col := case_when(!!col == "Accuratezza_FlessibilitàCognitiva" ~ "Cognitive Flexibility Acc.",
                                  !!col == "Accuratezza_Inibizione" ~ "Inhibition Acc.",
                                  !!col == "Accuratezza_Memoria_Lavoro" ~ "Working Memory Acc.",
                                  !!col == "Accuratezza_Planning" ~ "Planning Acc.",
                                  !!col == "Problem_solving" ~ "Problem Solving",
                                  !!col == "Tempo_Pianificazione_Planning" ~ "Planning Time",
                                  TRUE ~ !!col))
}

# add_paper_cond ----------------------------------------------------------

add_paper_cond <- function(data,
                           m_cor = 0.7, 
                           b_cor = 0.5, 
                           mu_cor = 0.5,
                           model_t = "Multivariate",
                           meta = "Fixed"){
    
    cond_to_eval <- expr(morris_cor == m_cor & boren_cor == b_cor & 
                             multi_cor == mu_cor & model_type == model_t & 
                             meta_type == meta)
    
    data %>% 
        mutate(paper_cor = ifelse(eval(cond_to_eval),
                                  "yes", "no"))
    
}