# fit_multi_random --------------------------------------------------------

# This function fit a multivariate model using the cov_matrix argument

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

# fit_multi_fixed ---------------------------------------------------------

# This function fit a multivariate fixed effect meta-analysis

fit_multi_fixed <- function(data, cov_matrix){
    rma.mv(
        yi = eff_size,
        V = cov_matrix,
        mods = ~ 0 + outcome2,
        data = data)
}

# fit_uni_fixed -----------------------------------------------------------

# This function fit a fixed-effect model. In the project is used to fit
# a specific moderator level

fit_uni_fixed <- function(data){
    rma(yi = eff_size, 
        vi = eff_size_var,
        method = "FE",
        data = data)
}


# fit_uni_random ----------------------------------------------------------

# This function fit a random-effect model. In the project is used to fit
# a specific moderator level

fit_uni_random <- function(data){
    rma(yi = eff_size, 
        vi = eff_size_var,
        method = "REML",
        data = data)
}