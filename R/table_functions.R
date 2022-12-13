# multivariate_table ------------------------------------------------------

#' Create a table using \code{flextable} for a fitted multivariate model
#' 
#' @param data_row A tibble with a row of the nested tibble with all fitted models indicating
#' a specific combination of correlations values.
#' @param digits An integer indicating how many digits to include after the comma. Default to \code{3}.
#' @return A \code{flextable}object
#' @export

multivariate_table <- function(data, col, digits = 3, fit){
    data %>% 
        filter(type == "summary") %>% 
        p_value_to_string(p.value) %>% 
        round_numeric_cols(digits = digits) %>% 
        # create 95% CI column
        mutate("95% CI" = paste0("[", conf.low, ", ", conf.high, "]")) %>% 
        select(term, estimate, std.error, "95% CI", statistic, p.value) %>% 
        flextable::flextable() %>% 
        flextable::set_header_labels(term = "Outcome",
                                     estimate = "d~ppc2~",
                                     std.error = "SE",
                                     conf.low = "95% CI",
                                     conf.high = "95% CI",
                                     statistic = fit$test,
                                     p.value = "p") %>% 
        flextable::theme_vanilla() %>% 
        flextable::autofit() %>% 
        flextable::merge_h(part = "header") %>% 
        flextable::align(part = "header", align = "center") %>% 
        flextable::align(j = c(2:6), part = "body", align = "center") %>% 
        flextable::italic(part = "header", j = 6) %>% 
        add_qm_stat(fit) %>% 
        add_cor_table(data) %>%
        fontsize(size = 10, part = "footer") %>% 
        ftExtra::colformat_md(part = "all")
}

# get_data -------------------------------------------------------------

#' Pull out from a list with a single element.
#' 
#' @param data the main dataframe with all variables
#' @param col the name of the columncolV to select
#' @param m_cor the pre-post (morris) correlation
#' @param b_cor the aggregation (borenstein) correlation
#' @param mu_cor the multivariate correlation
#' @return the dataframe given the input parameters

get_data <- function(data, col, 
                     m_cor = 0.7, 
                     b_cor = 0.5, 
                     mu_cor = 0.5){
    
    col <- rlang::enexpr(col)
    data %>% 
        filter(morris_cor == m_cor, boren_cor == b_cor, multi_cor == mu_cor) %>% 
        select(morris_cor, boren_cor, multi_cor, !!col) %>% 
        unnest(!!col)
}

get_model <- function(data, col, 
                     m_cor = 0.7, 
                     b_cor = 0.5, 
                     mu_cor = 0.5){
    col <- rlang::enexpr(col)
    data %>% 
        filter(morris_cor == m_cor, boren_cor == b_cor, multi_cor == mu_cor) %>% 
        pull(!!col) %>% 
        `[[`(1)
}

# p_value_to_string -------------------------------------------------------

#' Convert the columns with p-values into a string column depending on the
#' value. If the p-value is lower than 0.001 return the string \code{< 0.001}
#' otherwise return the value as string.
#' 
#' @param data A tibble
#' @param p_col Name of the column with the p-value information
#' @param digits An integer indicating how many digits to include
#' @return A tibble

p_value_to_string <- function(data, p_col, digits = 3){
    p_col <- rlang::enexpr(p_col) # use lazy evaluation
    data %>% 
        mutate(!!p_col := round(!!p_col, digits),
               !!p_col := ifelse(!!p_col < 0.001, "< 0.001", !!p_col))
}

# round_numeric_cols------------------------------------------------------

#' Round all numeric columns of a dataframe/tibble with a specified rounding
#' parameter
#' @param data A tibble/dataframe
#' @param digits An integer indicating the amount of digits
#' @return A tibble

round_numeric_cols <- function(data, digits){
    data %>% 
        mutate(across(where(is.numeric), ~ round(.x, digits = digits)))
}


# add_cor_table -----------------------------------------------------------

#' Add a flextable::footnote to a \code{flextable} object with information about correlations
#' for the specific model.
#' @param tab A \code{flextable} object
#' @param data_row A tibble with a row of the nested tibble with all fitted models indicating
#' a specific combination of correlations values.
#' @return A \code{flextable} object

add_cor_table <- function(tab, data){
    tab %>% 
        add_footer_lines(value = as_paragraph("\u03c1", as_sub("pre-post"), " = ", data$morris_cor[1], ", ",
                                      "\u03c1", as_sub("agg"), " = ", data$boren_cor[1], ", ",
                                      "\u03c1", as_sub("multi"), " = ", data$multi_cor[1]))
}

get_qm_stat <- function(fit){
    qm_stat <- list(
        QM = round(fit$QM, 3),
        QMdf = fit$QMdf[!is.na(fit$QMdf)],
        QMp  = ifelse(fit$QMp < 0.001, "p < 0.001", paste("p = ", round(fit$QMp, 3)))
    )
    qm_stat$type <- ifelse(length(qm_stat$QMdf) < 2, "Chisq", "F")
    return(qm_stat)
}

add_qm_stat <- function(tab, fit){
    qm_stat <- get_qm_stat(fit)

    if(qm_stat$type == "Chisq"){
        tab %>% 
            add_footer_lines(value = as_paragraph(as_i("Omnibus Test"), "  ", 
                                                  "\u03C7", as_sub(qm_stat$QMdf), "=", 
                                                  "  ", qm_stat$QM, "  ", qm_stat$QMp))
    }else{
        tab %>% 
            add_footer_lines(values = as_paragraph(as_i("Omnibus Test"), "  ", 
                                                   "F", as_sub(paste(qm_stat$QMdf[1], qm_stat$QMdf[2])), 
                                                   "=", "  ", qm_stat$QM, "  ", qm_stat$QMp))
    }

}

add_par_header <- function(tab, par_list){
    for(i in seq_along(par_list)){
        tab <- tab %>% 
            compose(j = i, part = "header", value = par_list[[i]])
    }
    return(tab)
}

get_qe_stat <- function(fit){
    qe_stat <- list(
        QE = round(fit$QE, 3),
        QEdf = fit$k - fit$p,
        QEp  = ifelse(fit$QEp < 0.001, "p < 0.001", paste("p = ", round(fit$QEp, 3)))
    )
    return(qe_stat)
}

add_qe_stat <- function(tab, fit){
    qe_stat <- get_qe_stat(fit)
    tab %>% 
        flextable::add_footer_lines(values = as_paragraph(as_i("Residual Heterogeneity"), "  ", 
                                                          "Q", as_sub(qe_stat$QEdf), "=", 
                                                          "  ", qe_stat$QE, "  ", qe_stat$QEp))
}