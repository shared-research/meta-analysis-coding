library(magrittr)
library(dplyr)

# Testing Dppc2 using the example from Metafor
# https://www.metafor-project.org/doku.php/analyses:morris2008

## Dataframe for testing

datT <- data.frame(
    m_pre   = c(30.6, 23.5, 0.5, 53.4, 35.6),
    m_post  = c(38.5, 26.8, 0.7, 75.9, 36.0),
    sd_pre  = c(15.0, 3.1, 0.1, 14.5, 4.7),
    sd_post = c(11.6, 4.1, 0.1, 4.4, 4.6),
    ni      = c(20, 50, 9, 10, 14),
    ri      = c(0.47, 0.64, 0.77, 0.89, 0.44))

datC <- data.frame(
    m_pre   = c(23.1, 24.9, 0.6, 55.7, 34.8),
    m_post  = c(19.7, 25.3, 0.6, 60.7, 33.4),
    sd_pre  = c(13.8, 4.1, 0.2, 17.3, 3.1),
    sd_post = c(14.8, 3.3, 0.2, 17.9, 6.9),
    ni      = c(20, 42, 9, 11, 14),
    ri      = c(0.47, 0.64, 0.77, 0.89, 0.44))

colnames(datT) <- paste0(colnames(datT), "_t")
colnames(datC) <- paste0(colnames(datC), "_c")

dat <- cbind(datT, datC)

expect_res <- data.frame(yi = c(0.77, 0.80, 1.20, 1.05, 0.44),
                         vi = c(0.11, 0.04, 0.14, 0.07, 0.16))

my_fun_res <- dat %>% 
    mutate(yi = get_dpcc2(mt_pre = m_pre_t, mt_post = m_post_t, mc_pre = m_pre_c, mc_post = m_post_c,
                          st_pre = sd_pre_t, st_post = sd_post_t, sc_pre = sd_pre_c, sc_post = sd_post_c,
                          nt = ni_t, nc = ni_c),
           vi = get_dpcc2_var(yi, nt = ni_t, nc = ni_c, rho = ri_c)) %>% 
    select(yi, vi) %>% 
    round(., 2)

# Testing

test_that("The dpcc2 is calculated as Metafor", {
    expect_equal(my_fun_res$yi, expect_res$yi)
})

test_that("The dpcc2 variance is calculated as Metafor", {
    expect_equal(my_fun_res$vi, expect_res$vi)
})
