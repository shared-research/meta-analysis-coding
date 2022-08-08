# Run Analysis

system("Rscript --encoding=UTF-8 scripts/0_pre_processing.R")
system("Rscript --encoding=UTF-8 scripts/1_analysis.R")
system("Rscript --encoding=UTF-8 scripts/2_tables_figures.R")
system("Rscript --encoding=UTF-8 scripts/3_supplementary.R")

# Update Method

compile_doc("docs/method/method_results.Rmd")

# Update Supplementary

compile_doc("docs/supplementary/supplementary.Rmd")


