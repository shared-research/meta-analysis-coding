# success -----------------------------------------------------------------

success <- function(msg){
    cli::cli_alert_success(msg)
}

# compile_doc -------------------------------------------------------------

compile_doc <- function(file){
    suppressWarnings(suppressMessages(
        rmarkdown::render(file, 
                          output_format = c("bookdown::html_document2",
                                            "bookdown::pdf_document2"),
                          quiet = TRUE)))
    success(paste(basename(file), "compiled!"))
}

# get_rmd_format ----------------------------------------------------------

get_rmd_format <- function(){
    output <- knitr::opts_knit$get("rmarkdown.pandoc.to")
    format <- ifelse(grepl("html", output), "html", "pdf")
    return(format)
}

# conditionally -----------------------------------------------------------

# thanks to https://community.rstudio.com/t/conditional-pipelines/6076/2
# use the function factories approach https://adv-r.hadley.nz/function-factories.html
# create a function that can be conditionally evaluated. ..1 return the first element of
# the ... list

conditionally <- function(fun){
    function(..., execute) {
        if (execute) fun(...) else ..1
    }
}

# message -----------------------------------------------------------------

test <- function(res, name){
    if(res){
        cli::cli_alert_success(name)
    }else{
        cli::cli_alert_danger(name)
    }
}


# get_all_packages --------------------------------------------------------

get_all_packages <- function(exclude = "metaCoding"){
    exclude <- paste0(exclude, collapse = "|")
    pkgs <- invisible(renv::dependencies())
    pkgs <- unique(pkgs$Package)
    pkgs[!grepl(exclude, pkgs)]
}