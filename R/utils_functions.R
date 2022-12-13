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

# get_rmd_comments --------------------------------------------------------

get_rmd_comments <- function(file){
    file_lines <- readLines(file, warn = FALSE)
    has_comment <- stringr::str_detect(file_lines, "<!--")
    file_lines_comments <- file_lines[has_comment]
    comments <- stringr::str_match_all(file_lines_comments, "<!--\\s*(.*?)\\s*-->")
    comments <- unlist(sapply(comments, function(x) x[,2]))
    return(comments)
}

# write_rmd_comments ------------------------------------------------------

write_rmd_comments <- function(file){
    comments <- get_rmd_comments(file)
    cat(paste("-", comments), sep = "\n")
}

# get_rmd_format ----------------------------------------------------------

get_rmd_format <- function(){
    output <- knitr::opts_knit$get("rmarkdown.pandoc.to")
    format <- ifelse(grepl("html", output), "html", "pdf")
    return(format)
}

# Find Files --------------------------------------------------------------

get_all_files <- function(dir){
    all_files <- list.files(here::here(dir), recursive = T, no.. = FALSE)
    all_files[endsWith(all_files, ".R") | endsWith(all_files, ".Rmd")]
}

get_all_files_content <- function(dir){
    files <- get_all_files(dir)
    files_content <- lapply(files, readLines) %>% suppressWarnings()
    names(files_content) <- files
    return(files_content)
}

find_fun_within_file <- function(file, fun){
    sum(grepl(fun, file))
}

get_files_fun <- function(fun, dir = "."){
    files_content <- get_all_files_content(dir)
    files_id <- sapply(files_content, 
                       find_fun_within_file, 
                       fun)
    
    out <- names(files_content)[files_id > 0]
    out <- paste(out, "---", files_id[files_id > 0], "times")
    out_nice <- ifelse(grepl("R/", out), cli::col_blue(out), out)
    msg <- map(out_nice, cli::cli_alert_success)
}

get_bib_database <- function(destination, overwrite = FALSE){
    link <- "https://raw.githubusercontent.com/filippogambarota/bib-database/main/references.bib"
    if(!file.exists(destination) | overwrite){
        download.file(url = link, destfile = destination, quiet = TRUE)
        success("Bibtex database downloaded!")
    }else{
        cli::cli_alert_warning("The bibtex file already exist! Set overwrite = TRUE to download again!")
    }
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