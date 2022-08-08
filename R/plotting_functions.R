# remove_x_text -----------------------------------------------------------

remove_x_text <- function(){
    theme(axis.title.x = element_blank())
}

# remove_y_text -----------------------------------------------------------

remove_y_text <- function(){
    theme(axis.title.y = element_blank())
}

# prep_data_forest_multi --------------------------------------------------

prep_data_forest_multi <- function(dat_agg, dat_model){
    
    dat_forest <- dat_agg %>% 
        select(paper_id, outcome2, eff_size, eff_size_se)
    
    dat_forest_fit <- dat_model %>% 
        filter(type == "summary") %>% 
        mutate(paper_id = "Average") %>% 
        rename("outcome2" = term,
               "eff_size" = estimate,
               "eff_size_se" = std.error) %>% 
        select(paper_id, outcome2, eff_size, eff_size_se)
    
    dat_forest_multi <- bind_rows(dat_forest, dat_forest_fit)
    
    dat_forest_multi %>% 
        mutate(.lower = eff_size - eff_size_se*1.96,
               .upper = eff_size + eff_size_se*1.96,
               paper_id = factor(paper_id),
               paper_id = fct_relevel(paper_id, "Average")) %>% 
        clean_outcome_names(outcome2)
}

# forest_plot_multi -------------------------------------------------------

forest_plot_multi <- function(dat_forest_multi){
    
    # TODO try to draw diamond
    
    dat_forest_multi %>% 
        ggplot() +
        geom_hline(yintercept = 0.95, alpha = 0.2) +
        
        # 0 Effect
        
        geom_vline(xintercept = 0, 
                   alpha = 0.5,
                   size = 0.5, 
                   col = "black", 
                   linetype = "dashed") +
        
        # Average Effect
        
        # geom_rect(data = dat_forest_multi %>% filter(paper_id == "Average"),
        #           aes(xmin = .lower,
        #               xmax = .upper,
        #               ymax = 1.01,
        #               ymin = 0.90),
        #           fill = "firebrick1") +
        
        # All papers
        
        geom_pointrange(data = dat_forest_multi %>% filter(paper_id != "Average"),
                        aes(x = eff_size, 
                            y = factor(paper_id),
                            xmin = .lower,
                            xmax = .upper),
                        show.legend = FALSE,
                        size = 0.70,
                        shape = 15) +
        
        # Outcome
        
        facet_grid(~outcome2, labeller = label_wrap_gen(width = 4, multi_line = TRUE)) +
        
        # Tweaks
        
        scale_y_discrete(breaks = levels(dat_forest_multi$paper_id),
                         limits = c(levels(dat_forest_multi$paper_id)[1],
                                    levels(dat_forest_multi$paper_id)[-1])) +
        cowplot::theme_map() +
        remove_y_text() +
        remove_x_text() +
        theme(
            axis.text.y = element_text(
                face = ifelse(levels(dat_forest_multi$paper_id) == "Average", "bold", "plain"),
                size = ifelse(levels(dat_forest_multi$paper_id) == "Average", 15, 12),
                margin = margin(t = 0, r = 10, b = 0, l = 0),
                hjust = 0
            ),
            strip.text = element_text(
                face = "bold"
            ),
            axis.text.x = element_text(size = 13)
        ) +
        cowplot::panel_border(size = 0.5, linetype = 1, color = "grey")
}

# create_cor_labels -------------------------------------------------------

create_cor_labels <- function(data, col, what = c("pre-post", "agg", "multi")){
    col <- rlang::enexpr(col)
    col_s <- as.character(col)
    what <- match.arg(what)
    latex2exp::TeX(sprintf("$\\rho_{%s} = %s$",
                what,
                unique(data[[col_s]])))
}


# remove_element_ggplot ---------------------------------------------------

# thanks to https://stackoverflow.com/a/52028670

remove_element_ggplot <- function(plot, to_remove){
    p_tab <- ggplotGrob(plot)
    keep <- !p_tab$layout$name %in% to_remove
    p_tab$layout <- p_tab$layout[keep, , drop = FALSE]
    p_tab$grobs <- p_tab$grobs[keep]
    grid::grid.newpage()
    grid::grid.draw(p_tab)
}

# ggplotgrob_save ---------------------------------------------------------

ggplotgrob_save <- function(plot, 
                            file, 
                            device = c("png", "svg", "tiff", "pdf"), 
                            width, 
                            height){
    
    device <- match.arg(device)
    
    dev <- switch(device,
                  "png" = png,
                  "svg" = svg,
                  "tiff" = tiff,
                  "pdf" = pdf
    )
    
    dev(file, width = width, height = height)
    plot
    dev.off()
}

# sensitivity_plot --------------------------------------------------------

sensitivity_plot <- function(data){
    data %>% 
        prep_data_sens() %>% 
        ggplot(aes(x = factor(morris_cor))) +
        ggh4x::facet_nested(model_type + meta_type ~ multi_cor_parse + boren_cor_parse,
                            labeller = label_parsed) +
        cowplot::theme_minimal_grid(font_size = 12) +
        cowplot::panel_border() +
        ggtitle(data$outcome[1]) +
        xlab(latex2exp::TeX("$\\rho_{pre-post}$"))
    
}

# prep_data_sens ----------------------------------------------------------

prep_data_sens <- function(data){
    data %>% 
        add_paper_cond() %>% 
        mutate(sign = ifelse(p.value < 0.05, "yes", "no"))
}

# sensitivity_plot_estimate -----------------------------------------------

sensitivity_plot_estimate <- function(data){
    # symbol for plotting
    sign <- 24
    not_sign <- 21
    
    data %>% 
        sensitivity_plot() +
        geom_hline(yintercept = 0,
                   #linetype = "dashed",
                   col = "red",
                   size = 0.5,
                   alpha = 0.3) +
        geom_segment(aes(y = conf.low,
                         yend = conf.high,
                         xend = factor(morris_cor)),
                     size = 0.5) +
        geom_point(aes(y = estimate,
                       fill = paper_cor,
                       shape = sign),
                   color = "black",
                   size = 2,
                   show.legend = FALSE) +
        ylab(latex2exp::TeX("$d_{ppc2}$")) +
        {
           if(length(unique(data$sign)) > 1){
               scale_shape_manual(values = c(not_sign,sign))
           }else{
               scale_shape_manual(values = sign)
           }
        } +
        #scale_shape_manual(values = c(21,24)) +
        scale_fill_manual(values = c("darkgrey", "firebrick1"))
    
    
}

# sensitivity_plot_pvalue -------------------------------------------------

sensitivity_plot_pvalue <- function(data){
    data %>% 
        sensitivity_plot() +
        geom_line(aes(y = p.value, group = 0),
                  size = 0.5) +
        geom_hline(yintercept = 0.05,
                   #linetype = "dashed",
                   col = "red",
                   size = 0.5,
                   alpha = 0.3) +
        geom_point(aes(y = p.value,
                       size = 2,
                       color = interaction(sign, paper_cor),
                       shape = paper_cor),
                   show.legend = FALSE) +
        ylab("p value") +
        scale_color_manual(values = c("black", "lightseagreen", "firebrick1"))
    
}