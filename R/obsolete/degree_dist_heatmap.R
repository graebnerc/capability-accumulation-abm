#' Creates a heatmap
make_heatmap <- function(data_used, x_variable, y_variable, 
                         color_variable, title_used){
  ggplot(
    data = data_used, 
    aes(
      x=!!as.name(x_variable), 
      y=!!as.name(y_variable), 
      fill=!!as.name(color_variable)
    )
  ) +
    geom_tile() +
    labs(title = title_used) +
    scale_y_discrete(expand = expansion()) +
    scale_x_discrete(expand = expansion()) +
    scale_fill_viridis() +
    theme_icae() +
    theme(
      axis.title = element_blank(), 
      panel.border = element_blank(), 
      axis.line = element_blank(),
      legend.title = element_blank())
}

x_var <- "prod_sp_val_alct"
y_var <- "prod_space_struc"
color_vars <- c(
  "share_produced_products", 
  "price_produced_prod_mean", 
  "comp_produced_prod_mean"
)

agg_tmax_data <- data.table::fread(
  here("output/topdist/topdist_tmax.csv")) %>%
  group_by(!!as.name(x_var), !!as.name(y_var)) %>%
  summarise(
    across(.cols = all_of(color_vars), 
           ~mean(.x, na.rm = F)), .groups="drop"
  ) %>%
  mutate(
    prod_sp_val_alct = ifelse(
      prod_sp_val_alct=="degree_centrality", 
      "Degree centrality", "Random allocation"),
    prod_space_struc = ifelse(
      prod_space_struc=="BA_4", "Barabasi-Albert", ifelse(
        prod_space_struc=="random_0.25", 
        "Random network",ifelse(
          prod_space_struc=="full",  "Complete network", NA))),
    prod_space_struc = factor(
      prod_space_struc, 
      levels=c("Random network", "Barabasi-Albert", "Complete network"))
  )

title_plot <- "Impact on share of produced products\n"
share <- make_heatmap(agg_tmax_data, x_var, y_var, color_vars[1], title_plot)

title_plot <- "Impact on prices of produced products\n"
price <- make_heatmap(agg_tmax_data, x_var, y_var, color_vars[2], title_plot)

title_plot <- "Impact on complexity of produced products\n"
comp <- make_heatmap(agg_tmax_data, x_var, y_var, color_vars[3], title_plot)

full_map <- ggarrange(
  share, price, comp, ncol = 3, labels = paste0(LETTERS[1:3], ")"))

full_map <- ggpubr::annotate_figure(
  full_map, top = "Joint effect of topology and allocation of complexity")

file_name <- "text/figures/topdist_heatmap.pdf"

ggsave(filename = here(file_name), plot = full_map, width = 12, height = 5)

paste0("Figure saved in: ", file_name)


# Create alternative with violin plots-----------------------------------------

#' Computes median and IQR for a data set
median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}
