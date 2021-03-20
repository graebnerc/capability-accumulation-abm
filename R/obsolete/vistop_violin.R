
# Create alternative with violin plots-----------------------------------------

#' Computes median and IQR for a data set
median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}

#' Computes mean and standard deviation for a data set
mean_std <- function(x) {
  data.frame(y = mean(x), # Median
             ymin = mean(x) - 0.5 * sd(x) ,
             ymax = mean(x) + 0.5 * sd(x)
  )  
}

make_group_violins <- function(base_data,x_variable,y_variable, 
                               fill_variable, title, y_title, x_title,
                               point_function="mean", order_y_var=FALSE, 
                               violin_width=1.0){
  if (point_function=="median"){
    sum_fun <- "median_IQR"
  } else if (point_function=="mean"){
    sum_fun <- "mean_std"
  } else{
    stop("Wrong point function provided!")
  }
  ggplot(data = base_data, 
         mapping= aes(x=UQ(as.name(x_variable)), 
                      y=UQ(as.name(y_variable)),
                      fill=UQ(as.name(fill_variable)),
                      color=UQ(as.name(fill_variable))
         )
  ) +
    geom_violin(width = violin_width, color="black", 
                draw_quantiles = c(0.25, 0.5, 0.75)
    ) + 
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE, direction = -1) +
    labs(title = title, y=y_title, x=x_title) +
    coord_flip() + 
    theme_icae() +
    theme(plot.title = element_text(size=13),
          axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 0, 
                                     hjust = 0.5, 
                                     vjust = 0.6),
          axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 12), 
          panel.grid.major.x = element_blank())
}


agg_tmax_data <- data.table::fread(
  here("output/vistop/vistop_tmax.csv")) %>%
  mutate(max_info = factor(max_info)) %>%
  mutate(
    prod_space_struc = ifelse(
      prod_space_struc=="BA_4", "Barabasi-Albert", ifelse(
        prod_space_struc=="random_0.25", 
        "Random network",ifelse(
          prod_space_struc=="full",  "Complete network", NA)))
  )

fill_variable <- "prod_space_struc"
x_title <- TeX("Maximum range of vision $\\Upsilon$")
x_variable <- "max_info"

# The effect on complexity of produced products-------------------------------------

title <- "Effect on average complexity"
y_variable <- "comp_produced_prod_mean"
y_title <- "Average complexity of produced products"

topdist_comp <- make_group_violins(
  agg_tmax_data, x_variable, y_variable, fill_variable,
  title, y_title, x_title, violin_width=1.1) +
  theme(axis.text.x = element_text(angle = 0))
topdist_comp

# The effect on prices of produced products-------------------------------------

title <- "Effect on prices"
y_variable <- "price_produced_prod_mean"
y_title <- "Average prices of produced products"

topdist_prices <- make_group_violins(
  agg_tmax_data, x_variable, y_variable, fill_variable,
  title, y_title, x_title, violin_width=1.1) +
  theme(axis.text.x = element_text(angle = 0))
topdist_prices

# The effect on share of produced products-------------------------------------

title <- "Effect on produced products"
y_title <- "Share of produced products"
y_variable <- "share_produced_products"

topdist_share <- make_group_violins(
  agg_tmax_data, x_variable, y_variable, fill_variable,
  title, y_title, x_title, violin_width=1.1) +
  theme(axis.text.x = element_text(angle = 0)) +
  scale_y_continuous(limits = c(0, 0.8), labels = scales::percent_format()) 
topdist_share

# The full violin plot -------------------------------------

full_violin <- ggarrange(
  topdist_share, topdist_prices, topdist_comp, ncol = 3, 
  labels = paste0(LETTERS[1:3], ")"), common.legend = T, legend = "bottom")

full_violin <- ggpubr::annotate_figure(
  full_violin, 
  top = "Joint effect of topology and allocation of complexity\n", 
  fig.lab.size = 16)

file_name <- "text/figures/vistop_dist.pdf"

ggsave(filename = here(file_name), 
       plot = full_violin, width = 12, height = 5)

paste0("Figure saved in: ", file_name)
