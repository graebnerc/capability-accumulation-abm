options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(icaeDesign))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(latex2exp))

#' Computes mean and standard deviation for a data set
mean_std <- function(x) {
  data.frame(y = mean(x), # Median
             ymin = mean(x) - 0.5 * sd(x) ,
             ymax = mean(x) + 0.5 * sd(x)
  )  
}

make_group_violins <- function(
  base_data,x_variable,y_variable, fill_variable, title, y_title, 
  point_function="mean", order_y_var=FALSE, violin_width=1.0){
  if (point_function=="median"){
    sum_fun <- "median_IQR"
  } else if (point_function=="mean"){
    sum_fun <- "mean_std"
  } else{
    stop("Wrong point function provided!")
  }
  ggplot(data = base_data, 
         mapping= aes(
           x=UQ(as.name(x_variable)), 
           y=UQ(as.name(y_variable)),
           fill=UQ(as.name(fill_variable)),
           color=UQ(as.name(fill_variable))
         )
  ) +
    geom_violin(
      width = violin_width, color="grey", 
      draw_quantiles = c(0.25, 0.5, 0.75)
      ) + 
    scale_fill_viridis(discrete = TRUE) +
    scale_color_viridis(discrete = TRUE, direction = -1) +
    ggtitle(title) +
    xlab("") + ylab(y_title) +
    theme_icae() +
    theme(
      plot.title = element_text(size=14),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(
        angle = 0, hjust = 0.5, vjust = 0.6),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 13), 
      panel.grid.major.x = element_blank(), 
      legend.text = element_text(size = 14)
      )
}


agg_tmax_data <- data.table::fread(
  here("output/topdist/topdist_tmax.csv")) %>%
  mutate(
    prod_sp_val_alct = ifelse(
      prod_sp_val_alct=="degree_centrality", 
      "Degree centrality", "Random allocation"),
    prod_space_struc = ifelse(
      prod_space_struc=="BA_4", "Barabasi-Albert", ifelse(
        prod_space_struc=="random_0.25", 
        "Random",ifelse(
          prod_space_struc=="full",  "Complete", NA))),
    prod_space_struc = factor(
      prod_space_struc, 
      levels=c("Random", "Barabasi-Albert", "Complete"))
  )

title <- "The effect on average complexity"
y_title <- "Average complexity of produced products"
x_variable <- "prod_space_struc"
fill_variable <- "prod_sp_val_alct"
y_variable <- "comp_produced_prod_mean"

topdist_comp <- make_group_violins(
  agg_tmax_data, x_variable, y_variable, fill_variable,
  title, y_title, violin_width=1.1) +
  theme(axis.text.x = element_text(angle = 0))
topdist_comp

title <- "The effect on prices"
y_title <- "Average price of produced products"
x_variable <- "prod_space_struc"
fill_variable <- "prod_sp_val_alct"
y_variable <- "price_produced_prod_mean"

topdist_prices <- make_group_violins(
  agg_tmax_data, x_variable, y_variable, fill_variable,
  title, y_title, violin_width=1.1) +
  theme(axis.text.x = element_text(angle = 0))
topdist_prices

title <- "The effect on produced products"
y_title <- "Share of produced products"
x_variable <- "prod_space_struc"
fill_variable <- "prod_sp_val_alct"
y_variable <- "share_produced_products"

topdist_share <- make_group_violins(
  agg_tmax_data, x_variable, y_variable, fill_variable,
  title, y_title, violin_width=1.1) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 100, accuracy = 1)
    ) +
  theme(axis.text.x = element_text(angle = 0))
topdist_share

full_violin <- ggarrange(
  topdist_share, topdist_comp, topdist_prices, ncol = 3, 
  labels = paste0(LETTERS[1:3], ")"), common.legend = T, legend = "bottom")

full_violin <- ggpubr::annotate_figure(
  full_violin, 
  top = ggpubr::text_grob(
    "Joint effect of topology and allocation of complexity", size = 16), 
  fig.lab.size = 18)

file_name <- "text/figures/topdist_dist.pdf"

ggsave(filename = here(file_name), plot = full_violin, width = 12, height = 5)

paste0("Figure saved in: ", file_name)
