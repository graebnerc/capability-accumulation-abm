options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(icaeDesign))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(latex2exp))
source(here("R/visualize_funs.R"))
stat_used <- "mean"

fill_colors <- viridis(6)
color_colors <- viridis(6, direction = -1)
col_mapping_fill <- c(
  "Random"=fill_colors[1], "PLC"=fill_colors[2], 
  "Ring"=fill_colors[3], "Barabasi-Albert"=fill_colors[4],
  "Regular"=fill_colors[5], "Complete"=fill_colors[6]
)
col_mapping_col <- c(
  "Random"=color_colors[1], "PLC"=color_colors[2], 
  "Ring"=color_colors[3], "Barabasi-Albert"=color_colors[4],
  "Regular"=color_colors[5], "Complete"=color_colors[6]
)

agg_tmax_data <- data.table::fread(
  here("output/topology/topology_tmax.csv"))

# The effect of topology on average complexity --------------------------------

title <- "Topology and average complexity"
y_title <- "Average complexity of produced products"
x_variable <- "prod_space_struc"
y_variable <- "comp_produced_prod_mean"
topology_complexity <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, 
  plot_kind = "violin", violin_width=1.1) +
  theme(
    axis.text.x = element_text(angle = 15)) + 
  scale_fill_manual(values = col_mapping_fill) +
  scale_color_manual(values = col_mapping_col) +
  guides(color = guide_legend(ncol = 6), fill = guide_legend(ncol = 6))

# The effect of topology on produced products ---------------------------------

x_variable <- "prod_space_struc"
y_variable <- "share_produced_products"
title <- "Topology and produced products"
y_title <- "Share of of produced products"

topology_share_products <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used) +
  theme(axis.text.x = element_text(angle = 15)) + 
  scale_y_continuous(
    labels = scales::percent_format(scale = 100, accuracy = 1)
    ) +
  scale_fill_manual(values = col_mapping_fill) +
  scale_color_manual(values = col_mapping_col) +
  guides(color = guide_legend(ncol = 6), fill = guide_legend(ncol = 6))

# The effect of topology on prices --------------------------------------------

x_variable <- "prod_space_struc"
y_variable <- "price_produced_prod_mean"
title <- "Topology and prices"
y_title <- "Mean of prices of produced products"

topology_price_mean <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used) +
  theme(axis.text.x = element_text(angle = 15)) + 
  scale_fill_manual(values = col_mapping_fill) +
  scale_color_manual(values = col_mapping_col) +
  guides(color = guide_legend(ncol = 6), fill = guide_legend(ncol = 6))


# Full plot -------------------------------------------------------------------

file_name <- paste0("text/figures/topology_dist_", stat_used, ".pdf")

topology_total <- ggpubr::ggarrange(
  topology_share_products, 
  topology_price_mean, 
  topology_complexity,
  ncol = 3, nrow = 1, 
  labels = paste0(LETTERS[1:3], ")"), 
  font.label = list(size=17), 
  common.legend = T, legend = "bottom")

ggsave(plot = topology_total, 
       filename = here(file_name), 
       height = 7, width = 16) 

paste0("Saved figure in: ", file_name)
