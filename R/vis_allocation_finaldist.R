options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(icaeDesign))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(latex2exp))
source(here("R/visualize_funs.R"))
stat_used <- "mean"

fill_colors <- viridis(4)
color_colors <- viridis(4, direction = -1)
col_mapping_fill <- c(
  "Eigenvector"=fill_colors[1], "Closeness"=fill_colors[2], 
  "Degree"=fill_colors[3], "Random"=fill_colors[4]
)
col_mapping_col <- c(
  "Eigenvector"=color_colors[1], "Closeness"=color_colors[2], 
  "Degree"=color_colors[3], "Random"=color_colors[4]
)

# The effect of complexity allocation on average complexity -------------------

agg_tmax_data <- data.table::fread(
  here("output/allocation/allocation_tmax.csv"))
x_variable <- "prod_sp_val_alct"

title <- "The effect on average complexity"
y_title <- "Average complexity of produced products"
y_variable <- "comp_produced_prod_mean"
topology_complexity <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used) + 
  scale_fill_manual(values = col_mapping_fill) +
  scale_color_manual(values = col_mapping_col) +
  scale_x_discrete(
    labels=c("Eigenvector" = "Eigenvect.", "Closeness" = "Closen.")
    )

# The effect of complexity allocation on produced products --------------------

y_variable <- "share_produced_products"
title <- "The effect on produced products"
y_title <- "Share of of produced products"

topology_share_products <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 100, accuracy = 1)
    ) + 
  scale_fill_manual(values = col_mapping_fill) +
  scale_color_manual(values = col_mapping_col) +
  scale_x_discrete(
    labels=c("Eigenvector" = "Eigenvect.", "Closeness" = "Closen.")
  )

# The effect of complexity allocation on prices -------------------------------

y_variable <- "price_produced_prod_mean"
title <- "The effect on prices"
y_title <- "Mean of prices of produced products"

topology_price_mean <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used) + 
  scale_y_continuous(
    expand = expansion(mult = c(0,0), add = c(0, 0)), 
    limits = c(0, 320)
    ) +
  scale_fill_manual(values = col_mapping_fill) +
  scale_color_manual(values = col_mapping_col) +
  scale_x_discrete(
    labels=c("Eigenvector" = "Eigenvect.", "Closeness" = "Closen.")
  )

# Full plot- ------------------------------------------------------------------

file_name <- paste0("figures/allocation_dist_", stat_used,".pdf")

topology_total <- ggpubr::ggarrange(
  topology_share_products, 
  topology_price_mean, 
  topology_complexity,
  ncol = 3, nrow = 1, labels = paste0(LETTERS[1:3], ")"),
  font.label = list(size=18) , common.legend = T, legend = "bottom")

topology_total <- ggpubr::annotate_figure(
  topology_total, 
  top = ggpubr::text_grob(
    "Effect of the allocation of product complexity values", size = 24), 
  fig.lab.size = 18)

ggsave(plot = topology_total, 
       filename = here(file_name), 
       height = 7, width = 16)

paste0("Saved figure in: ", file_name)
