options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(icaeDesign))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(latex2exp))
source(here("R/visualize_funs.R"))
stat_used <- "mean"

legend_title_size <- 17

# The effect of alpha on average complexity -----------------------------------

agg_tmax_data <- data.table::fread(
  here("output/delta_small/delta_small_tmax.csv")) %>%
  mutate(delta_coefficient=as.factor(delta_coefficient))
x_variable <- "delta_coefficient"

title <- latex2exp::TeX("$\\alpha$ and average complexity")
y_title <- "Average complexity of produced products"
y_variable <- "comp_produced_prod_mean"

alpha_complexity <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F) +
  labs(fill=TeX("Values for $\\alpha$:")) + 
  theme(legend.title = element_text(size = legend_title_size))

# The effect of alpha on produced products ------------------------------------

title <- latex2exp::TeX("$\\alpha$ and produced products")
y_title <- "Share of of produced products"
y_variable <- "share_produced_products"

alpha_share_products <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F) +
  labs(fill=TeX("Values for $\\alpha$:")) + 
  scale_y_continuous(
    labels = scales::percent_format(scale = 100, accuracy = 1),
    limits = c(0, 0.305)) +
  theme(legend.title = element_text(size = legend_title_size))

# The range of vision and average complexity-----------------------------------

agg_tmax_data <- data.table::fread(
  here("output/range_vis/range_vis_tmax.csv")) %>%
  mutate(max_info=as.factor(max_info))
x_variable <- "max_info"

title <- latex2exp::TeX("Vision and average complexity")
y_title <- "Average complexity of produced products"
y_variable <- "comp_produced_prod_mean"

rangevis_complexity <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F) +
  labs(fill=TeX("Values for $\\Upsilon$:")) + 
  theme(legend.title = element_text(size = legend_title_size))

# The range of vision and produced products------------------------------------

title <- latex2exp::TeX("Vision and produced products")
y_variable <- "share_produced_products"
y_title <- "Share of of produced products"

rangevis_share_products <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F) +
  labs(fill=TeX("Values for $\\Upsilon$:")) + 
  scale_y_continuous(
    labels = scales::percent_format(scale = 100, accuracy = 1)
  ) +
  theme(legend.title = element_text(size = legend_title_size))

# The effect of the saturation threshold on average complexity-----------------

agg_tmax_data <- data.table::fread(
  here("output/saturation/saturation_tmax.csv")) %>%
  mutate(nominal_demand=as.factor(nominal_demand))
x_variable <- "nominal_demand"

title <- latex2exp::TeX("Saturation threshold and average complexity")
y_title <- "Average complexity of produced products"
y_variable <- "comp_produced_prod_mean"

saturation_complexity <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F) +
  labs(fill=TeX("Values for $q^{max}$:")) + 
  theme(legend.title = element_text(size = legend_title_size))

# The effect of the saturation threshold on produced products------------------

title <- latex2exp::TeX("Saturation threshold and produced products")
y_variable <- "share_produced_products"
y_title <- "Share of of produced products"

saturation_share_products <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F) +
  labs(fill=TeX("Values for $q^{max}$:")) + 
  scale_y_continuous(
    labels = scales::percent_format(scale = 100, accuracy = 1)) +
  theme(legend.title = element_text(size = legend_title_size))

# The effect of the saturation threshold on prices-----------------------------

title <- latex2exp::TeX("Saturation threshold and prices")
y_variable <- "price_produced_prod_mean"
y_title <- "Average price of produced products"

saturation_prices <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F) +
  labs(fill=TeX("Values for $q^{max}$:")) + 
  theme(legend.title = element_text(size = legend_title_size))

# The full plot----------------------------------------------------------------

file_name <- paste0("text/figures/sensitivity_", stat_used, ".pdf")

sensitivity_total <- ggpubr::ggarrange(
  ggpubr::ggarrange(
    saturation_share_products, saturation_complexity, saturation_prices, 
    ncol = 3, nrow = 1, labels = paste0(LETTERS[1:3], ")"), 
    font.label = list(size=17),
    common.legend = T, legend = "bottom"),
  ggpubr::ggarrange(
    ggpubr::ggarrange(
      alpha_share_products, alpha_complexity, 
      ncol = 2, nrow = 1, labels = paste0(LETTERS[4:5], ")"), 
      font.label = list(size=17),
      common.legend = T, legend = "bottom"
    ),
    ggpubr::ggarrange(
      rangevis_share_products + guides(fill=guide_legend(ncol = 6)), 
      rangevis_complexity + guides(fill=guide_legend(ncol = 6)), 
      ncol = 2, nrow = 1, labels = paste0(LETTERS[6:7], ")"), 
      font.label = list(size=17),
      common.legend = T, legend = "bottom"
    ),
    ncol = 2, nrow = 1),
  ncol = 1, nrow = 2)

ggsave(plot = sensitivity_total, 
       filename = here(file_name), 
       height = 12, width = 18) 

paste0("Saved figure in: ", file_name)
