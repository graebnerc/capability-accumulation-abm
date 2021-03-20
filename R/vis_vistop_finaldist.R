options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(icaeDesign))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(latex2exp))

#' Creates a heatmap
make_heatmap <- function(data_used, x_variable, y_variable, 
                         color_variable, title_used){
  ggplot(
    data = data_used, 
    aes(x=!!as.name(x_variable), y=!!as.name(y_variable), 
        fill=!!as.name(color_variable))
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
      axis.text = element_text(size=12),
      plot.title = element_text(size=14),
      legend.title = element_blank())
}

x_var <- "max_info"
y_var <- "prod_space_struc"
color_vars <- c("share_produced_products", "price_produced_prod_mean", 
                "comp_produced_prod_mean")

agg_tmax_data <- data.table::fread(
  here("output/vistop/vistop_tmax.csv")) %>%
  filter(max_info>0.01) %>%
  group_by(!!as.name(x_var), !!as.name(y_var)) %>%
  summarise(across(.cols = all_of(color_vars), 
                   ~mean(.x, na.rm = F)), .groups="drop") %>%
  mutate(
    max_info = factor(max_info),
    prod_space_struc = ifelse(
      prod_space_struc=="BA_4", "B-A", ifelse(
        prod_space_struc=="random_0.25", 
        "Random",ifelse(
          prod_space_struc=="full",  "Complete", NA))),
    prod_space_struc = factor(
      prod_space_struc, 
      levels=c("Random", "B-A", "Complete"))
  )

title_plot <- "Share of produced products\n"
share <- make_heatmap(agg_tmax_data, x_var, y_var, color_vars[1], title_plot)

title_plot <- "Prices of produced products\n"
price <- make_heatmap(agg_tmax_data, x_var, y_var, color_vars[2], title_plot)

title_plot <- "Complexity of produced products\n"
comp <- make_heatmap(agg_tmax_data, x_var, y_var, color_vars[3], title_plot)

full_map <- ggarrange(
  share, price, comp, ncol = 3, labels = paste0(LETTERS[1:3], ")"))

full_map <- ggpubr::annotate_figure(
  full_map, 
  top = ggpubr::text_grob(TeX(
    "Joint effect of topology and the maximum range of vision $\\Upsilon_{max}$"
    ), size = 16)
  )

file_name <- "text/figures/vistop_heatmap.pdf"

ggsave(filename = here(file_name), plot = full_map, width = 12, height = 5)

paste0("Figure saved in: ", file_name)
