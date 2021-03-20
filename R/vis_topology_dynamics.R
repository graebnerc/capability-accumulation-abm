options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(icaeDesign))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(latex2exp))

# The dynamics plot for network topology--------------------------------------

update_dyn_plot <- function(x){
  x +
    scale_x_continuous(expand = c(0, 0)) +
    guides(color=guide_legend(ncol = 6)) +
    theme_icae() + 
    theme(
      axis.title.x = element_blank(),
      plot.title = element_text(size=20),
      axis.text.x = element_text(size = 16), 
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16), 
      legend.text = element_text(size=18)
      )
}

agg_dynamics_data <- data.table::fread(
  here("output/topology/topology_agg.csv"))%>%
  filter(t>1)

errorbar_width <- 8
errorbar_alpha <- 0.5
errorbar_steps <- 50

dyn_plot_shares <- ggplot(
  agg_dynamics_data, 
  aes(
    x=t, 
    y=share_produced_products_mean, 
    color=prod_space_struc)
  ) +
  ggtitle("Share of produced products") +
  geom_line(key_glyph = "rect", size=1.1) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  geom_errorbar(
    data = dplyr::filter(
      agg_dynamics_data, t %in% seq(1, 500, by = errorbar_steps)
      ),
    stat = "identity", aes(
      ymin=share_produced_products_mean-0.5*share_produced_products_std,
      ymax=share_produced_products_mean+0.5*share_produced_products_std),
    alpha=errorbar_alpha, width=errorbar_width) +
  ylab("Average share of produced products") +
  viridis::scale_color_viridis(
    discrete = T, option = "D",
    labels = c(
      "regular_4"=unname(TeX("Regular ($m=4$)")),
      "BA_4"=unname(TeX("Scale-free ($m=4$)")),
      "powerlaw_cluster_4_0.7"=unname(TeX("PLC ($m=4,\ p=.7$)")),
      "ring_2"=unname(TeX("Ring")),
      "random_0.25"=unname(TeX("Random ($p=.25$)")),
      "full"=unname(TeX("Complete"))
    )
  )
dyn_plot_shares <- update_dyn_plot(dyn_plot_shares)

dyn_plot_prices <- ggplot(
  agg_dynamics_data, aes(
    x=t, 
    y=price_produced_prod_mean_mean, 
    color=prod_space_struc)
  ) +
  ggtitle("Prices of produced products") +
  geom_line(key_glyph = "rect", size=1.1) + 
  geom_errorbar(
    data = dplyr::filter(
      agg_dynamics_data, t %in% seq(1, 500, by = errorbar_steps)),
    stat = "identity", aes(
      ymin=price_produced_prod_mean_mean-0.5*price_produced_prod_mean_std,
      ymax=price_produced_prod_mean_mean+0.5*price_produced_prod_mean_std),
    alpha=errorbar_alpha, width=errorbar_width) +
  ylab("Average prices of produced products") +
  viridis::scale_color_viridis(
    discrete = T, option = "D",
    labels = c(
      "regular_4"=unname(TeX("Regular ($m=4$)")),
      "BA_4"=unname(TeX("Scale-free ($m=4$)")),
      "powerlaw_cluster_4_0.7"=unname(TeX("PLC ($m=4,\ p=.7$)")),
      "ring_2"=unname(TeX("Ring")),
      "random_0.25"=unname(TeX("Random ($p=.25$)")),
      "full"=unname(TeX("Complete"))
    )
  )
dyn_plot_prices <- update_dyn_plot(dyn_plot_prices)

dyn_plot_complexity <- ggplot(
  agg_dynamics_data, 
  aes(
    x=t,
    y=comp_produced_prod_mean_mean, 
    color=prod_space_struc)
  ) +
  ggtitle("Complexity of produced products") +
  geom_line(key_glyph = "rect", size=1.1) + 
  geom_errorbar(
    data = dplyr::filter(
      agg_dynamics_data, t %in% seq(1, 500, by = errorbar_steps)),
    stat = "identity", aes(
      ymin=comp_produced_prod_mean_mean-0.5*comp_produced_prod_mean_std,
      ymax=comp_produced_prod_mean_mean+0.5*comp_produced_prod_mean_std),
    alpha=errorbar_alpha, width=errorbar_width) +
  ylab("Average complexity of produced products")  +
  viridis::scale_color_viridis(
    discrete = T, option = "D",
    labels = c(
      "regular_4"=unname(TeX("Regular ($m=4$)")),
      "BA_4"=unname(TeX("Scale-free ($m=4$)")),
      "powerlaw_cluster_4_0.7"=unname(TeX("PLC ($m=4,\ p=.7$)")),
      "ring_2"=unname(TeX("Ring")),
      "random_0.25"=unname(TeX("Random ($p=.25$)")),
      "full"=unname(TeX("Complete"))
    )
  )
dyn_plot_complexity <- update_dyn_plot(dyn_plot_complexity) +
  scale_x_continuous(limits = c(1, 512), expand = expansion())

# Final plot-------------------------------------------------------------------

dyn_plot <- ggarrange(
  dyn_plot_shares, dyn_plot_prices, dyn_plot_complexity, 
  ncol = 3, labels = paste0(LETTERS, ")"), font.label = list(size=18),
  legend = "bottom", common.legend = T)

file_name <- "text/figures/topology_dynamics.pdf"

ggsave(plot = dyn_plot, 
       filename = here(file_name), 
       height = 6, width = 16) 

paste0("Saved figure in: ", file_name)
