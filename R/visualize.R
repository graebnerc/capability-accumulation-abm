library(tidyverse)
library(feather)
library(here)
library(icaeDesign)
library(ggpubr)
library(viridis)
library(latex2exp)
options(warn = -1)
stat_used <- "mean"

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

#' Create a point plot 
#' 
#' Plots a point summarizing the data using `point_function` and add a 
#'  pointrange for the IQR
make_point_plot <- function(base_data,x_variable,y_variable, title, y_title, 
                            point_function="median", order_y_var=TRUE, 
                            plot_kind="violin", violin_width=1.0){
  if (point_function=="median"){
    sum_fun <- "median_IQR"
  } else if (point_function=="mean"){
    sum_fun <- "mean_std"
  } else{
    stop("Wrong point function provided!")
  }
  if (order_y_var){
    base_data <- base_data %>%
      mutate(UQ(as.name(x_variable)) := fct_reorder(
        UQ(as.name(x_variable)), UQ(as.name(y_variable)), 
        .fun=match.fun(point_function) )
      )
  }
  if (plot_kind=="line"){
    final_plot <- base_data %>%
      ggplot(data = ., 
             mapping= aes(x=UQ(as.name(x_variable)), 
                          y=UQ(as.name(y_variable)))) +
      geom_linerange(stat = "summary",size=2, color="grey", #"#f0f0f0",# alpha=0.5,
                     fun.data = match.fun(sum_fun)
      ) +
      geom_point(stat = "summary", fun.data = match.fun(sum_fun), 
                 size=5, color="#0e1c53")
  } else if(plot_kind=="violin"){
    final_plot <- base_data %>%
      ggplot(data = ., 
             mapping= aes(x=UQ(as.name(x_variable)), 
                          y=UQ(as.name(y_variable)),
                          fill=UQ(as.name(x_variable)),
                          color=UQ(as.name(x_variable))
                          )
             ) +
      geom_violin(width = violin_width, color=NA) + 
      geom_boxplot(width=0.1, color="grey", alpha=0.2, show.legend = FALSE) +
      geom_point(stat = "summary", fun.data = match.fun(sum_fun), 
                 size=5, show.legend = FALSE) +
      scale_fill_viridis(discrete = TRUE) +
      scale_color_viridis(discrete = TRUE, direction = -1)
  }
  final_plot +
    ggtitle(title) +
    xlab("") + ylab(y_title) +
    theme_icae() +
    theme(plot.title = element_text(size=16),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 0, 
                                     hjust = 0.5, 
                                     vjust = 0.6),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14), 
          panel.grid.major.x = element_blank())
}

# Figure 2: the dynamics plot for network topology-----------------------------
agg_dynamics_data <- data.table::fread(
  here("output/topology/topology_agg.csv"))%>%
  filter(t>1)

errorbar_width <- 8
errorbar_alpha <- 0.5
errorbar_steps <- 50

dyn_plot_shares <- ggplot(agg_dynamics_data, 
                   aes(x=t,
                       y=share_produced_products_mean, 
                       color=prod_space_struc)) +
  ggtitle("Share of produced products") +
  geom_line(key_glyph = "rect", size=1.1) +
  geom_errorbar(
    data = dplyr::filter(agg_dynamics_data, t %in% seq(1, 500, by = errorbar_steps)),
    stat = "identity", aes(
    ymin=share_produced_products_mean-0.5*share_produced_products_std,
    ymax=share_produced_products_mean+0.5*share_produced_products_std),
    alpha=errorbar_alpha, width=errorbar_width) +
  ylab("Average share of produced products") +
  scale_x_continuous(expand = c(0, 0)) +
  guides(color=guide_legend(ncol = 6)) +
  theme_icae() + viridis::scale_color_viridis(
    discrete = T, option = "D",
      labels = c(
      "regular_4"=unname(TeX("Regular ($m=4$)")),
      "BA_4"=unname(TeX("Scale-free ($m=4$)")),
      "powerlaw_cluster_4_0.7"=unname(TeX("PLC ($m=4,\ p=.7$)")),
      "ring_2"=unname(TeX("Ring")),
      "random_0.25"=unname(TeX("Random ($p=.25$)")),
      "full"=unname(TeX("Complete"))
      )
  ) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=17),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), 
        legend.text = element_text(size=18))
dyn_plot_shares

dyn_plot_prices <- ggplot(agg_dynamics_data, 
                          aes(x=t,
                              y=price_produced_prod_mean_mean, 
                              color=prod_space_struc)) +
  ggtitle("Prices of produced products") +
  geom_line(key_glyph = "rect", size=1.1) + 
  geom_errorbar(
    data = dplyr::filter(agg_dynamics_data, t %in% seq(1, 500, by = errorbar_steps)),
    stat = "identity", aes(
      ymin=price_produced_prod_mean_mean-0.5*price_produced_prod_mean_std,
      ymax=price_produced_prod_mean_mean+0.5*price_produced_prod_mean_std),
    alpha=errorbar_alpha, width=errorbar_width) +
  ylab("Average prices of produced products") +
  scale_x_continuous(expand = c(0, 0)) +
  guides(color=guide_legend(ncol = 6)) +
  theme_icae() + viridis::scale_color_viridis(
    discrete = T, option = "D",
    labels = c(
      "regular_4"=unname(TeX("Regular ($m=4$)")),
      "BA_4"=unname(TeX("Scale-free ($m=4$)")),
      "powerlaw_cluster_4_0.7"=unname(TeX("PLC ($m=4,\ p=.7$)")),
      "ring_2"=unname(TeX("Ring")),
      "random_0.25"=unname(TeX("Random ($p=.25$)")),
      "full"=unname(TeX("Complete"))
    )
  ) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=17),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), 
        legend.text = element_text(size=18))

dyn_plot_complexity <- ggplot(agg_dynamics_data, 
                          aes(x=t,
                              y=comp_produced_prod_mean_mean, 
                              color=prod_space_struc)) +
  ggtitle("Complexity of produced products") +
  geom_line(key_glyph = "rect", size=1.1) + 
  geom_errorbar(
    data = dplyr::filter(agg_dynamics_data, t %in% seq(1, 500, by = errorbar_steps)),
    stat = "identity", aes(
      ymin=comp_produced_prod_mean_mean-0.5*comp_produced_prod_mean_std,
      ymax=comp_produced_prod_mean_mean+0.5*comp_produced_prod_mean_std),
    alpha=errorbar_alpha, width=errorbar_width) +
  ylab("Average complexity of produced products")  +
  scale_x_continuous(expand = c(0, 0)) +
  guides(color=guide_legend(ncol = 6)) +
  theme_icae() + viridis::scale_color_viridis(
    discrete = T, option = "D",
    labels = c(
      "regular_4"=unname(TeX("Regular ($m=4$)")),
      "BA_4"=unname(TeX("Scale-free ($m=4$)")),
      "powerlaw_cluster_4_0.7"=unname(TeX("PLC ($m=4,\ p=.7$)")),
      "ring_2"=unname(TeX("Ring")),
      "random_0.25"=unname(TeX("Random ($p=.25$)")),
      "full"=unname(TeX("Complete"))
    )
  ) +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=17),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14), 
        legend.text = element_text(size=18))

dyn_plot <- ggarrange(dyn_plot_shares, dyn_plot_prices, dyn_plot_complexity, 
                      ncol = 3, labels = paste0(LETTERS, ")"),
                      legend = "bottom", common.legend = T)
dyn_plot

ggsave(plot = dyn_plot, 
       filename = here("text/figures/topology_dynamics_R-rev.pdf"), 
       height = 7, width = 18) 

# The violin plots-------------------------------------------------------------

agg_tmax_data <- data.table::fread(
  here("output/topology/topology_tmax.csv"))

# Figure 3: Topology - Single plots -------------------------------------------

title <- "The effect of topology on average complexity"
y_title <- "Average complexity of produced products"
x_variable <- "prod_space_struc"
y_variable <- "comp_produced_prod_mean"
topology_complexity <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                       title, y_title, point_function=stat_used, 
                                       plot_kind = "violin", violin_width=1.1) +
  theme(axis.text.x = element_text(angle = 15))
topology_complexity

x_variable <- "prod_space_struc"
y_variable <- "output_all_prdcts"
title <- "The effect of topology on total output"
y_title <- "Total output"

topology_total_output <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                         title, y_title, point_function=stat_used) +
  theme(axis.text.x = element_text(angle = 15))


x_variable <- "prod_space_struc"
y_variable <- "share_produced_products"
title <- "The effect of topology on produced products"
y_title <- "Share of of produced products"

topology_share_products <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                           title, y_title, point_function=stat_used) +
  theme(axis.text.x = element_text(angle = 15))


x_variable <- "prod_space_struc"
y_variable <- "price_produced_prod_mean"
title <- "The effect of topology on prices"
y_title <- "Mean of prices of produced products"

topology_price_mean <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                       title, y_title, point_function=stat_used) +
  theme(axis.text.x = element_text(angle = 15))
topology_price_mean

x_variable <- "prod_space_struc"
y_variable <- "price_produced_prod_std"
title <- "The effect of topology on price dispersion"
y_title <- "Std of prices of produced products"

topology_price_std <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                      title, y_title, point_function=stat_used) +
  theme(axis.text.x = element_text(angle = 15))

# Figure 3: Topology - Full plot -------------------------------------------

topology_total <- ggpubr::ggarrange(#topology_total_output, 
  topology_share_products, 
  topology_price_mean, 
  topology_price_std, 
  topology_complexity,
  ncol = 2, nrow = 2, labels = "AUTO")
file_name <- paste0("text/figures/topology_total_", stat_used, "-rev.pdf")
ggsave(plot = topology_total, 
       filename = here(file_name), 
       height = 14, width = 16) # 14

file_name <- paste0("text/figures/topology_total_red_", stat_used, "-rev.pdf")
topology_total <- ggpubr::ggarrange(#topology_total_output, 
  topology_share_products, 
  topology_price_mean, 
  topology_complexity,
  ncol = 3, nrow = 1, labels = paste0(LETTERS[1:3], ")"))
ggsave(plot = topology_total, 
       filename = here(file_name), 
       height = 7, width = 16) # 14

# Figure X: Distribution - Single plots ---------------------------------------

agg_tmax_data <- data.table::fread(
  here("output/comp_dist/comp_dist_tmax.csv"))
x_variable <- "prod_sp_val_dist"

title <- "The effect of complexity distribution on average complexity"
y_title <- "Average complexity of produced products"
y_variable <- "comp_produced_prod_mean"
topology_complexity <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                       title, y_title, point_function=stat_used)

y_variable <- "output_all_prdcts"
title <- "The effect of complexity distribution on total output"
y_title <- "Total output"

topology_total_output <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                         title, y_title, point_function=stat_used)

y_variable <- "share_produced_products"
title <- "The effect of complexity distribution on produced products"
y_title <- "Share of of produced products"

topology_share_products <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                           title, y_title, point_function=stat_used)

y_variable <- "price_produced_prod_mean"
title <- "The effect of complexity distribution on prices"
y_title <- "Mean of prices of produced products"

topology_price_mean <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                       title, y_title, point_function=stat_used)

y_variable <- "price_produced_prod_std"
title <- "The effect of complexity distribution on price dispersion"
y_title <- "Std of prices of produced products"

topology_price_std <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                      title, y_title, point_function=stat_used)

# Figure X: Distribution - Full plot ------------------------------------------
file_name <- paste0("text/figures/distribution_total_", stat_used, "-rev.pdf")
topology_total <- ggpubr::ggarrange(
  topology_share_products, topology_complexity,
  ncol = 2, nrow = 1, labels = c("A)", "B)"))
ggsave(plot = topology_total, filename = here(file_name), 
       height = 7, width = 16) 

# Figure 5: Allocation - single plots -----------------------------------------

agg_tmax_data <- data.table::fread(
  here("output/allocation/allocation_tmax.csv"))
x_variable <- "prod_sp_val_alct"

title <- "The effect on average complexity"
y_title <- "Average complexity of produced products"
y_variable <- "comp_produced_prod_mean"
topology_complexity <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                       title, y_title, point_function=stat_used)

y_variable <- "output_all_prdcts"
title <- "The effect of complexity allocation on total output"
y_title <- "Total output"

topology_total_output <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                         title, y_title, point_function=stat_used)

y_variable <- "share_produced_products"
title <- "The effect on produced products"
y_title <- "Share of of produced products"

topology_share_products <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                           title, y_title, point_function=stat_used)

y_variable <- "price_produced_prod_mean"
title <- "The effect on prices"
y_title <- "Mean of prices of produced products"

topology_price_mean <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                       title, y_title, point_function=stat_used)

y_variable <- "price_produced_prod_std"
title <- "The effect on price dispersion"
y_title <- "Std of prices of produced products"

topology_price_std <- make_point_plot(agg_tmax_data, x_variable, y_variable, 
                                      title, y_title, point_function=stat_used)

# Figure 5: Allocation - full plot- -----------------------------------------
file_name <- paste0("text/figures/allocation_total_", stat_used, "-rev.pdf")
topology_total <- ggpubr::ggarrange(#topology_total_output, 
                                    topology_share_products, 
                                    topology_price_mean, 
                                    topology_price_std, 
                                    topology_complexity,
                                    ncol = 2, nrow = 2, labels = "AUTO")
ggsave(plot = topology_total, 
       filename = here(file_name), 
       height = 14, width = 16) # 14

file_name <- paste0("text/figures/allocation_total_red_", stat_used,"-rev.pdf")
topology_total <- ggpubr::ggarrange(#topology_total_output, 
  topology_share_products, 
  topology_price_mean, 
  topology_complexity,
  ncol = 3, nrow = 1, labels = paste0(LETTERS[1:3], ")"))
ggsave(plot = topology_total, 
       filename = here(file_name), 
       height = 7, width = 16) # 14

# Figure 6: the relevance of delta/alpha - single plots------------------------

agg_tmax_data <- data.table::fread(
  here("output/delta_small/delta_small_tmax.csv")) %>%
  mutate(delta_coefficient=as.factor(delta_coefficient))
x_variable <- "delta_coefficient"

title <- latex2exp::TeX("The effect of $\\alpha$ on average complexity")
y_title <- "Average complexity of produced products"
y_variable <- "comp_produced_prod_mean"
topology_complexity <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F)

y_variable <- "output_all_prdcts"
title <- latex2exp::TeX("The effect of  $\\alpha$ on total output")
y_title <- "Total output"

topology_total_output <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F)

y_variable <- "share_produced_products"
title <- latex2exp::TeX("The effect of  $\\alpha$ on produced products")
y_title <- "Share of of produced products"

topology_share_products <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F)

y_variable <- "price_produced_prod_mean"
title <- latex2exp::TeX("The effect of  $\\alpha$ on prices")
y_title <- "Mean of prices of produced products"

topology_price_mean <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F)

y_variable <- "price_produced_prod_std"
title <- latex2exp::TeX("The effect of  $\\alpha$ on price dispersion")
y_title <- "Std of prices of produced products"

topology_price_std <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F)

# Figure 6: the relevance of delta/alpha - full plot---------------------------
file_name <- paste0("text/figures/delta_total_", stat_used,"-rev.pdf")
topology_total <- ggpubr::ggarrange(#topology_total_output, 
                                    topology_share_products, 
                                    topology_price_mean, 
                                    topology_price_std, 
                                    topology_complexity,
                                    ncol = 2, nrow = 2, labels = "AUTO")
ggsave(plot = topology_total, 
       filename = here(file_name), 
       height = 14, width = 16) # 14

file_name <- paste0("text/figures/delta_total_red_", stat_used, "-rev.pdf")
topology_total <- ggpubr::ggarrange(#topology_total_output, 
  topology_share_products, 
  #topology_price_mean, 
  topology_complexity,
  ncol = 2, nrow = 1, labels = paste0(LETTERS[1:3], ")"))
ggsave(plot = topology_total, 
       filename = here(file_name), 
       height = 7, width = 16) # 14

# Figure Y: the role of the range of vision - single plots---------------------
agg_tmax_data <- data.table::fread(
  here("output/range_vis/range_vis_tmax.csv")) %>%
  mutate(max_info=as.factor(max_info))
x_variable <- "max_info"

title <- latex2exp::TeX("Vision and average complexity")
y_title <- "Average complexity of produced products"
y_variable <- "comp_produced_prod_mean"
rangevis_complexity <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F)

y_variable <- "output_all_prdcts"
title <- latex2exp::TeX("Vision and total output")
y_title <- "Total output"

rangevis_total_output <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F)

y_variable <- "share_produced_products"
title <- latex2exp::TeX("Vision and produced products")
y_title <- "Share of of produced products"

rangevis_share_products <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F)

y_variable <- "price_produced_prod_mean"
title <- latex2exp::TeX("Vision and prices")
y_title <- "Mean of prices of produced products"

rangevis_price_mean <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F)

y_variable <- "price_produced_prod_std"
title <- latex2exp::TeX("Vision and price dispersion")
y_title <- "Std of prices of produced products"

rangevis_price_std <- make_point_plot(
  agg_tmax_data, x_variable, y_variable, 
  title, y_title, point_function=stat_used, order_y_var = F)

# Figure Y: the role of the range of vision - full plot------------------------
file_name <- paste0("text/figures/range_vis", stat_used,"-rev.pdf")
rangevis_total <- ggpubr::ggarrange(#rangevis_total_output, 
  rangevis_share_products, 
  rangevis_price_mean, 
  rangevis_price_std, 
  rangevis_complexity,
  ncol = 2, nrow = 2, labels = "AUTO")
ggsave(plot = rangevis_total, 
       filename = here(file_name), 
       height = 14, width = 16) # 14

file_name <- paste0("text/figures/range_vis_red_", stat_used, "-rev.pdf")
rangevis_total <- ggpubr::ggarrange(#rangevis_total_output, 
  rangevis_share_products, 
  #rangevis_price_mean, 
  rangevis_complexity,
  ncol = 2, nrow = 1, labels = paste0(LETTERS[1:3], ")"))
ggsave(plot = rangevis_total, 
       filename = here(file_name), 
       height = 7, width = 16) # 14

# Complexity and range of vision combined--------------------------------------

file_name <- paste0("text/figures/alpha_range_vis_red_", stat_used, "-rev.pdf")
topology_total <- ggpubr::ggarrange(
  ggpubr::ggarrange(
    topology_share_products, topology_complexity, 
    ncol = 2, nrow = 1, labels = paste0(LETTERS[1:2], ")"), 
    common.legend = T, legend = "bottom"
  ),
  ggpubr::ggarrange(
    rangevis_share_products + guides(fill=guide_legend(ncol = 6)), 
    rangevis_complexity + guides(fill=guide_legend(ncol = 6)), 
    ncol = 2, nrow = 1, labels = paste0(LETTERS[3:4], ")"), 
    common.legend = T, legend = "bottom"
    ),
  ncol = 2, nrow = 1)
ggsave(plot = topology_total, 
       filename = here(file_name), 
       height = 7, width = 18) # 14
