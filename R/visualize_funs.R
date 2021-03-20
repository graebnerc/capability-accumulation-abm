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
make_point_plot <- function(
  base_data,x_variable, y_variable, title, y_title, 
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
      mutate(
        UQ(as.name(x_variable)) := fct_reorder(
          UQ(as.name(x_variable)), UQ(as.name(y_variable)), 
          .fun=match.fun(point_function) )
        )
  }
  if (plot_kind=="line"){
    final_plot <- base_data %>%
      ggplot(data = ., 
             mapping= aes(
               x=UQ(as.name(x_variable)), 
               y=UQ(as.name(y_variable)))
             ) +
      geom_linerange(
        stat = "summary",size=2, color="grey",
        fun.data = match.fun(sum_fun)
      ) +
      geom_point(stat = "summary", fun.data = match.fun(sum_fun), 
                 size=5, color="#0e1c53")
  } else if(plot_kind=="violin"){
    final_plot <- base_data %>%
      ggplot(data = ., 
             mapping= aes(
               x=UQ(as.name(x_variable)), 
               y=UQ(as.name(y_variable)),
               fill=UQ(as.name(x_variable)),
               color=UQ(as.name(x_variable))
               )
      ) +
      geom_violin(width = violin_width, color=NA) + 
      geom_boxplot(
        width=0.1, color="grey", alpha=0.2, show.legend = FALSE
        ) +
      geom_point(
        stat = "summary", fun.data = match.fun(sum_fun), 
        size=5, show.legend = FALSE) +
      scale_fill_viridis(discrete = TRUE) +
      scale_color_viridis(discrete = TRUE, direction = -1, option = "A")
  }
  final_plot +
    ggtitle(title) +
    xlab("") + ylab(y_title) +
    theme_icae() +
    theme(
      plot.title = element_text(size=20),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(
        angle = 0, hjust = 0.5, vjust = 0.6),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 18), 
      axis.title.x = element_text(size = 16), 
      legend.text = element_text(size = 19),
      panel.grid.major.x = element_blank())
}