###Grapher and set duration 

source(here("Code","Helper_Function_Code", "Grapher_vert_hor.R"))

grapher_SA <- function(x_df){
  
  tmp <- split(x_df, x_df$change)
  
  
  
  ###This creates the boundary lines 
  grapher_lines_1 <- cbind(
                      grapher_mortality_boundary(tmp [[1]]),
                      change = unique(tmp[[1]]$change))
  grapher_lines_2<- cbind(
                      grapher_mortality_boundary(tmp [[2]]), 
                      change=  unique(tmp[[2]]$change))
    
  grapher_f <- rbind.data.frame(grapher_lines_1, grapher_lines_2)
  
  ###This finds the optimal points for the SA 
  optimal_points <- do.call(rbind,
                            lapply(tmp, function(x) 
                             x[which.max(x$end_fitness),
                               c("B_V","C_V", 'change')]))
  
  ###This finds the optimal points for the SA 
  longest_points <- do.call(rbind,
                            lapply(tmp, function(x) 
                            x[which.max(x$endtime),
                            c("B_V","C_V", 'change')]))
  
  ### The figure maker
SA_GG_Fitness <-
  ggplot(
    data = subset(
      x_df,
      x_df$status != "Fail"
    ),
    aes(x = B_V, y = C_V)
  ) +
  geom_raster(aes(fill = end_fitness)) +
  geom_raster(
    data = subset(
      x_df,
      x_df$status ==
        "Fail"
    ), aes(x = B_V, y = C_V),
    fill = "#d1dbe8"
  ) +
  geom_point(data = optimal_points, aes(
    x = B_V, y = C_V,
  ), 
  col = "#0037ff",
  size = 2, 
  shape = 5, stroke = 1.1) +
  geom_segment(
    data = grapher_f,
    aes(
      x = x, xend = xend,
      y = y, yend = yend,
      color = id
    ), size = 1.1,
    lineend = "round"
  ) +
  facet_wrap(~change) +
  scale_fill_viridis(name = "Cumulative \ntransmission potential",
  option = "magma") +
  scale_color_manual(values = c("fail" = 'black', 
                                'mort' = '#b8fcd5'),
                     guide = 'none')+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Transmission investment (c)") +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "right",
    strip.background = element_blank()
  ) +
  annotate("text",
    x = 8, y = 0.93, label = "Unestablished \ninfection",
    size = 5
  ) +
  annotate("text",
    x = 45, y = 0.1, label = "Host \nMortality",
    size = 5, color = "#b8fcd5"
  ) +
  annotate("point", x = 8, y = 0.54, size= 2, shape = 21,
           stroke = 1.1,
           col = "#0037ff" )


### The figure maker
SA_GG_Duration <-
  ggplot(
    data = subset(
      x_df,
      x_df$status != "Fail"
    ),
    aes(x = B_V, y = C_V)
  ) +
  geom_raster(aes(fill = endtime)) +
  geom_raster(
    data = subset(
      x_df,
      x_df$status ==
        "Fail"
    ), aes(x = B_V, y = C_V),
    fill = "#d1dbe8"
  ) +
  geom_point(data = longest_points, aes(
    x = B_V, y = C_V,
  ), 
  col = "#ff52a3",
  size = 2, 
  shape = 5, stroke = 1.1) +
  geom_segment(
    data = grapher_f,
    aes(
      x = x, xend = xend,
      y = y, yend = yend,
      color = id
    ), size = 1.1,
    lineend = "round"
  ) +
  facet_wrap(~change) +
  scale_fill_distiller(
    name = "Acute phase \n duration (days)",
    type = "seq",
    direction = -1,
    palette = "Greys",
    na.value = "black") +
  scale_color_manual(values = c("fail" = 'black', 
                                'mort' = '#b8fcd5'),
                     guide = 'none')+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Transmission investment (c)") +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "right",
    strip.background = element_blank()
  ) +
  annotate("text",
           x = 8, y = 0.93, label = "Unestablished \ninfection",
           size = 5
  ) +
  annotate("text",
           x = 45, y = 0.1, label = "Host \nMortality",
           size = 5, color = "#b8fcd5"
  ) +
  annotate("point", x = 8, y = 0.59, size= 2, shape = 21,
           stroke = 1.1,
           col = "#ff52a3" )

return(SA_GG_Fitness / SA_GG_Duration )

}
