###Grapher and set duration 

source(here("Code","Helper_Function_Code", "Grapher_vert_hor.R"))

grapher_SA <- function(x_df, opt_df,long_df){
  
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
     subset(
       x_df,
       x_df$status !=
         "Fail"
     ),
     aes( x = B_V, y = C_V)
   ) +
   geom_raster(aes(fill = end_fitness)) +
   geom_raster(
     data = subset(
       x_df,
       x_df$status ==
         "Fail"
     ),
     aes(x = B_V, y = C_V), fill = "#d1dbe8"
   ) +
   geom_segment(
     data = grapher_f,
     aes(
       x = x, xend = xend,
       y = y, yend = yend,
       color = id
     ), size = 1.1,
     lineend = "round"
   ) +
    facet_wrap(~change)+
   scale_color_manual(
     values = c(
       "fail" = "black",
       "mort" = "#b8fcd5"
     ),
     guide = "none"
   ) +
   scale_fill_viridis(
     name = "Cumulative transmission \npotential",
     option = "magma"
   ) +
   scale_x_continuous(expand = c(0, 0)) +
   scale_y_continuous(expand = c(0, 0)) +
   xlab(expression(paste("Burst size", "( ", beta, ")"))) +
   ylab("Transmission investment (c)") +
   theme(
     text = element_text(size = 14),
     axis.text = element_text(size = 14, color = "black"),
     axis.title = element_text(size = 14, color = "black"),
     legend.position = "top"
   ) +
   geom_point(data = optimal_points ,
              aes( x = B_V, y= C_V),
              color = '#FF116B', size = 2)+
   geom_point(data = opt_df,
              aes(x = B_V, y = C_V),
              color = '#FF116B', size = 2, shape = 5)+

   annotate("text",
            x = 8, y = 0.93, label = "Unestablished \ninfection",
            size = 5
   ) +
   annotate("text",
            x = 45, y = 0.1, label = "Host \nMortality",
            size = 5, color = "#b8fcd5"
   )
 
### The figure maker
SA_GG_Duration <-
  ggplot(
    subset(
      x_df,
      x_df$status !=
        "Fail"
    ),
    aes( x = B_V, y = C_V)
  ) +
  geom_raster(aes(fill = endtime)) +
  geom_raster(
    data = subset(
      x_df,
      x_df$status ==
        "Fail"
    ),
    aes(x = B_V, y = C_V), fill = "#d1dbe8"
  ) +
  geom_segment(
    data = grapher_f,
    aes(
      x = x, xend = xend,
      y = y, yend = yend,
      color = id
    ), size = 1.1,
    lineend = "round"
  ) +
  facet_wrap(~change)+
  scale_color_manual(
    values = c(
      "fail" = "black",
      "mort" = "#b8fcd5"
    ),
    guide = "none"
  ) +
  scale_fill_viridis(
    name = "Acute \n phase (days)",
    option = "viridis"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Transmission investment (c)") +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "top"
  ) +
  geom_point(data =   longest_points ,
             aes( x = B_V, y= C_V),
             color = 'white', size = 2)+
  geom_point(data =  long_df,
             aes(x = B_V, y = C_V),
             color = 'white', size = 2, shape = 5)+

  annotate("text",
           x = 8, y = 0.93, label = "Unestablished \ninfection",
           size = 5
  ) +
  annotate("text",
           x = 45, y = 0.1, label = "Host \nMortality",
           size = 5, color = "#b8fcd5"
  )

return(SA_GG_Fitness / SA_GG_Duration )

}
