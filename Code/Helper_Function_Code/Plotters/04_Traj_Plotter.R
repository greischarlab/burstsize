
TRAJ_PLOT_GG <- function(x,y){


###Produces the mortality and non-establishing infection lines
horizontal_vert_df <- grapher_mortality_boundary(x)


 traj_plot_GG<-  ggplot(y, 
         aes(x = B_V, y = C_V, 
             color = as.factor(lifespan))) +
  geom_text(aes(x= B_V, y= C_V, label = lifespan))+
  geom_point(size = 3.0)+
  geom_tile(
    data = subset(
      x,
      x$status ==
        "Fail"
    ),
    aes(x = B_V, y = C_V), fill = "#d1dbe8",
    inherit.aes = FALSE
    
  ) +
  geom_tile(
    data = subset(
      x,
      x$status ==
        "mort"
    ),
    aes(x = B_V, y = C_V), fill = "black",
    inherit.aes = FALSE
    
  )+
  geom_path(size=1,lineend = "round", group =1) +
  
  scale_color_viridis(option='plasma',discrete = TRUE,)+
  new_scale_color() +
  
  geom_segment(
    data = horizontal_vert_df,
    aes(
      x = x, xend = xend,
      y = y, yend = yend,
      color = id), size = 0.5,
    
    lineend = "round",
    inherit.aes = FALSE
  )+
  scale_color_manual(
    values = c(
      "fail" = "black",
      "mort" = "#b8fcd5"
    ),
    guide = "none"
  )+
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Transmission investment (c)") +
  theme(panel.background=element_rect(fill ='#e3e0da'),
        panel.border = element_rect(color = 'black', fill = NA, size =1),
        panel.grid  = element_blank(),
        
        text = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        legend.position = "none"
  ) +
  coord_cartesian(xlim =c(0,50), ylim = c(0,1))




}
