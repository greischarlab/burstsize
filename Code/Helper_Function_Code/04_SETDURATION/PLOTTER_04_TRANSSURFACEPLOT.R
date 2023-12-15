###This makes figure 6


Set_Duration_Plotter <- function(x,y,y_2){
  
  ggplot(subset(x,x$id !='Fail' & x$B_V < 51)) + 
  geom_raster(aes( x= time, y= B_V, fill = (cutprob))) +
  geom_line(data = y ,aes( x= lifespan, y=B_V), size = 1, 
             color ="orange")+
  geom_point(data = y_2 ,aes( x= lifespan, y=B_V), size = 3, 
              fill  ="orange", color = 'black', shape = 21)+
  scale_fill_viridis(option = 'mako', discrete = 'TRUE',
                     name = "Transmission \n Probability (f)",
                     guide = guide_coloursteps(even.steps = FALSE))+
  xlab("Days post-infection") +  
  ylab(expression(paste("Burst size", "( ", beta, ")"))) + 
  scale_x_continuous(expand=c(0,0), limits =c(0,51))+
  scale_y_continuous(expand=c(0,0))+ 
  geom_hline(yintercept = 25.3, color = 'red')+
  theme_classic()+
  theme(axis.text = element_text(size = 14,color = 'black'),
        axis.title = element_text(size = 15, color = 'black'),
        legend.position = 'right',
        legend.key.height = unit(1.5, "cm"))
}



