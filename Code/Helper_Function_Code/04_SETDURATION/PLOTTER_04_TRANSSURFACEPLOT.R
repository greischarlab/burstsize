
Set_Duration_Plotter <- function(x,y){
  
  ggplot(subset(x,x$id !='Fail' & x$B_V < 38.6)) + 
  geom_raster(aes( x= time, y= B_V, fill = (cutprob))) +
  geom_point(data = y,aes( x= lifespan, y=B_V), size = 3, 
             color ="#FF116B")+
  scale_fill_viridis(option = 'mako', discrete = 'TRUE',
                     name = "Transmission \n Probability",
                     guide = guide_coloursteps(even.steps = FALSE))+
  xlab("Days post infection") + ylab("Burst size") + 
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), limits= c(14,39))+ 
  geom_vline(xintercept = c(5,30,50), color = 'black',
             linetype = 2, alpha = 0.8)+
  geom_hline(yintercept = 24.2, color = 'red')+
  theme_classic()+
  theme(axis.text = element_text(size = 14,color = 'black'),
        axis.title = element_text(size = 15, color = 'black'),
        legend.position = 'right',
        legend.key.height = unit(1.5, "cm"))
}



