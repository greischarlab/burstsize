
###A helper function for the main script 04_Set_Duration_SurfacePlot.R

###This is the function that shows how the optimal burst 
###size changes with the infection length, assuming that the 
###transmission investment is fixed.

###Input:
###x

###Output: 
###A ggplot

Optimalburstsize_SetDuration_Lineplotter <- function(x){
  gg_plot <-  ggplot(x , aes(x = lifespan, y= B_V)) + 
    geom_line() + geom_point()+
    xlab("Total infection length (days)") + 
    ylab("Optimal burst size") +
    geom_hline(yintercept = 24, color = '#39c49f', linetype = 2, size = 1.1)+
    theme_classic() + 
    theme(axis.text = element_text(size =14, color = 'black'),
          axis.title = element_text(size = 15, color = 'black'))
  
  plot(gg_plot)
}
