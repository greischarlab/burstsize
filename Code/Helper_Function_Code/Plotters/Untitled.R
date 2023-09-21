RC_GG_Plotter <- function(x){
RC_GG <- 
  ggplot(data = x[[1]],
         aes(x = lifespan, y= RC, color = group, group = group))+
  geom_line(size = 1) + 
  geom_point(data = subset(x[[2]],
                           x[[2]]$group != 'min'),
             aes( x= lifespan, y = RC),
             size = 3 ) +
  scale_color_manual(values = c('fit' = 'black', 'max'='green'))+
  scale_x_continuous(limits=c(0,100)) +
  xlab("Total infection length") + 
  ylab(expression(paste("Replicative capacity (","(1-c)", beta,")")))+ 
  theme_classic() + 
  theme(
    axis.text = element_text(size = 14,color = 'black'),
    axis.title = element_text(size = 14),
    strip.background =  element_blank(),
    strip.text = element_text(size = 16),
    strip.text.y.right = element_text(angle = 0),
    panel.spacing = unit(1.5, "lines"), 
    legend.position = 'none')

print(RC_GG)

}
