##################
###Set durations##
##################

###Burst Size Versus Transmission Investment ###
B_V = seq(1, 50, 0.5) #Burst size
C_V = seq(.01, 1, 0.01) #Transmission investment 
B_V_C_V <- expand.grid(B_V = B_V, C_V = C_V) #Different combination 

###I am able to figure out which of the burst and transmission 
###are unable to establish the infection
p_val = 2.5e-6
mu_M = 48
R_val = 8500000

RM_limit_1 <- ((p_val * R_val) + mu_M)/(p_val * R_val)

B_V_C_V$Establish <- ifelse((1-B_V_C_V$C_V) * B_V_C_V$B_V >= RM_limit_1,
                            "Establish","Fail")

###The time point of interests 
time_point_of_interest <- c(10,20,30,40,50, 60)


Set_Duration_DF <- Calculate_Set_Duration_CP(FULL_MODEL_PC,
                          B_V_C_V,
                          time_point_of_interest)
        
Grapher_Mortality_DF <- NULL
for (k in seq(1,length(Set_Duration_DF))){
  Grapher_Mortality_DF[[k]]<-  
    cbind.data.frame(grapher_mortality_boundary(Set_Duration_DF[[k]]),
                                     end_time = time_point_of_interest[[k]])
}



Grapher_Mortality_DF<- do.call(rbind, Grapher_Mortality_DF)
Set_Duration_DF <- do.call(rbind, Set_Duration_DF)


###Finding the maximum fitness points for each of the time point of interest

Max_Point_DF <- do.call(rbind, by(Set_Duration_DF, Set_Duration_DF$end_time, function (x)
                x[which.max(x$fitness),], simplify = FALSE))

###Day Labeler
Day_Names <- c(
  `10` = "A. Day 10",
  `20` = "B. Day 20",
  `30` = "C. Day 30",
  `40` = "D. Day 40",
  `50` = "E. Day 50",
  `60` = "F. Day 60"
)

GG_Fitness_Set_Duration <- 
  ggplot(Set_Duration_DF,
         aes(x=B_V, y= C_V))+
  geom_raster(aes(fill = fitness))+
  geom_raster(data = subset(Set_Duration_DF, 
                            Set_Duration_DF$status == 'Fail'),
              aes(x = B_V, y = C_V), fill = '#d1dbe8')+
  geom_segment(data =
                 Grapher_Mortality_DF , 
               aes(x=x, xend=xend,
                   y=y, yend=yend,color=id),
                   size=1.2,lineend='round')+
  geom_point(data = Max_Point_DF, aes(x = B_V, y = C_V),
             size = 2, shape = 21, col = "#ffb5df",
             stroke = 1) + 
  facet_wrap(~end_time,  labeller = as_labeller(Day_Names)) + 
  scale_fill_viridis(option='inferno', name = "Cumulative \ntransmisison potential")+
  scale_color_manual(values = c("fail" = 'black', 
                                'mort' = '#b8fcd5'),
                     guide = 'none') +   
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  xlab(expression("Burst size" ~beta)) +
  ylab("Transmission investment (c)") + 
  annotate('text', x = 12.5, y = 0.92, label = 'Unestablished \ninfection',
             size = 5)+
  annotate('text', x = 43, y = 0.1, label = 'Host \nMortality',
             size = 5, color = '#b8fcd5') + 
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15), 
    legend.position = 'top')
  

GG_Fitness_Set_Duration 

ggsave(here("Figures","Raw",
            "Fitness_Surface_Plot_SD.pdf"), 
       width = 12, height= 9)
