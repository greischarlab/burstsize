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

###

time_point_of_interest <- c(10,30,50,70,90)


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


GG_Fitness_Set_Duration <- 
  ggplot(Set_Duration_DF , 
         aes(x=B_V, y= C_V))+
  geom_raster(aes(fill = fitness))+
  geom_segment(data=
                 Grapher_Mortality_DF, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='#31f5a3',size=1.2,
               lineend='round')+
  facet_wrap(~end_time) + 
  scale_fill_viridis(option='inferno')+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab(expression("Burst size" ~beta))+
  ylab("Transmission investment (c)") 



ggsave(here("Figures","Manuscript_Figures",
            "Fitness_Surface_Plot_SD.pdf"), 
       width = 10.5, height= 8)
