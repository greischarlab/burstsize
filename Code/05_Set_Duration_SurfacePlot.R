##################
###Set durations##
##################

###The set duration is a supplementary analysis to show that
###our way of choosing the end of the acute phase works.


############
###DAY 10###
############

Fitness_MODEL_PC_10<- lapply(FULL_MODEL_PC, function(x) subset(as.data.frame(x),
                                                               as.data.frame(x)$time <= 10))

Fitness_Cut_PC_10_cp <- data.frame(fitness=do.call(rbind,
                                                   lapply(Fitness_MODEL_PC_10 ,
                                                          Gametocyte_Fitness,"PC")))

###The maximum fitness is Burst size =50, transmission investment
###90%.


Fitness_Cut_PC_10_cp$B_V <- R_C_V$B_V
Fitness_Cut_PC_10_cp$C_V <- R_C_V$C_V 
Fitness_Cut_PC_10_cp$Death <- FULL_T_F_PC$TF.V2



GG_Fitness_Cut_PC_10_cp <- 
  ggplot(Fitness_Cut_PC_10_cp, aes(x=B_V, y= C_V))+
  geom_raster(aes(fill = fitness))+
  geom_raster(data=unestablished_index ,aes(x=B_V, y=C_V),
              fill='#DBDEE8')+
  geom_raster(data=mort_fitness,aes(x= B_V, y=C_V,
                                    fill= endfitness))+
  geom_segment(data=vert_seg_mort, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='#31f5a3',size=1.2,
               lineend='round')+
  geom_segment(data=hor_seg_mort, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='#31f5a3',size=1.2,
               lineend='round')+
  geom_segment(data=vert_seg, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='black',size=0.8,lineend='round')+
  geom_segment(data=hor_seg, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='black',size=0.8,
               lineend='round')+
  scale_fill_viridis(name="Cumulative fitness",option='inferno',limits=c(0,45))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab(expression("Burst size" ~beta))+
  ylab("Transmission investment (c)")+  ggtitle("A.Day 10")+
  theme(axis.text= element_text(size=14,color='black'),
        axis.title= element_text(size=14,color='black'))+
  annotate("point",x=50, y = 0.90,size=3,col='red')



###
###############
###DAY 20###
############

Fitness_MODEL_PC_20<- lapply(FULL_MODEL_PC, function(x) 
  subset(as.data.frame(x),
         as.data.frame(x)$time <= 20))

Fitness_Cut_PC_20_cp <- data.frame(fitness=do.call(rbind,
                                                   lapply(Fitness_MODEL_PC_20 ,
                                                          Gametocyte_Fitness,"PC")))


Fitness_Cut_PC_20_cp$B_V <- R_C_V$B_V
Fitness_Cut_PC_20_cp$C_V <- R_C_V$C_V 
Fitness_Cut_PC_20_cp$Death <- FULL_T_F_PC$TF.V2

###A burst size of 21 and transmission investment of 80% is the most
###optimal

GG_Fitness_Cut_PC_20_cp <-  ggplot(
  Fitness_Cut_PC_20_cp, 
  aes(x=B_V, y= C_V))+
  geom_raster(aes(fill =fitness))+
  geom_raster(data=unestablished_index ,aes(x=B_V, y=C_V),
              fill='#DBDEE8')+
  geom_raster(data=mort_fitness,aes(x= B_V, y=C_V, fill=
                                      endfitness))+
  geom_segment(data=vert_seg_mort, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='#31f5a3',size=1.2,
               lineend='round')+
  geom_segment(data=hor_seg_mort, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='#31f5a3',size=1.2,
               lineend='round')+
  geom_segment(data=vert_seg, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='black',size=0.8,lineend='round')+
  geom_segment(data=hor_seg, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='black',size=0.8,
               lineend='round')+
  scale_fill_viridis(name="Cumulative fitness",option='inferno',
                     limits=c(0,45))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab(expression("Burst size" ~beta))+
  ylab("Transmission investment (c)")+
  ggtitle("B. Day 20")+
  theme(axis.text= element_text(size=14,color='black'),
        axis.title= element_text(size=14,color='black'))+
  annotate("point",x=21, y = 0.80,size=3,col='red')


#########################################################
###Day 30- Not Accounting for the post-infection gametocytes#
#########################################################

Fitness_MODEL_PC_30<- lapply(FULL_MODEL_PC, function(x) subset(as.data.frame(x),
                                                               as.data.frame(x)$time <= 30))

Fitness_Cut_PC_30_cp <- data.frame(fitness=do.call(rbind,
                                                   lapply(Fitness_MODEL_PC_30 ,
                                                          Gametocyte_Fitness,"PC")))
###Optimal burst size is burst size 13 and and a transmission investmnet of 70%

Fitness_Cut_PC_30_cp$B_V <- R_C_V$B_V
Fitness_Cut_PC_30_cp$C_V <- R_C_V$C_V 
Fitness_Cut_PC_30_cp$Death <- FULL_T_F_PC$TF.V2

GG_Fitness_Cut_PC_30_cp<- 
  ggplot(
    Fitness_Cut_PC_30_cp, 
    aes(x=B_V, y= C_V))+
  geom_raster(aes(fill =fitness))+
  geom_raster(data=unestablished_index ,aes(x=B_V, y=C_V),
              fill='#DBDEE8')+
  geom_raster(data=mort_fitness,aes(x= B_V, y=C_V, fill=
                                      endfitness))+
  geom_segment(data=vert_seg_mort, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='#31f5a3',size=1.2,
               lineend='round')+
  geom_segment(data=hor_seg_mort, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='#31f5a3',size=1.2,
               lineend='round')+
  geom_segment(data=vert_seg, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='black',size=0.8,lineend='round')+
  geom_segment(data=hor_seg, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='black',size=0.8,
               lineend='round')+
  scale_fill_viridis(name="Cumulative fitness",option='inferno',
                     limits=c(0,45))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab(expression("Burst size" ~beta))+
  ylab("Transmission investment (c)")+
  ggtitle("C. Day 30")+
  theme(axis.text= element_text(size=14,color='black'),
        axis.title= element_text(size=14,color='black'))+
  annotate("point",x=13, y = 0.70,size=3,col='red')


#########################################################
###Day 50- Not Accounting for the post-infection gametocytes#
#########################################################

Fitness_MODEL_PC_50<- lapply(FULL_MODEL_PC, function(x) subset(as.data.frame(x),
                                                               as.data.frame(x)$time <= 50))

Fitness_Cut_PC_50_cp <- data.frame(fitness=do.call(rbind,
                                                   lapply(Fitness_MODEL_PC_50 ,
                                                          Gametocyte_Fitness,"PC")))


Fitness_Cut_PC_50_cp$B_V <- R_C_V$B_V
Fitness_Cut_PC_50_cp$C_V <- R_C_V$C_V 
Fitness_Cut_PC_50_cp$Death <- FULL_T_F_PC$TF.V2

###burst size =44 and transmission investment = 80%

GG_Fitness_Cut_PC_50_cp<- 
  ggplot(Fitness_Cut_PC_50_cp, 
         aes(x=B_V, y= C_V))+
  geom_raster(aes(fill =fitness))+
  geom_raster(data=unestablished_index ,aes(x=B_V, y=C_V),
              fill='#DBDEE8')+
  geom_raster(data=mort_fitness,aes(x= B_V, y=C_V, fill=
                                      endfitness))+
  scale_fill_viridis(name="Cumulative fitness",option='inferno',
                     limits=c(0,45))+
  geom_segment(data=vert_seg_mort, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='#31f5a3',size=1.2,
               lineend='round')+
  geom_segment(data=hor_seg_mort, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='#31f5a3',size=1.2,
               lineend='round')+
  geom_segment(data=vert_seg, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='black',size=0.8,lineend='round')+
  geom_segment(data=hor_seg, 
               aes(x=x, xend=xend,
                   y=y, yend=yend),color='black',size=0.8,
               lineend='round')+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab(expression("Burst size" ~beta))+
  ylab("Transmission investment (c)")+
  ggtitle("D. Day 50")+
  theme(axis.text= element_text(size=14,color='black'),
        axis.title= element_text(size=14,color='black'))+
  annotate("point",x=44, y = 0.8, size =3, col='red')

(GG_Fitness_Cut_PC_10_cp+GG_Fitness_Cut_PC_20_cp) /
  (GG_Fitness_Cut_PC_30_cp+GG_Fitness_Cut_PC_50_cp) +
  plot_layout(guides='collect')

ggsave(here("Figures","Manuscript_Figures",
            "Fitness_Surface_Plot_SD.pdf"), 
       width = 10.5, height= 8)