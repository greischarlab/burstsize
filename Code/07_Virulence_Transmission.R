########################################
###RED BLOOD CELLS/GAMETOCYTES/CUMULATIVE ###
#############################################
########################################################
###Understanding how the viable burst size combinations#
###influence virulence, transmission, and such.     ####
########################################################

###This is for the 


RBC_Amount <-    do.call(rbind,
                         lapply(FULL_MODEL_PC, 
                                function(x) 
                                  min(as.data.frame(x[,"R"]))))


G_Amount <-    do.call(rbind,
                       lapply(FULL_MODEL_PC, 
                              function(x) 
                                max(as.data.frame(x[,"G"]))))

RBC_Amount_F <- as.data.frame(RBC_Amount)
G_Amount_F <- as.data.frame(G_Amount )

RBC_Amount_F$B_V <- R_C_V$B_V
RBC_Amount_F$C_V <- R_C_V$C_V
RBC_Amount_F$Death <- FULL_T_F_PC$TF.V2

G_Amount_F$B_V <- R_C_V$B_V
G_Amount_F$C_V <- R_C_V$C_V
G_Amount_F$Death <- FULL_T_F_PC$TF.V2
Fitness_Cut_PC$Death <- FULL_T_F_PC$TF.V2

RBC_Amount_F_35 <- subset(RBC_Amount_F, 
                          RBC_Amount_F$C_V==0.35 &
                            RBC_Amount_F$Death == 1 &
                            RBC_Amount_F$B_V >= 5.5)

G_Amount_F_35 <- subset(G_Amount_F, 
                        G_Amount_F$C_V==0.35 &
                          G_Amount_F$Death==1 &
                          G_Amount_F $B_V >= 5.5)

Fitness_Cut_PC_35<- subset(Fitness_Cut_PC, 
                           Fitness_Cut_PC$C_V==0.35 &
                             Fitness_Cut_PC$Death==1 &
                             Fitness_Cut_PC$B_V >=5.5)

Duration_Initial_PC_2_35 <- 
  subset(Duration_Initial_PC_2, 
         Duration_Initial_PC_2$C_V==0.35
         &
           Duration_Initial_PC_2$Death==1 &
           Duration_Initial_PC_2$B_V >=5.5)


Fitness_Cut_PC_35$RBC <- RBC_Amount_F_35$V1
Fitness_Cut_PC_35$G <- G_Amount_F_35 $V1
Fitness_Cut_PC_35$Opt <- RBC_Amount_F_35$Opt
Fitness_Cut_PC_35$Duration <- Duration_Initial_PC_2_35$endtime


###Proportion of red blood cells compared too burst size
rbc_burst <- ggplot(Fitness_Cut_PC_35,
                    aes (x=B_V,y= (1- RBC/8500000)))+
  geom_line()+
  geom_point(size =3)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  xlab("Burst size ")+
  ylab("Proportion of red blood cells lost")



###Proportion of red blood cells lost compared
###to the proportion of maximum gametocytes 
rbc_gam <- ggplot(Fitness_Cut_PC_35,
                  aes (x= (1- RBC/8500000), y= G/231044.2))+
  geom_line()+
  geom_point(size =3)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  xlab("Proportion of red blood cells lost")+
  ylab("Proportion of maximum gametocytes")


###################################################
###Gametocyes to cumulative transmission potential#
###################################################

gam_fit <- 
  ggplot(
    Fitness_Cut_PC_35,
    aes(x= G/231044.2,
        y = fitness/23.89940))+
  geom_line()+
  geom_point(size=3, 
             lwd=2)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  xlab("Proportion of maximum gametocytes")+
  ylab("Proportion of maximum fitness")+
  theme(legend.position = "none")


rbc_duration <- ggplot(
  Fitness_Cut_PC_35,
  aes (x= 1- RBC/8500000,
       y = Duration))+
  geom_line()+
  geom_point(size=3,lwd=2)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Duration (days)")+
  xlab("Proportion of red blood cells lost")+
  theme(legend.position = "none")

gam_fit /(rbc_burst + rbc_gam + rbc_duration )
