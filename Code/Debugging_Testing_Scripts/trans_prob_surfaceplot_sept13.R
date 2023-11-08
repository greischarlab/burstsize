

PRIPC_MED_90 <-Trans_SurfacePlot_Maker(x=FULL_MODEL_100_MED,
                                       y=Fitness_MODEL_PC_FULL_MED, CV_vec = c(0.90),4358.965)
PRIPC_MED_MID <-Trans_SurfacePlot_Maker(x=FULL_MODEL_100_MED,
                                        y=Fitness_MODEL_PC_FULL_MED, CV_vec = c(0.50),4358.965)
PRIPC_MED_70 <-Trans_SurfacePlot_Maker(x=FULL_MODEL_100_MED,
                                       y=Fitness_MODEL_PC_FULL_MED, CV_vec = c(0.75),4358.965)

PRIPC_MED_90$time <- round(PRIPC_MED_90$time,3)
PRIPC_MED_MID$time <- round(PRIPC_MED_MID$time,3)
PRIPC_MED_70$time <- round(PRIPC_MED_70$time,3)

PRIPC_MED_MID$pripc <- round(PRIPC_MED_MID$pripc,2)
PRIPC_MED_90$pripc <- round(PRIPC_MED_90$pripc,2)
PRIPC_MED_70$pripc <- round(PRIPC_MED_70$pripc,2)

PRIPC_MED_90$cutprob <- cut(PRIPC_MED_90$pripc, breaks=c(-1,0,0.05,0.1,.15,.2,.25,.3,.35,.4,
                                                         .45,.50,.55,.60,.65,.70,.75,.80,
                                                         .85,.90,.95,1))
PRIPC_MED_70$cutprob <- cut(PRIPC_MED_70$pripc, breaks=c(-1,0,0.05,0.1,.15,.2,.25,.3,.35,.4,
                                                         .45,.50,.55,.60,.65,.70,.75,.80,
                                                         .85,.90,.95,1))
PRIPC_MED_MID$cutprob <- cut(PRIPC_MED_MID$pripc, breaks=c(-1,0,0.05,0.1,.15,.2,.25,.3,.35,.4,
                                                           .45,.50,.55,.60,.65,.70,.75,.80,
                                                           .85,.90,.95,1))

a<-ggplot(PRIPC_MED_90,
          aes( x= time, y= B_V, fill = (cutprob), z = pripc )) + 
  geom_raster() +
  scale_fill_viridis(option = 'mako', discrete = 'TRUE')+
  xlab("Days post infection") + ylab("Burst size") + 
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  ggtitle('TRANSMISSION INVESTMENT = 90%') 

b<-  ggplot(PRIPC_MED_MID,
            
            aes( x= time, y= B_V, fill = (cutprob), z = pripc )) + 
  geom_raster() +
  scale_fill_viridis(option = 'mako', discrete = 'TRUE')+
  xlab("Days post infection") + ylab("Burst size") + 
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  ggtitle('TRANSMISSION INVESTMENT = 50%') 

c<- ggplot(PRIPC_MED_70, 
           aes( x= time, y= B_V, fill = (cutprob), z = pripc )) + 
  geom_raster() +
  scale_fill_viridis(option = 'mako', discrete = 'TRUE')+
  xlab("Days post infection") + ylab("Burst size") + 
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  ggtitle('TRANSMISSION INVESTMENT = 10%') 


(a/b/c) +plot_layout(guides = 'collect')