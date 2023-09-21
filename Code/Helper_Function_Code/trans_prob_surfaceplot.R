
Trans_SurfacePlot_Maker <- function(x, y, CV_vec,initial_value) {

###Look at the fitness data.frame (y)
  y_subsetted <- subset(y, round(y$C_V,4) %in% CV_vec)

  ###################################
  ###For the unestablished infection#
  ###################################
  
  y_subsetted_fail <- subset(y_subsetted, y_subsetted$status == 'Fail')
  
  y_subsetted_fail_list <- NULL
  
  for(k in seq(1,nrow(y_subsetted_fail))){
    
    tmp <-  y_subsetted_fail[k,]
    
    y_subsetted_fail_list[[k]] =
                data.frame(time =  seq(0,100,0.1),
               B_V = tmp$B_V,
               C_V = tmp$C_V,
               pripc = 0,
               endtime = 0,
               peaktime = 0
               )
    
    
  }
  y_df_fail <- do.call(rbind, y_subsetted_fail_list )
  y_df_fail$id <- 'Fail'
  ###For the infections that induces mortality
  
  
  y_subsetted_mort <- subset(y_subsetted, y_subsetted$status == 'mort')
  
  if(nrow(y_subsetted_mort) != 0){
  FULL_MODEL_PC_100_MORT <- mcmapply(Simulator_Malaria_BC,
                                c( y_subsetted_mort$B_V),
                                c( y_subsetted_mort$C_V),
                                c(initial_value),
                                mc.cores = 4,
                                SIMPLIFY = FALSE)
  
y_subsetted_mort_list = NULL

for (k in seq(1,length(FULL_MODEL_PC_100_MORT))){
  
   tmp = FULL_MODEL_PC_100_MORT[[k]]
   tmp2 =   y_subsetted_mort[k,]
   tmp$time <- round(tmp$time, 3)
   time_df <-data.frame(time =  round(seq(0,100,0.1),4)
                        )
   full_df<-  left_join(time_df, tmp, by = 'time')
   
   full_df$G[is.na(full_df$G) == TRUE] <- 0
   
   full_df <- full_df[order(full_df$time),]
  
   
   y_subsetted_mort_list[[k]] =
     data.frame(time =  full_df$time,
                B_V = unique(tmp$B_V),
                C_V = unique(tmp$C_V),
                pripc = PrI_PC(full_df$G),
                endtime = tmp2$endtime,
                peaktime = 0
     )
   
}

   y_df_mort <- do.call(rbind, y_subsetted_mort_list )
  }else{
  y_df_mort <-      data.frame(time = NA,
                               B_V = NA,
                               C_V = NA,
                               pripc = NA,
                               endtime =NA,
                               peaktime = NA
  )
  
  }
y_df_mort$id <- 'mort'
  

####SUCCESSFUL INFECTIONS
   
   
   FULLMODEL_Subsetted <- subset(x,  round(x$C_V,4) %in% CV_vec &
                                   x$time <=100)
   
   y_subsetted_success<- subset(y_subsetted, y_subsetted$status == 'success')
   
   CV_Subsetted_Status_SUCCESS <- merge(FULLMODEL_Subsetted, y_subsetted_success, by= c("B_V", "C_V"))
   
   CV_Subsetted_Status_Split <- split(CV_Subsetted_Status_SUCCESS,CV_Subsetted_Status_SUCCESS$B_V)
   
   y_subsetted_success_list <- NULL
   
   for (k in seq(1,length(CV_Subsetted_Status_Split))){
     
     tmp <-    CV_Subsetted_Status_Split [[k]]
     
     y_subsetted_success_list [[k]] =
       data.frame(time =   tmp$time,
                  B_V = unique(tmp$B_V),
                  C_V = unique(tmp$C_V),
                  pripc = PrI_PC(  tmp $G),
                  endtime = unique(tmp$endtime),
                  peaktime = 0
       )
   }
   y_df_success <- do.call(rbind,   y_subsetted_success_list )
  y_df_success$id <- 'success'

   
   y_df_ALL <- rbind(y_df_mort, y_df_fail, y_df_success) 
   


return(y_df_ALL)
}



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
