###A helper function for the main script 04_Set_Duration_SurfacePlot.R

###This is the plotting function that shows how the transmission
### probability changes with the burst size over the infection- 

###Input:
#1) x - Full Model Data Frame
#2) y - Fitness Data frame
#3) CV_vec - what transmission investments are you interested in
#4) initial_value = what is the initial value of the infected RBC you're
### interested in
###Output: 
###A list specifically for the trans_surfaceplot_maker 

Trans_SurfacePlot_Maker <- function(x, y, CV_vec, initial_value) {

     ###Look at the fitness data.frame (y) and subset the 
     ### transmission investment of interest
     y_subsetted <- subset(y, round(y$C_V,4) %in% CV_vec)

     ###################################
     ###For the unestablished infection#
     ###################################
  
      ###Failed infections
      y_subsetted_fail <- subset(y_subsetted, y_subsetted$status == 'Fail')
    
      y_subsetted_fail_list <- NULL
  
        for(k in seq(1,nrow(y_subsetted_fail))){
    
          tmp <-  y_subsetted_fail[k,] #Iterate over each strain
    
          y_subsetted_fail_list[[k]] =
                data.frame(time =  seq(0,100,0.1),
                B_V = tmp$B_V,
                C_V = tmp$C_V,
                pripc = 0, # Failed infections are always 0
                endtime = 0, #No infection so 0
                peaktime = 0 #No infection so 0
               )
           }
    
      y_df_fail <- do.call(rbind, y_subsetted_fail_list)
      y_df_fail$id <- 'Fail'
    
      ##############################################
      ###For the infections that induces mortality##
      ##############################################
  
      y_subsetted_mort <- subset(y_subsetted, y_subsetted$status == 'mort')
  
      ###If there are strains that induces host mortality
        if(nrow(y_subsetted_mort) != 0){
  
  ###Run simulations
      FULL_MODEL_PC_100_MORT <- mcmapply(Simulator_Malaria_BC,
                                c(y_subsetted_mort$B_V),
                                c(y_subsetted_mort$C_V),
                                c(initial_value),
                                mc.cores = 4,
                                SIMPLIFY = FALSE)
  
      y_subsetted_mort_list = NULL

      ###The entire point of this function is that some strains kill
      ###The infection at day 4 or day 30, I want to ensure that 
      ###that all days are included even if the host is dead 
      
      for (k in seq(1,length(FULL_MODEL_PC_100_MORT))){
  
       tmp = FULL_MODEL_PC_100_MORT[[k]] #For the full model of interest
       tmp2 =   y_subsetted_mort[k,] ### get the information
       tmp$time <- round(tmp$time, 3)
       time_df <-data.frame(time =  round(seq(0,100,0.1),4))
 
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
  }
        else{
         y_df_mort <-      data.frame(time = NA,
                               B_V = NA,
                               C_V = NA,
                               pripc = NA,
                               endtime =NA,
                               peaktime = NA)
  
  }
       y_df_mort$id <- 'mort'
  
       ########################################################
       ####SUCCESSFUL INFECTIONS THAT DO NOT KILL THE HOST####
      #######################################################
   
       FULLMODEL_Subsetted <- subset(x,  
                                   round(x$C_V,4) %in% CV_vec &
                                   x$time <=100)
   
    
        y_subsetted_success<- subset(y_subsetted, y_subsetted$status == 'success')
   
        CV_Subsetted_Status_SUCCESS <- merge(FULLMODEL_Subsetted, y_subsetted_success, by= c("B_V", "C_V"))
   
        CV_Subsetted_Status_Split <- split(CV_Subsetted_Status_SUCCESS,CV_Subsetted_Status_SUCCESS$B_V)
   
        y_subsetted_success_list <- NULL
   
        for (k in seq(1,length(CV_Subsetted_Status_Split))){
     
               tmp <- CV_Subsetted_Status_Split [[k]]
     
               y_subsetted_success_list [[k]] =
                  data.frame(time =   tmp$time,
                   B_V = unique(tmp$B_V),
                   C_V = unique(tmp$C_V),
                   pripc = PrI_PC(tmp $G),
                   endtime = unique(tmp$endtime),
                   peaktime = 0)
   }
        y_df_success <- do.call(rbind, y_subsetted_success_list )
        y_df_success$id <- 'success'
        y_df_ALL <- rbind(y_df_mort, y_df_fail, y_df_success) 
   
         return(y_df_ALL)
}

