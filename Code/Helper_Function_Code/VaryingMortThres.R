###SUPP 

Finder_RM_DeathLimit <- function(x_list, mu_M_c, DeathLimit){
  
  ###Shape parameter (to ensure we're getting the right columns)
  n1 = 100
  ###Check first if the infection induces mortality
  ###Look at the red blood cells time series only  
  RBC_TS <- data.frame(time = x_list$time,
                       RBC = x_list[,"R"])
  
  ### Look at the red blood cells to see when it is under 6.5 * 10^5
  index_time_mort <- RBC_TS[which(RBC_TS$RBC <= DeathLimit),]
  
  ###If the index_time_mort is empty, then the host survived
  ###and we note that
  
  if (nrow(index_time_mort) == 0){
    mort_time <- cbind.data.frame(surv = 'surv',
                                  time = NA)
  } else {
    
    ###Use approxfun to calculate exaclty when the curve hits the mortality
    ###threshold
    
    RBC_TS_Func <- approxfun(RBC_TS[,'time'], RBC_TS[,"RBC"] - DeathLimit )
    death_time <- uniroot(RBC_TS_Func , c(0,index_time_mort$time[1]))
    
    mort_time <- cbind.data.frame(surv = 'mort',
                                  time = death_time$root)
  }
  ###This now checks the Acute Phase Duration for the Mice 
  ###That Survived
  
  if(mort_time$surv == "surv"){
    p =  2.5e-6
    
    rate = (1 - unique(x_list$C_V)) * 
      unique(x_list$B_V)*((x_list[,"R"] * p)/((p * x_list[,"R"]) + 
                                                mu_M_c))
    
    RM_time_df <-  cbind.data.frame(time = x_list[,'time'],
                                    rate = rate, 
                                    B_V = unique(x_list$B_V),
                                    C_V = unique(x_list$C_V))
    
    min_RM <- RM_time_df[which.min(RM_time_df$rate),]
    
    end_time <- subset(RM_time_df, 
                       RM_time_df$time >= min_RM$time & 
                         RM_time_df$rate >= 1)[1,'time']
    
    
    
    return(data.frame(
      endtime =  end_time,
      up_down =  'up',
      end_fitness = NA,
      status = 'success'))
    
    
  }else{
    
    ###Look at the infected red blood cells time series only  
    ###Calculating the fitness
    G_TS <- data.frame(time = x_list[,'time'],
                       G = x_list[,"G"])
    
    #Prevent negative gamaetocytes 
    G_TS$G[G_TS$G < 0 ] <- 0 
    
    truncate_G_TS <- subset(G_TS, G_TS$time <= mort_time$time)
    
    
    end_fitness_mort <-  sum(PrI_PC(truncate_G_TS$G) * 1/10)
    
    return(data.frame(
      endtime = mort_time$time,
      up_down = "up" ,
      end_fitness = end_fitness_mort, 
      status = 'mort'))
  }
}



grapher_mortality_boundary_Supp <- function (x){
  
  tmp_df <- data.frame(B_V = x$B_V,
                       C_V = x$C_V,
                       Status = x$status)
  
  tmp_df$Status_Num <- ifelse(tmp_df$Status == 'mort', 1, 0)
  
  tmp_mat_cast_mort <- acast(tmp_df,
                             C_V ~ B_V, 
                             value.var = c("Status_Num"))
  
  vertical_df_mort <- NULL
  for (n in seq(1,nrow(tmp_mat_cast_mort))){
    
    row_numbers = as.numeric(rownames(tmp_mat_cast_mort)[n])
    
    tmp_vec_mort <- tmp_mat_cast_mort[n, 2:ncol(tmp_mat_cast_mort)] - 
      tmp_mat_cast_mort[n, 1:ncol(tmp_mat_cast_mort) - 1] 
    
    vertical_mort <- as.numeric(names(tmp_vec_mort[tmp_vec_mort == 1 ]))
    
    vertical_F_mort <- ifelse(length(vertical_mort)==0, NA, vertical_mort)
    
    row_F_mort <- ifelse(length(vertical_mort) == 0, NA, row_numbers)
    
    vertical_df_mort[[n]] <- data.frame(
      x = vertical_F_mort - 0.5 , 
      xend = vertical_F_mort - 0.5 , 
      y = row_F_mort -  0.05/2, 
      yend = row_F_mort +  0.05/2)
    
  }
  
  vert_seg_mort <- na.omit(do.call(rbind,vertical_df_mort))
  vert_seg_mort$id = "mort"
  
  ####################################
  ###Horizontal lines for mortality###
  ####################################
  horizontal_df_mort <- NULL
  for (n in seq(1,ncol(tmp_mat_cast_mort))){
    col_mort = as.numeric(colnames(tmp_mat_cast_mort)[n])
    
    tmp_vec_mort <- tmp_mat_cast_mort[2:nrow(tmp_mat_cast_mort), n] -
      tmp_mat_cast_mort[1:nrow(tmp_mat_cast_mort) - 1,n] 
    
    horizontal_mort <- as.numeric(names(tmp_vec_mort[tmp_vec_mort == -1]))
    
    horizontal_F_mort <- ifelse(length(horizontal_mort) == 0, NA, 
                                horizontal_mort)
    
    col_F_mort <- ifelse(length(horizontal_mort) == 0, NA, col_mort)
    
    horizontal_df_mort [[n]] <- data.frame(
      x = col_F_mort - 0.5 ,
      xend = col_F_mort + 0.5 , 
      y = horizontal_F_mort -  0.05/2, 
      yend = horizontal_F_mort  - 0.05/2)
    
  }
  hor_seg_mort <- na.omit(do.call(rbind,horizontal_df_mort))
  hor_seg_mort$id <- 'mort'
  ################################
  ###Does not establish infection#
  ################################
  
  tmp_df <-    data.frame(B_V = as.factor(x$B_V),
                          C_V = as.factor(x$C_V),
                          Status = x$status)
  
  tmp_df$Status_Num_fail <- ifelse(tmp_df$Status == 'Fail', 0, 1)
  
  tmp_mat_cast_fail <- acast(tmp_df,
                             C_V ~ B_V, 
                             value.var = c("Status_Num_fail"))
  
  
  vertical_df_fail <- NULL
  for (n in seq(1,nrow(tmp_mat_cast_fail))){
    row = as.numeric(rownames(tmp_mat_cast_fail)[n])
    
    tmp_vec <- tmp_mat_cast_fail[n,2:ncol(tmp_mat_cast_fail)] -
      tmp_mat_cast_fail[n,1:ncol(tmp_mat_cast_fail) - 1] 
    
    vertical <- as.numeric(names(tmp_vec[tmp_vec == 1]))
    
    vertical_F <- ifelse(length(vertical) == 0, NA, vertical)
    row_F <- ifelse(length(vertical) == 0, NA, row)
    
    vertical_df_fail[[n]] <- data.frame(
      x = vertical_F  - 0.5 , 
      xend = vertical_F - 0.5 , 
      y = row_F -  0.05/2, 
      yend = row_F +  0.05/2)
    
  }
  
  vert_seg_fail <- na.omit(do.call(rbind, vertical_df_fail))
  vert_seg_fail$id <-  "fail"
  
  horizontal_df_fail <- NULL
  
  for (n in seq(1,ncol(tmp_mat_cast_fail))){
    
    col = as.numeric(colnames(tmp_mat_cast_fail)[n])
    
    tmp_vec <- tmp_mat_cast_fail[2:nrow(tmp_mat_cast_fail),n] - 
      tmp_mat_cast_fail[1:nrow(tmp_mat_cast_fail)-1,n] 
    
    horizontal <- as.numeric(names(tmp_vec[tmp_vec == -1]))
    
    horizontal_F <- ifelse(length(horizontal) == 0, NA, horizontal)
    col_F <- ifelse(length(horizontal) == 0, NA, col)
    
    horizontal_df_fail [[n]] <- data.frame(
      x = col_F - 0.5 ,
      xend = col_F + 0.5 , 
      y = horizontal_F - 0.05/2, 
      yend = horizontal_F  -  0.05/2)
    
  }
  hor_seg_fail <- na.omit(do.call(rbind,horizontal_df_fail))
  hor_seg_fail$id <-  "fail"
  
  full_df_fail_mort <- rbind.data.frame(vert_seg_mort,
                                        hor_seg_mort,
                                        vert_seg_fail,
                                        hor_seg_fail)
  
  return(full_df_fail_mort)
}

