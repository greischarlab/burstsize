
###This is for making the multi-panel figures with the 
###min rbc lost to burst size, cumulative gametocyte number to
###burst size, acute phase (length in days) to burst size.

###x_list <- is the FULL_MODEL OUTPUT (Includes RBC and G)
###y_list <- is the FITNESS MODEL
###C_V_Value is the transmission investment of interest(usually, it's the 
###optimal one)

Virulence_Transmission_DF_Calculator<- function(x_list,
                                                y_list,
                                                C_V_value ){
  
  ###Info of the strains that lead to successful infection (
  ###both surviving and leading to mortality)
  Fitness_Subset <- subset(y_list,
                           y_list$C_V == C_V_value &
                           y_list$status != 'Fail')
  
  ###I order the burst size just to make sure the combining
  ###goes smoothly
  Fitness_Subset <-  Fitness_Subset[order( Fitness_Subset$B_V),]

  ###Info of the strains that will induce host mortality
  Fitness_Subset_M <- subset(y_list,
                             y_list$C_V == C_V_value &
                             y_list$status == 'mort')
  
  ###I order the burst size just to make sure the combining
  ###goes smoothly
  Fitness_Subset_M <-  Fitness_Subset_M[order( Fitness_Subset_M$B_V),]
  
  ###Burst sizes of interest
  Burst_size_Interest <-  Fitness_Subset$B_V
  
  ###Look at the full model list
  Full_Subset <- subset(x_list,
                        round(x_list$C_V,3) == C_V_value &
                        x_list$B_V %in% Burst_size_Interest)
  
  ###Pull out the RBC min
  RBC_Min <- do.call(rbind, by(Full_Subset, Full_Subset$B_V, 
                                  function(x) 
                                  x[which.min(x[,"R"]),],
                                  simplify = FALSE))
  
  RBC_Min <- RBC_Min[order(RBC_Min$B_V),]
  
  #####################################################
  ####Gametocytes require more work due to some burst##
  ###size leading to host mortality ###################
  ##################################
  
  ###This looks at the full model that leads to mortality
  Full_Subset_M <- subset(Full_Subset, 
                          Full_Subset$B_V %in% 
                          c(Fitness_Subset_M$B_V))
  
  split_Full_Subset_M <- split(Full_Subset_M, Full_Subset_M$B_V)
  
  G_Cum_Mortality = NULL
  for (k in seq(1, length(split_Full_Subset_M))){
    tmp_endtime = Fitness_Subset_M[k,'endtime']
    gam_approx_fun <- approxfun( split_Full_Subset_M[[k]][,c('time','G')])
    
    
    total_G<- sum(diff(gam_approx_fun(seq(1,tmp_endtime, length = 1000))))
    
    

    G_Cum_Mortality[[k]] = cbind.data.frame(G = total_G, 
                                  B_V = unique(split_Full_Subset_M[[k]]$B_V))
  }
  
  G_Cum_Mortality <- do.call(rbind, G_Cum_Mortality)
  
  
  ###Calculating total gametocytes for the acute phase
  
  ###This looks for the burst size of successful infections that don't lead
  ###host mortality
  Fit_Subset_S <- subset(Fitness_Subset, Fitness_Subset$status == 'success')
  
  Full_Subset_S <- subset(Full_Subset, Full_Subset$B_V %in% Fit_Subset_S$B_V)
  
  split_Full_Subset_S <- split(Full_Subset_S, Full_Subset_S$B_V)
  

  G_Cum_Success = NULL
  for (k in seq(1, length(split_Full_Subset_S))){
    tmp_endtime = Fit_Subset_S [k,'endtime']
    gam_approx_fun <- approxfun(split_Full_Subset_S[[k]][,c('time','G')])
    
    
    total_G<- sum(diff(gam_approx_fun(seq(1,tmp_endtime, length = 1000))))
    
    
    
    G_Cum_Success[[k]] = cbind.data.frame(G = total_G, 
                                            B_V = unique(split_Full_Subset_S[[k]]$B_V))
  }
  G_Cum_Success <- do.call(rbind, G_Cum_Success)
  
  colnames(G_Cum_Success) <- c('G','B_V')
  G_Cum_ALL <- rbind(G_Cum_Mortality, G_Cum_Success)
  
  G_Cum <- data.frame(G_Cum_ALL[order(G_Cum_ALL[,'B_V']),])
  

  ###Has the burst size/transmission investment strains with differences
  ###in min RBC, max Gametocyte, cumulative transmission potential (CTP),
  ###length of the acute phase, and the status of the infection (should 
  ###alive or dead no FAIL)
  Full_DF <- cbind.data.frame(RMin = log10(RBC_Min$R), 
                              GCum =  log10(G_Cum$G),
                              B_V = G_Cum$B_V, 
                              CTP = Fitness_Subset$end_fitness,
                              End_Time = Fitness_Subset$endtime,
                              Status = Fitness_Subset$status)
  
  return(Full_DF)
  
}
