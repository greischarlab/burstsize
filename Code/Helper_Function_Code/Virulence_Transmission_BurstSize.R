
###This is for making the multi-panel figures with 
###Cumulative transmission potential versus max. gametocyte number,
###min rbc lost to burst size, max gametocyte number to burst size,
###acute phase to burst size 




Virulence_Transmission_DF_Calculator<- function(x_list,y_list,C_V_value){
  
  ###Info of the strains that lead to successful infection
  Fitness_Subset<- subset(y_list,
                           y_list$C_V == C_V_value &
                           y_list$status != 'Fail')
  
  Fitness_Subset <-  Fitness_Subset[order( Fitness_Subset$B_V),]

  ###Info of the strains that will induce host mortality
  Fitness_Subset_M <- subset(y_list,
                             y_list$C_V == C_V_value &
                             y_list$status == 'mort')
  
  Fitness_Subset_M <-  Fitness_Subset_M[order( Fitness_Subset_M$B_V),]
  
  ###Burst size that does not induce mortality AND establish infection
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
  
  ##################################
  ####Gametocytes require more work#
  ##################################
  Full_Subset_M <- subset(Full_Subset, 
                          Full_Subset$B_V %in% 
                          c(41,42))
  
  split_Full_Subset_M <- split(Full_Subset_M, Full_Subset_M$B_V)
  
  G_Max_Mortality = NULL
  for (k in seq(1, length(split_Full_Subset_M))){
    tmp_endtime = Fitness_Subset_M[k,'endtime']
    tmp = split_Full_Subset_M[[k]][split_Full_Subset_M[[k]]$time <= 3.0,]
    
    G_Max_Mortality[[k]] = cbind.data.frame(G = max(tmp$G), 
                                  B_V = unique(split_Full_Subset_M[[k]]$B_V))
  }
  
  G_Max_Mortality <- do.call(rbind, G_Max_Mortality)
  
  
  Fit_Subset_S <- subset(Fitness_Subset, Fitness_Subset$status == 'success')
  
  Full_Subset_S <- subset(Full_Subset, Full_Subset$B_V %in% Fit_Subset_S$B_V)
  
  G_Max_Success <-    do.call(rbind, 
                                by(Full_Subset_S, Full_Subset_S$B_V, 
                                function(x) 
                                x[which.max(x[,"G"]),c("G","B_V")],
                                simplify = FALSE))
  
  G_Max_ALL <- rbind(G_Max_Mortality, G_Max_Success)
  
  G_Max <- G_Max_ALL[order(G_Max_ALL$B_V),]
  
  ###Has the burst size/transmission investment strains with differences
  ###in min RBC, max Gametocyte, cumulative transmission potential (CTP),
  ###length of the acute phase, and the status of the infection (should 
  ###alive or dead no FAIL)
  Full_DF <- cbind.data.frame(RMin = RBC_Min$R, 
                              GMax = G_Max$G,
                              B_V = RBC_Min$B_V, 
                              C_V = RBC_Min$C_V,
                              CTP = Fitness_Subset$end_fitness,
                              End_Time = Fitness_Subset$endtime,
                              Status = Fitness_Subset$status)
  
  return(Full_DF)
  
}
