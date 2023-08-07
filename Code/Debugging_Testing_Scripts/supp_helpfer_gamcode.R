###supp gam code


Gametocyte_Fitness_time <- function(x,mean_span){
  
  gam_dat <- data.frame(time = x['time'],
                        gam = x[,"G"])
  
  ###If it's negative then, make it 0 just in case
  gam_dat$gam[gam_dat$gam < 0] <- 0            
  
  end_fitness <- (cumsum(exp(-(1/mean_span)*gam_dat $time)* PrI_PC(gam_dat$gam) * 1/10))
  
  return(end_fitness)
}

#########################################################################
###Input: Data.frame of all the life-stages from the de-solve simulator #
###Output: The maximum cumulative transmission potential at the end of the 
#infection length of interest  
#########################################################################  
Gametocyte_Fitness_LXMX <- function(x, mean_span){
  
  ###All of the list elements should have the same infection 
  ###length- this means that it is a 100 day infection
  gam_dat <- data.frame(time = x['time'],
                        gam = x[,"G"])
  
  ###If it's negative then, make it 0 just in case
  gam_dat $gam[gam_dat $gam < 0] <- 0            
  
  
  ###Calculate the max transmission potential 
  end_fitness <- max(cumsum(exp(-(1/mean_span)*gam_dat$time)*(PrI_PC(gam_dat $gam) * 1/10)))
  
  return(end_fitness)
}


Gametocyte_Fitness_EqualWeighting <- function(x, life_span){
  
  gam_dat <- data.frame(time = x['time'],
                        gam = x[,"G"])
  
  ###If it's negative then, make it 0 just in case
  gam_dat$gam[gam_dat$gam < 0] <- 0        
  
  gam_dat_func <- approxfun(gam_dat)
  
  gam_dat_trunc <-gam_dat_func(seq(1,life_span,1/10))
  
  end_fitness <- max(cumsum(1*(PrI_PC(gam_dat_trunc) * 1/10)))
  
  return(end_fitness)
}


###################################################
###INPUT: The element of the full deSolve output###
### as well as the mu_M                         ###
###################################################
Finder_RM_lxmx <- function(x_list, mu_M_c,
                           mean_span, life_time){
  
  ###Shape parameter (to ensure we're getting the right columns)
  n1 = 100
  ###Check first if the infection induces mortality
  ###Look at the red blood cells time series only  
  RBC_TS <- data.frame(time = x_list$time,
                       RBC = x_list[,"R"])
  
  ### Look at the red blood cells to see when it is under 6.5 * 10^5
  index_time_mort <- RBC_TS[which(RBC_TS$RBC <= 6.5 * 10^5),]
  
  ###If the index_time_mort is empty, then the host survived
  ###and we note that it 'surv'
  if (nrow(index_time_mort) == 0){
    mort_time <- cbind.data.frame(surv = 'surv',
                                  time = NA)
  } else {  ###If the index_time_mort is not empty, then host RBC 
    ###dips
    
    ###Use approxfun to calculate exaclty when the curve hits the mortality
    ###threshold
    RBC_TS_Func <- approxfun(RBC_TS[,'time'], RBC_TS[,"RBC"] - 6.5 * 10^5 )
    
    ###We use uniroot to figure out when it hits 0
    death_time <- uniroot(RBC_TS_Func , c(0,index_time_mort$time[1]))
    
    mort_time <- cbind.data.frame(surv = 'mort',
                                  time = death_time$root)
  }
  ###This now checks the Acute Phase Duration for the Mice 
  ###That Survived
  
  if(mort_time$surv == "surv"){
    p <-   2.5e-6
    B_V_Interest <- unique(x_list$B_V)
    C_V_Interest <- unique(x_list$C_V)
    
    rate = (1 -  C_V_Interest) * 
      B_V_Interest*((p * x_list[,"R"] )/((p * x_list[,"R"]) + 
                                           mu_M_c))
    
    RM_time_df <-  cbind.data.frame(time = x_list[,'time'],
                                    rate = rate, 
                                    B_V =  B_V_Interest,
                                    C_V = C_V_Interest)
    
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
    
    
    
    
    end_fitness_mort <-  max(cumsum(exp(-(1/mean_span)* truncate_G_TS$time)*
                                      (PrI_PC(truncate_G_TS$G) * 1/10)))
    
    return(data.frame(
      endtime = mort_time$time,
      up_down = "up" ,
      end_fitness = end_fitness_mort, 
      status = 'mort'))
  }
}


