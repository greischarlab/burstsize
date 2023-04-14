########################################
###These are the functions that you use#
###to quantify fitness for all analyses#
########################################
#######################################################################
###Function: Daily transmission calculator based on gametocyte density# 
###Input: Gametocyte density                                          #
###Output: Probability of transmitting to a mosquito                  #
#######################################################################
PrI_PC <- function(G){
           exp(-12.69 + 3.6 * log10(G))/(1+exp(-12.69 + 3.6 * log10(G)))
}
#########################################################################
###Input: Data.frame of all the life-stages from the de-solve simulator #
###Output: Cumulative daily transmission across time                    #
#########################################################################    
TranP_Func <- function(x) {
  gam = as.data.frame(x)$G
  gam[gam < 0 ] <- 0  #ensures that any gametocytes that are negative is 0
  df <- data.frame(
    time = as.data.frame(x)$time,
    tp = cumsum(PrI_PC(gam) * (1/10)))
  return(df)
}

#########################################################################
###Input: Data.frame of all the life-stages from the de-solve simulator #
###Output: The maximum cumulative transmission potential                #
#########################################################################  

Gametocyte_Fitness <- function(x){
  
  gam_dat <- data.frame(time = x['time'],
                        gam = x[,"G"])
  
  ###If it's negative then, make it 0 just in case
  gam_dat$gam[gam_dat$gam < 0] <- 0            
  
  end_fitness <- max(cumsum(PrI_PC(gam_dat) * 1/10))
  
  return(end_fitness)
}

###################################################
###INPUT: The element of the full deSolve output###
### as well as the mu_M                         ###
###################################################
Finder_RM <- function(x_list, mu_M_c){
    
    ###Shape parameter (to ensure we're getting the right columns)
    n1 = 100
    ###Check first if the infection induces mortality
    ###Look at the red blood cells time series only  
    RBC_TS <- data.frame(time = x_list$time,
                         RBC = x_list[,"R"])
  
     ### Look at the red blood cells to see when it is under 6.5 * 10^5
    index_time_mort <- RBC_TS[which(RBC_TS$RBC <= 6.5 * 10^5),]
  
     ###If the index_time_mort is empty, then the host survived
     ###and we note that
    
    if (nrow(index_time_mort) == 0){
        mort_time <- cbind.data.frame(surv = 'surv',
                                  time = NA)
    } else {
      
      ###Use approxfun to calculate exaclty when the curve hits the mortality
      ###threshold
      
      RBC_TS_Func <- approxfun(RBC_TS[,'time'], RBC_TS[,"RBC"] - 6.5 * 10^5 )
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

  
  
  