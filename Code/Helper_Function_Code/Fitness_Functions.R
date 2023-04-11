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


###DEFUNCT: But keeping this in
###This is the main function used to figure out when to figure out
###if the infection leads to either mortality or failed/successful
### infections. If the infection kills the host, I calculate the 
###fitness up to the time of death. If the infection does not kill
###the host, I check to see if it's a failed infection. If it's a
###failed infection, the initial IRBC decreases and thus the fitness
###is automatically set to 0. If it's a successful infection, then
###the end fitness is set to NA and I note the duration of the acute 
###phase 

###################################################
###INPUT: The element of the full deSolve output###
###################################################
#Finder_Peak_Duration_Mortality <- function(x){
#        
#        ###Shape parameter (to ensure we're getting the right columns)
#        n1=100
#        ###Check first if it induces mortality
#        ###Look at the red blood cells time series only  
#        
#        RBC_TS <- data.frame(time = x$time,
#                             RBC = x[,"R"])
#        
#        ### Look at the red blood cells to see when it is under 6.5 * 10^5
#        index_time_mort <- RBC_TS[which(RBC_TS$RBC <= 6.5 * 10^5),]
#    
#        ###If the index_time_mort is empty, then the host survived
#        ###and we note that
#        if (nrow(index_time_mort) == 0){
#          mort_time <- cbind.data.frame(surv = 'surv',
#                                        time = NA)
#          } else {
#          mort_time <- cbind.data.frame(surv = 'mort',
#                                        time = index_time_mort$time[1])
#          }
#        
#        ###This now checks the Acute Phase Duration for the Mice 
#        ###That Survived
#        if(mort_time$surv == "surv"){
#  
#        ###Look at the infected red blood cells time series only  
#        IRBC_TS <- data.frame(time = x$time,
#                              IRBC = rowSums(x[,3:(n1+2)]))
#  
#        ### Look at the maximum IRBC peak and pull out the time     
#        ### if time = 0 for the index time, then this means that IRBC
#        ### goes down- we should keep a note of it                   
#        index_time <-  IRBC_TS[which.max(IRBC_TS$IRBC),] # find the maximum
#  
#        #If the peak is 0 that suggests that the infection never establishes
#        up_down_factor <- ifelse(index_time$time == 0, 'down', 'up')
#  
#        ###The end of the acute phase
#        end_duration <- (index_time$time * 2.0)
#        
#        end_fitness_status <-  ifelse(end_duration == 0, 0, NA)
#        
#        infection_status <- ifelse(end_duration ==0, 'failed', 'success')  
#        
#        return(data.frame(
#          peaktime = index_time$time,
#          endtime = end_duration,
#          up_down =  up_down_factor,
#          end_fitness = end_fitness_status,
#          status = infection_status))
#        
#        
#        }else{
#          
#        ###Look at the infected red blood cells time series only  
#        ###Calculating the fitness
#        G_TS <- data.frame(time = x$time,
#                           G = x[,"G"])
#        
#        #Prevent negative gamaetocytes 
#        G_TS$G[G_TS$G < 0 ] <- 0 
#        
#        truncate_G_TS <- subset(G_TS, G_TS$time <= mort_time)
#        end_fitness_mort <-  sum(PrI_PC(truncate_G_TS$G) * 1/10)
#        
#        return(data.frame(
#               peaktime = NA,
#               endtime = mort_time$time,
#               up_down = "up" ,
#               end_fitness = end_fitness_mort, 
#               status = 'mort'))
#        }
#}



###This is the main function used to figure out when to figure out
###if the infection leads to either mortality or failed/successful
### infections. If the infection kills the host, I calculate the 
###fitness up to the time of death. If the infection does not kill
###the host, I check to see if it's a failed infection. If it's a
###failed infection, the initial IRBC decreases and thus the fitness
###is automatically set to 0. If it's a successful infection, then
###the end fitness is set to NA and I note the duration of the acute 
###phase 

###################################################
###INPUT: The element of the full deSolve output###
###################################################

Finder_King <- function(x){
  
  ###Shape parameter (to ensure we're getting the right columns)
  n1 = 100
  ###Check first if it induces mortality
  ###Look at the red blood cells time series only  
  RBC_TS <- data.frame(time = x$time,
                       RBC = x[,"R"])
  
  ### Look at the red blood cells to see when it is under 6.5 * 10^5
  index_time_mort <- RBC_TS[which(RBC_TS$RBC <= 6.5 * 10^5),]
  
  ###If the index_time_mort is empty, then the host survived
  ###and we note that
  if (nrow(index_time_mort) == 0){
    mort_time <- cbind.data.frame(surv = 'surv',
                                  time = NA)
  } else {
    mort_time <- cbind.data.frame(surv = 'mort',
                                  time = index_time_mort$time[1])
  }
  
  ###This now checks the Acute Phase Duration for the Mice 
  ###That Survived
  if(mort_time$surv == "surv"){
    
    Initial_IRBC <- rowSums(x[,3:(n1+2)])[1]
    
    ###Look at the infected red blood cells time series only  
    IRBC_TS <- data.frame(time = x$time,
                          IRBC = rowSums(x[,3:(n1+2)]))
    
    IRBC_TS_Peak <- IRBC_TS[which.max(IRBC_TS$IRBC),]
    
    ### Look at the maximum IRBC peak and pull out the time     
    ### if time = 0 for the index time, then this means that IRBC
    ### goes down- we should keep a note of it 
    
    if (IRBC_TS_Peak$time == 0){  #this means the infection fails to establish
      end_time <- 0 
    }
    else{
    
    post_peak_ts <- subset(IRBC_TS,IRBC_TS$time > IRBC_TS_Peak$time)
      
    ###This simply means that the IRBC abundance is now less than the initial
    ###IRBC and I get the first time 
    end_time <-  post_peak_ts[post_peak_ts$IRBC < Initial_IRBC ,]$time[1]
    }
    
    #If the peak is 0 that suggests that the infection never establishes
    up_down_factor <- ifelse(end_time == 0, 'down', 'up')
    
    ###If the end time is 0, that means that the infection never establishes
    ###and the fitness is set to 0. If it's not 0, the host survives
    ###and we must run a different function to get the fitness
    end_fitness_status <-  ifelse(end_time == 0, 0, NA)
    
    infection_status <- ifelse(end_time == 0, 'failed', 'success')  
    
    return(data.frame(
      peaktime = IRBC_TS_Peak$time,
      endtime =  end_time ,
      up_down =  up_down_factor,
      end_fitness = end_fitness_status,
      status = infection_status))
    
    
  }else{
    
    ###Look at the infected red blood cells time series only  
    ###Calculating the fitness
    G_TS <- data.frame(time = x$time,
                       G = x[,"G"])
    
    #Prevent negative gamaetocytes 
    G_TS$G[G_TS$G < 0 ] <- 0 
    
    truncate_G_TS <- subset(G_TS, G_TS$time <= mort_time$time)
    
  
    end_fitness_mort <-  sum(PrI_PC(truncate_G_TS$G) * 1/10)
    
    return(data.frame(
      peaktime = NA,
      endtime = mort_time$time,
      up_down = "up" ,
      end_fitness = end_fitness_mort, 
      status = 'mort'))
  }
}


###This is the main function used to figure out when to figure out
###if the infection leads to either mortality or successful infection
### infections. If the infection kills the host, I calculate the 
###fitness up to the time of death. If the infection does not kill
###the host, it is a successful infection

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
        mort_time <- cbind.data.frame(surv = 'mort',
                                  time = index_time_mort$time[1])
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
      
      end_time <- subset(RM_time_df, RM_time_df$time >= min_RM$time & 
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
  
  
  
  
  
  