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
###Output: The maximum cumulative transmission potential                #
#########################################################################  
Gametocyte_Fitness <- function(x){
  
  gam_dat <- data.frame(time = x['time'],
                        gam = x[,"G"])
  
  ###If it's negative then, make it 0 just in case
  gam_dat$gam[gam_dat$gam < 0] <- 0            
  
  end_fitness <- max(cumsum(PrI_PC(gam_dat$gam) * 1/10))
  
  return(end_fitness)
}
############################################################
###Truncator################################################
Gametocyte_Fitness_EqualWeighting <- function(x, life_span){
  
  gam_dat <- data.frame(time = x['time'],
                        gam = x[,"G"])
  
  ###If it's negative then, make it 0 just in case
  gam_dat$gam[gam_dat$gam < 0] <- 0        
  
  gam_dat_trunc <- subset(gam_dat, gam_dat$time <= life_span)
  
  end_fitness <- sum(PrI_PC(gam_dat_trunc) * 1/10)
  
  return(end_fitness)
}
###


###################################################
###INPUT: The element of the full deSolve output###
### as well as the mu_M                         ###
###################################################
Finder_RM <- function(x_list, mu_M_c) {
  
    ###If it's NA that means that it never kills the host,
    ###If it's a numeric value, that means it kills the host.
    infection_length <- unique(x_list$infection_length)
  
    ###If the infection_length is NA...
    if(is.null(infection_length) == FALSE){
    
      ###Look at the gametocyte time series 
      ###Calculating the fitness
      G_TS <- data.frame(time = x_list[,'time'],
                         G = x_list[,"G"])
      
      #Prevent negative gametocytes if there are any
      G_TS$G[G_TS$G < 0 ] <- 0 
      
      ###Let's approximate when the gametocyte time series hit
      ###The time...
      G_TS_Function <- approxfun(G_TS)
      
      G_TS_Function(seq(0,infection_length, 1/10))
      
      truncate_G_TS <- G_TS_Function(seq(0,infection_length, 1/10))
      
      end_fitness_mort <-  sum(PrI_PC(truncate_G_TS) * 1/10)
      
      return(data.frame(
        endtime = infection_length,
        up_down = "up" ,
        end_fitness = end_fitness_mort, 
        status = 'mort'))
      
    }
    else
      {
    
     p =  4.0e-6
     unique_B_V <- unique(x_list$B_V)
     unique_C_V <- unique(x_list$C_V)
     
     rate = (1 - unique_C_V) * unique_B_V *
                ((x_list[,"R"] * p)/((p * x_list[,"R"]) + mu_M_c))
    
      RM_time_df <-  cbind.data.frame(time = x_list[,'time'],
                                      rate = rate, 
                                      B_V = unique_B_V ,
                                      C_V = unique_C_V)
      
      min_RM <- RM_time_df[which.min(RM_time_df$rate),]
      
      end_time <- subset(RM_time_df, 
                            RM_time_df$time >= min_RM$time & 
                            RM_time_df$rate >= 1)[1,'time']
      
    
    return(data.frame(
      endtime =  end_time,
      up_down =  'up',
      end_fitness = NA,
      status = 'success'))
    
    
    }
}
 

###For finding the durations across different lists

Duration_Finder <- function(x_list, mu_M_c){
  
  tmp <- do.call(
    rbind,
    mclapply(x_list,
             Finder_RM,
             mu_M_c = mu_M_c,
             mc.cores = 2))
  
  return(tmp)
}


 



###For the sensitivity analysis finding the Finder
Finder_RM_SA <- function(x_list, sens_var) {
  
  ###Seperate out the list element by the minus or plus change
  tmp_x_list <- split(x_list, x_list$change)
  
  if (sens_var == "UM") { ###If the sens_var is UM then we need to account
    ###for the difference in mortality rates
    mu_M_c <- c(36,60)
  } else {                ###If the sens_var is not UM then we need to account
    mu_M_c <- 48
  }
  
  ### Shape parameter (to ensure we're getting the right columns)
  n1 <- 100
  
  ###################################################
  ### Check first if the infection induces mortality#
  ### Look at the red blood cells time series only  #
  ###################################################
  RBC_TS <- lapply(tmp_x_list, function(x) {
    data.frame(
      time = x$time,
      RBC = x[, "R"]
    )
  })
  ###################################################################
  ### Look at the red blood cells to see when it is under 6.5 * 10^5#
  ###################################################################
  index_time_mort <- lapply(
    RBC_TS,
    function(x) {
      x[which(x$RBC   <= 6.5 * 10^5), ]
    }
  )
  
  ### If the index_time_mort is empty, then the host survived
  ### and we note that
  mort_list <- NULL
  
  for (k in seq(1, 2)) {
    tmp <- index_time_mort[[k]]
    
    if (nrow(tmp) == 0) {
      mort_time <- cbind.data.frame(
        surv = "surv",
        time = NA
      )
      mort_list[[k]] <- mort_time
    } else {
      RBC_TS_Func <- approxfun(RBC_TS[[k]][, "time"], 
                               RBC_TS[[k]][, "RBC"] - 6.5 * 10^5)
      death_time <- uniroot(RBC_TS_Func, c(0, tmp$time[1]))
      mort_time <- cbind.data.frame(
        surv = "mort",
        time = death_time$root
      )
      mort_list[[k]] <- mort_time
    }
  }
  ##########################################################
  ### This now checks the Acute Phase Duration for the Mice#
  ### That Survived                                        #
  ##########################################################
  
  end_time_list <- NULL
  
  for (k in seq(1, 2)){
    
    tmp <- mort_list[[k]]
    tmp_x_list_df <- tmp_x_list[[k]]
    
    ###If the tmp leads to survival
    if (tmp$surv == "surv") {
      p <- 2.5e-6 #this is a parameter value that does not change 
      
      rate <- (1 - unique(tmp_x_list_df$C_V)) *
        unique(tmp_x_list_df$B_V) * ((tmp_x_list_df[, "R"] * p) / ((p * tmp_x_list_df[, "R"]) +
                                                                     mu_M_c))
      
      RM_time_df <- cbind.data.frame(
        time = tmp_x_list_df[, "time"],
        rate = rate,
        B_V = unique(tmp_x_list_df$B_V),
        C_V = unique(tmp_x_list_df$C_V)
      )
      
      min_RM <- RM_time_df[which.min(RM_time_df$rate), ]
      
      end_time <- subset(
        RM_time_df,
        RM_time_df$time >= min_RM$time &
          RM_time_df$rate >= 1
      )[1, "time"] ###Endtime is the first time when the 
      
      end_time_list[[k]] <-
        data.frame(
          endtime = end_time,
          up_down = "up",
          end_fitness = NA,
          status = "success",
          B_V = unique(tmp_x_list_df$B_V),
          C_V = unique(tmp_x_list_df$C_V),
          change = unique(tmp_x_list_df$change)
        )
    } else {
      ### Look at the infected red blood cells time series only
      ### Calculating the fitness
      G_TS <- data.frame(
        time = tmp_x_list_df[, "time"],
        G = tmp_x_list_df[, "G"]
      )
      
      # Prevent negative gamaetocytes
      G_TS$G[G_TS$G < 0] <- 0
      
      truncate_G_TS <- subset(G_TS, G_TS$time <=  mort_list[[k]]$time)
      
      
      end_fitness_mort <- sum(PrI_PC(truncate_G_TS$G) * 1 / 10)
      
      end_time_list[[k]] <- 
        data.frame(
          endtime = tmp$time,
          up_down = "up",
          end_fitness = end_fitness_mort,
          status = "mort",
          B_V = unique(tmp_x_list_df$B_V),
          C_V = unique(tmp_x_list_df$C_V),
          change = unique(tmp_x_list_df$change)
          
        )
    }
  }
  return(do.call(rbind, end_time_list))
  
}
