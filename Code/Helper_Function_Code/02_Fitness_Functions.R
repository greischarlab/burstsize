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
  
  end_fitness <- sum(PrI_PC(gam_dat_trunc$gam) * 1/10)
  
  return(cbind.data.frame(end_fitness, life_span, 
                          B_V = unique(x$B_V), 
                          C_V = unique(x$C_V)))
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
  
    ###If the infection_length is NA.
    if(is.na(infection_length) == FALSE){
    
      ###Look at the gametocyte time series 
      ###Calculating the fitness
      G_TS <- data.frame(time = x_list[,'time'],
                         G = x_list[,"G"])
      
      #Prevent negative gametocytes if there are any
      G_TS$G[G_TS$G < 0 ] <- 0 
      
      ###Let's approximate when the gametocyte time series hit
      ###The time...
      G_TS_Function <- approxfun(G_TS)
    
      truncate_G_TS <- G_TS_Function(seq(0,infection_length, 1/10))
      
      end_fitness_mort <-  sum(PrI_PC(truncate_G_TS) * 1/10)
      
      df<- data.frame(
        endtime = infection_length,
        up_down = "up" ,
        end_fitness = end_fitness_mort, 
        status = 'mort')
      
      }else{
    
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
      
    
    df <- data.frame(
      endtime =  end_time,
      up_down =  'up',
      end_fitness = NA,
      status = 'success')

      }
    
    return(df)
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


