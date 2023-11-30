# HEADER --------------------------------------------
#
# Author: Damie Pak
# Copyright (c) Damie Pak, 2023
# Email:  dp495@cornell.edu
#
# Date: 2023-04-21
#
# Script Name: 06_Full_Model_simulator_Gflux.R
#
# Script Description: This is the script
# that looks at the optimal transmission investment and
# look at all the burst sizes and spit out the necessary
# data file to do more analysis on 
#
# Notes:
#
###################
### Set durations##
###################

FULL_MODEL_SIMULATOR_GFLUX <- function(x,initial_value){
  
  ###find the optimal transmission investment
  optimum_strategy <- x[which.max(x$end_fitness),]
  
  optimum_CV <- optimum_strategy$C_V
  optimum_BV <- optimum_strategy$B_V
### burst sizes that successfully establish infection (regardless of
###rather it kills the host or not)
  BV_OptimalCV <- subset(x,
                  x$C_V ==  optimum_CV &
                  x$status != 'Fail')

Burstsize_OptimalCV <- BV_OptimalCV$B_V
Endtime_OptimalCV <- BV_OptimalCV$endtime

### We need to rerun the Fitness_Model for the Gflux
Fitness_MODEL_PC_GFlux <- mcmapply(Simulator_MalariaPC_DDE_BC_GFLUX,
                                   c(Burstsize_OptimalCV),
                                   c(optimum_CV),
                                   c(initial_value),
                                   c(Endtime_OptimalCV),
                                   mc.cores = 1,
                                   SIMPLIFY = FALSE)



for (strain in seq(1,nrow(BV_OptimalCV))){
  Fitness_MODEL_PC_GFlux[[strain]]$status <- BV_OptimalCV$status[[strain]]
}

Gflux_DF <- do.call(rbind,
                    lapply(Fitness_MODEL_PC_GFlux, Virulence_Gam_Finder))


Gflux_DF$optimal <- ifelse(Gflux_DF$B_V == optimum_BV , 1,NA)


return(Gflux_DF)
}



FULL_MODEL_SIMULATOR_GFLUX_ALL <- function(x, initial_value){

  BV_CV_Sucess <- subset(x,x$status != 'Fail')
  
    
  ### We need to rerun the Fitness_Model for the Gflux
  Fitness_MODEL_PC_GFlux <- mcmapply(Simulator_MalariaPC_DDE_BC_GFLUX,
                                     c(BV_CV_Sucess$B_V),
                                     c(BV_CV_Sucess$C_V),
                                     c(initial_value),
                                     c(BV_CV_Sucess$endtime),
                                      mc.cores = 3,
                                      SIMPLIFY = FALSE)
  
  
  
  for (strain in seq(1,nrow(BV_CV_Sucess))){
    Fitness_MODEL_PC_GFlux[[strain]]$status <- BV_CV_Sucess$status[[strain]]
  }
  
  Gflux_DF <- do.call(rbind,
                      lapply(Fitness_MODEL_PC_GFlux, Virulence_Gam_Finder))
  
  
  
  
  return(Gflux_DF)
}
