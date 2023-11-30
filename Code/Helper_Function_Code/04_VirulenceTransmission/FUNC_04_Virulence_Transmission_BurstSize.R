# HEADER --------------------------------------------
#
# Author: Damie Pak
# Copyright (c) Damie Pak, 2023
# Email:  dp495@cornell.edu
#
# Date: 2023-04-21
#
# Script Name: Virulence_Transmission_BurstSize.R
#
# Script Description: This is the script that looks at the fitness
# output and investigate how different burst size affects the 
# virulence, the total cumulative gametocyctes, and the 
#
# Notes: 
#

###Minimum RBC finder

Virulence_Gam_Finder <- function (fitness_list){

 RMin <-  (fitness_list[which.min(fitness_list$R),]$R)

 if(unique(fitness_list$status) == 'success'){
   
   ###If the infection is successful, just find the maximum 
   
  GMax <- (fitness_list[which.max(fitness_list$Gflux),]$Gflux)

 }else{
   
   
   ###If the infection induces host mortality, truncate timeseries
   ###to before the death time and find the cumulative maximum gametocyte
   ###flux
   
   subsetted_Gflux_func <- approxfun(fitness_list$time,
                                     fitness_list$Gflux)
  
   
   GMax <- (subsetted_Gflux_func(unique(fitness_list$endtime)))
 }


 full_df <- data.frame(RMin, GMax, 
                       B_V = unique(fitness_list$B_V),
                       C_V = unique(fitness_list$C_V),
                       endtime = unique(fitness_list$endtime),
                       status = unique(fitness_list$status))
 
 return(full_df)
}





