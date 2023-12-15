
# HEADER --------------------------------------------
#
# Author: Damie Pak
# Copyright (c) Damie Pak, 2023
# Email:  dp495@cornell.edu
#
# Date: 2023-11-29
#
# Script Name: FUNC_03_Threepointstrategy .R
#
# Script Description: This is the script that identifies 
#three strategies on the varying duration surface plot.
#The most virulent strategy (lowest RBC), the greatest
#amount of gametocyte, and the longest acute phase
#
# Notes: 
#

library(here)
### Packages to load
source(here("Code", "Helper_Function_Code", "FUNC_00_Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "FUNC_00_Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "04_VirulenceTransmission", "FUNC_04_Virulence_Transmission_BurstSize.R"))
source(here("Code","Simulator_Code", "Simulator_Main_Gflux.R"))
source(here("Code","Helper_Function_Code","04_VirulenceTransmission", "FUNC_04_Glux_Simulator.R"))
###
Fitness_MODEL_PC_FULL_MED <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_med.csv"
))

###THESE ARE FOR THE POINTS THAT I NEED###
FULL_MODEL_SIMULATOR_GFLUX_MED_ALL <- FULL_MODEL_SIMULATOR_GFLUX_ALL(
                                      Fitness_MODEL_PC_FULL_MED,4385.96491)


###Lowest RBC
surviving_fit <- subset(FULL_MODEL_SIMULATOR_GFLUX_MED_ALL ,FULL_MODEL_SIMULATOR_GFLUX_MED_ALL $status != 'mort')
LowestR <-surviving_fit[which.min(surviving_fit $R),]
###Greatest cumulative gametocytes
GreatestGmax <- FULL_MODEL_SIMULATOR_GFLUX_MED_ALL[which.max(FULL_MODEL_SIMULATOR_GFLUX_MED_ALL $GMax),]
###Longest acute phase
LongestAcute <- FULL_MODEL_SIMULATOR_GFLUX_MED_ALL[which.max(FULL_MODEL_SIMULATOR_GFLUX_MED_ALL $endtime),]


Three_Point_Strategies_Surfaceplot <- rbind(LowestR,
                                            GreatestGmax, 
                                            LongestAcute)
Three_Point_Strategies_Surfaceplot$id <- c("RBC", "Gmax","Longest")

