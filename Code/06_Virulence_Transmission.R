# HEADER --------------------------------------------
#
# Author: Damie Pak
# Copyright (c) Damie Pak, 2023
# Email:  dp495@cornell.edu
#
# Date: 2023-04-21
#
# Script Name: 06_Virulence_Transmission .R
#
# Script Description: This is the script that looks at the fitness
# output and investigate how different burst size affects the 
# virulence, the total cumulative gametocyctes, and the acute phase
# duration
#
# Notes: 
#

########################################################
### Understanding how the viable burst size combinations#
### influence virulence, transmission, and such.     ####
########################################################


# HEADER --------------------------------------------
#
# Author: Damie Pak
# Copyright (c) Damie Pak, 2023
# Email:  dp495@cornell.edu
#
# Date: 2023-04-21
#
# Script Name: 06_Virulence_Transmission .R
#
# Script Description: This is the script that looks at the fitness
# output and investigate how different burst size affects the 
# virulence, the total cumulative gametocyctes, and the acute phase
# duration
#
# Notes: 
#

########################################################
### Understanding how the viable burst size combinations#
### influence virulence, transmission, and such.     ####
########################################################

library(here)

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "02_Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "06_Virulence_Transmission_BurstSize.R"))
source(here("Code","Simulator_Code", "Simulator_Main_Gflux.R"))
source(here("Code","Helper_Function_Code", "06_Glux_Simulator.R"))
source(here("Code", "Helper_Function_Code", "Plotters", "06_Triplot_Plotter.R"))
###
Fitness_MODEL_PC_FULL_MED <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_med.csv"
))


###Run a specific simulating function that gives me the total cumualtive gametocytes
###over the course of the infection
FULL_MODEL_SIMULATOR_GFLUX_MED <- FULL_MODEL_SIMULATOR_GFLUX(Fitness_MODEL_PC_FULL_MED,
                                                              4385.96491)
Triplot_Plotter(FULL_MODEL_SIMULATOR_GFLUX_MED, 4385.96491) 


###THESE ARE FOR THE POINTS THAT I NEED###
FULL_MODEL_SIMULATOR_GFLUX_MED_ALL <- FULL_MODEL_SIMULATOR_GFLUX_ALL(Fitness_MODEL_PC_FULL_MED,4385.96491)

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

