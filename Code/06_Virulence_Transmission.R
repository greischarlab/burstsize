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
source(here("Code", "Helper_Function_Code", "Virulence_Transmission_BurstSize.R"))
source(here("Code","Simulator_Code", "Simulator_Main_Gflux.R"))
###

Fitness_MODEL_PC_FULL_LOW_SUPP <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_43.58965.csv"
))

Fitness_MODEL_PC_FULL_MED <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_4358.965.csv"
))

Fitness_MODEL_PC_HIGH_SUPP <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_43589.65.csv"
))

FULL_MODEL_SIMULATOR_GFLUX_LOW <- FULL_MODEL_SIMULATOR_GFLUX (Fitness_MODEL_PC_FULL_LOW_SUPP,
                                                              43.58965)
Triplot_Plotter(FULL_MODEL_SIMULATOR_GFLUX_LOW , 43.58965)

FULL_MODEL_SIMULATOR_GFLUX_MED <- FULL_MODEL_SIMULATOR_GFLUX (Fitness_MODEL_PC_FULL_MED,
                                                              4358.965)
Triplot_Plotter(FULL_MODEL_SIMULATOR_GFLUX_MED,      4358.965)

FULL_MODEL_SIMULATOR_GFLUX_HIGH <- FULL_MODEL_SIMULATOR_GFLUX (Fitness_MODEL_PC_HIGH_SUPP,
                                                              43589.6)
Triplot_Plotter(FULL_MODEL_SIMULATOR_GFLUX_HIGH , 43589.6)
