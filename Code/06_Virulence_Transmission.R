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
source(here("Code","Helper_Function_Code", "06_Glux_Simulator.R"))
source(here("Code", "Helper_Function_Code", "Plotters", "06_Triplot_Plotter.R"))
###

Fitness_MODEL_PC_FULL_LOW_SUPP <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_low.csv"
))

Fitness_MODEL_PC_FULL_MED <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_med.csv"
))

Fitness_MODEL_PC_HIGH_SUPP <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_high.csv"
))

FULL_MODEL_SIMULATOR_GFLUX_LOW <- FULL_MODEL_SIMULATOR_GFLUX(Fitness_MODEL_PC_FULL_LOW_SUPP,
                                                              43.85965)
Triplot_Plotter(FULL_MODEL_SIMULATOR_GFLUX_LOW , 43.58965) + ggtitle("Low inoculation")

FULL_MODEL_SIMULATOR_GFLUX_MED <- FULL_MODEL_SIMULATOR_GFLUX (Fitness_MODEL_PC_FULL_MED,
                                                              4385.96491)
Triplot_Plotter(FULL_MODEL_SIMULATOR_GFLUX_MED, 4385.96491) 

FULL_MODEL_SIMULATOR_GFLUX_HIGH <- FULL_MODEL_SIMULATOR_GFLUX (Fitness_MODEL_PC_HIGH_SUPP,
                                                               438596.49123)
Triplot_Plotter(FULL_MODEL_SIMULATOR_GFLUX_HIGH , 438596.49123)


###

FULL_MODEL_SIMULATOR_GFLUX_MED <- FULL_MODEL_SIMULATOR_GFLUX (Fitness_MODEL_PC_FULL_MED,
                                                              4385.96491)
Triplot_Plotter(FULL_MODEL_SIMULATOR_GFLUX_MED, 4385.96491) 

