# HEADER --------------------------------------------
#
# Author: Damie Pak
# Copyright (c) Damie Pak, 2023
# Email:  dp495@cornell.edu
#
# Date: 2023-04-21
#
# Script Name: 04_Set_Duration_SurfacePlot.R
#
# Script Description: This is a script that looks at the cumulative
# transmission potential at specific end points.
#
# Notes:
#
###################
### Set durations##
###################
library(here)
library(magrittr)
### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "02_Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "average_life_span_RC.R"))
source(here("Code", "Helper_Function_Code", "Grapher_vert_hor.R"))

ifelse(file.exists(here("Output",
                        "Full_Model",
                        "FULL_MODEL_100_PC_DT.csv"))== TRUE,
                        print("Available and ready to rumble!"),
                        source(here("Code","Helper_Function_Code",
                                    "SD_100Day_Infection_Creator.R")))
  
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


FULL_MODEL_100_LOW <- as.data.frame(fread(here("Output",
                                           "Full_Model",
                                           "FULL_MODEL_100_PC_DTLow.csv")))

FULL_MODEL_100_MED <- as.data.frame(fread(here("Output",
                                               "Full_Model",
                                               "FULL_MODEL_100_PC_DTMed.csv")))

FULL_MODEL_100_HIGH <- as.data.frame(fread(here("Output",
                                               "Full_Model",
                                               "FULL_MODEL_100_PC_DTHigh.csv")))

###This creates the necessary data.frame to then make the data-visualization
Opt_Burst_MED<- Optimalburstsize_SetDuration_Finder(FULL_MODEL_100_MED ,Fitness_MODEL_PC_FULL_MED, 0.75)

###

Opt_Burst_MED_GG <- Optimalburstsize_SetDuration_Lineplotter(Opt_Burst_MED[[3]])



TRAJ_PLOT_GG(Fitness_MODEL_PC_FULL_MED,RC_MED[[3]]) +ggtitle('hi') coord_cartesian(xlim =c(25,50),
                                                                      ylim =c(0.75,1))
=
