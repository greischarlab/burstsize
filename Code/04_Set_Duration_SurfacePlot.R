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
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
### Packages to load
source(here("Code", "Simulator_Code", "Simulator_Main_2.R"))
source(here("Code", "Helper_Function_Code", "02_Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "average_life_span_RC.R"))
source(here("Code", "Helper_Function_Code","FUNC_04_SetDuration_Plotter.R"))
source(here("Code", "Helper_Function_Code","PLOTTER_04_TRANSSURFACEPLOT.R"))

source(here("Code", "Helper_Function_Code","04_SET_DURATION_HELPER_FUNCTION"))

ifelse(file.exists(here("Output",
                        "Full_Model",
                        "FULL_MODEL_100_PC_DT.csv"))== TRUE,
                        print("Available and ready to rumble!"),
                        source(here("Code","Helper_Function_Code",
                                    "SD_100Day_Infection_Creator.R")))


###WE ARE INTERESTED IN THE MEDIUM INOCULUM

Fitness_MODEL_PC_FULL_MED <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_med.csv"
))

###USE THIS FOR THE SURFACE PLOT OF THE TRANSMISSION PROBABILITY
FULL_MODEL_100_MED <- as.data.frame(fread(here("Output",
                                               "Full_Model",
                                               "FULL_MODEL_100_PC_DTMed.csv")))

###USE THIS FOR Optimal burst size information that is overlaid on the 
###Surface plot
FULL_MODEL_100_MED_NO_DEATH <- as.data.frame(fread(here("Output",
                                                        "Full_Model",
                                                        "FULL_MODEL_100_PC_DTMed0.75.csv")))
#############################################################################
###This creates the necessary data.frame to then make the data-visualization#
#############################################################################
Opt_Burst_MED <- Optimalburstsize_SetDuration_Finder(FULL_MODEL_100_MED_NO_DEATH,
                                                     Fitness_MODEL_PC_FULL_MED, 
                                                     0.75)



PRIPC_MED_75 <-Trans_SurfacePlot_Maker(x = FULL_MODEL_100_MED,
                                       y = Fitness_MODEL_PC_FULL_MED, 
                                       CV_vec = c(0.75),
                                       4385.96491)

PRIPC_MED_75 $time <- round(PRIPC_MED_75 $time,3)
PRIPC_MED_75 $pripc <- round(PRIPC_MED_75 $pripc,2)


PRIPC_MED_75$cutprob <- cut(PRIPC_MED_75$pripc, breaks=c(-1,0,0.05,0.1,.15,.2,.25,.3,.35,.4,
                                                         .45,.50,.55,.60,.65,.70,.75,.80,
                                                         .85,.90,.95,1))



Set_Duration_Plotter(PRIPC_MED_75,Opt_Burst_MED[[3]])

ggsave(here("Figures", "Raw", "TransProb_SurfacePlot.pdf"),
       width = 7, height =4.5, unit = "in")

