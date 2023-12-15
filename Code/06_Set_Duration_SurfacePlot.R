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
source(here("Code", "Helper_Function_Code", "FUNC_00_Packages_Loader.R"))
### Packages to load
source(here("Code", "Simulator_Code", "Simulator_Main_2.R"))
source(here("Code", "Helper_Function_Code", "FUNC_00_Fitness_Functions.R"))
###Specific to 04
source(here("Code", "Helper_Function_Code", "04_SETDURATION","FUNC_04_average_life_span_RC.R"))
source(here("Code", "Helper_Function_Code","04_SETDURATION","PLOTTER_04_TRANSSURFACEPLOT.R"))
source(here("Code", "Helper_Function_Code","04_SETDURATION","FUNC_04_average_life_span_RC.R"))
source(here("Code", "Helper_Function_Code","04_SETDURATION","FUNC_04_OptB_Finder_SetDuration.R"))

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

###USE THIS FOR Optimal burst size information that is overlaid on the 
###Surface plot
FULL_MODEL_100_MED_NO_DEATH <- as.data.frame(fread(here("Output",
                                                        "Full_Model",
                                                        "FULL_MODEL_100_PC_DTMed0.76.csv")))
#############################################################################
###This creates the necessary data.frame to then make the data-visualization#
#############################################################################
Opt_Burst_MED <- Optimalburstsize_SetDuration_Finder(FULL_MODEL_100_MED_NO_DEATH,
                                                     Fitness_MODEL_PC_FULL_MED, 
                                                     0.76)



PRIPC_MED <-Trans_SurfacePlot_Maker(x = FULL_MODEL_100_MED_NO_DEATH,
                                       y = Fitness_MODEL_PC_FULL_MED, 
                                       CV_vec = c(0.76),
                                       4385.96491)

PRIPC_MED $time <- round(PRIPC_MED$time,3)
PRIPC_MED $pripc <- round(PRIPC_MED $pripc,2)


PRIPC_MED$cutprob <- cut(PRIPC_MED$pripc, breaks=c(-1,0,0.05,0.1,.15,.2,.25,.3,.35,.4,
                                                         .45,.50,.55,.60,.65,.70,.75,.80,
                                                         .85,.90,.95,1))



Set_Duration_Plotter(PRIPC_MED,Opt_Burst_MED[[3]],Opt_Burst_MED[[2]])

ggsave(here("Figures", "TransProb_SurfacePlot.pdf"),
       width = 7, height =4.5, unit = "in")

