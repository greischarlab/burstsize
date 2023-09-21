# HEADER --------------------------------------------
#
# Author: Damie Pak
# Copyright (c) Damie Pak, 2023
# Email:  dp495@cornell.edu
#
# Date: 2023-04-21
#
# Script Name: 03_Main_Analysis_SurfacePlot.R
#
# Script Description: This is the script that takes the output
#from the previous analysis and produce the first surface plot 
#
# Notes: Originally, had the fitness and duration plot together,
# but changed it so that only the fitness is in main and the
# duration is in the supp.
#

### Packages to load
library(here)
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "Grapher_vert_hor.R"))
source(here("Code", "Helper_Function_Code","Plotters", "03_SurfacePlot_Maker.R"))

######################################################
### Makes Figure 2, which is the surface plot of the##
### acute phase duration and the fitness
####################################################
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


LOW_GG_SurfacePlot <- MAIN_SURFACEPLOT_GG_GRAPHER_FIT(Fitness_MODEL_PC_FULL_LOW_SUPP,"Low")
MED_GG_SurfacePlot <- MAIN_SURFACEPLOT_GG_GRAPHER_FIT(Fitness_MODEL_PC_FULL_MED ,"Med")
HIGH_GG_SurfacePlot <- MAIN_SURFACEPLOT_GG_GRAPHER_FIT(Fitness_MODEL_PC_HIGH_SUPP  ,"High")

###IF you want to see the 

Best_Strategy_Finder(Fitness_MODEL_PC_FULL_LOW_SUPP)
Best_Strategy_Finder(Fitness_MODEL_PC_FULL_MED )
Best_Strategy_Finder(Fitness_MODEL_PC_HIGH_SUPP)

###I think I want the
MED_GG_SurfacePlot

ggsave(here("Figures", "Raw", "MED_SurfacePlots.pdf"),
       width = 8.5, height =8.5, unit = "in")

LOW_GG_SurfacePlot  + HIGH_GG_SurfacePlot 

ggsave(here("Figures", "Raw", "Supp_Low_High_SurfacePlots.pdf"),
       width = 11, height =5.5, unit = "in")


###DURATION SUPP PLOTS
LOW_GG_SurfacePlot_D<- MAIN_SURFACEPLOT_GG_GRAPHER_DURATION(Fitness_MODEL_PC_FULL_LOW_SUPP,"Low")
MED_GG_SurfacePlot_D<- MAIN_SURFACEPLOT_GG_GRAPHER_DURATION(Fitness_MODEL_PC_FULL_MED ,"Med")
HIGH_GG_SurfacePlot_D<- MAIN_SURFACEPLOT_GG_GRAPHER_DURATION(Fitness_MODEL_PC_HIGH_SUPP  ,"High")


Longest_Finder(Fitness_MODEL_PC_FULL_LOW_SUPP)
Longest_Finder(Fitness_MODEL_PC_FULL_MED )
Longest_Finder(Fitness_MODEL_PC_HIGH_SUPP)

LOW_GG_SurfacePlot_D + MED_GG_SurfacePlot_D + HIGH_GG_SurfacePlot_D

ggsave(here("Figures", "Raw", "Supp_Duration_SurfacePlots.pdf"),
       width = 12, height =5.5, unit = "in")
