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
source(here("Code", "Helper_Function_Code", "FUNC_00_Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "FUNC_00_Grapher_vert_hor.R"))
source(here("Code", "Helper_Function_Code", "FUNC_00_best_long_strategyfinder.R"))
source(here("Code", "Helper_Function_Code", "FUNC_00_Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code","03_VD_SURFACEPLOT","FUNC_03_StrategyLocators.R"))
source(here("Code", "Helper_Function_Code","03_VD_SURFACEPLOT","PLOTTER_03_SurfacePlot_Maker.R"))

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
HIGH_GG_SurfacePlot <- MAIN_SURFACEPLOT_GG_GRAPHER_FIT(Fitness_MODEL_PC_HIGH_SUPP,"High")


###IF you want to see the optimal strategy (B_V/C_V)
Best_Strategy_Finder(Fitness_MODEL_PC_FULL_LOW_SUPP); #46.5/0.92
Best_Strategy_Finder(Fitness_MODEL_PC_FULL_MED ) #15.5, 0.76
Best_Strategy_Finder(Fitness_MODEL_PC_HIGH_SUPP) #7, 0.47


###THIS IS THE MAIN FIGURE 03!
MED_GG_SurfacePlot + geom_point(data = Three_Point_Strategies_Surfaceplot,
                                aes(x =B_V, y=C_V), size = 2,color = 'white')

ggsave(here("Figures", "Raw", "03_MED_2_SurfacePlots.pdf"),
       width = 8.5, height =8.5, unit = "in")

###This is the Suppelmentary Figure 03!
LOW_GG_SurfacePlot  + HIGH_GG_SurfacePlot 
ggsave(here("Figures", "Raw", "SUPP_Low_High_SurfacePlots.pdf"),
       width = 11, height =5.5, unit = "in")


##########################
###DURATION SUPP PLOTS####
##########################

###These surface plots are for how long the acute phase is for 
###the different innocula 

LOW_GG_SurfacePlot_D<- MAIN_SURFACEPLOT_GG_GRAPHER_DURATION(Fitness_MODEL_PC_FULL_LOW_SUPP,"Low")
MED_GG_SurfacePlot_D<- MAIN_SURFACEPLOT_GG_GRAPHER_DURATION(Fitness_MODEL_PC_FULL_MED ,"Med")
HIGH_GG_SurfacePlot_D<- MAIN_SURFACEPLOT_GG_GRAPHER_DURATION(Fitness_MODEL_PC_HIGH_SUPP  ,"High")

###IF you want to see the longest strategy (B_V/C_V)
Longest_Finder(Fitness_MODEL_PC_FULL_LOW_SUPP) #4/0.07
Longest_Finder(Fitness_MODEL_PC_FULL_MED ) #4/0.07
Longest_Finder(Fitness_MODEL_PC_HIGH_SUPP) #4/0.07

LOW_GG_SurfacePlot_D + MED_GG_SurfacePlot_D + HIGH_GG_SurfacePlot_D 

ggsave(here("Figures", "Raw", "SUPP_Duration_SurfacePlots.pdf"),
       width = 16, height =5.5, unit = "in")
