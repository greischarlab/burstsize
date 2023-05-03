########################################################
### Understanding how the viable burst size combinations#
### influence virulence, transmission, and such.     ####
########################################################

library(here)

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "Surfaceplot_SA_Grapher.R"))
source(here("Code", "Helper_Function_Code", "Grapher_vert_hor.R"))
###We need the data here 

FITNESS_L <- as.data.frame(fread(here("Output","Fitness_Model",
                                           "Fitness_L.csv")))
FITNESS_RI <- as.data.frame(fread(here("Output","Fitness_Model",
                                          "Fitness_RI.csv")))                     
FITNESS_UM <- as.data.frame(fread(here("Output","Fitness_Model",
                                           "Fitness_UM.csv")))       

###Lambda
grapher_SA (FITNESS_L)
ggsave(here("Figures", "Raw", "SA_L_Duration_Fitness_SurfacePlot.pdf"),
       width = 14, height =10, unit = 'in')

###Merozoite Mortality
grapher_SA (FITNESS_UM)
ggsave(here("Figures", "Raw", "SA_UM_Duration_Fitness_SurfacePlot.pdf"),
       width = 14, height =10, unit = 'in')


grapher_SA (FITNESS_RI)
ggsave(here("Figures", "Raw", "SA_RI_Duration_Fitness_SurfacePlot.pdf"),
       width = 14, height =10, unit = 'in')



