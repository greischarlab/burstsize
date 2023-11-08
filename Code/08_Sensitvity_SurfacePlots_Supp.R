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


Opt_Med_Strategy <- data.frame(B_V = 14.5, C_V =  0.75)
Long_Med_Strategy <- data.frame(B_V = 4, C_V =  0.09)

###Lambda



grapher_SA (FITNESS_L,Opt_Med_Strategy,Long_Med_Strategy)
by(FITNESS_L, FITNESS_L$change,  Best_Strategy_Finder)
by(FITNESS_L, FITNESS_L$change,  Longest_Finder )

ggsave(here("Figures", "Raw", "SA_L_Duration_Fitness_SurfacePlot.pdf"),
       width = 10, height =10, unit = 'in')

###Merozoite Mortality
grapher_SA (FITNESS_UM,Opt_Med_Strategy,Long_Med_Strategy)
by(FITNESS_UM, FITNESS_UM$change,  Best_Strategy_Finder)
by(FITNESS_UM, FITNESS_UM$change,  Longest_Finder )

ggsave(here("Figures", "Raw", "SA_UM_Duration_Fitness_SurfacePlot.pdf"),
       width = 10, height =10, unit = 'in')


grapher_SA (FITNESS_RI,Opt_Med_Strategy,Long_Med_Strategy)
by(FITNESS_RI, FITNESS_RI$change,  Best_Strategy_Finder)
by(FITNESS_UM, FITNESS_RI$change,  Longest_Finder )
ggsave(here("Figures", "Raw", "SA_RI_Duration_Fitness_SurfacePlot.pdf"),
       width = 10, height =10, unit = 'in')



