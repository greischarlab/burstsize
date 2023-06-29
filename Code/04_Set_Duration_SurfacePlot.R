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
source(here("Code", "Helper_Function_Code", "Fitness_Function.R"))
source(here("Code", "Helper_Function_Code", "average_life_span_RC.R"))


ifelse(file.exists(here("Output",
                        "Full_Model",
                        "FULL_MODEL_100_PC_DT.csv"))== TRUE,
                        print("Available and ready to rumble!"),
                        source(here("Code","Helper_Function_Code",
                                    "SD_100Day_Infection_Creator.R")))
       

FITNESS_MODEL <- 
  as.data.frame(fread(here(
  "Output", "Fitness_Model",
  "Fitness_MODEL_PC_FULL.csv"
)))

FULL_MODEL_100 <- as.data.frame(fread(here("Output",
                                           "Full_Model",
                                           "FULL_MODEL_100_PC_DT.csv")))

###1) MERGE the full_model and the fitness model,
###2) filter only for successful infections and 
###3) split them 
FULL_MODEL_100_F <- merge(FULL_MODEL_100, FITNESS_MODEL, by =c("B_V","C_V")) %>%
                    filter(., status == 'success') %>%
                    split(., list(.$B_V,.$C_V),drop = TRUE)


###Different average lifespans of the equation: DAY 5 to DAY 100, by 5 days
average_lifespan <- c(seq(5, 100, 0.1))

###This is for the exponential weighting
exponential_RC <- average_life_span_RC(FULL_MODEL_100_F, 'exponential')
###This is for the equal weighting
equal_RC <- average_life_span_RC(FULL_MODEL_100_F, 'equal')
###Combine them all together
full_RC <- rbind.data.frame(exponential_RC, equal_RC)

###Days of interest
doi <-  c(5,25,50,75)
doi_points <- filter(full_RC, full_RC$lifespan %in% c(5,25,50,75))



###This is the figure that makes the replicative capacity versus the 
###average infection length 

New_Label_Weighting <- c("equal" = "Equal\nweighting",
                         "exponential" = "Exponential\nweighting")

RC_GG <- 
  ggplot(data = full_RC, 
         aes(x = lifespan, y= RC)) + 
  geom_line(size = 0.8) + 
  geom_vline(xintercept = c(5,25,50,75), 
             linetype = 2) +
  geom_point(data = doi_points, aes( x= lifespan, y = RC),
             size = 3 ) +
  geom_text(data = doi_points,
             aes(x = lifespan , y = RC , 
                 label = paste(expression(beta), "=",B_V,"\n","c =", C_V)),
            nudge_x = 5, nudge_y = - 0.5) + 
  facet_wrap(~id, ncol = 1, strip.position = 'right',
             labeller = labeller(id = New_Label_Weighting),
             scales = 'free_x') + 
  scale_x_continuous(limits=c(0,100)) +
  xlab("Average infection length") + 
  ylab(expression(paste("Replicative capacity (","(1-c)", beta,")")))+ 
  theme_classic() + 
  theme(
        axis.text = element_text(size = 14,color = 'black'),
        axis.title = element_text(size = 14),
        strip.background =  element_blank(),
        strip.text = element_text(size = 16),
        strip.text.y.right = element_text(angle = 0),
        panel.spacing = unit(1.5, "lines"))


ggsave(here("Figures", "Raw", "Standard_Duration_Lineplots.pdf"),
       width = 10, height = 8, units= 'in')
