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
       
Fitness_MODEL_PC_FULL <- 
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
FULL_MODEL_100_F <- merge(FULL_MODEL_100, 
                          Fitness_MODEL_PC_FULL, by =c("B_V","C_V")) %>%
                           filter(., status == 'success') %>%
                          split(., list(.$B_V,.$C_V),drop = TRUE)
  

###Different average lifespans of the equation: DAY 5 to DAY 100, by 5 days
average_lifespan <- c(seq(5, 100,1))

###This is for the equal weighting
equal_RC <- average_life_span_RC(FULL_MODEL_100_F, 
                                 average_lifespan )

###Days of interest
doi <-  c(5,20,25,50,75)
doi_points <- filter(equal_RC, equal_RC$lifespan %in% c(5,20,25,50,75))

equal_RC_nomin <- subset(equal_RC,equal_RC$group == 'fit' &
                           equal_RC$lifespan %in% c(5,20,25,30,40,50))

###This is the figure that makes the replicative capacity versus the 
###average infection length 

RC_GG <- 
  ggplot(data = equal_RC,
         aes(x = lifespan, y= RC, color = group, group = group))+
  geom_line(size = 1) + 
  geom_point(data = subset(doi_points,
                           doi_points$group != 'min'),
             aes( x= lifespan, y = RC),
             size = 3 ) +
  scale_color_manual(values = c('fit' = 'black', 'max'='green'))+
  scale_x_continuous(limits=c(0,100)) +
  xlab("Total infection length") + 
  ylab(expression(paste("Replicative capacity (","(1-c)", beta,")")))+ 
  theme_classic() + 
  theme(
        axis.text = element_text(size = 14,color = 'black'),
        axis.title = element_text(size = 14),
        strip.background =  element_blank(),
        strip.text = element_text(size = 16),
        strip.text.y.right = element_text(angle = 0),
        panel.spacing = unit(1.5, "lines"), 
        legend.position = 'none')



###
######################################################
### Makes Figure 2, which is the surface plot of the##
### acute phase duration and the fitness
####################################################


###Produces the mortality and non-establishing infection lines
horizontal_vert_df <- grapher_mortality_boundary(Fitness_MODEL_PC_FULL)


traj_plot_GG <- 
  ggplot(equal_RC_nomin, 
         aes(x = B_V, y = C_V, 
        color = as.factor(lifespan))) +
  geom_text(aes(x= B_V, y= C_V, label = lifespan))+
  geom_point(size = 3.0)+
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL,
      Fitness_MODEL_PC_FULL$status ==
        "Fail"
    ),
    aes(x = B_V, y = C_V), fill = "#d1dbe8",
    inherit.aes = FALSE
    
  ) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL,
      Fitness_MODEL_PC_FULL$status ==
        "mort"
    ),
    aes(x = B_V, y = C_V), fill = "black",
    inherit.aes = FALSE
    
  )+
  geom_path(size=1,lineend = "round", group =1) +
  
  scale_color_viridis(option='plasma',discrete = TRUE,)+
  new_scale_color() +
 
  geom_segment(
    data = horizontal_vert_df,
    aes(
      x = x, xend = xend,
      y = y, yend = yend,
    color = id), size = 0.5,
    
    lineend = "round",
    inherit.aes = FALSE
  )+
  scale_color_manual(
    values = c(
      "fail" = "black",
      "mort" = "#b8fcd5"
    ),
    guide = "none"
  )+
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Transmission investment (c)") +
 theme(panel.background=element_rect(fill ='#e3e0da'),
       panel.border = element_rect(color = 'black', fill = NA, size =1),
       panel.grid  = element_blank(),
     
         text = element_text(size = 14, color = "black"),
         axis.text = element_text(size = 14, color = "black"),
         axis.title = element_text(size = 15, color = "black"),
         legend.position = "none"
       ) +
  coord_cartesian(xlim =c(20,50), ylim = c(0.7,1))


traj_plot_GG +
  RC_GG +  plot_annotation(tag_levels = "A")


ggsave(here("Figures", "Raw", "Standard_Duration_Plots_Alt.pdf"),
       height = 4, width = 10,
       units = 'in')
