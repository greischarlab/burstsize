###################
### Set durations##
###################
library(here)
library(magrittr)
### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "average_life_span_RC.R"))
source(here("Code", "Helper_Function_Code", "Grapher_vert_hor.R"))

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
FULL_MODEL_100_F <- merge(FULL_MODEL_100, 
                          FITNESS_MODEL, by =c("B_V","C_V")) %>%
  filter(., status == 'success') %>%
  split(., list(.$B_V,.$C_V),drop = TRUE)


Surface_Plot_Standard_Duration_20 <- Surface_Plot_Standard_Duration(FULL_MODEL_100_F ,20)

optimal <- Surface_Plot_Standard_Duration_20[which.max(Surface_Plot_Standard_Duration_20$end_fitness),]


horizontal_vert_df <- grapher_mortality_boundary(FITNESS_MODEL)



ggplot(Surface_Plot_Standard_Duration_20, aes(x = B_V, y= C_V, fill = end_fitness))+
  geom_raster()+
  geom_point(data = optimal, aes(x= B_V, y= C_V))+ 
  annotate("text",  x=23, y= 0.79, label = 17.82)+
  annotate("point", x = 8,y =0.42,col='red')+
  annotate("text", x = 8, y=0.44, label =  17.20)+
  annotate("point", x= 10, y = .42, col='blue')+
  annotate("text", x = 13, y=0.42, label =  15.81)+
    geom_raster(
    data = subset(
      FITNESS_MODEL,
      FITNESS_MODEL$status ==
        "Fail"
    ),
    aes(x = B_V, y = C_V), fill = "#d1dbe8",
    inherit.aes = FALSE
    
  ) +
  geom_raster(
    data = subset(
      FITNESS_MODEL,
      FITNESS_MODEL$status ==
        "mort"
    ),
    aes(x = B_V, y = C_V), fill = "black",
    inherit.aes = FALSE
    
  )+
  scale_fill_viridis(option='inferno')+
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
  ggtitle("Day 20")+
  scale_x_continuous(expand=c(0,0))+ 
  scale_y_continuous(expand=c(0,0))+
  theme(panel.background=element_rect(fill ='#e3e0da'),
        panel.border = element_rect(color = 'black', fill = NA, size =1),
        panel.grid  = element_blank(),
        
        text = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15, color = "black"))+
  annotate("text",
           x = 8, y = 0.93, label = "Unestablished \ninfection",
           size = 5
  ) +
  annotate("text",
           x = 45, y = 0.1, label = "Host \nMortality",
           size = 5, color = "#b8fcd5"
  )


ggsave(here("Figures", "Raw", "SUPP_Fitness_SurfacePlot_Day20.pdf"),
       width = 10, height = 8, unit = "in"
)
