###Makes Figure 2, which is the surface plot of the 
###acute phase duration and the fitness

Fitness_MODEL_PC_FULL <- read.csv(here("Output", "Fitness_Model",
                              "Fitness_MODEL_PC_FULL.csv"))

source(here("Code", "Helper_Function_Code", "Grapher_vert_hor.R"))

horizontal_vert_df <- grapher_mortality_boundary(Fitness_Model_PC_FULL)



#####################################
###The duration of the acute phases##
#####################################
GG_Duration_Main <-
  ggplot(subset(Fitness_MODEL_PC_FULL,
                Fitness_MODEL_PC_FULL$status != 
                  'Fail'),
      aes(x = B_V, y = C_V)) +
  geom_raster(aes(fill = endtime)) +
  geom_raster(data = subset(Fitness_MODEL_PC_FULL,
                            Fitness_MODEL_PC_FULL$status == 
                              'Fail'),
              aes(x=B_V, y = C_V), fill = '#d1dbe8')+
  geom_segment(
       data = horizontal_vert_df,
             aes(x = x, 
                 xend = xend,
                 y = y, 
                 yend = yend,
                 color = id),
             size = 1.2, 
      lineend = "round") +
  scale_color_manual(values = c("fail" = 'black', 
                                'mort' = '#b8fcd5'),
                     guide = 'none') + 
  scale_fill_distiller(
      name = "Acute phase \n duration (days)",
      type = "seq",
      direction = -1,
      palette = "Greys",
      na.value = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression("Burst size" ~ beta)) +
  ylab("Transmission investment (c)") +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    legend.position="top") +
  annotate("point", x = 8, y = 0.59, shape =21, size = 2, col = "#0037ff",
           stroke = 1) + 
  annotate("point", x = 8, y = 0.54, shape = 21, size = 2, col = "#ff52a3",
           stroke = 1) + 
  annotate('text', x = 8, y = 0.93, label = 'Unestablished \ninfection',
           size = 5)+
  annotate('text', x = 45, y = 0.1, label = 'Host \nMortality',
         size = 5, color = '#b8fcd5')

  


###########################
###Plots with the fitness##
###########################
GG_Fitness_Cut_PC <- 
  ggplot(subset(Fitness_MODEL_PC_FULL,
                Fitness_MODEL_PC_FULL$status != 
                  'Fail'),
         aes(x = B_V, y = C_V)) +
  geom_raster(aes(fill = end_fitness)) +
  geom_raster(data = subset(Fitness_MODEL_PC_FULL,
                            Fitness_MODEL_PC_FULL$status == 
                              'Fail'),
              aes(x=B_V, y = C_V), fill = '#d1dbe8')+
  geom_segment(
      data = horizontal_vert_df,
      aes(x = x, xend = xend,
          y = y, yend = yend,
          color = id), size = 1.1, 
          lineend = "round")+
  scale_color_manual(values = c("fail" = 'black', 
                                'mort' = '#b8fcd5'),
                     guide = 'none')+
  scale_fill_viridis(
    name = "Cumulative transmission \npotential",
    option = "magma")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression("Burst size" ~ beta)) +
  ylab("Transmission investment (c)") +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position="top") +
  annotate("point", x = 8, y = 0.59, shape =21, size = 2, col = "#0037ff",
           stroke = 1) + 
  annotate("point", x = 8, y = 0.54, shape = 21, size = 2, col = "#ff52a3",
           stroke = 1) + 
  annotate('text', x = 8, y = 0.93, label = 'Unestablished \ninfection',
           size = 5)+
  annotate('text', x = 45, y = 0.1, label = 'Host \nMortality',
           size = 5, color = '#b8fcd5')


GG_Duration_Main + GG_Fitness_Cut_PC + plot_annotation(tag_levels = "A")

ggsave(here("Figures", "Raw", "02_Duration_Fitness_SurfacePlot.pdf"),
       width = 12, height =6, unit = 'in')
