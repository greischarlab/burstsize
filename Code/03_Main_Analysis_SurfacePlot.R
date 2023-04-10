###Makes Figure 2, which is the surface plot of the 
###acute phase duration and the fitness

Fitness_MODEL_PC_FULL <- read.csv(here("Output", "Fitness_Full",
                              "Fitness_MODEL_PC_FULL.csv"))

source(here("Code", "Helper_Function_Code", "Grapher_vert_hor.R"))

horizontal_vert_df <- grapher_mortality_boundary(Full_Data_PC)


#####################################
###The duration of the acute phases##
#####################################
GG_Duration_Main <-
  ggplot(Fitness_MODEL_PC_FULL,
    aes(x = B_V, y = C_V)) +
  geom_raster(aes(fill = endtime)) +
  geom_segment(
    data = horizontal_vert_df,
    aes(
      x = x, xend = xend,
      y = y, yend = yend,
      color = id),
      size = 1.2, 
      lineend = "round") +
  scale_fill_distiller(
      name = "Duration \n (days)",
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
    axis.title = element_text(size = 15, color = "black")) +
  annotate("point", x = 6, y = 0.35, size = 3, col = "red")

###########################
###Plots with the fitness##
###########################
GG_Fitness_Cut_PC <- 
  ggplot(Fitness_MODEL_PC_FULL, aes(x = B_V, y = C_V)) +
  geom_raster(aes(fill = end_fitness/max(end_fitness)))+
  geom_segment(
      data = horizontal_vert_df,
      aes(x = x, xend = xend,
          y = y, yend = yend,
          color = id), size = 1.1, 
          lineend = "round")+
  scale_fill_viridis(
    name = "Proportion of\n max. cumulative \n transmission potential",
    option = "magma", limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression("Burst size" ~ beta)) +
  ylab("Transmission investment (c)") +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, color = "black")) +
  annotate("point", x = 6, y = 0.35, size = 2, col = "red")

GG_Duration_Main + GG_Fitness_Cut_PC 
  