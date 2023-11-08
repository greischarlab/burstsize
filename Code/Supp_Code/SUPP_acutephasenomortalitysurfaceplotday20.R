### Burst Size and Transmission Investment ###
B_V <- seq(1, 50, 0.5) # Burst size
C_V <- seq(.01, 1, 0.01) # Transmission investment
B_V_C_V <- expand.grid(B_V = B_V, C_V = C_V, initialvalue = 4385.96491) # Different combinations

###We can already take out parameter combinations that would not
###lead to the establishment of infection
p_val <- 4.0e-6 #maximum infection rate 
mu_M <- 48  #the mu_mortality 
R_val <- 8500000 #initial red blood cell
initial_RM_modifier <- 1.5 #

RM_limit_1 <- initial_RM_modifier * ((p_val * R_val) + mu_M) / (p_val * R_val)

B_V_C_V$Establish <- ifelse((1 - B_V_C_V$C_V) * B_V_C_V$B_V >= RM_limit_1,
                            "Establish", "Fail")

### Simulate infections that are successful (may kill host!)
B_V_C_V_F <- subset(B_V_C_V, B_V_C_V$Establish == "Establish")


Failed_B_V_C_V <- subset(B_V_C_V, B_V_C_V$Establish == "Fail")

Duration_Initial_PC_FAIL <-
  data.frame(
    endtime = 0,
    up_down = 0,
    end_fitness = 0,
    status = "Fail",
    B_V = Failed_B_V_C_V$B_V,
    C_V = Failed_B_V_C_V$C_V
  )

Duration_Initial_PC_SUCCESS_20 <-
  data.frame(
    endtime = 20,
    up_down = "up",
    end_fitness = 0,
    status = "Success",
    B_V = B_V_C_V_F $B_V,
    C_V = B_V_C_V_F $C_V
  )

Fitness_MODEL_PC_NODEATH_20 <- mcmapply(Simulator_MalariaPC_DDE_BC_Cut,
                                     c(B_V_C_V_F  $B_V),
                                     c(B_V_C_V_F $C_V),
                                     c( 4385.965),
                                     20,
                                     mc.cores = 2,
                                     SIMPLIFY = FALSE
)

### Now we know what the fitness is
Duration_Initial_PC_SUCCESS_20 $end_fitness <-
  unlist(lapply(Fitness_MODEL_PC_NODEATH_20, Gametocyte_Fitness))

### These are the fitness model data.frame that should work
Fitness_MODEL_PC_FULL_NODEATH_20 <- rbind.data.frame(
  Duration_Initial_PC_SUCCESS_20 ,
  Duration_Initial_PC_FAIL
)


### Write into a CSV TO BE SAVED 
write.csv(Fitness_MODEL_PC_FULL_NODEATH_20, file = here(
  "Output", "Fitness_Model",'Fitness_MODEL_PC_FULL_NODEATH_20.csv'))


optimum_strategy <- Fitness_MODEL_PC_FULL_NODEATH_20[which.max(Fitness_MODEL_PC_FULL_NODEATH_20 $end_fitness),]


GG_Fitness_Cut_PC <-
  ggplot(
    subset(
      Fitness_MODEL_PC_FULL_NODEATH_20 ,
      Fitness_MODEL_PC_FULL_NODEATH_20 $status !=
        "Fail"
    ),
    aes(x = B_V, y = C_V)
  ) +
  geom_raster(aes(fill = end_fitness)) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL_NODEATH_20 ,
      Fitness_MODEL_PC_FULL_NODEATH_20 $status ==
        "Fail"
    ),
    aes(x = B_V, y = C_V), fill = "#d1dbe8"
  ) +
  scale_fill_viridis(
    name = "Cumulative transmission \npotential",
    option = "magma"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Transmission investment (c)") +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "top"
  ) +
  geom_point(data = optimum_strategy ,
             aes( x = B_V, y= C_V),
             color = '#FF116B', size = 2)+
  geom_segment(
    data = optimum_strategy,
    aes(
      x = B_V, xend = B_V, y = 0.01, yend = C_V)
    , col = "#FF116B",
    linetype = 2) +
  geom_segment(data = optimum_strategy ,
               aes(
                 x =1, xend = B_V, y = C_V, yend = C_V), col = "#FF116B",
               linetype = 2) +
  
  annotate("text",
           x = 8, y = 0.93, label = "Unestablished \ninfection",
           size = 5
  ) 



ggsave(here("Figures", "Raw", "Supp_nomort_Surfaceplot_20.pdf"),
       width = 8, height =8, unit = "in")
