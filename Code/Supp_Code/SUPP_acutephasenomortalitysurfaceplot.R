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

### These infections are successful OR lead to host mortality
FULL_MODEL_PC_NODEATH <- mcmapply(Simulator_Malaria_BC_NODEATH,
                          c(B_V_C_V_F$B_V),
                          c(B_V_C_V_F$C_V),
                          c(B_V_C_V_F$initialvalue),
                          mc.cores = 3,
                          SIMPLIFY = FALSE)

### Combine all list elements to make it easier to save/read
FULL_MODEL_PC_DT_NODEATH <- do.call(rbind, FULL_MODEL_PC_NODEATH)

### Write into a CSV TO BE SAVED
write.csv(FULL_MODEL_PC_DT_NODEATH, file = here(
  "Output", "Full_Model","FULL_MODEL_PC_DT_NODEATH.csv"))

###REMOVE right now to save space
remove(FULL_MODEL_PC_DT)




Duration_Initial_PC_NODEATH <- 
  do.call(
    rbind,
    mclapply(FULL_MODEL_PC_NODEATH ,
             Finder_RM_Nodeath,
             mu_M_c = 48,
             mc.cores = 2
    )
  )

Duration_Initial_PC_NODEATH$B_V <- B_V_C_V_F$B_V
Duration_Initial_PC_NODEATH$C_V <- B_V_C_V_F$C_V

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

### These are the B_V/C_V that would lead to mortality,
### This means that we know the end time (point of death) and the
### cumulative transmission potential


Duration_Initial_PC_SUCCESS_NODEATH <- subset(
  Duration_Initial_PC_NODEATH,
  Duration_Initial_PC_NODEATH$status == "success"
)


Fitness_MODEL_PC_NODEATH <- mcmapply(Simulator_MalariaPC_DDE_BC_Cut,
                             c(Duration_Initial_PC_SUCCESS_NODEATH $B_V),
                             c(Duration_Initial_PC_SUCCESS_NODEATH $C_V),
                             c( 4385.965),
                             c(Duration_Initial_PC_SUCCESS_NODEATH $endtime),
                             mc.cores = 2,
                             SIMPLIFY = FALSE
)


### Now we know what the fitness is
Duration_Initial_PC_SUCCESS_NODEATH$end_fitness <-
  unlist(lapply(Fitness_MODEL_PC_NODEATH, Gametocyte_Fitness))

### These are the fitness model data.frame that should work
Fitness_MODEL_PC_FULL_NODEATH <- rbind.data.frame(
  Duration_Initial_PC_SUCCESS_NODEATH,
  Duration_Initial_PC_FAIL
)


### Write into a CSV TO BE SAVED 
write.csv(Fitness_MODEL_PC_FULL_NODEATH , file = here(
  "Output", "Fitness_Model",'Fitness_MODEL_PC_FULL_NODEATH.csv'))


optimum_strategy <- Fitness_MODEL_PC_FULL_NODEATH[which.max(Fitness_MODEL_PC_FULL_NODEATH$end_fitness),]


GG_Fitness_Cut_PC <-
  ggplot(
    subset(
      Fitness_MODEL_PC_FULL_NODEATH,
      Fitness_MODEL_PC_FULL_NODEATH$status !=
        "Fail"
    ),
    aes(x = B_V, y = C_V)
  ) +
  geom_raster(aes(fill = end_fitness)) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL_NODEATH,
      Fitness_MODEL_PC_FULL_NODEATH$status ==
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
  ) +
  annotate("text",
           x = 45, y = 0.1, label = "Host \nMortality",
           size = 5, color = "#b8fcd5"
  )

ggsave(here("Figures", "Raw", "Supp_nomort_Surfaceplot.pdf"),
       width = 8, height =8, unit = "in")
