
library(here)
# This is the main simulation assuming nothing has changed
# about the initial red blood cell density or the replenishment rate.
# This is one of the longer code to source

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))

### Main modeling code
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_CUT.cpp"))
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_UNCUT.cpp"))
source(here("Code", "Simulator_Code", "Simulator_Main.R"))
source(here("Code", "Helper_Function_Code", "Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "VaryingMortThres.R"))


### Burst Size Versus Transmission Investment ###
B_V <- seq(1, 50, 1) # Burst size
C_V <- seq(0.01, 1, 0.05) # Transmission investment
B_V_C_V <- expand.grid(B_V = B_V, C_V = C_V) # Different combinations


### I am able to figure out which of the burst and transmission
### are unable to establish the infection so that I can exclude from
### simulation (the acute time is set to 0 and transmission potential is set to 0)
p_val <- 2.5e-6
mu_M <- 48
R_val <- 8500000

RM_limit_1 <- ((p_val * R_val) + mu_M) / (p_val * R_val)

B_V_C_V$Establish <- ifelse((1 - B_V_C_V$C_V) * B_V_C_V$B_V >= RM_limit_1,
                            "Establish", "Fail"
)

### Simulate infections that are successful
B_V_C_V_F <- subset(B_V_C_V, B_V_C_V$Establish == "Establish")

### These infections are successful OR lead to host mortality
FULL_MODEL_PC_SUPP <- mcmapply(Simulator_Malaria_BC,
                          c(B_V_C_V_F$B_V),
                          c(B_V_C_V_F$C_V),
                          mc.cores = 5,
                          SIMPLIFY = FALSE
)


Duration_Initial_PC_SUPP <- 
  do.call(
    rbind,
    mclapply(FULL_MODEL_PC_SUPP,
             Finder_RM_DeathLimit ,
             mu_M_c = 48,
             DeathLimit = 2.33* 10^6,
             mc.cores = 2
    )
  )

Duration_Initial_PC_SUPP2 <- 
  do.call(
    rbind,
    mclapply(FULL_MODEL_PC_SUPP,
             Finder_RM_DeathLimit ,
             mu_M_c = 48,
             DeathLimit = 0.1* 10^6,
             mc.cores = 2
    )
  )

### Assign the burst size and transmission investment
Duration_Initial_PC_SUPP$B_V <- B_V_C_V_F$B_V
Duration_Initial_PC_SUPP$C_V <- B_V_C_V_F$C_V

Duration_Initial_PC_SUPP2$B_V <- B_V_C_V_F$B_V
Duration_Initial_PC_SUPP2$C_V <- B_V_C_V_F$C_V




### I know which B_V/C_V combos would lead to failed infections so,
### I create a data.frame for them
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
Duration_Initial_PC_MORT <- subset(
  
  Duration_Initial_PC_SUPP,
  
  Duration_Initial_PC_SUPP$status != "success"
)

Duration_Initial_PC_MORT2 <- subset(
  
  Duration_Initial_PC_SUPP2,
  
  Duration_Initial_PC_SUPP2$status != "success"
)

Duration_Initial_PC_SUCCESS <- subset(
  
  Duration_Initial_PC_SUPP,
  
  Duration_Initial_PC_SUPP$status == "success"
)

Duration_Initial_PC_SUCCESS2 <- subset(
  
  Duration_Initial_PC_SUPP2,
  
  Duration_Initial_PC_SUPP2$status == "success"
)


### Only run the function for the burst size and transmission investment
### combination that leads to successful infection that does not induce
### host mortality

Fitness_MODEL_PC <- mcmapply(Simulator_MalariaPC_DDE_BC_Cut,
                             c(Duration_Initial_PC_SUCCESS$B_V),
                             c(Duration_Initial_PC_SUCCESS$C_V),
                             c(Duration_Initial_PC_SUCCESS$endtime),
                             mc.cores = 1,
                             SIMPLIFY = FALSE
)


Fitness_MODEL_PC2 <- mcmapply(Simulator_MalariaPC_DDE_BC_Cut,
                             c(Duration_Initial_PC_SUCCESS2$B_V),
                             c(Duration_Initial_PC_SUCCESS2$C_V),
                             c(Duration_Initial_PC_SUCCESS2$endtime),
                             mc.cores = 5,
                             SIMPLIFY = FALSE
)

### Now we know what the fitness is

Duration_Initial_PC_SUCCESS$end_fitness <-
  unlist(lapply(Fitness_MODEL_PC, Gametocyte_Fitness))

Duration_Initial_PC_SUCCESS2$end_fitness <-
  unlist(lapply(Fitness_MODEL_PC2, Gametocyte_Fitness))

### These are the fitness model data.frame that should work
Fitness_MODEL_PC_FULL <- rbind.data.frame(
  Duration_Initial_PC_SUCCESS,
  Duration_Initial_PC_FAIL,
  Duration_Initial_PC_MORT
)

Fitness_MODEL_PC_FULL2 <- rbind.data.frame(
  Duration_Initial_PC_SUCCESS2,
  Duration_Initial_PC_FAIL,
  Duration_Initial_PC_MORT2
)

horizontal_vert_df <- grapher_mortality_boundary_Supp(Fitness_MODEL_PC_FULL)
horizontal_vert_df2 <- grapher_mortality_boundary_Supp(Fitness_MODEL_PC_FULL2)


FITNESS1_GG <- ggplot(subset(Fitness_MODEL_PC_FULL,
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
  xlab(expression(paste("Burst size", "( ",beta,")"))) +
  ylab("Transmission investment (c)") + 
theme(
  text = element_text(size = 14),
  axis.text = element_text(size = 14, color = "black"),
  axis.title = element_text(size = 14, color = "black"),
  legend.position="top") +
  annotate("point", x = 7, y = 0.46, shape = 21, size = 2, col = "#00AE97",
           stroke = 2) + 
  annotate('text', 
           x = 3.5, 
           y = 0.93, 
           label = 'Unestablished \ninfection',
           size = 5)+
  annotate('text', 
           x = 45, 
           y = 0.1, 
           label = 'Host \nMortality',
           size = 5, 
           color = '#b8fcd5')


FITNESS2_GG <- ggplot(subset(Fitness_MODEL_PC_FULL2,
              Fitness_MODEL_PC_FULL2$status != 
                'Fail'),
       aes(x = B_V, y = C_V)) +
  geom_raster(aes(fill = end_fitness)) +
  geom_raster(data = subset(Fitness_MODEL_PC_FULL2,
                            Fitness_MODEL_PC_FULL2$status == 
                              'Fail'),
              aes(x=B_V, y = C_V), fill = '#d1dbe8')+
  geom_segment(
    data = horizontal_vert_df2,
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
  xlab(expression(paste("Burst size", "( ",beta,")"))) +
  ylab("Transmission investment (c)") + 
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position="top") +
  annotate("point", x = 7, y = 0.46, shape = 21, size = 2, col = "#00AE97",
           stroke = 2) + 
  annotate('text', 
           x = 3.5, 
           y = 0.93, 
           label = 'Unestablished \ninfection',
           size = 5)+
  annotate('text', 
           x = 45, 
           y = 0.1, 
           label = 'Host \nMortality',
           size = 5, 
           color = '#b8fcd5')

FITNESS1_GG + FITNESS2_GG 


ggsave(here("Figures", "Raw", "Supp_Mort_Thres.pdf"),
       width = 12, height =5, units = 'in')
