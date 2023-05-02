########################################################
### Understanding how the viable burst size combinations#
### influence virulence, transmission, and such.     ####
########################################################

library(here)

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "sensitivity_analysis_functions.R"))
source(here("Code", "Simulator_Code", "Simulator_Main.R"))

sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_CUT.cpp"))
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_UNCUT.cpp"))




### Burst Size Versus Transmission Investment ###
B_V <- seq(1, 50, 0.5) # Burst size
C_V <- seq(.01, 1, 0.01) # Transmission investment
mu_M_vec <- c(48 * 0.75, 48 * 1.25)
RI_vec <-  c(8500000 * 0.75, 8500000 * 1.25)
L_vec <- c(277500,462500)

B_V_C_V_UM <- expand.grid(B_V = B_V, C_V = C_V, mu_M = mu_M_vec) # Different combinations
B_V_C_V_RI <- expand.grid(B_V = B_V, C_V = C_V, RI = RI_vec)
B_V_C_V_L<- expand.grid(B_V = B_V, C_V = C_V, L = L_vec)

### I am able to figure out which of the burst and transmission
### are unable to establish the infection so that I can exclude from
### simulation (the acute time is set to 0 and transmission potential is set to 0)

p_val <- 2.5e-6
mu_M <- 48
R_val <- 8500000

RM_limit_UM <- ((p_val * R_val) + mu_M_vec) / (p_val * R_val)
RM_limit_RI <- ((p_val * RI_vec) + mu_M ) / (p_val * RI_vec)
RM_limit <- ((p_val * R_val) + mu_M ) / (p_val * R_val)

B_V_C_V_UM_F <- B_V_C_V_Establisher("UM")
B_V_C_V_RI_F <-  B_V_C_V_Establisher("RI")
B_V_C_V_L_F <-  B_V_C_V_Establisher("L")

### Simulate infections that are successful
B_V_C_V_UM_F_E <- subset(B_V_C_V_UM_F, B_V_C_V_UM_F$Establish == "Establish")
B_V_C_V_RI_F_E <- subset(B_V_C_V_RI_F, B_V_C_V_RI_F$Establish == "Establish")
B_V_C_V_L_F_E <-  subset(B_V_C_V_L_F , B_V_C_V_L_F $Establish == "Establish")

#save(B_V_C_V_UM_F_E, file = here("Cluster", "B_V_C_V_UM_F_E.RData"))
#save(B_V_C_V_RI_F_E, file = here("Cluster", "B_V_C_V_RI_F_E.RData"))
#save(B_V_C_V_L_F_E, file = here("Cluster", "B_V_C_V_L_F_E.RData"))

#############################
###UM: Merozoite mortality###
#############################

### These infections are successful OR lead to host mortality
FULL_MODEL_PC_UM <- mcmapply(Sens_Analysis,
                          c(B_V_C_V_UM_F_E$B_V),
                          c(B_V_C_V_UM_F_E$C_V),
                          c(0.25),
                          "UM",
                          mc.cores = 3,
                          SIMPLIFY = FALSE)

### Combine
FULL_MODEL_PC_UM_DT <- do.call(rbind, FULL_MODEL_PC_UM)

### Write into a CSV
write.csv(FULL_MODEL_PC_UM_DT, file = here(
  "Output", "Full_Model",
  "SA_FULL_MODEL_PC_UM_DT.csv"
))
remove(FULL_MODEL_PC_UM_DT)




Duration_UM_PC <- do.call(
  rbind,
  mclapply(FULL_MODEL_PC_UM,
           Finder_RM_SA,
           sens_var = 'UM',
           mc.cores = 4))


FITNESS_UM <- Fitness_Finder_SA(
                  B_V_C_V_UM_F,
                  "UM",
                  Duration_UM_PC)

write.csv(FITNESS_UM , file = here(
  "Output", "Fitness_Model",
  "FITNESS_UM.csv"
))
remove(FULL_MODEL_PC_UM )
remove(Duration_UM_PC)
remove(FITNESS_UM )
gc()


#############################
###RI: Initial RBC        ###
#############################

### These infections are successful OR lead to host mortality
FULL_MODEL_PC_RI <- mcmapply(Sens_Analysis,
                             c(B_V_C_V_RI_F_E$B_V),
                             c(B_V_C_V_RI_F_E$C_V),
                             c(0.25),
                             "RI",
                             mc.cores = 3,
                             SIMPLIFY = FALSE)

FULL_MODEL_PC_PI <- lapply(FULL_MODEL_PC_RI, function (x) 
                            do.call(rbind,x))


### Combine
FULL_MODEL_PC_RI_DT <- do.call(rbind, FULL_MODEL_PC_RI)


### Write into a CSV
write.csv(FULL_MODEL_PC_RI_DT, file = here(
  "Output", "Full_Model",
  "SA_FULL_MODEL_PC_RI_DT.csv"
))
remove(FULL_MODEL_PC_RI_DT)

Duration_RI_PC <- do.call(
  rbind,
  mclapply(FULL_MODEL_PC_PI_F,
           Finder_RM_SA,
           sens_var = 'RI',
           mc.cores = 1))


FITNESS_RI <- Fitness_Finder_SA(
  B_V_C_V_RI_F,
  "RI",
  Duration_RI_PC)

write.csv(FITNESS_RI , file = here(
  "Output", "Fitness_Model",
  "FITNESS_RI.csv"
))

remove(FULL_MODEL_PC_RI )
remove(Duration_RI_PC)
remove(FITNESS_RI )
gc()




### These infections are successful OR lead to host mortality
FULL_MODEL_PC_L <- mcmapply(Sens_Analysis,
                             c(B_V_C_V_L_F_E$B_V),
                             c(B_V_C_V_L_F_E$C_V),
                             c(0.25),
                             "lambda",
                             mc.cores = 3,
                             SIMPLIFY = FALSE)

### Combine
FULL_MODEL_PC_L_DT <- do.call(rbind, FULL_MODEL_PC_L)

### Write into a CSV
write.csv(FULL_MODEL_PC_L_DT, file = here(
  "Output", "Full_Model",
  "SA_FULL_MODEL_PC_lambda_DT.csv"
))
remove(FULL_MODEL_PC_L_DT)

Duration_L_PC <- do.call(
  rbind,
  mclapply(FULL_MODEL_PC_L,
           Finder_RM_SA,
           sens_var = 'L',
           mc.cores = 4))


FITNESS_L <- Fitness_Finder_SA(
  B_V_C_V_L_F,
  "L",
  Duration_L_PC)

write.csv(FITNESS_L , file = here(
  "Output", "Fitness_Model",
  "FITNESS_L.csv"
))

remove(FULL_MODEL_PC_L )
remove(Duration_L_PC)
remove(FITNESS_L)
gc()




