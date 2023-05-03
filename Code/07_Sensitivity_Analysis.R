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


B_V_C_V_UM_F <- B_V_C_V_Establisher("UM",0.25)
B_V_C_V_RI_F <-  B_V_C_V_Establisher("RI",0.25)
B_V_C_V_L_F <-  B_V_C_V_Establisher("L", 0.25)

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
                          c(B_V_C_V_UM_F_E[,3]),
                          "UM",
                          mc.cores = 5,
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
                             c(B_V_C_V_RI_F_E[,3]),
                             "RI",
                             mc.cores = 5,
                             SIMPLIFY = FALSE)

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
  mclapply(FULL_MODEL_PC_RI,
           Finder_RM_SA,
           sens_var = 'RI',
           mc.cores = 5))


FITNESS_RI <- Fitness_Finder_SA(
  B_V_C_V_RI_F,
  Duration_RI_PC)

write.csv(FITNESS_RI , file = here(
  "Output", "Fitness_Model",
  "FITNESS_RI.csv"
))
gc()

remove(FULL_MODEL_PC_RI )
remove(Duration_RI_PC)
remove(FITNESS_RI )
gc()




### These infections are successful OR lead to host mortality
FULL_MODEL_PC_L <- mcmapply(Sens_Analysis,
                             c(B_V_C_V_L_F_E$B_V),
                             c(B_V_C_V_L_F_E$C_V),
                             c(B_V_C_V_L_F_E[,3]),
                            "lambda",
                             mc.cores = 5,
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
  Duration_L_PC)

write.csv(FITNESS_L , file = here(
  "Output", "Fitness_Model",
  "FITNESS_L.csv"
))

remove(FULL_MODEL_PC_L )
remove(Duration_L_PC)
remove(FITNESS_L)
gc()




