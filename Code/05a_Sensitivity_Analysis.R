########################################################
### Understanding how the viable burst size combinations#
### influence virulence, transmission, and such.     ####
########################################################

library(here)

### Packages to load
source(here("Code", "Helper_Function_Code", "FUNC_00_Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "00_Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "07_sensitivity_analysis_functions.R"))
source(here("Code", "Simulator_Code", "Simulator_Main_2.R"))

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
#B_V, C_V, initial_value, change, id
FULL_MODEL_PC_UM <- mcmapply(Sens_Analysis,
                          c(B_V_C_V_UM_F_E$B_V),
                          c(B_V_C_V_UM_F_E$C_V),
                          c(4385.96491),
                          c(B_V_C_V_UM_F_E[,3]),
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


###Because the mu_M changes... We need to this some data-wrangling to make sure
###that I'm getting the duration correctly 

FULL_MODEL_PC_UM_Changes <- split(FULL_MODEL_PC_UM_DT, FULL_MODEL_PC_UM_DT$change)
FULL_MODEL_PC_UM_36<- split(FULL_MODEL_PC_UM_Changes[[1]], list(FULL_MODEL_PC_UM_Changes[[1]]$B_V,
                                                                FULL_MODEL_PC_UM_Changes[[1]]$C_V),
                            drop = TRUE)
FULL_MODEL_PC_UM_60<- split(FULL_MODEL_PC_UM_Changes[[2]], list(FULL_MODEL_PC_UM_Changes[[2]]$B_V,
                                                           FULL_MODEL_PC_UM_Changes[[2]]$C_V),
                            drop = TRUE)

FULL_36_MODEL_BVCV <- do.call(rbind, lapply(FULL_MODEL_PC_UM_36, function(x) unique(x[,c("B_V", "C_V")])))
FULL_60_MODEL_BVCV <- do.call(rbind, lapply(FULL_MODEL_PC_UM_60, function(x) unique(x[,c("B_V", "C_V")])))

###Get how long the acute phase ends for a mortality rate of 36
Duration_UM_PC_36 <- do.call(
  rbind,
  mclapply(FULL_MODEL_PC_UM_36,
           Finder_RM,
           mu_M_c = 36,
           mc.cores = 2))

Duration_UM_PC_36$B_V <-FULL_36_MODEL_BVCV$B_V
Duration_UM_PC_36$C_V <-FULL_36_MODEL_BVCV$C_V
Duration_UM_PC_36$change <-36

###Get how long the acute phase ends for a mortality rate of 60

Duration_UM_PC_60 <- do.call(
  rbind,
  mclapply(FULL_MODEL_PC_UM_60,
           Finder_RM,
           mu_M_c = 60,
           mc.cores = 2))

Duration_UM_PC_60$B_V <-FULL_60_MODEL_BVCV$B_V
Duration_UM_PC_60$C_V <-FULL_60_MODEL_BVCV$C_V
Duration_UM_PC_60$change<-60

Duration_UM_PC <- rbind(Duration_UM_PC_36,Duration_UM_PC_60)

FITNESS_UM <- Fitness_Finder_SA(
                  B_V_C_V_UM_F,
                  Duration_UM_PC,
                   4385.96491, "UM") ###This is the fitness function

write.csv(FITNESS_UM , file = here(
  "Output", "Fitness_Model",
  "FITNESS_UM.csv"
))


#############################
###RI: Initial RBC        ###
#############################

### These infections are successful OR lead to host mortality
FULL_MODEL_PC_RI <- mcmapply(Sens_Analysis,
                             c(B_V_C_V_RI_F_E$B_V),
                             c(B_V_C_V_RI_F_E$C_V),
                             c(4385.96491),
                             c(B_V_C_V_RI_F_E[,3]),
                             "RI",
                             mc.cores = 2,
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
           Finder_RM,
           mu_M_c = 48,
           mc.cores = 5))

Duration_RI_PC$B_V <- B_V_C_V_RI_F_E$B_V
Duration_RI_PC$C_V <- B_V_C_V_RI_F_E$C_V
Duration_RI_PC$change<- B_V_C_V_RI_F_E$RI


FITNESS_RI <- Fitness_Finder_SA(
  B_V_C_V_RI_F,
  Duration_RI_PC,
  4385.96491,"RI")

write.csv(FITNESS_RI , file = here(
  "Output", "Fitness_Model",
  "FITNESS_RI.csv"
))


### These infections are successful OR lead to host mortality
FULL_MODEL_PC_L <- mcmapply(Sens_Analysis,
                             c(B_V_C_V_L_F_E$B_V),
                             c(B_V_C_V_L_F_E$C_V),
                            c(4385.96491),
                             c(B_V_C_V_L_F_E[,3]),
                            "lambda",
                             mc.cores = 2,
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
  mclapply(FULL_MODEL_PC_L ,
           Finder_RM,
           mu_M_c = 48,
           mc.cores = 4))

Duration_L_PC$B_V <- B_V_C_V_L_F_E $B_V
Duration_L_PC$C_V <- B_V_C_V_L_F_E $C_V
Duration_L_PC$change<- B_V_C_V_L_F_E$L



FITNESS_L <- Fitness_Finder_SA(
  B_V_C_V_L_F,
  Duration_L_PC,
  4385.96491,"lambda")



write.csv(FITNESS_L , file = here(
  "Output", "Fitness_Model",
  "FITNESS_L.csv"
))

###After you have all the files saved you can proceed to code 
###08_Sensitivity_SurfacePlots_Supp.R and 09_Sensitivity_LinePlot 
