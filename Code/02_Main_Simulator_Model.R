##############################################################
###This produces the necessary data files that is needed for ##
###further analysis                                         ##
##############################################################

library(here)
#This is the main simulation assuming nothing has changed
#about the initial red blood cell density or the replenishment rate

###Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))

###Main modeling code 
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_CUT.cpp"))
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_UNCUT.cpp"))
source(here("Code", "Simulator_Code", "Simulator_Main.R"))
source(here("Code", "Helper_Function_Code", "Fitness_Functions.R"))

###Burst Size Versus Transmission Investment ###
B_V = seq(1, 50, 0.5) #Burst size
C_V = seq(.01, 1, 0.01) #Transmission investment 
B_V_C_V <- expand.grid(B_V = B_V, C_V = C_V) #Different combination 

###I am able to figure out which of the burst and transmission

p_val =   2.5e-6
mu_M = 48
R_val = 8500000

RM_limit_1 <- ((p_val * R_val) + mu_M)/(p_val*R_val)

B_V_C_V$Establish <- ifelse((1-B_V_C_V$C_V)* B_V_C_V$B_V > RM_limit_1,
                            "Establish","Fail")

###Simulate infections that are successful

B_V_C_V_F <- subset(B_V_C_V, B_V_C_V$Establish == "Establish")

FULL_MODEL_PC <- mcmapply(Simulator_Malaria_BC, 
                          c(B_V_C_V_F$B_V),
                          c(B_V_C_V_F$C_V), 
                            mc.cores = 7,
                            SIMPLIFY = FALSE)    

FULL_MODEL_PC_DT <- do.call(rbind, FULL_MODEL_PC)

write.csv(FULL_MODEL_PC_DT, file = here("Output", "Full_Model",
                                        "FULL_MODEL_PC_DT.csv"))

######################################################
###This is the duration of the acute phase (RM WAY)#
######################################################
Duration_Initial_PC_RM <- do.call(rbind, 
                                 mclapply(FULL_MODEL_PC,
                                         Finder_RM, 
                                         mu_M_c = 48, 
                                         mc.cores = 5))


###Assign the burst size and transmission investment 
Duration_Initial_PC$B_V <- B_V_C_V_F$B_V
Duration_Initial_PC$C_V <- B_V_C_V_F$C_V

###Only run the function for the burst size and transmission investment
###combination that leads to successful infection that does not induce
###host mortality

DURATION_INITIAL_PC_FAIL_MORT <-  subset(Duration_Initial_PC, 
                                   Duration_Initial_PC$status != "success")


Duration_Initial_PC_SUCCESS <- subset(Duration_Initial_PC, 
                                      Duration_Initial_PC$status == "success")

Fitness_MODEL_PC <- mcmapply(Simulator_MalariaPC_DDE_BC_Cut,
                    c(Duration_Initial_PC_SUCCESS$B_V),
                    c(Duration_Initial_PC_SUCCESS$C_V),
                    c(Duration_Initial_PC_SUCCESS$endtime),
                    mc.cores = 9,
                    SIMPLIFY = FALSE)

Duration_Initial_PC_SUCCESS$end_fitness<- 
  unlist(lapply(Fitness_MODEL_PC, Gametocyte_Fitness))


Fitness_Model_PC_FULL <- rbind.data.frame(Duration_Initial_PC_SUCCESS, 
                                          DURATION_INITIAL_PC_FAIL_MORT)


write.csv(Fitness_Model_PC_FULL , file = here("Output",
                                              "Fitness_Full",
                                              "Fitness_MODEL_PC_FULL.csv"))




