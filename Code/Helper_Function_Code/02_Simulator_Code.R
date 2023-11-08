###########################################################
###THIS IS THE FUNCTION THAT SIMULATES INFECTIONS         #
###AND TELLS YOU THE FITNESS AT THE END OF THE ACUTE PHASE#
###########################################################

#Input: 
#1)initial_value - The different intial values for the infected RBC inoculum
#2) mu_M_value - the mortality rate of merozoite, crucial for figuring out
# which infections are unestablished
#3 id - give it a low, med, high id


FULL_MODEL_SIMULATING_Duration <- function(initial_value, mu_M_value, id){

  ##########################################################
  ###THE FIRST PART OF THE MODEL, trying to figure out what#
  ###parameter combinations we want to be looking at       #
  ##########################################################
  
  ### Burst Size and Transmission Investment ###
  B_V <- seq(1, 50, 0.5) # Burst size
  C_V <- seq(.01, 1, 0.01) # Transmission investment
  B_V_C_V <- expand.grid(B_V = B_V, C_V = C_V, initialvalue = initial_value) # Different combinations
  
  ###We can already take out parameter combinations that would not
  ###lead to the establishment of infection
  p_val <- 4.0e-6 #maximum infection rate 
  mu_M <- mu_M_value #the mu_mortality 
  R_val <- 8500000 #initial red blood cell
  initial_RM_modifier <- 1.5 #

  RM_limit_1 <- initial_RM_modifier * ((p_val * R_val) + mu_M) / (p_val * R_val)

  B_V_C_V$Establish <- ifelse((1 - B_V_C_V$C_V) * B_V_C_V$B_V >= RM_limit_1,
                            "Establish", "Fail")

  ### Simulate infections that are successful (may kill host!)
  B_V_C_V_F <- subset(B_V_C_V, B_V_C_V$Establish == "Establish")

  ### These infections are successful OR lead to host mortality
  FULL_MODEL_PC <- mcmapply(Simulator_Malaria_BC,
                          c(B_V_C_V_F$B_V),
                          c(B_V_C_V_F$C_V),
                          c(B_V_C_V_F$initialvalue),
                          mc.cores = 5,
                          SIMPLIFY = FALSE)
  
  ### Combine all list elements to make it easier to save/read
  FULL_MODEL_PC_DT <- do.call(rbind, FULL_MODEL_PC)
  
  ### Write into a CSV TO BE SAVED
  write.csv(FULL_MODEL_PC_DT, file = here(
    "Output", "Full_Model",
    paste("FULL_MODEL_PC_DT_",id,".csv")))
  
  ###REMOVE right now to save space
  remove(FULL_MODEL_PC_DT)
  
  #########################################################
  ###THIS THEN FIGURES OUT THE DURATION OF THE ACUTE PHASE#
  #########################################################
  Duration_Initial_PC <- 
    do.call(
      rbind,
      mclapply(FULL_MODEL_PC,
               Finder_RM,
               mu_M_c = 48,
               mc.cores = 2
      )
    )

  Duration_Initial_PC$B_V <- B_V_C_V_F$B_V
  Duration_Initial_PC$C_V <- B_V_C_V_F$C_V
  
  
  ###Now we can find the parameter combinations that lead to 
  ###unestablished infection- we then set the values 
  ###manually to 9
  
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
    Duration_Initial_PC,
    Duration_Initial_PC$status != "success"
  )
  
  Duration_Initial_PC_SUCCESS <- subset(
    Duration_Initial_PC,
    Duration_Initial_PC$status == "success"
  )
  
  
  Fitness_MODEL_PC <- mcmapply(Simulator_MalariaPC_DDE_BC_Cut,
                               c(Duration_Initial_PC_SUCCESS$B_V),
                               c(Duration_Initial_PC_SUCCESS$C_V),
                               c(initial_value),
                               c(Duration_Initial_PC_SUCCESS$endtime),
                               mc.cores = 5,
                               SIMPLIFY = FALSE
  )
  
  Duration_Initial_PC_SUCCESS$end_fitness <-
    unlist(lapply(Fitness_MODEL_PC, Gametocyte_Fitness))
  
  ### These are the fitness model data.frame that should work
  Fitness_MODEL_PC_FULL <- rbind.data.frame(
    Duration_Initial_PC_SUCCESS,
    Duration_Initial_PC_FAIL,
    Duration_Initial_PC_MORT
  )
  
  ### Write into a CSV TO BE SAVED 
  write.csv(Fitness_MODEL_PC_FULL, file = here(
    "Output", "Fitness_Model",
    paste("FITNESS_MODEL_PC_",id,".csv", sep = "")))
  
  
  return(Fitness_MODEL_PC_FULL)
  
  
}



