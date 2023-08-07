###This is the data.frame I need to create
###for the set_duration_surface plot
###Load in the fitness data.frame
Fitness_MODEL_PC_FULL <- as.data.frame(fread(here(
  "Output", "Fitness_Model",
  "Fitness_MODEL_PC_FULL.csv"
)))


###Only get the ones that do not kill the host

SUCCESS_BV_CV<- subset(Fitness_MODEL_PC_FULL,
                       Fitness_MODEL_PC_FULL$status == 'success')


FULL_MODEL_PC_100 <- mcmapply(Simulator_MalariaPC_DDE_BC_Cut,
                             c(SUCCESS_BV_CV$B_V),
                             c(SUCCESS_BV_CV$C_V),
                             100,
                             mc.cores = 4,
                             SIMPLIFY = FALSE)

 
 FULL_MODEL_PC_100_DT <- do.call(rbind, FULL_MODEL_PC_100)


 write.csv(FULL_MODEL_PC_100_DT, file = here(
  "Output", "Full_Model",
  "FULL_MODEL_100_PC_DT.csv"
 ))
 
 
 