
FITNESS_MODEL <- as.data.frame(fread(here(
  "Output", "Fitness_Model",
  "Fitness_MODEL_PC_FULL.csv"
)))


FULL_MODEL_PC_100 <- mcmapply(Simulator_MalariaPC_DDE_BC_Cut,
                             c(SUCCESS_BV_CV$B_V),
                             c(SUCCESS_BV_CV$C_V),
                             100,
                             mc.cores = 4,
                             SIMPLIFY = FALSE)

 Combine
 FULL_MODEL_PC_100_DT <- do.call(rbind, FULL_MODEL_PC_100)


 write.csv(FULL_MODEL_PC_100_DT, file = here(
  "Output", "Full_Model",
  "FULL_MODEL_100_PC_DT.csv"
 ))
 
