
C_V <- seq(.01, 1, 0.01)

FULL_MODEL_PC <- mcmapply(Simulator_Malaria_BC_ALTERNATIVE,
                          10,
                          C_V,
                          mc.cores = 5,
                          SIMPLIFY = FALSE
)

tryit <- Surface_Plot_Standard_Duration(FULL_MODEL_PC,20)

plot(tryit$C_V, tryit$end_fitness, xlab = "Transmission investment",
ylab="End fitness at Day 20", main = 'p = 4.0e-6, I = 43859.65')
abline(v = 0.42)

Duration_Initial_PC <- 
  do.call(
    rbind,
    mclapply(FULL_MODEL_PC,
             Finder_RM,
             mu_M_c = 48,
             mc.cores = 2
    )
  )

Duration_Initial_PC$C_V <- C_V


success <- subset(Duration_Initial_PC , Duration_Initial_PC$status == 'success')
