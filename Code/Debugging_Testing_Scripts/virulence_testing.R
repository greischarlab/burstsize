###TEsting script for the virulence-transmission script analysis
source(here("Code","Helper_Function_Code", "Virulence_Transmission_BurstSize.R"))
source(here("Code", "Simulator_Code", "Simulator_Main_Gflux.R"))
tmp_Death <- Simulator_MalariaPC_DDE_BC_GFLUX(20,0.54,5.863140)
plot(tmp_Death$time, tmp_Death$R,type='l')
abline(v = 5.863140)
abline(h = 6.5*10^5)


plot(tmp_Death$time,tmp_Death$Gflux)

max(tmp_Death$Gflux)
abline(h = 4626761, col = 'blue')
abline(v = 5.863140)
abline(v = 7.863140)
