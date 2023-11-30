
sourceCpp(here("Code","RCPP_Code", "rcpp_malaria_dynamics_CUT_gflux.cpp"))


Simulator_MalariaPC_DDE_BC_GFLUX<- function(B_V, C_V, initial_value, endtime) {
  
  parameters_n <-c(
    lambda = 370000, # Replenishment rate of RBC (#SimulatedTimeSeries.R)
    K = 19968254, # Carrying capacity of RBC population in the absence of mortality
      pmax = 4.0e-6, # Rate of infection (From American Naturalist- Greischar et al. 2014)
      muR = 0.025, # Daily mortality rate of red blood cells
      muI = 0.025, # Daily mortality rate of infected red blood cells
      c = C_V, # Transmission investment (THE VARYING FACTOR)
      B = B_V, # The burst size (THE VARYING FACTOR)
      alpha1 = 1, # The rate of development of parasite in iRBC
      alpha2 = 1 / 2, # The rate of development
      muM = 48, # Background mortality of the merozoite
      muG = 4, # background mortality of the immature/mature gametocytes
      n1 = 100, # shape parameter controlling the variability in asexual devleopment
      n2 = 100, # shape parameter controlling the variability in sexual development
      endtime = endtime # When the asexual replication end
      )

    ### The number of subcompartments for infected rbc(n1)
    n1 <- parameters_n["n1"]
    ### The number of subcompartments for immature gametocytes
    n2 <- parameters_n["n2"]
  
  ### The initial numbers
  inits_n <- c(
    R = 8500000,
    I = rep(initial_value/n1, n1),
    M = 0,
    IG = rep(0, n2),
    G = 0,
    IGflux = 0,
    Gflux = 0
  )
  
  times <- seq(0, 100, by = 1 / 10)
  
  out_DDE <- ode(
    y = inits_n, times = times,
    func = Erlang_Malaria_Cut_GFLUX,
    parms = parameters_n
  )
  
  return(data.frame(out_DDE[, c("time","R","IGflux","Gflux")], 
                    B_V = B_V, C_V = C_V,endtime = endtime))
}
