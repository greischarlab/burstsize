

###These are the three functions for simulating the malaria infections
###with changes to the parameters (muM (merozoite mortality ), lambda
###(red blood cell replenishment), and RI ###(initla red blood cell)).
###These are not run directly but ran through the Sens_Analysis function

###Input: B_V And C_V Are the vector of burst size and transmission
###investment specifically and the percent modifier is how much the original
###parameter is varied by (default is 25%). The output from these models 
###are the RBC/gam abundance over time.
rootfun <- function (t, y, parms) {return(y['R'] - 6.5*10^5)}


Simulator_Malaria_BC_MuM <- function(B_V, C_V, initial_value, change) {

  params_muM <-
    c(lambda = 370000, # Replenishment rate of RBC (#SimulatedTimeSeries.R)
      K = 19968254, # Carrying capacity of RBC population in the absence of mortality
      pmax = 4.0e-6, # Rate of infection (From American Naturalist- Greischar et al. 2014)
      muR = 0.025, # Daily mortality rate of red blood cells
      muI = 0.025, # Daily mortality rate of infected red blood cells
      c = C_V, # Transmission investment (THE VARYING FACTOR)
      B = B_V, # The burst size (THE VARYING FACTOR)
      alpha1 = 1, # The rate of development of parasite in iRBC
      alpha2 = 1 / 2, # The rate of development
      muM = change, # Background mortality of the merozoite
      muG = 4, # background mortality of the immature/mature gametocytes
      n1 = 100, # shape parameter controlling the variability in asexual devleopment
      n2 = 100 # shape parameter controlling the variability in sexual development
    )


  n1 <-  params_muM[['n1']] 
  ### The number of subcompartments for immature gametocytes
  n2 <-  params_muM[['n2']]
  
  ### The initial numbers
  inits_n <- c(
    R = 8500000, # RBC
    I = rep(initial_value/ n1, n1), # Infected RBC- note the uniform distribution
    M = 0, # merozoite
    IG = rep(0, n2), # immature gametocytes
    G = 0
  ) # gametocytes

  times <- seq(0, 100, by = 1 / 10)

  out_DDE_mu_M <-
    ode(
      y = inits_n,
      times = times,
      func = Erlang_Malaria,
      parms = params_muM,
      rootfun = rootfun)



  df_mu_M <-  data.frame(out_DDE_mu_M[, c("time", "R", "G")], 
                                B_V = B_V, 
                                C_V = C_V,
                                initialvalue = initial_value,
                                change = change,
                                infection_length = 
                                  ifelse(!is.null(attributes(out_DDE_mu_M )$troot),
                                         attributes(out_DDE_mu_M )$troot,
                                         NA))
  



  return(df_mu_M)
}

Simulator_Malaria_BC_Lambda <- function(B_V, C_V, initial_value, change) {
 
  params_lambda <-
    c(
      lambda = change, # Replenishment rate of RBC (#SimulatedTimeSeries.R)
      K = (change* 8500000) / (change - (0.025 * 8500000)), # Carrying capacity of RBC population in the absence of mortality
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
      n2 = 100 # shape parameter controlling the variability in sexual development
    )

  
  n1 <- params_lambda[['n1']] 
  ### The number of subcompartments for immature gametocytes
  n2 <- params_lambda[['n2']]
  
  
   ### The initial numbers
  inits_n <- c(
    R = 8500000, # RBC
    I = rep(initial_value / n1, n1), # Infected RBC- note the uniform distribution
    M = 0, # merozoite
    IG = rep(0, n2), # immature gametocytes
    G = 0
  ) # gametocytes
  
  
  times <- seq(0, 100, by = 1 / 10)
  
  
  out_DDE_L <-
    ode(
      y = inits_n,
      times = times,
      func = Erlang_Malaria,
      parms = params_lambda,
      rootfun = rootfun)
  
  
  
  df_L<- data.frame(out_DDE_L[, c("time", "R", "G")], 
                    B_V = B_V, 
                    C_V = C_V,
                    initialvalue = initial_value,
                    change = change,
                    infection_length = 
                      ifelse(!is.null(attributes(out_DDE_L)$troot),
                             attributes(out_DDE_L)$troot,
                             NA))
  
  
  
  return(df_L)
}
Simulator_Malaria_BC_RI <- function(B_V, C_V, initial_value, change) {

  params_RI <-
    c(
      lambda = 370000, # Replenishment rate of RBC (#SimulatedTimeSeries.R)
      K = (370000 * change) / (370000 - (0.025 * change)), # Carrying capacity of RBC population in the absence of mortality
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
      n2 = 100 # shape parameter controlling the variability in sexual development
    )
  
  n1 <- params_RI[['n1']] 
  ### The number of subcompartments for immature gametocytes
  n2 <- params_RI[['n2']]
  

  ### The initial numbers
  inits_RI <- c(
    R = change,
    I = rep(initial_value / n1 , n1),
    M = 0,
    IG = rep(0,n2),
    G = 0
  )

  times <- seq(0, 100, by = 1 / 10)
  

  out_DDE_RI <-
    ode(
      y = inits_RI,
      times = times,
      func = Erlang_Malaria,
      parms = params_RI,
      rootfun = rootfun
      
    )

 

  df_RI <- data.frame(out_DDE_RI[, c("time", "R", "G")], 
                      B_V = B_V, 
                      C_V = C_V,
                      initialvalue = initial_value,
                      change= change,
                      infection_length = 
                        ifelse(!is.null(attributes(out_DDE_RI)$troot),
                               attributes(out_DDE_RI)$troot,
                               NA))
  return(df_RI)
}



###These are the three functions for simulating the malaria infections
###with changes to the parameters (muM (merozoite mortality ), lambda
###(red blood cell replenishment), and RI ###(initla red blood cell)).
###These are not run directly but ran through the Sens_Analysis function

###Input: B_V And C_V Are the vector of burst size and transmission
###investment specifically and the percent modifier is how much the original
###parameter is varied by (default is 25%). The output from these models 
###are the RBC/gam abundance over time.


Simulator_MalariaPC_DDE_BC_MuM_Cut <- function(B_V, C_V, change,initial_value, endtime) {
  parameters_muM <-
    c(
      lambda = 370000, # Replenishment rate of RBC (#SimulatedTimeSeries.R)
      K = 19968254, # Carrying capacity of RBC population in the absence of mortality
      pmax = 4.0e-6, # Rate of infection (From American Naturalist- Greischar et al. 2014)
      muR = 0.025, # Daily mortality rate of red blood cells
      muI = 0.025, # Daily mortality rate of infected red blood cells
      c = C_V, # Transmission investment (THE VARYING FACTOR)
      B = B_V, # The burst size (THE VARYING FACTOR)
      alpha1 = 1, # The rate of development of parasite in iRBC
      alpha2 = 1 / 2, # The rate of development
      muM = change, # Background mortality of the merozoite
      muG = 4, # background mortality of the immature/mature gametocytes
      n1 = 100, # shape parameter controlling the variability in asexual devleopment
      n2 = 100, # shape parameter controlling the variability in sexual development
      endtime = endtime # When the asexual replication end
    )
  
  ### The number of subcompartments for infected rbc(n1)
  n1 <- parameters_muM["n1"]
  ### The number of subcompartments for immature gametocytes
  n2 <- parameters_muM["n2"]
  
  ### The initial numbers
  inits_n <- c(
    R = 8500000,
    I = rep(initial_value / n1, n1),
    M = 0,
    IG = rep(0, n2),
    G = 0
  )
  
  times <- seq(0, 100, by = 1 / 10)
  
  out_DDE <- ode(
    y = inits_n, times = times,
    func = Erlang_Malaria_Cut,
    parms = parameters_muM 
  )
  
  return(data.frame(out_DDE[, c("time", "R", "G")], B_V = B_V, C_V = C_V,
                    change =  change))
}

Simulator_MalariaPC_DDE_BC_Lambda_Cut <- function(B_V, C_V, change,initial_value, endtime) {
  parameters_Lambda <-
    c(
      lambda = change, # Replenishment rate of RBC (#SimulatedTimeSeries.R)
      K = (change * 8500000) / (change - (0.025 * 8500000)), # Carrying capacity of RBC population in the absence of mortality
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
  n1 <- parameters_Lambda["n1"]
  ### The number of subcompartments for immature gametocytes
  n2 <- parameters_Lambda["n2"]
  
  ### The initial numbers
  inits_n <- c(
    R = 8500000,
    I = rep(initial_value / n1, n1),
    M = 0,
    IG = rep(0, n2),
    G = 0
  )
  
  times <- seq(0, 100, by = 1 / 10)
  
  out_DDE <- ode(
    y = inits_n, times = times,
    func = Erlang_Malaria_Cut,
    parms = parameters_Lambda
  )
  
  return(data.frame(out_DDE[, c("time", "R", "G")], B_V = B_V, C_V = C_V,
                    change =  change))
}
Simulator_MalariaPC_DDE_BC_RI_Cut <- function(B_V, C_V, change, initial_value, endtime) {
 
  parameters_RI <-
    c(
      lambda = 370000, # Replenishment rate of RBC (#SimulatedTimeSeries.R)
      K = (370000 * change) / (370000 - (0.025 * change)), # Carrying capacity of RBC population in the absence of mortality
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
      n2 = 100, # shape parameter controlling the variability in sexual development,
      endtime = endtime
    )
  
  
  ### The number of subcompartments for infected rbc(n1)
  n1 <- parameters_RI["n1"]
  ### The number of subcompartments for immature gametocytes
  n2 <- parameters_RI["n2"]
  
  ### The initial numbers
  inits_RI <- c(
    R = change,
    I = rep(initial_value / n1 , n1),
    M = 0,
    IG = rep(0,n2),
    G = 0
  )
  
  times <- seq(0, 100, by = 1 / 10)
  
  out_DDE <- ode(
    y = inits_RI, times = times,
    func = Erlang_Malaria_Cut,
    parms =  parameters_RI
  )
  
  return(data.frame(out_DDE[, c("time", "R", "G")], B_V = B_V, C_V = C_V,
                    change =  change))
}


#####################
### Runs everything###
#####################
Sens_Analysis <- function(B_V, C_V,initalvalue,change, sens_var) {
  function_sens_var <-
    switch(sens_var,
      "UM" = Simulator_Malaria_BC_MuM,
      "lambda" = Simulator_Malaria_BC_Lambda,
      "RI" = Simulator_Malaria_BC_RI
    )


  return(function_sens_var(B_V, C_V, initalvalue,change))
}



###This is the function for figuring out which burst size and transmission
###investment can lead to establishing of infections depending on what the
###percent modifier 
###Input: B_V And C_V Are the vector of burst size and transmission
###investment specifically and the percent modifier is how much the original
###parameter is varied by (default is 25%). The output from these models 
###are the RBC/gam abundance over time.

B_V_C_V_Establisher <- function(sens_var,percent_modifier) {
 
  B_V <- seq(1, 50, 0.5) # Burst size
  C_V <- seq(.01, 1, 0.01) # Transmission investment
  mu_M_vec <- c(48 * (1 - percent_modifier), 48 * (1 + percent_modifier)) ###The different merozoite mortality rate
  RI_vec <-  c(8500000 * (1 - percent_modifier), 8500000 *  (1 + percent_modifier)) ###The different RI rate 
  L_vec <- c(277500,462500)
  
  B_V_C_V_UM <- expand.grid(B_V = B_V, C_V = C_V, mu_M = mu_M_vec) # Different combinations
  B_V_C_V_RI <- expand.grid(B_V = B_V, C_V = C_V, RI = RI_vec)
  B_V_C_V_L <- expand.grid(B_V = B_V, C_V = C_V, L = L_vec)
  
  ### I am able to figure out which of the burst and transmission
  ### are unable to establish the infection so that I can exclude from
  ### simulation (the acute time is set to 0 and transmission potential is set to 0)
  
  p_val <- 4.0e-6 #Same for everyone
  mu_M <- 48 #Only for RI and L 
  R_val <- 8500000 #Only for L and mu M
  initial_RM_modifier <- 1.5 #
  
  
  RM_limit_UM <- (1/((100* (1))/(100 + 0.025))^100)* initial_RM_modifier*((p_val * R_val) + mu_M_vec) / (p_val * R_val)
  RM_limit_RI <- (1/((100* (1))/(100 + 0.025))^100)*initial_RM_modifier*((p_val * RI_vec) + mu_M ) / (p_val * RI_vec)
  RM_limit <- (1/((100* (1))/(100 + 0.025))^100)*initial_RM_modifier*((p_val * R_val) + mu_M ) / (p_val * R_val)
  
  
  limit <- switch(sens_var,
    "UM" = RM_limit_UM, #two value
    "RI" = RM_limit_RI, #two value
    "L" = RM_limit #one value
  )
  
  BVCV <- switch(sens_var,
    "UM" = B_V_C_V_UM,
    "RI" = B_V_C_V_RI,
    "L" = B_V_C_V_L
  )

  num_index <- switch(sens_var,
    "UM" = 2,
    "RI" = 2,
    "L" = 1
  )

  BVCV_temp <- split(BVCV, BVCV[, 3])
  
  B_V_interest <- BVCV_temp[[1]]$B_V
  C_V_interest <- BVCV_temp[[1]]$C_V
  
  
  BVCV_temp[[1]]$Establish <- ifelse((1 - C_V_interest ) *   B_V_interest >=
    limit[[1]], "Establish", "Fail")

  if (num_index == 2) {
    BVCV_temp[[2]]$Establish <- ifelse((1 -   C_V_interest ) *  B_V_interest  >=
      limit[[2]],"Establish", "Fail"
    )

    return(rbind.data.frame(BVCV_temp[[1]], BVCV_temp[[2]]))
  
  } 
  else{
    
    BVCV_temp[[2]]$Establish <- ifelse((1 - C_V_interest) *   B_V_interest >=
                                         limit[[1]],
                                       "Establish", "Fail"
    )
    
    
    return(rbind.data.frame(BVCV_temp[[1]],  BVCV_temp[[2]]))
  }
}


###THis is the function that looks at the successful infections (aka infections
###That don't fail to establish) and see if it kills the host. If it kills
###the host, we calculate when the host dies (when the RBC dip below the 
###mortality threshold and calculate the cumulative transmission to the time of
###death). However, if the infection does not kill, then we keep the list empty
###to be populated by the next function- we run this function to figure out
###what the length of the acute phase is



###For the sensitivity analysis finding the Finder
###This function is for finding the fitness of the successful infections
###that do not establish infection  
Fitness_Finder_SA <- function(BC_vec,
                              Duration,
                              initial_value,
                              sens_var) {
  sim <-
    switch(sens_var,
      "UM" = Simulator_MalariaPC_DDE_BC_MuM_Cut,
      "lambda" = Simulator_MalariaPC_DDE_BC_Lambda_Cut,
      "RI" = Simulator_MalariaPC_DDE_BC_RI_Cut
    )


  ### I know which B_V/C_V combos would lead to failed infections so,
  ### I create a data.frame for them
  Failed_B_V_C_V <- subset(BC_vec, BC_vec$Establish == "Fail")

  Duration_FAIL <-
    data.frame(
      endtime = 0,
      up_down = 0,
      end_fitness = 0,
      status = "Fail",
      B_V = Failed_B_V_C_V$B_V,
      C_V = Failed_B_V_C_V$C_V,
      change = Failed_B_V_C_V[, 3]
    )

  ### These are the B_V/C_V that would lead to mortality,
  ### This means that we know the end time (point of death) and the
  ### cumulative transmission potential

  ### If it's not a success, then it's death
  Duration_MORT <- subset(
    Duration,
    Duration$status != "success"
  )


  ### If it's a success, then the host survived
  Duration_SUCCESS <- subset(
    Duration,
    Duration$status == "success"
  )


  ### Only run the function for the burst size and transmission investment
  ### combination that leads to successful infection that does not induce
  ### host mortality

  Fitness_MODEL <-
    mcmapply(
      sim,
      c(Duration_SUCCESS$B_V),
      c(Duration_SUCCESS$C_V),
      c(Duration_SUCCESS$change),
      c(initial_value),
      c(Duration_SUCCESS$endtime),
      mc.cores = 1,
      SIMPLIFY = FALSE
    )

  ### Now we know what the fitness is
  Duration_SUCCESS$end_fitness <-
    unlist(lapply(
      Fitness_MODEL,
      Gametocyte_Fitness
    ))

  ### These are the fitness model data.frame that should work
  Fitness_MODEL_FULL <-
    rbind.data.frame(
      Duration_SUCCESS,
      Duration_FAIL,
      Duration_MORT
    )

  return(Fitness_MODEL_FULL)
}
