

###These are the three functions for simulating the malaria infections
###with changes to the parameters (muM (merozoite mortality ), lambda
###(red blood cell replenishment), and RI ###(initla red blood cell)).
###These are not run directly but ran through the Sens_Analysis function

###Input: B_V And C_V Are the vector of burst size and transmission
###investment specifically and the percent modifier is how much the original
###parameter is varied by (default is 25%). The output from these models 
###are the RBC/gam abundance over time.
Simulator_Malaria_BC_MuM <- function(B_V, C_V, percent_modifier) {

  params_muM_1 <-
    c(lambda = 370000, # Replenishment rate of RBC (#SimulatedTimeSeries.R)
      K = 19968254, # Carrying capacity of RBC population in the absence of mortality
      pmax = 2.5e-6, # Rate of infection (From American Naturalist- Greischar et al. 2014)
      muR = 0.025, # Daily mortality rate of red blood cells
      muI = 0.025, # Daily mortality rate of infected red blood cells
      c = C_V, # Transmission investment (THE VARYING FACTOR)
      B = B_V, # The burst size (THE VARYING FACTOR)
      alpha1 = 1, # The rate of development of parasite in iRBC
      alpha2 = 1 / 2, # The rate of development
      muM = mu_M_vec[[1]], # Background mortality of the merozoite
      muG = 4, # background mortality of the immature/mature gametocytes
      n1 = 100, # shape parameter controlling the variability in asexual devleopment
      n2 = 100 # shape parameter controlling the variability in sexual development
    )

  params_muM_2 <-
    c(lambda = 370000, # Replenishment rate of RBC (#SimulatedTimeSeries.R)
      K = 19968254, # Carrying capacity of RBC population in the absence of mortality
      pmax = 2.5e-6, # Rate of infection (From American Naturalist- Greischar et al. 2014)
      muR = 0.025, # Daily mortality rate of red blood cells
      muI = 0.025, # Daily mortality rate of infected red blood cells
      c = C_V, # Transmission investment (THE VARYING FACTOR)
      B = B_V, # The burst size (THE VARYING FACTOR)
      alpha1 = 1, # The rate of development of parasite in iRBC
      alpha2 = 1 / 2, # The rate of development
      muM = mu_M_vec[[2]], # Background mortality of the merozoite
      muG = 4, # background mortality of the immature/mature gametocytes
      n1 = 100, # shape parameter controlling the variability in asexual devleopment
      n2 = 100 # shape parameter controlling the variability in sexual development
    )

  ### The number of subcompartments for infected rbc(n1)
  n1 <-  params_muM_1[['n1']] 
  ### The number of subcompartments for immature gametocytes
  n2 <-   params_muM_1[['n2']]
  
  ### The initial numbers
  inits_n <- c(
    R = 8500000, # RBC
    I = rep(43859.65 / n1, n1), # Infected RBC- note the uniform distribution
    M = 0, # merozoite
    IG = rep(0, n2), # immature gametocytes
    G = 0
  ) # gametocytes

  times <- seq(0, 100, by = 1 / 10)

  out_DDE_1_mu_M <-
    ode(
      y = inits_n,
      times = times,
      func = Erlang_Malaria,
      parms = params_muM_1
    )

  out_DDE_2_mu_M <-
    ode(
      y = inits_n,
      times = times,
      func = Erlang_Malaria,
      parms = params_muM_2
    )


  df_1_mu_M <- data.frame(out_DDE_1_mu_M[, c("time", "R", "G")],
    B_V = B_V, C_V = C_V, change =   mu_M_vec[[1]]
  )
  df_2_mu_M <- data.frame(out_DDE_2_mu_M[, c("time", "R", "G")],
    B_V = B_V, C_V = C_V, change =   mu_M_vec[[2]]
  )

  out_list_mu_M <- rbind(df_1_mu_M, df_2_mu_M)


  return(out_list_mu_M)
}

Simulator_Malaria_BC_Lambda <- function(B_V, C_V, percent_modifier) {
  lambda_control <- 370000
  lambda_vec <- c(lambda_control * (1 - percent_modifier), lambda_control * (1 + percent_modifier))

  params_lambda_1 <-
    c(
      lambda = lambda_vec[[1]], # Replenishment rate of RBC (#SimulatedTimeSeries.R)
      K = (lambda_vec[[1]] * 8500000) / (lambda_vec[[1]] - (0.025 * 8500000)), # Carrying capacity of RBC population in the absence of mortality
      pmax = 2.5e-6, # Rate of infection (From American Naturalist- Greischar et al. 2014)
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

  params_lambda_2 <-
    c(
      lambda = lambda_vec[[2]], # Replenishment rate of RBC (#SimulatedTimeSeries.R)
      K = (lambda_vec[[2]] * 8500000) / (lambda_vec[[2]] - 0.025 * 8500000), # Carrying capacity of RBC population in the absence of mortality
      pmax = 2.5e-6, # Rate of infection (From American Naturalist- Greischar et al. 2014)
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

  inits_n <- c(
    R = 8500000,
    I = rep(
      (43859.65 / params_lambda_1[["n1"]]),
      params_lambda_1[["n1"]]
    ),
    M = 0,
    IG = rep(0, params_lambda_1[["n2"]]),
    G = 0
  )

  times <- seq(0, 100, by = 1 / 10)

  out_DDE_1_lambda_1 <-
    ode(
      y = inits_n,
      times = times,
      func = Erlang_Malaria,
      parms = params_lambda_1
    )

  out_DDE_2_lambda_2 <-
    ode(
      y = inits_n,
      times = times,
      func = Erlang_Malaria,
      parms = params_lambda_2
    )


  df_1_lambda <- data.frame(out_DDE_1_lambda_1[, c("time", "R", "G")], B_V = B_V, C_V = C_V, change = lambda_vec[[1]])
  df_2_lambda <- data.frame(out_DDE_2_lambda_2[, c("time", "R", "G")], B_V = B_V, C_V = C_V, change = lambda_vec[[2]])

  out_list_lambda <- rbind(df_1_lambda, df_2_lambda)


  return(out_list_lambda)
}
Simulator_Malaria_BC_RI <- function(B_V, C_V, percent_modifier) {
  RI_control <- 8500000
  RI_vec <- c(RI_control * (1 - percent_modifier), RI_control * (1 + percent_modifier))

  params_RI_1 <-
    c(
      lambda = 370000, # Replenishment rate of RBC (#SimulatedTimeSeries.R)
      K = (370000 * RI_vec[[1]]) / (370000 - (0.025 * RI_vec[[1]])), # Carrying capacity of RBC population in the absence of mortality
      pmax = 2.5e-6, # Rate of infection (From American Naturalist- Greischar et al. 2014)
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

  params_RI_2 <-
    c(
      lambda = 370000, # Replenishment rate of RBC (#SimulatedTimeSeries.R)
      K = (370000 * RI_vec[[2]]) / (370000 - (0.025 * RI_vec[[2]])), # Carrying capacity of RBC population in the absence of mortality
      pmax = 2.5e-6, # Rate of infection (From American Naturalist- Greischar et al. 2014)
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


  ### The initial numbers
  inits_RI_1 <- c(
    R = RI_vec[[1]],
    I = rep((43859.65 / params_RI_1[["n1"]]), params_RI_1[["n1"]]),
    M = 0,
    IG = rep(0, params_RI_1[["n2"]]),
    G = 0
  )

  inits_RI_2 <- c(
    R = RI_vec[[2]],
    I = rep((43859.65 / params_RI_1[["n1"]]), params_RI_1[["n1"]]),
    M = 0,
    IG = rep(0, params_RI_1[["n2"]]),
    G = 0
  )


  out_DDE_RI_1 <-
    ode(
      y = inits_RI_1,
      times = times,
      func = Erlang_Malaria,
      parms = params_RI_1
    )

  out_DDE_RI_2 <-
    ode(
      y = inits_RI_2,
      times = times,
      func = Erlang_Malaria,
      parms = params_RI_2
    )


  df_1_RI <- data.frame(out_DDE_RI_1[, c("time", "R", "G")], B_V = B_V, C_V = C_V, change = RI_vec[[1]])
  df_2_RI <- data.frame(out_DDE_RI_2[, c("time", "R", "G")], B_V = B_V, C_V = C_V, change = RI_vec[[2]])

  out_list_RI <- rbind(df_1_RI, df_2_RI)

  return(out_list_RI)
}



#####################
### Runs everything###
#####################
Sens_Analysis <- function(B_V, C_V, percent_modifier, sens_var) {
  function_sens_var <-
    switch(sens_var,
      "UM" = Simulator_Malaria_BC_MuM,
      "lambda" = Simulator_Malaria_BC_Lambda,
      "RI" = Simulator_Malaria_BC_RI
    )


  return(function_sens_var(B_V, C_V, percent_modifier))
}


B_V_C_V_Establisher <- function(sens_var, percent_modifier) {
 
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
  
  p_val <- 2.5e-6 #Same for everyone
  mu_M <- 48 #Only for RI and L 
  R_val <- 8500000 #Only for L Aand mu M
  
  RM_limit_UM <- ((p_val * R_val) + mu_M_vec) / (p_val * R_val)
  RM_limit_RI <- ((p_val * RI_vec) + mu_M ) / (p_val * RI_vec)
  RM_limit <- ((p_val * R_val) + mu_M ) / (p_val * R_val)
  
  
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


Finder_RM_SA <- function(x_list, sens_var) {
  
  ###Seperate out the list element by the minus or plus change
  tmp_x_list <- split(x_list, x_list$change)
  
  if (sens_var == "UM") { ###If the sens_var is UM then we need to account
    ###for the difference in mortality rates
    mu_M_c <- c(36,60)
  } else {                ###If the sens_var is not UM then we need to account
    mu_M_c <- 48
  }
  
  ### Shape parameter (to ensure we're getting the right columns)
  n1 <- 100
  
  ### Check first if the infection induces mortality
  ### Look at the red blood cells time series only
  
  RBC_TS <- lapply(tmp_x_list, function(x) {
    data.frame(
      time = x$time,
      RBC = x[, "R"]
    )
  })
  ### Look at the red blood cells to see when it is under 6.5 * 10^5
  
  index_time_mort <- lapply(
    RBC_TS,
    function(x) {
      x[which(x$RBC
              <= 6.5 * 10^5), ]
    }
  )
  
  ### If the index_time_mort is empty, then the host survived
  ### and we note that
  mort_list <- NULL ###The empty list that we use to figure out 
                    ### rather the host survive or not
  
  for (k in seq(1, 2)) {
    
    tmp <- index_time_mort[[k]]

    if (nrow(tmp) == 0) { #If it is empty that means that there is no 
                          #RBC that dips below the mortality threshold
      mort_time <- cbind.data.frame(
        surv = "surv",
        time = NA
      )
      
      mort_list[[k]] <- mort_time
    
   } else {
      RBC_TS_Func <- approxfun(
        RBC_TS[[k]][, "time"],
        RBC_TS[[k]][, "RBC"] - 6.5 * 10^5
      )
      
      death_time <- uniroot(RBC_TS_Func, c(0, tmp$time[1]))
      mort_time <- cbind.data.frame(
        surv = "mort",
        time = death_time$root
      )
      
      mort_list[[k]] <- mort_time
    }
  }
  ### This now checks the Acute Phase Duration for the Mice
  ### That Survived
  end_time_list <- NULL
  
  for (k in seq(1, 2)){
    tmp <- mort_list[[k]]
    
    tmp_x_list_df <- tmp_x_list[[k]]
    
    if (tmp$surv == "surv") {
      
      p <- 2.5e-6  #this is a parameter value that does not change 
      
      if (length(mu_M_c == 2)) {
     
      rate <- (1 - unique(tmp_x_list_df$C_V)) *
               unique(tmp_x_list_df$B_V) * ((tmp_x_list_df[, "R"] * p) / 
                                              ((p * tmp_x_list_df[, "R"]) + mu_M_c[[k]]))
      }
      else{
        
        rate <- (1 - unique(tmp_x_list_df$C_V)) *
          unique(tmp_x_list_df$B_V) * ((tmp_x_list_df[, "R"] * p) / 
                                         ((p * tmp_x_list_df[, "R"]) + mu_M_c))
      }
      
      RM_time_df <- cbind.data.frame(
        time = tmp_x_list_df[, "time"],
        rate = rate,
        B_V = unique(tmp_x_list_df$B_V),
        C_V = unique(tmp_x_list_df$C_V)
      )
      
      min_RM <- RM_time_df[which.min(RM_time_df$rate), ]
      
      end_time <- subset(
        RM_time_df,
        RM_time_df$time >= min_RM$time & RM_time_df$rate >= 1
      )[1, "time"] ###Endtime is the first time when the 
      
      end_time_list[[k]] <-
        data.frame(
          endtime = end_time,
          up_down = "up",
          end_fitness = NA,
          status = "success",
          B_V = unique(tmp_x_list_df$B_V),
          C_V = unique(tmp_x_list_df$C_V),
          change = unique(
            tmp_x_list_df$change)
        )
    } else {
      ### Look at the infected red blood cells time series only
      ### Calculating the fitness
      G_TS <- data.frame(
        time = tmp_x_list_df[, "time"],
        G = tmp_x_list_df[, "G"]
      )
      
      # Prevent negative gamaetocytes
      G_TS$G[G_TS$G < 0] <- 0
      
      truncate_G_TS <- subset(G_TS, G_TS$time <= mort_time$time)
      
      
      end_fitness_mort <- sum(PrI_PC(truncate_G_TS$G) * 1 / 10)
      
      end_time_list[[k]] <- 
        data.frame(
          endtime = tmp$time,
          up_down = "up",
          end_fitness = end_fitness_mort,
          status = "mort",
          B_V = unique(x_list$B_V),
          C_V = unique(x_list$C_V),
          change = unique( tmp_x_list_df$change)
          
        )
    }
  }
  return(do.call(rbind, end_time_list))
  
}


###For the sensitivity analysis finding the Finder


###This function is for finding the fitness of the successful infections
###that do not establish infection  

Fitness_Finder_SA <- function(BC_vec, 
                              sens_var, 
                              Duration) {

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
            change = Failed_B_V_C_V[,3]
          )
      
        ### These are the B_V/C_V that would lead to mortality,
        ### This means that we know the end time (point of death) and the
        ### cumulative transmission potential
      
        ###If it's not a success, then it's death
        Duration_MORT <- subset(
            Duration,
            Duration$status != "success"
          )
      
      
        ###If it's a success, then the host survived
        Duration_SUCCESS <- subset(
            Duration,
            Duration$status == "success"
          )
      
      
        ### Only run the function for the burst size and transmission investment
        ### combination that leads to successful infection that does not induce
        ### host mortality
      
        Fitness_MODEL <- 
          mcmapply(
          Simulator_MalariaPC_DDE_BC_Cut,
          c(Duration_SUCCESS$B_V),
          c(Duration_SUCCESS$C_V),
          c(Duration_SUCCESS$endtime),
          mc.cores = 4,
          SIMPLIFY = FALSE
        )
      
        ### Now we know what the fitness is
        Duration_SUCCESS$end_fitness <-
                 unlist(lapply(Fitness_MODEL, 
                               Gametocyte_Fitness))
      
        ### These are the fitness model data.frame that should work
        Fitness_MODEL_FULL <- 
          rbind.data.frame(
          Duration_SUCCESS,
          Duration_FAIL,
          Duration_MORT
        )
      
        return(Fitness_MODEL_FULL)
}
