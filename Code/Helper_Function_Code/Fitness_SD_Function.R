
###Function that allows one to 

###x_list is the full model list 
Calculate_Set_Duration_CP <- function(x_list,
                                      params, 
                                      vec_select, 
                                      lxmx) {
  
  
  
  lxmx_yes = switch(lxmx,
    'yes' = Gametocyte_Fitness_LXMX,
    'no' = Gametocyte_Fitness)
  
  lxmx_yes_RM = switch(lxmx,
                    'yes' = Finder_RM_lxmx,
                    'no' =Finder_RM)
  
  
  ### Take the full model and truncate it based on the maximum
  ### point of interest
  max_time_interest <- lapply(
    x_list,
    function(x) {
      subset(
        as.data.frame(x),
        as.data.frame(x)$time <= max(vec_select)
      )
    }
  )
  ### For the different time point of interest, we're going to
  ### calculate the cumulative transmission potential for the
  ### infections that don't kill their hosts

  CP_LIST <- NULL

  for (k in seq(1, length(vec_select))) {
    ### truncates the time point to the time of interest
    list_time_interest <-
       lapply(max_time_interest,
        function(x) {
        subset(
          as.data.frame(x),
          as.data.frame(x)$time <=
            vec_select[[k]]
        )
      }
    )
    
    
    BV_CV <- do.call(rbind,
            lapply(list_time_interest , function (x) unique(x[,c("B_V", "C_V")])))
    

    ######################################################
    ### look to make sure which leads to mortality or not#
    ######################################################
  
    duration_interest <- do.call(
      rbind,
      mclapply(list_time_interest,
        lxmx_yes_RM ,
        mu_M_c = 48,
        mc.cores = 2
      )
    )

    duration_interest$B_V <- BV_CV$B_V
    duration_interest$C_V <- BV_CV$C_V
    

    ### If it does not lead to mortality and is NA, make it success
    duration_interest$status[is.na(duration_interest$status) == TRUE] <- "success"

    ### Figure out what rows lead to successful
    interest_success <- which(duration_interest$status == "success")

    ### Get the B_V AND C_V
    burst_size_transmission_investment <- subset(
      duration_interest,
      duration_interest$status ==
        "success"
    )[, c("B_V", "C_V")]
    
    ##################################
    ### Pull out from the full model##
    ###################################
    list_success_interest <- list_time_interest[interest_success]

    ##############
    ####SUCCESS###
    ##############
    Fitness.TimePoint_DF <- data.frame(
      burst_size_transmission_investment,
      status = "success",
      end_time = vec_select[[k]],
      fitness =
        do.call(
          rbind,
          lapply(
            list_success_interest, lxmx_yes
          )
        )
    )
    ### Let's get information for burst size and transmission investment
    ### that lead to host mortality
    Fitness.Death <- subset(
      duration_interest,
      duration_interest$status == "mort"
    )

    Fitness.Death_DF <- data.frame(
      B_V = Fitness.Death$B_V,
      C_V = Fitness.Death$C_V,
      status = "mort",
      end_time = vec_select[[k]],
      fitness = Fitness.Death$end_fitness
    )

    ### Let's get information for burst size and transmission investment
    ### that lead to no infection
    
    
    ########################
    ###FAILED INFECTIONS ###
    ########################
    Fitness.Fail <- subset(
      params,
      params$Establish == "Fail"
    )

    Fitness.Fail_DF <- data.frame(
      B_V = Fitness.Fail$B_V,
      C_V = Fitness.Fail$C_V,
      status = "Fail",
      end_time = vec_select[[k]],
      fitness = 0
    )


    CP_LIST[[k]] <- rbind.data.frame(
      Fitness.TimePoint_DF,
      Fitness.Death_DF,
      Fitness.Fail_DF
    )
  }

  return(CP_LIST)
}


