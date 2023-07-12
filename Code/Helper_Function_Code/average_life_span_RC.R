average_life_span_RC <- function(full_split, weighting_var,lifespan) {
  
  Weighting_Function <-
    switch(weighting_var,
      "equal" = Gametocyte_Fitness_EqualWeighting,
      "exponential" = Gametocyte_Fitness_LXMX)

  full_lifespan <- NULL

  for (k in seq(1, length(lifespan))) {
    lifespan_interest <- lifespan[[k]]

    gam_lxmx <- do.call(
      rbind.data.frame,
      lapply(
        full_split,
        Weighting_Function,
        lifespan_interest
      )
    )

    colnames(gam_lxmx) <- "end_fitness"

    gam_lxmx$B_V <- do.call(rbind, lapply(full_split, function(x) unique(x$B_V)))
    gam_lxmx$C_V <- do.call(rbind, lapply(full_split, function(x) unique(x$C_V)))
    
    gam_lxmx$id <- do.call(
      rbind,
      lapply(
        full_split,
        function(x) {
          paste(unique(x$B_V), "_", unique(x$C_V))
        }
      )
    )

    tmp_fit <- gam_lxmx[which.max(round(gam_lxmx$end_fitness,4)), ]
    tmp_fit$RC <- tmp_fit$B_V * (1 - tmp_fit$C_V)
    tmp_fit$lifespan <- lifespan_interest
    tmp_fit$group <- 'fit'
   
    
    fixed_optimal_CV_df <- subset(gam_lxmx, gam_lxmx$C_V %in%   tmp_fit $C_V)
    fixed_optimal_CV_df[which.max(fixed_optimal_CV_df$B_V),]
   
    tmp_max <-  fixed_optimal_CV_df[which.max(fixed_optimal_CV_df$B_V),]
    tmp_max$RC <-  tmp_max$B_V * (1 - tmp_max$C_V)
    tmp_max$lifespan <- lifespan_interest
    tmp_max$group <- 'max'
     
    tmp_min <-  fixed_optimal_CV_df[which.min(fixed_optimal_CV_df$B_V),]
    tmp_min$RC <-  tmp_min$B_V * (1 - tmp_min$C_V)
    tmp_min$lifespan <- lifespan_interest
    tmp_min$group <- 'min'
      
     
    tmp_all <- rbind(tmp_fit, tmp_max, tmp_min)
    
   
    
  full_lifespan[[k]] <- tmp_all 
    
  }

  full_lifespan_f <- do.call(rbind, full_lifespan)

  return(full_lifespan_f)
}

