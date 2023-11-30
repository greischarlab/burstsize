average_life_span_RC <- function(full_split, lifespan) {
   
  full_lifespan <- NULL
    for (k in seq(1, length(lifespan))) {
    lifespan_interest <- lifespan[[k]]
    
    ###
    gam_lxmx <- do.call(
      rbind.data.frame,
      lapply(
        full_split,
        Gametocyte_Fitness_EqualWeighting ,
        lifespan_interest
      )
    )

    colnames(gam_lxmx)[1] <- "end_fitness"

    ###Replicative capacity for the best strategy
    tmp_fit <- gam_lxmx[which.max(round(gam_lxmx$end_fitness,3)), ]
    tmp_fit$lifespan <- lifespan_interest
    tmp_fit$group <- 'fit'
   
    tmp_all <- tmp_fit
    
   
    
  full_lifespan[[k]] <- tmp_all 
    
  }

  full_lifespan_f <- do.call(rbind, full_lifespan)

  return(full_lifespan_f)
}

average_life_span_ALL <- function(full_x, lifespan) {
  
 splitted_full <- split(full_x, list(full_x$B_V, full_x$C_V), drop = TRUE)
  
  full_lifespan <- NULL
  
  for (k in seq(1, length(lifespan))) {
    lifespan_interest <- lifespan[[k]]
    
    ###
    gam_lxmx <- do.call(
      rbind.data.frame,
      lapply(
        splitted_full ,
        Gametocyte_Fitness_EqualWeighting ,
        lifespan_interest
      )
    )
    
    
    full_lifespan[[k]] <- gam_lxmx
    
  }
  
  full_lifespan_f <- do.call(rbind, full_lifespan)
  
  return(full_lifespan_f)
}
