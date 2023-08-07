average_life_span_RC <- function(full_split, lifespan) {
   
  full_lifespan <- NULL
    for (k in seq(1, length(lifespan))) {
    lifespan_interest <- lifespan[[k]]
    
    ###
    gam_lxmx <- do.call(
      rbind.data.frame,
      lapply(
        full_split,
        Gametocyte_Fitness_EqualWeighting  ,
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

    ###Replicative capacity for the best strategy
    tmp_fit <- gam_lxmx[which.max(round(gam_lxmx$end_fitness,3)), ]
    tmp_fit$RC <- tmp_fit$B_V * (1 - tmp_fit$C_V)
    tmp_fit$lifespan <- lifespan_interest
    tmp_fit$group <- 'fit'
   
    
    fixed_optimal_CV_df <- subset(gam_lxmx, gam_lxmx$C_V %in%   tmp_fit $C_V)
    fixed_optimal_CV_df[which.max(fixed_optimal_CV_df$B_V),]
   
    tmp_max <-  fixed_optimal_CV_df[which.max(fixed_optimal_CV_df$B_V),]
    tmp_max$RC <-  tmp_max$B_V * (1 - tmp_max$C_V)
    tmp_max$lifespan <- lifespan_interest
    tmp_max$group <- 'max'
     
     
    tmp_all <- rbind(tmp_fit, tmp_max)
    
   
    
  full_lifespan[[k]] <- tmp_all 
    
  }

  full_lifespan_f <- do.call(rbind, full_lifespan)

  return(full_lifespan_f)
}

