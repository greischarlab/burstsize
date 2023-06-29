average_life_span_RC <- function(full_split, weighting_var) {
  
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

    tmp_max <- gam_lxmx[which.max(gam_lxmx$end_fitness), ]
    tmp_max$RC <- tmp_max$B_V * (1 - tmp_max$C_V)
    tmp_max$lifespan <- lifespan_interest
    tmp_max$id <- weighting_var
    full_lifespan[[k]] <- tmp_max
  }

  full_lifespan_f <- do.call(rbind, full_lifespan)

  return(full_lifespan_f)
}
