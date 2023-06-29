

Best_Replicative_Capacity_Finder <- function(x){

full_lifespan <- NULL

for(k in seq(1, length(lifespan))){
  lifespan_interest <- lifespan[[k]]
  
  gam_lxmx <- do.call(rbind.data.frame,
                      lapply(FULL_MODEL_100,
                             Gametocyte_Fitness_LXMX, 
                             lifespan_interest))
  
  colnames(gam_lxmx) <- "end_fitness"
  
  gam_lxmx$B_V <- do.call(rbind,lapply(FULL_MODEL_100, function (x) unique(x$B_V)))
  gam_lxmx$C_V <-do.call(rbind,lapply(FULL_MODEL_100, function (x) unique(x$C_V)))
  
  tmp_max <- gam_lxmx[which.max(gam_lxmx$end_fitness),]
  
  tmp_max$mer <- tmp_max$B_V * tmp_max$C_V
  
  tmp_max$lifespan <-lifespan_interest 
  
  full_lifespan[[k]] <- tmp_max
  
}
full_lifespan_f <- do.call(rbind, full_lifespan)

}

PEAK_IDENTIFIER <- function(){
  tmp <- do.call(rbind,lapply(FULL_MODEL_100_Split, function(x) x[which.min(x$R),]))
  tmp$id <-  paste(tmp$B_V,"_", tmp$C_V)
  return(tmp)
  
  
  
}
