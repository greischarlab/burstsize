Surface_Plot_Standard_Duration<- function(x_list,time){
  
  list_20 <- lapply(x_list, function(x) subset(x, x$time<= 20))
  
  gam_Fit <- do.call(rbind.data.frame,lapply(list_20, Gametocyte_Fitness))
 colnames(gam_Fit) <- 'end_fitness'
  gam_Fit$B_V <- do.call(rbind, lapply(x_list, function(x) unique(x$B_V)))
  gam_Fit$C_V <- do.call(rbind, lapply(x_list, function(x) unique(x$C_V)))
  
  return( gam_Fit)
}
