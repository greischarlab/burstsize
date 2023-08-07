####################################################
###Testing script to see what the initial RMs are###
####################################################
p_val <- 4.0e-6 
mu_M <- 48
R_val <- 8500000


initial_RM_calc <- function(B_V, C_V){
 return((1 - C_V) * B_V * ((p_val*R_val)/((p_val*R_val) + mu_M)))
}


### Burst Size Versus Transmission Investment ###
B_V <- seq(1, 50, 0.5) # Burst size
C_V <- seq(.01, 1, 0.01) # Transmission investment
B_V_C_V <- expand.grid(B_V = B_V, C_V = C_V) 

initialrate <- data.frame(mapply(initial_RM_calc, 
       B_V = c(B_V_C_V$B_V),
       C_V = c(B_V_C_V$C_V)))

colnames(initialrate)[1] <- "initial_RM"
initialrate$B_V <- B_V_C_V$B_V
initialrate$C_V <- B_V_C_V$C_V


ggplot(initialrate, aes(x = B_V, y = C_V, fill = initial_RM))+
  geom_raster()+
  scale_fill_viridis(name = "Initial RM", option= 'turbo')+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab("Burst size")+
  ylab("Transmission investment")


ggsave(here("Figures","Miscellaneous Figures","initialrate_surfaceplot.pdf"), 
       width = 7, height = 6, units = 'in')
