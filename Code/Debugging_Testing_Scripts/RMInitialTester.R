###########################################
###Trying different initial RM modifier###
##########################################
library(here)


### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))

### Main modeling code
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_CUT.cpp"))
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_UNCUT.cpp"))

initialrate_modifier_BVCV <- function(modifier){
p_val <- 4.0e-6 
mu_M <- 48
R_val <- 8500000
RM_limit_1 <- modifier * ((p_val * R_val) + mu_M) / (p_val * R_val)

B_V <- seq(1, 50, 0.5) # Burst size
C_V <- seq(.01, 1, 0.01) # Transmission investment
B_V_C_V <- expand.grid(B_V = B_V, C_V = C_V) # Different combinations


B_V_C_V$Establish <- ifelse((1 - B_V_C_V$C_V) * B_V_C_V$B_V >= RM_limit_1,
                            "Establish", "Fail"
)

return(B_V_C_V)
}

BVCV_1.1<- initialrate_modifier_BVCV(1.1)
BVCV_1.5<- initialrate_modifier_BVCV(1.5)
BVCV_2.0<- initialrate_modifier_BVCV(2)
BVCV_3.0<- initialrate_modifier_BVCV(3)
BVCV_4.0<- initialrate_modifier_BVCV(4)

###Successful infection 
BVCV_1.1_establish <- subset(BVCV_1.1,BVCV_1.1$Establish != 'Establish')
BVCV_1.5_establish <- subset(BVCV_1.5,BVCV_1.5$Establish != 'Establish')
BVCV_2.0_establish <- subset(BVCV_2.0,BVCV_2.0$Establish != 'Establish')
BVCV_3.0_establish <- subset(BVCV_3.0,BVCV_3.0$Establish != 'Establish')
BVCV_4.0_establish <- subset(BVCV_4.0,BVCV_4.0$Establish != 'Establish')


############The original model ########################################
p_val <- 4.0e-6 
mu_M <- 48
R_val <- 8500000
initial_RM_modifier <- 1
RM_limit_1 <- initial_RM_modifier * ((p_val * R_val) + mu_M) / (p_val * R_val)

B_V_C_V$Establish <- ifelse((1 - B_V_C_V$C_V) * B_V_C_V$B_V >= RM_limit_1,
                            "Establish", "Fail"
)

### Simulate infections that are successful
B_V_C_V_F <- subset(B_V_C_V, B_V_C_V$Establish == "Establish")

### These infections are successful OR lead to host mortality
FULL_MODEL_PC_100Day <-   mcmapply(Simulator_Malaria_BC_ALTERNATIVE,
                          c(B_V_C_V_F$B_V),
                          c(B_V_C_V_F$C_V),
                          mc.cores = 5,
                          SIMPLIFY = FALSE
)
### Combine
FULL_MODEL_PC_DT <- do.call(rbind,FULL_MODEL_PC_100Day)

### Write into a CSV
write.csv(FULL_MODEL_PC_DT, file = here(
  "Output", "Full_Model",
  "FULL_MODEL_PC100_2_DT.csv"
))


Duration_Initial_PC_100 <- 
  do.call(
    rbind,
    mclapply(FULL_MODEL_PC_100Day,
             Finder_RM,
             mu_M_c = 48,
             mc.cores = 2
    )
  )


### Assign the burst size and transmission investment
Duration_Initial_PC_100 $B_V <- B_V_C_V_F$B_V
Duration_Initial_PC_100$C_V <- B_V_C_V_F$C_V


### I know which B_V/C_V combos would lead to failed infections so,
### I create a data.frame for them
Failed_B_V_C_V <- subset(B_V_C_V, B_V_C_V$Establish == "Fail")

Duration_Initial_PC_FAIL <-
  data.frame(
    endtime = 0,
    up_down = 0,
    end_fitness = 0,
    status = "Fail",
    B_V = Failed_B_V_C_V$B_V,
    C_V = Failed_B_V_C_V$C_V
  )


### These are the B_V/C_V that would lead to mortality,
### This means that we know the end time (point of death) and the
### cumulative transmission potential
Duration_Initial_PC_MORT <- subset(
  Duration_Initial_PC_100 ,
  Duration_Initial_PC_100 $status != "success"
)


Duration_Initial_PC_SUCCESS <- subset(
  Duration_Initial_PC_100 ,
  Duration_Initial_PC_100 $status == "success"
)


### Only run the function for the burst size and transmission investment
### combination that leads to successful infection that does not induce
### host mortality

Fitness_MODEL_PC <- mcmapply(Simulator_MalariaPC_DDE_BC_Cut,
                             c(Duration_Initial_PC_SUCCESS$B_V),
                             c(Duration_Initial_PC_SUCCESS$C_V),
                             c(Duration_Initial_PC_SUCCESS$endtime),
                             mc.cores = 5,
                             SIMPLIFY = FALSE
)

### Now we know what the fitness is
Duration_Initial_PC_SUCCESS$end_fitness <-
  unlist(lapply(Fitness_MODEL_PC, Gametocyte_Fitness))

### These are the fitness model data.frame that should work
Fitness_MODEL_PC_FULL_100 <- rbind.data.frame(
  Duration_Initial_PC_SUCCESS,
  Duration_Initial_PC_FAIL,
  Duration_Initial_PC_MORT
)

### This is saved
write.csv(Fitness_MODEL_PC_FULL_100 , file = here(
  "Output",
  "Fitness_Model",
  "Fitness_MODEL_PC_FULL100.csv"
))




###1) MERGE the full_model and the fitness model,
###2) filter only for successful infections and 
###3) split them 
FULL_MODEL_100_F <- merge( FULL_MODEL_PC_DT, 
  Fitness_MODEL_PC_FULL_100, by =c("B_V","C_V")) %>%
  filter(., status == 'success') %>%
  split(., list(.$B_V,.$C_V),drop = TRUE)


###Different average lifespans of the equation: DAY 5 to DAY 100, by 5 days
average_lifespan <- c(seq(5, 100,1))

###This is for the equal weighting
equal_RC <- average_life_span_RC(FULL_MODEL_100_F, 'equal',
                                 average_lifespan )

###Days of interest
doi <-  c(5,20,25,50,75)
doi_points <- filter(equal_RC, equal_RC$lifespan %in% c(5,20,25,50,75))

equal_RC_nomin <- subset(equal_RC,equal_RC$group == 'fit' &
                           equal_RC$lifespan %in%  seq(0,50,5))



###

RC_GG <- 
  ggplot(data = subset(equal_RC,equal_RC$group != 'min'), 
         aes(x = lifespan, y= RC, color = group, group = group))+
  geom_line(size = 1) + 
  geom_point(data = subset(doi_points,
                           doi_points$group != 'min'),
             aes( x= lifespan, y = RC),
             size = 3 ) +
  scale_color_manual(values = c('fit' = 'black', 'max'='#50de66'))+
  scale_x_continuous(limits=c(0,100)) +
  xlab("Total infection length") + 
  ylab(expression(paste("Replicative capacity (","(1-c)", beta,")")))+ 
  theme_classic() + 
  theme(
    axis.text = element_text(size = 14,color = 'black'),
    axis.title = element_text(size = 14),
    strip.background =  element_blank(),
    strip.text = element_text(size = 16),
    strip.text.y.right = element_text(angle = 0),
    panel.spacing = unit(1.5, "lines"), 
    legend.position = 'none')





###ZOOMED IN


original2 <- 
  ggplot(equal_RC_nomin, 
         aes(x = B_V, y = C_V, 
             color = as.factor(lifespan))) +
  geom_text(aes(x= B_V+1.5, y= C_V, label = lifespan))+
  geom_point(size = 3.0)+
  geom_path(size=1,lineend = "round", group =1) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL_100,
      Fitness_MODEL_PC_FULL_100$status ==
        "Fail"
    ),
    aes(x = B_V, y = C_V), fill = "#d1dbe8",
    inherit.aes = FALSE
    
  ) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL_100,
      Fitness_MODEL_PC_FULL_100$status ==
        "mort"
    ),
    aes(x = B_V, y = C_V), fill = "black",
    inherit.aes = FALSE
    
  )+
  scale_color_viridis(option='rocket',discrete = TRUE,)+
  new_scale_color() +
  
  geom_segment(
    data = horizontal_vert_df,
    aes(
      x = x, xend = xend,
      y = y, yend = yend,
      color = id), size = 0.5,
    
    lineend = "round",
    inherit.aes = FALSE
  )+
  scale_color_manual(
    values = c(
      "fail" = "black",
      "mort" = "#b8fcd5"
    ),
    guide = "none"
  )+
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Transmission investment (c)") +
  
  scale_x_continuous(limits=c(45,51),expand=c(0,0))+ 
  scale_y_continuous(limits=c(.75, 1),expand=c(0,0))+
  theme(panel.background=element_rect(fill ='#e3e0da'),
        panel.border = element_rect(color = 'black', fill = NA, size =1),
        panel.grid  = element_blank(),
        
        text = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        legend.position = "none"
  ) +
  ggtitle("Initial RM Must be greater or equal to 1")







###Produces the mortality and non-establishing infection lines
horizontal_vert_df <- grapher_mortality_boundary(Fitness_MODEL_PC_FULL_100)

original <- 
  ggplot(equal_RC_nomin, 
         aes(x = B_V, y = C_V, 
             color = as.factor(lifespan))) +
  geom_text(aes(x= B_V+1.5, y= C_V, label = lifespan))+
  geom_point(size = 3.0)+
  geom_path(size=1,lineend = "round", group =1) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL_100,
      Fitness_MODEL_PC_FULL_100$status ==
        "Fail"
    ),
    aes(x = B_V, y = C_V), fill = "#d1dbe8",
    inherit.aes = FALSE
    
  ) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL_100,
      Fitness_MODEL_PC_FULL_100$status ==
        "mort"
    ),
    aes(x = B_V, y = C_V), fill = "black",
    inherit.aes = FALSE
    
  )+
  scale_color_viridis(option='rocket',discrete = TRUE,)+
  new_scale_color() +
  
  geom_segment(
    data = horizontal_vert_df,
    aes(
      x = x, xend = xend,
      y = y, yend = yend,
      color = id), size = 0.5,
    
    lineend = "round",
    inherit.aes = FALSE
  )+
  scale_color_manual(
    values = c(
      "fail" = "black",
      "mort" = "#b8fcd5"
    ),
    guide = "none"
  )+
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Transmission investment (c)") +
  
  scale_x_continuous(limits=c(1,51),expand=c(0,0))+ 
  scale_y_continuous(limits=c(.01, 1),expand=c(0,0))+
  theme(panel.background=element_rect(fill ='#e3e0da'),
        panel.border = element_rect(color = 'black', fill = NA, size =1),
        panel.grid  = element_blank(),
        
        text = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        legend.position = "none"
  ) +
  ggtitle("Initial RM Must be greater or equal to 1")



###



RM1.5 <- ggplot(equal_RC_nomin, 
      aes(x = B_V, y = C_V, 
          color = as.factor(lifespan))) +
  geom_text(aes(x= B_V+1.5, y= C_V, label = lifespan))+
  geom_point(size = 3.0)+
  geom_path(size=1,lineend = "round", group =1) +
  geom_raster(
    data = BVCV_1.5_establish,
    aes(x = B_V, y = C_V), fill = "#d1dbe8",
    inherit.aes = FALSE
    
  ) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL_100,
      Fitness_MODEL_PC_FULL_100$status ==
        "mort"
    ),
    aes(x = B_V, y = C_V), fill = "black",
    inherit.aes = FALSE
    
  ) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL_100,
      Fitness_MODEL_PC_FULL_100$status ==
        "mort"
    ),
    aes(x = B_V, y = C_V), fill = "black",
    inherit.aes = FALSE
    
  )+
  scale_color_viridis(option='rocket',discrete = TRUE,)+
  new_scale_color() +
  
  geom_segment(
    data = horizontal_vert_df,
    aes(
      x = x, xend = xend,
      y = y, yend = yend,
      color = id), size = 0.5,
    
    lineend = "round",
    inherit.aes = FALSE
  )+
  scale_color_manual(
    values = c(
      "fail" = "black",
      "mort" = "#b8fcd5"
    ),
    guide = "none"
  )+
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Transmission investment (c)") +
  scale_x_continuous(limits=c(1,51),expand=c(0,0))+ 
  scale_y_continuous(limits=c(.01, 1),expand=c(0,0))+
  theme(panel.background=element_rect(fill ='#e3e0da'),
        panel.border = element_rect(color = 'black', fill = NA, size =1),
        panel.grid  = element_blank(),
        text = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        legend.position = "none") +
  ggtitle("Initial RM Must be greater or equal to 1.5")




RM2 <- ggplot(equal_RC_nomin, 
                aes(x = B_V, y = C_V, 
                    color = as.factor(lifespan))) +
  geom_text(aes(x= B_V+1.5, y= C_V, label = lifespan))+
  geom_point(size = 3.0)+
  geom_path(size=1,lineend = "round", group =1) +
  geom_raster(
    data = BVCV_2.0_establish,
    aes(x = B_V, y = C_V), fill = "#d1dbe8",
    inherit.aes = FALSE
    
  ) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL_100,
      Fitness_MODEL_PC_FULL_100$status ==
        "mort"
    ),
    aes(x = B_V, y = C_V), fill = "black",
    inherit.aes = FALSE
    
  ) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL_100,
      Fitness_MODEL_PC_FULL_100$status ==
        "mort"
    ),
    aes(x = B_V, y = C_V), fill = "black",
    inherit.aes = FALSE
    
  )+
  scale_color_viridis(option='rocket',discrete = TRUE,)+
  new_scale_color() +
  
  geom_segment(
    data = horizontal_vert_df,
    aes(
      x = x, xend = xend,
      y = y, yend = yend,
      color = id), size = 0.5,
    
    lineend = "round",
    inherit.aes = FALSE
  )+
  scale_color_manual(
    values = c(
      "fail" = "black",
      "mort" = "#b8fcd5"
    ),
    guide = "none"
  )+
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Transmission investment (c)") +
  scale_x_continuous(limits=c(1,51),expand=c(0,0))+ 
  scale_y_continuous(limits=c(.01, 1),expand=c(0,0))+
  theme(panel.background=element_rect(fill ='#e3e0da'),
        panel.border = element_rect(color = 'black', fill = NA, size =1),
        panel.grid  = element_blank(),
        text = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        legend.position = "none") +
  ggtitle("Initial RM Must be greater or equal to 2")


RM4 <- ggplot(equal_RC_nomin, 
              aes(x = B_V, y = C_V, 
                  color = as.factor(lifespan))) +
  geom_text(aes(x= B_V+1.5, y= C_V, label = lifespan))+
  geom_point(size = 3.0)+
  geom_path(size=1,lineend = "round", group =1) +
  geom_raster(
    data = BVCV_4.0_establish,
    aes(x = B_V, y = C_V), fill = "#d1dbe8",
    inherit.aes = FALSE
    
  ) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL_100,
      Fitness_MODEL_PC_FULL_100$status ==
        "mort"
    ),
    aes(x = B_V, y = C_V), fill = "black",
    inherit.aes = FALSE
    
  ) +
  geom_raster(
    data = subset(
      Fitness_MODEL_PC_FULL_100,
      Fitness_MODEL_PC_FULL_100$status ==
        "mort"
    ),
    aes(x = B_V, y = C_V), fill = "black",
    inherit.aes = FALSE
    
  )+
  scale_color_viridis(option='rocket',discrete = TRUE,)+
  new_scale_color() +
  
  geom_segment(
    data = horizontal_vert_df,
    aes(
      x = x, xend = xend,
      y = y, yend = yend,
      color = id), size = 0.5,
    
    lineend = "round",
    inherit.aes = FALSE
  )+
  scale_color_manual(
    values = c(
      "fail" = "black",
      "mort" = "#b8fcd5"
    ),
    guide = "none"
  )+
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Transmission investment (c)") +
  scale_x_continuous(limits=c(1,51),expand=c(0,0))+ 
  scale_y_continuous(limits=c(.01, 1),expand=c(0,0))+
  theme(panel.background=element_rect(fill ='#e3e0da'),
        panel.border = element_rect(color = 'black', fill = NA, size =1),
        panel.grid  = element_blank(),
        text = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        legend.position = "none") +
  ggtitle("Initial RM Must be greater or equal to 4")


original + RM1.5 + RM2 + RM4
