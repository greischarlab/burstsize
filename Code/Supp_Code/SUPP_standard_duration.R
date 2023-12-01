source(here("Code", "Simulator_Code", "Simulator_Main_2.R"))
source(here("Code", "Helper_Function_Code", "FUNC_00_Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "FUNC_02_Simulator_Code.R"))
source(here("Code", "Helper_Function_Code", "FUNC_00_Grapher_vert_hor.R"))

###Only get the ones that do not kill the host
SUCCESS_BV_CV <- subset(Fitness_MODEL_PC_FULL_MED,
                        Fitness_MODEL_PC_FULL_MED$status == 'success')

###Run the Model to 100
FULL_MODEL_PC_100_Success <- mcmapply(Simulator_Malaria_BC ,
                                     c(SUCCESS_BV_CV$B_V),
                                     c(SUCCESS_BV_CV$C_V),
                                     c(4385.96491),
                                     mc.cores = 2,
                                     SIMPLIFY = FALSE)

for(k in seq(1, length(FULL_MODEL_PC_100_Success))){
  FULL_MODEL_PC_100_Success[[k]]$status = SUCCESS_BV_CV[k,]$status
}

FULL_MODEL_PC_100_DT_Success <- do.call(rbind, FULL_MODEL_PC_100_Success)


write.csv(FULL_MODEL_PC_100_DT_Success, 
    here(file = 'Output',"Full_Model","FULL_MODEL_PC_100_Success_SUPP.csv"))



FULL_MODEL_100_MED_SPLIT <- split(FULL_MODEL_PC_100_DT_Success, 
                            list(FULL_MODEL_PC_100_DT_Success$B_V, 
                                 FULL_MODEL_PC_100_DT_Success$C_V),
                                 drop = TRUE)


Split_Timer <- function(x,time_val){
  tmp <- lapply(x, function(x) 
    subset(as.data.frame(x),
           as.data.frame(x)$time <= time_val))
     return(tmp)
}

SPLIT_TIME_10 <- Split_Timer(FULL_MODEL_100_MED_SPLIT ,10)
SPLIT_TIME_15 <- Split_Timer(FULL_MODEL_100_MED_SPLIT ,15)
SPLIT_TIME_20 <- Split_Timer(FULL_MODEL_100_MED_SPLIT ,20)
SPLIT_TIME_25 <- Split_Timer(FULL_MODEL_100_MED_SPLIT ,25)
SPLIT_TIME_50 <- Split_Timer(FULL_MODEL_100_MED_SPLIT ,50)
SPLIT_TIME_100 <- Split_Timer(FULL_MODEL_100_MED_SPLIT ,100)


GAM_Adder <- function(x, day){
  tmp <- do.call(rbind,lapply(x, function(x) c(unique(x$B_V), unique(x$C_V))))
  tmp2 <- do.call(rbind,lapply(x, function(x) Gametocyte_Fitness(x)))
  
  a<- cbind.data.frame(tmp, tmp2, day)
  
colnames(a) <- c("B_V","C_V","end_fitness","ID")
return(a)
}

TIME_10<-GAM_Adder(SPLIT_TIME_10,10)
TIME_15<-GAM_Adder(SPLIT_TIME_15,15)
TIME_20<-GAM_Adder(SPLIT_TIME_20,20)
TIME_25<-GAM_Adder(SPLIT_TIME_25,25)
TIME_50<-GAM_Adder(SPLIT_TIME_50,50)
TIME_100<-GAM_Adder(SPLIT_TIME_100,100)

optimum_finder_function <- function(x,ID_val){
optimum_strategy <- x[which.max(x$end_fitness),]
return(optimum_strategy)
}

optimum_10 <- optimum_finder_function(TIME_10,10)
optimum_15 <- optimum_finder_function(TIME_15,15)
optimum_20 <- optimum_finder_function(TIME_20,20)
optimum_25 <- optimum_finder_function(TIME_25,25)
optimum_50 <-optimum_finder_function(TIME_50,50)
optimum_100 <- optimum_finder_function(TIME_100,100)



OPT_DF <-rbind(optimum_10,
               optimum_15,
               optimum_20,
      optimum_25 ,
      optimum_50 ,
      optimum_100)


optimum_finder_76_function <- function(x,ID_val){
  tmp <- subset(x, x$C_V == 0.76)
  optimum_strategy <-  tmp [which.max(tmp$end_fitness),]
return(optimum_strategy)
}

C_V_optimum_10 <- optimum_finder_76_function (TIME_10,10)
C_V_optimum_15 <- optimum_finder_76_function (TIME_15,15)
C_V_optimum_20 <- optimum_finder_76_function (TIME_20,20)
C_V_optimum_25 <- optimum_finder_76_function (TIME_25,25)
C_V_optimum_50 <-optimum_finder_76_function (TIME_50,50)
C_V_optimum_100 <- optimum_finder_76_function (TIME_100,100)


CV_OPT_DF <-rbind(C_V_optimum_10,
                  C_V_optimum_15,
                  C_V_optimum_20,
                  C_V_optimum_25 ,
                  C_V_optimum_50 ,
                  C_V_optimum_100)

TIME_FULL <- rbind.data.frame(TIME_10,TIME_15,TIME_20,TIME_25,TIME_50,TIME_100)


  ###Produces the mortality and non-establishing infection lines
  horizontal_vert_df <- grapher_mortality_boundary(Fitness_MODEL_PC_FULL_MED)

  
  ggplot(TIME_FULL)+
    geom_raster( aes(x = B_V, y = C_V, fill = end_fitness))+
    geom_raster(data = subset(Fitness_MODEL_PC_FULL_MED ,
                            Fitness_MODEL_PC_FULL_MED $status=='mort'),
              aes(x= B_V, y= C_V, fill = end_fitness))+
    geom_raster(data = subset(Fitness_MODEL_PC_FULL_MED ,
                            Fitness_MODEL_PC_FULL_MED $status=='Fail'),
              aes(x= B_V, y= C_V, fill = end_fitness), fill = "#d1dbe8")+
    geom_point(data= OPT_DF,  aes(x= B_V, y = C_V),
               color = 'pink')+
    geom_point(data= CV_OPT_DF,  aes(x= B_V, y = C_V),
               color = 'white', size = 2)+
    facet_wrap(~ID)+
  geom_segment(
    data = horizontal_vert_df,
    aes(
      x = x, xend = xend,
      y = y, yend = yend,
      color = id
    ), size = 0.9,
    lineend = "round"
  ) +
    scale_color_manual(
      values = c(
        "fail" = "black",
        "mort" = "#b8fcd5"
      ),
      guide = "none"
    ) +
    scale_fill_viridis(
      name = "Cumulative transmission \npotential",
      option = "magma"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(expression(paste("Burst size", "( ", beta, ")"))) +
    ylab("Transmission investment (c)") +
    theme(
      text = element_text(size = 14),
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      legend.position = "top",
      strip.background = element_blank(),
      strip.text = element_text(size = 14))+
    geom_point(data =  OPT_DF ,
               aes( x = B_V, y= C_V),
               color = '#FF116B', size = 2)+
    geom_segment(
      data = OPT_DF ,
      aes(
        x = B_V, xend = B_V, y = 0.01, yend = C_V)
      , col = "#FF116B",
      linetype = 2) +
    geom_segment(data =OPT_DF  ,
                 aes(
                   x =1, xend = B_V, y = C_V, yend = C_V), col = "#FF116B",
                 linetype = 2) +
    geom_segment(
      data = CV_OPT_DF  ,
      aes(
        x = B_V, xend = B_V, y = 0.01, yend = C_V)
      , col = "white",
      linetype = 2) +
    geom_segment(data =CV_OPT_DF  ,
                 aes(
                   x =1, xend = B_V, y = C_V, yend = C_V), col = "white",
                 linetype = 2) +
    annotate("text",
             x = 12, y = 0.93, label = "Unestablished \ninfection",
             size = 4
    ) +
    annotate("text",
             x = 45, y = 0.1, label = "Host \nMortality",
             size = 4, color = "#b8fcd5"
    )
    
ggsave(here("Figures","Raw", "SUPP_SD.pdf"), width = 11, height = 8, units= 'in')
  


