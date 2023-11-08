
x = Fitness_MODEL_PC_FULL_MED


Day100_Simulator_Standard <- function(x,initialvalue, c_value_interest, name){


###Only get the ones that do not kill the host
SUCCESS_BV_CV <- subset(x,
                       x$status == 'success' & 
                       x$C_V %in%  c_value_interest)

FULL_MODEL_PC_100_Success<- mcmapply(Simulator_MalariaPC_DDE_BC_Cut  ,
                             c(SUCCESS_BV_CV$B_V),
                             c_value_interest,
                             c(initialvalue),
                             100,
                             mc.cores = 2,
                             SIMPLIFY = FALSE)

 
FULL_MODEL_PC_100_DT_Success <- do.call(rbind, FULL_MODEL_PC_100_Success)
FULL_MODEL_PC_100_DT_Success$status <- 'succes' 



###MORTALITY 
MORT_BV_CV <- subset(x,x$C_V %in%  c_value_interest& 
                                    x$status == "mort" )

FULL_MODEL_PC_100_Mort<- mcmapply(Simulator_MalariaPC_DDE_BC_Cut  ,
                                  c(MORT_BV_CV$B_V),
                                  c(c_value_interest),
                                  c(initialvalue),
                                  100,
                                  mc.cores = 2,
                                  SIMPLIFY = FALSE)
FULL_MODEL_PC_100_Mort_DT <- do.call(rbind, FULL_MODEL_PC_100_Mort)

FULL_MODEL_PC_100_Mort_FITNESS <- merge(FULL_MODEL_PC_100_Mort_DT , Fitness_MODEL_PC_FULL_MED)



FULL_MODEL_PC_Mort_100T_FITNESS_split <- split(FULL_MODEL_PC_100_Mort_FITNESS  ,
                                             list(FULL_MODEL_PC_100_Mort_FITNESS  $B_V,
                                                  FULL_MODEL_PC_100_Mort_FITNESS $C_V))

for (k in seq(1,length(FULL_MODEL_PC_Mort_100T_FITNESS_split ))){
  tmp = FULL_MODEL_PC_Mort_100T_FITNESS_split[[k]]
  tmp[tmp$time > unique(tmp$endtime),][,c("R","G")] <- 0
  FULL_MODEL_PC_Mort_100T_FITNESS_split [[k]] <- tmp
}


FULL_MORT <- do.call(rbind,FULL_MODEL_PC_Mort_100T_FITNESS_split)


FULL <- rbind(FULL_MORT[,c("time","R","G","B_V","C_V","status")],  
              FULL_MODEL_PC_100_DT_Success )


 write.csv(FULL, file = here(
  "Output", "Full_Model",
  paste("FULL_MODEL_100_PC_DT",name,c_value_interest,".csv",sep="")
 ))
}



Fitness_MODEL_PC_FULL_LOW_SUPP <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_low.csv"
))

Fitness_MODEL_PC_FULL_MED <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_med.csv"
))

Fitness_MODEL_PC_HIGH_SUPP <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_high.csv"
))

Day100_Simulator_Standard(Fitness_MODEL_PC_FULL_LOW_SUPP, 43.58965,"Low")
Day100_Simulator_Standard(Fitness_MODEL_PC_FULL_MED,4358.965,0.75,"Med")
Day100_Simulator_Standard(Fitness_MODEL_PC_HIGH_SUPP,43589.65, "High")



