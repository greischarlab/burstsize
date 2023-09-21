
Day100_Simulator_Standard <- function(x,initialvalue, name){


###Only get the ones that do not kill the host
SUCCESS_BV_CV <- subset(x,
                       x$status != 'success')


FULL_MODEL_PC_100_MORT<- mcmapply(Simulator_Malaria_BC ,
                             c(SUCCESS_BV_CV$B_V),
                             c(SUCCESS_BV_CV$C_V),
                             c(initialvalue),
                             100,
                             mc.cores = 4,
                             SIMPLIFY = FALSE)

 
 FULL_MODEL_PC_100_DT <- do.call(rbind, FULL_MODEL_PC_100)


 write.csv(FULL_MODEL_PC_100_DT, file = here(
  "Output", "Full_Model",
  paste("FULL_MODEL_100_PC_DT",name,".csv",sep="")
 ))
}



Fitness_MODEL_PC_FULL_LOW_SUPP <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_43.58965.csv"
))

Fitness_MODEL_PC_FULL_MED <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_4358.965.csv"
))

Fitness_MODEL_PC_HIGH_SUPP <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_43589.65.csv"
))

Day100_Simulator_Standard(Fitness_MODEL_PC_FULL_LOW_SUPP, 43.58965,"Low")
Day100_Simulator_Standard(Fitness_MODEL_PC_FULL_MED,4358.965,"Med")
Day100_Simulator_Standard(Fitness_MODEL_PC_HIGH_SUPP,43589.65, "High")



