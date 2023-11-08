library(here)

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code","Helper_Function_Code","Grapher_vert_hor.R"))
source(here("Code", "Helper_Function_Code", "02_Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "RM_Function.R"))
source(here("Code", "Helper_Function_Code", "Plotters","05_RM_TriplotPlotter.R"))

###



Fitness_MODEL_PC_FULL_MED <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_med.csv"
))

MODEL_PC_FULL_MED <- read.csv(here(
  "Output", "Full_Model",
  "FULL_MODEL_PC_DT_med.csv"
))


###Because the original simulation does not contain the 
###

###What are the burst size of interests

B_Interest <- c(14.5, 18.0)
C_Interest <- c(0.75)

###Look 


rate_PMR_dataframe1 <- rate_PMR_data(
  MODEL_PC_FULL_MED ,
  Fitness_MODEL_PC_FULL_MED ,
  c(B_Interest),
  c(C_Interest),
  48
)

B_V_Interest_one <- 11

FULL_One_MED<- Simulator_Malaria_BC(B_V_Interest_one,0.75,4385.96491)

Duration_Initial_One_MED <-  Finder_RM(FULL_One_MED, 48)
Duration_Initial_One_MED$B_V = B_V_Interest_one
Duration_Initial_One_MED$C_V = 0.75

Fitness_MODEL_One_MED<- mcmapply(Simulator_MalariaPC_DDE_BC_Cut,
                                   c(B_V_Interest_one),
                                   0.75,
                                   4385.96491,
                                   Duration_Initial_One_MED$endtime,
                                   mc.cores = 5,
                                   SIMPLIFY = FALSE)


Duration_Initial_One_MED$end_fitness <-
  unlist(lapply(Fitness_MODEL_One_MED, Gametocyte_Fitness))





rate_PMR_dataframe2 <- rate_PMR_data(
  FULL_One_MED,
  Duration_Initial_One_MED ,
  c(B_V_Interest_one) ,
  0.75,
  48
)




rate_PMR_FULL_1 <- rbind(rate_PMR_dataframe1 [[1]],
                         rate_PMR_dataframe2 [[1]])

rate_PMR_FULL_2 <- rbind(rate_PMR_dataframe1 [[2]],
                         rate_PMR_dataframe2 [[2]])



RM_Plotter(rate_PMR_FULL_1  ,rate_PMR_FULL_2)

ggsave(file = here("Figures", "Raw", "RM_3_Plot_2_Axis.pdf"), width =9, height = 11, units = "in")

####Inclusion...


