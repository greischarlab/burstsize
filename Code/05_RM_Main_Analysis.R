library(here)

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "02_Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "04_RM_Function.R"))

###


Fitness_MODEL_PC_FULL_MED <- read.csv(here(
  "Output", "Fitness_Model",
  "FITNESS_MODEL_PC_4358.965.csv"
))

MODEL_PC_FULL_MED <- read.csv(here(
  "Output", "Full_Model",
  "FULL_MODEL_PC_DT_ 4358.965 .csv"
))


###Because the original simulation does not contain the 
###

###What are the burst size of interests

B_Interest <- c(14.5, 15,15.5)
C_Interest <- c(0.75)

rate_PMR_dataframe1 <- rate_PMR_data(
  MODEL_PC_FULL_MED ,
  Fitness_MODEL_PC_FULL_MED ,
  c(B_Interest ),
  c(C_Interest),
  48
)



rate_PMR_FULL_TS <- rbind(rate_PMR_dataframe1[[1]])
rate_PMR_FULL_points <- rbind(rate_PMR_dataframe1[[2]])

RM_Plotter(rate_PMR_FULL_TS ,rate_PMR_FULL_points )

ggsave(file = here("Figures", "Raw", "RM_3_Plot.pdf"), width = 8, height = 9, units = "in")
