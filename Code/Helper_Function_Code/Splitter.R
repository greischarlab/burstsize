
library(here)


functi

FULL_MODEL_PC_lambda_DT <- as.data.frame(fread(here(
  "Output", "Full_Model",
  "SA_FULL_MODEL_PC_lambda_DT.csv"
)))

FULL_MODEL_PC_lambda <- FULL_MODEL_PC_lambda_DT[,-1]

FULL_MODEL_PC_lambda  <- split(FULL_MODEL_PC_lambda  , 
                           list(FULL_MODEL_PC_lambda$B_V,
                                FULL_MODEL_PC_lambda$C_V,
                                FULL_MODEL_PC_lambda$change),
                           drop = TRUE)

remove(FULL_MODEL_PC_lambda_DT)
gc()

