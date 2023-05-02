


FULL_MODEL_PC_L_DT <- as.data.frame(fread(here(
  "Output", "Full_Model",
  "SA_FULL_MODEL_PC_lambda_DT.csv"
)))

FULL_MODEL_PC_L <- FULL_MODEL_PC_L_DT  [,-1]

FULL_MODEL_PC_L  <- split(FULL_MODEL_PC_L , list(FULL_MODEL_PC_L $B_V,
                                                   FULL_MODEL_PC_L $C_V), 
                        drop = TRUE)

remove(FULL_MODEL_PC_L_DT)
gc()
