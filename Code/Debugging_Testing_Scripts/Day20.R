library(here)

source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
### Main modeling code
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_CUT.cpp"))
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_UNCUT.cpp"))
source(here("Code", "Simulator_Code", "Simulator_Main_2.R"))
source(here("Code", "Helper_Function_Code", "02_Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "02_Simulator_Code.R"))

LOW_Fitness <- read.csv(here("Output", "Fitness_Model",
                             "FITNESS_MODEL_PC_low.csv"))

C_V_DF <- subset(LOW_Fitness, 
       LOW_Fitness$B_V==10)

C_V_NoDeath <- subset(C_V_DF, C_V_DF$status != 'mort')


###
p_val <- 4.0e-6 #maximum infection rate 
mu_M <- 48 #the mu_mortality 
R_val <- 8500000 #initial red blood cell

RM_limit_1 <- 1 * ((p_val * R_val) + mu_M) / (p_val * R_val)

C_V_Estab  <- ifelse((1 - C_V_NoDeath$C_V ) * B_V >= RM_limit_1,
                            "Establish", "Fail")
C_V_NoDeath$status <- C_V_Estab 

FULL_CV <- subset(C_V_NoDeath, C_V_NoDeath$status!= "Fail")

FULL_MODEL_PC_20 <- mcmapply(Simulator_Malaria_BC,
                          10,
                          c(FULL_CV$C_V),
                          43.85965,
                          mc.cores = 5,
                          SIMPLIFY = FALSE
)

MODEL_20<- Surface_Plot_Standard_Duration(
  FULL_MODEL_PC_20,20)

ggplot(MODEL_20, aes( x= C_V, y= end_fitness)) + geom_line() + geom_vline(xintercept = 0.42)+
  theme_classic()+
  xlab("Transmission investment (c)") + 
  ylab("Cumulative transmission potential (f)")+
  theme(axis.text = element_text(size = 14,color = 'black'),
        axis.title = element_text(size = 15, color = 'black'),
        title = element_text(size = 16))+
  ggtitle("Day 20")


ggsave(here("Figures", "Raw", "Supp_Day20.pdf"),
       width = 7.5, height =7, unit = "in")
        