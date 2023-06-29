###Test case to make sure things work

library(here)

# This is the main simulation assuming nothing has changed
# about the initial red blood cell density or the replenishment rate.
# This is one of the longer code to source

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))

### Main modeling code
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_CUT.cpp"))
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_UNCUT.cpp"))
source(here("Code", "Simulator_Code", "Simulator_Main.R"))
source(here("Code", "Helper_Function_Code", "Fitness_Functions.R"))

###Definetly kills the host (burst size of 50 and transmission investment of 25%)
killing_strain_test <- Simulator_Malaria_BC(50, 0.25)
endtime_killing_strain <- Finder_RM(killing_strain_test,48)   

###The blue line and red line should intersect
plot(killing_strain_test[1:100,c('time',"R")],type = 'l')
abline(h = 6.5 * 10^5,col = 'red')
abline(v = endtime_killing_strain$endtime,col= 'blue')

###Eyeball to see if the end-fitness time and the cumulative fitness line 
###are the same
plot(killing_strain_test$time[1:100], cumsum(PrI_PC(killing_strain_test$G)[1:100] * 1/10))
abline(v = endtime_killing_strain$endtime,col= 'blue')

##############################
###Non establishing infection#
##############################
non_establishing_test <- Simulator_Malaria_BC_TESTING(5,0.5)
non_establishing_test_2 <- Simulator_Malaria_BC(5,0.5)

non_establishing_test_I <- rowSums(non_establishing_test[,3:102]) ###Total I

#Monotonically decreasing
all(diff(non_establishing_test_I) <0)

non_establishing_test <- Finder_RM(non_establishing_test_2, 48)   


##############################
###Successful infection.     #
##############################
establishing_test <- Simulator_Malaria_BC_TESTING(8,0.54)
establishing_test_2 <- Simulator_Malaria_BC(8,0.54)

establishing_test_I <- rowSums(establishing_test[,3:102]) ###Total I
endtime_establishing_test<- Finder_RM(establishing_test_2,48)   

plot(establishing_test_2[,'time'],establishing_test_I)
abline(v = 42.2,col = 'blue') 

###The blue line and red line should intersect
plot(establishing_test_2[1:500,c('time',"R")],type = 'l')
abline(h = 6.5 * 10^5,col = 'red') #Should not appear

Establishing_Cut <- Simulator_MalariaPC_DDE_BC_Cut_TESTING(8,0.54, 42.2)
##############################################################################
###Should see a drastic drop at where the acute infection ends and especially#
###two days after                                                            #
#############################################################################
plot(establishing_test_2$time, establishing_test_2$G,type='l')
lines(Establishing_Cut$time, Establishing_Cut$G, col = 'red')
abline(v = 42.2)
abline(v = 44.2)

