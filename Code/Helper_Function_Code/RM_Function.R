###R_M calculator 

###This is the function for calculating the 
###effective merozoite invading 
rate_PMR <- function(x, B_V, C_V, mu_M){
  p =   2.5e-6
  rate = (1 - C_V) * B_V*((x[,"R"] * p)/((p*x[,"R"]) + mu_M))
  return(cbind(time = x[,'time'],rate, B_V,C_V))
}


rate_PMR_data <- function(x, y, B_Vec, C_Value, mu_M){

  Data_Interest <- 
           subset(x, 
           round(x$B_V,4) %in% B_Vec &
           round(x$C_V,4) == C_Value)
  
  
  End_Interest <- 
    subset(y, 
           round(y$B_V,4) %in% B_Vec &
           round(y$C_V,4) == C_Value)
  
  
  Data_Interest_Split <- split(Data_Interest, Data_Interest$B_V)
  
  rate_PMR_DF <- NULL
 
  end_duration_DF <- NULL 
  ###Calculate the time changing PMR
  
  for (k in seq(1, length(Data_Interest_Split ))){
    
    
   PMR_B <- rate_PMR(Data_Interest_Split[[k]], B_Vec[[k]], C_Value, mu_M )
   which.time <- PMR_B[,'time'] <= End_Interest[k,]$endtime
   
   #Daily Transmission Probability
   Daily_Trans_Prob <- PrI_PC(Data_Interest_Split[[k]][,"G"])
   #Cumulative transmission probability
   Cum_Trans_Potential <-  cumsum(Daily_Trans_Prob * 1/10)
  
   ###To figure out when the point approximate
   PMR_B_Func<- approxfun(PMR_B[,"time"], PMR_B[,"rate"])
   PMR_B_Func_root<- approxfun(PMR_B[,"time"], PMR_B[,"rate"] - 1)
   
   Daily_Trans_Prob_Func <-  approxfun(Data_Interest_Split[[k]][,"time"], 
                                      Daily_Trans_Prob)
   Cum_Trans_Potential_Func<- approxfun(Data_Interest_Split[[k]][,"time"], 
                                        Cum_Trans_Potential)
   
   PMR_root <- uniroot(PMR_B_Func_root, c(0,20))
   
   PMR_END = PMR_B_Func(End_Interest[k,]$endtime)
   Daily_Trans_End = Daily_Trans_Prob_Func(End_Interest[k,]$endtime)
   Cum_Trans_End = Cum_Trans_Potential_Func(End_Interest[k,]$endtime)

   
   
   
   ###Helps with figuring out if the values are in the acute phae or not
    Acute_Cum_Trans = ifelse(Cum_Trans_Potential< End_Interest[k,]$endtime,
                          1, 0)
   
    end_duration_DF[[k]] <- cbind.data.frame(
                     B_V = End_Interest[k,]$B_V,
                     C_V = End_Interest[k,]$C_V,
                     endtime = End_Interest[k,]$endtime,
                     PMR_END = PMR_END,
                     Daily_Trans_End = Daily_Trans_End,
                     Cum_Trans_End =  Cum_Trans_End,
                     PMR_root =  PMR_root$root,
                     Cum_Trans_Root =  PMR_root$root + 2,
                     Daily_Trans_Root =  PMR_root$root + 2 )
   
   
   
   rate_PMR_DF[[k]] <- cbind.data.frame(PMR_B,
                                   Daily_Trans_Prob, 
                                   Cum_Trans_Potential, which.time)
  }
  return(list(do.call(rbind, rate_PMR_DF),(do.call(rbind,end_duration_DF))
            ))
   
}



  
  

