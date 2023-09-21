Replicative_Capacity_Finder <- function(Full, Fitness){

###1) MERGE the full_model and the fitness model,
###2) filter only for successful infections and 
###3) split them 
FULL_MODEL_100_F <- merge(Full, 
                          Fitness,
                          by =c("B_V","C_V")) %>%
                          filter(., status == 'success') %>%
                          split(., list(.$B_V,.$C_V),drop = TRUE)


###Different average lifespans of the equation: DAY 5 to DAY 100, by 5 days
average_lifespan <- c(seq(5, 100,1))

###This is for the equal weighting
equal_RC <- average_life_span_RC(FULL_MODEL_100_F, 
                                 average_lifespan)

###Days of interest
doi <-  c(5,20,25,50,75)
doi_points <- filter(equal_RC, equal_RC$lifespan %in% c(5,20,25,50,75))

equal_RC_nomin <- subset(equal_RC,equal_RC$group == 'fit' &
                           equal_RC$lifespan %in% c(5,20,25,30,40,50,70,100))

return(list(equal_RC, doi_points, equal_RC_nomin))

}

