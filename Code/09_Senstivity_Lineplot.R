########################################################
### Understanding how the viable burst size combinations#
### influence virulence, transmission, and such.     ####
########################################################

library(here)

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))

###We need the data here 

FITNESS_MODEL <- as.data.frame(fread(here("Output","Fitness_Model",
                                      "Fitness_MODEL_PC_FULL.csv")))
FITNESS_MODEL$change = "Original"
FITNESS_MODEL$grouping = "Original"

FITNESS_L <- as.data.frame(fread(here("Output","Fitness_Model",
                                      "Fitness_L.csv")))
FITNESS_L$grouping = "lambda"

FITNESS_RI <- as.data.frame(fread(here("Output","Fitness_Model",
                                       "Fitness_RI.csv")))   
FITNESS_RI$grouping = "RI"

FITNESS_UM <- as.data.frame(fread(here("Output","Fitness_Model",
                                       "Fitness_UM.csv")))       
FITNESS_UM$grouping = "UM"


FITNESS_ALL <- rbind(FITNESS_MODEL, FITNESS_L, FITNESS_RI, FITNESS_UM)
FITNESS_ALL_SPLIT <- split(FITNESS_ALL, FITNESS_ALL$change)
FITNESS_MAX <- do.call(rbind,by(FITNESS_ALL, FITNESS_ALL$change, function(x) x[which.max(x$end_fitness),],
                                simplify = FALSE))

FITNESS_ALL_ADDENDUM = NULL
for (k in seq(1,nrow(FITNESS_MAX))){
  tmp = FITNESS_MAX[k,]
  tmp_fitness = FITNESS_ALL_SPLIT[[k]]
  
  updown_CV <- subset(tmp_fitness,   tmp_fitness$C_V == tmp$C_V & 
                        tmp_fitness$status == 'success')
  leftright_BV <- subset(tmp_fitness,   tmp_fitness$B_V == tmp$B_V & 
                          tmp_fitness$status == 'success')
  
  min_CV = min(leftright_BV $C_V)
  max_CV = max(leftright_BV $C_V)
  
  min_BV = min(updown_CV$B_V)
  max_BV = max(updown_CV$B_V)
  
  
  FITNESS_ALL_ADDENDUM[[k]] <- cbind.data.frame(
    B_V = tmp$B_V,
    C_V= tmp$C_V,
    min_CV, max_CV, 
   min_BV,max_BV, change = tmp$change,
   grouping = tmp$grouping)
}




FITNESS_ALL_ADDENDUM <- do.call(rbind, FITNESS_ALL_ADDENDUM)

###SKETCHY WAY### - CHECK (WILL AUTOMATE IT LATEr)
FITNESS_ALL_ADDENDUM$minusplus <- c('p','m','m','p','p','m','NA')
FITNESS_ALL_MAX$minusplus <- c('p','m','m','p','p','m','NA')
FITNESS_ALL$minusplus
### Day Labeler
Parameter_Names<- c(
  'lambda' = "A. RBC replenishment",
  "RI" = "B. Initial RBC",
  'UM' = "C. Merozoite mortality rate"
)





ggplot(FITNESS_MAX[-7,],
       aes(x = B_V, y= C_V)) + 
  geom_point(aes(color = minusplus, 
                                                                size = end_fitness)
                                                           )+ 
  geom_segment(data =  FITNESS_ALL_ADDENDUM[-7,], aes(x= B_V, xend = B_V, 
                                                 y = min_CV, yend = max_CV, color = minusplus) )+
  geom_segment(data =  FITNESS_ALL_ADDENDUM[-7,], aes(x= min_BV, xend = max_BV, 
                                                 y = C_V, yend = C_V, color = minusplus) )+
  facet_wrap(~grouping,ncol = 1, labeller = as_labeller(Parameter_Names)) + 
  scale_size_continuous(range = c(1, 10))+
  scale_color_manual(values = c("#3DD8EE", '#FD7E6C')) + 
  xlab("Burst size") + 
  ylab("Transmission investment")+
  annotate("point", x = 8, y =0.54,size = 5, color = 'black') +
  annotate("segment", x = 8, xend = 8, y= 0.01, yend = 0.59 ) + 
  annotate("segment", x = 7.5, xend = 19.5, y= 0.54, yend = 0.54 ) + 
  theme_bw() + 
  theme(legend.position = 'none',
        strip.background = element_blank(), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 14)) 



FITNESS_54 <- subset(FITNESS_ALL, FITNESS_ALL$C_V == 0.54 &
                       FITNESS_ALL$change != 'Original')

ggplot(FITNESS_54, 
       aes( x = B_V, 
            y = end_fitness,
            group = change,
            color = change)) + 
  geom_line(size = 1.2) + 
  facet_wrap(~grouping) 
  



                                   
                                   
                                   