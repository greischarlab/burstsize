MED_LONG_OPTIMUM <- subset(FULL_MODEL_100_MED, 
                           FULL_MODEL_100_MED$B_V == 45.5 &
                            FULL_MODEL_100_MED$C_V == 0.89 &
                             FULL_MODEL_100_MED$time < 101)

MED_SHORT_OPTIMUM <- subset(FULL_MODEL_100_MED, 
                           FULL_MODEL_100_MED$B_V == 44 &
                             round(FULL_MODEL_100_MED$C_V,4) == 0.89 &
                           FULL_MODEL_100_MED$time < 101)


MED_SHORT_89<- subset(FULL_MODEL_100_MED, 
                              round(FULL_MODEL_100_MED$C_V,4) == 0.89&
                              FULL_MODEL_100_MED$time < 101)

MED_SHORT_89_Split <- split(MED_SHORT_89,MED_SHORT_89$B_V)
###HMM I think it's easier to figure out by fixing transmission investment


MED_ENDFITNESS_NULL <- NULL
lifespan_interest <- c(5,20,30,50,100)
for (k in seq(1, length(lifespan_interest ))){
  
tmp<-  do.call(rbind,lapply(MED_SHORT_89_Split,
               Gametocyte_Fitness_EqualWeighting,
                life_span =lifespan_interest[[k]]))
MED_ENDFITNESS_NULL[[k]] = tmp
}

MED_ENDFITNESS_FULL <- do.call(rbind, 
                               MED_ENDFITNESS_NULL )

MED_OPTIMUM_89 <- do.call(rbind,by(MED_ENDFITNESS_FULL, 
                     list(MED_ENDFITNESS_FULL$life_span),
                     function(x)
                     x[which.max(x$end_fitness),]))

ggplot(
  MED_OPTIMUM_89 , aes( x= life_span, y= B_V)) + 
  geom_line() + theme_bw() + 
  xlab("Infection length") + 
  ylab("Burst size") 

ggplot(MED_ENDFITNESS_FULL, 
       aes(x = B_V, y= end_fitness)) + geom_line() + facet_wrap(~life_span, scales = 'free_y')+
  xlab("Burst size") + ylab("End fitness at infection length (facet)")

############################
############################


MED_SHORT_75_LOWHIGH<- subset(FULL_MODEL_100_MED, 
                      round(FULL_MODEL_100_MED$C_V,4) == 0.89 &
                        round(FULL_MODEL_100_MED$B_V,4) %in% c(33,40,50) &
                        FULL_MODEL_100_MED$time < 101)

ggplot(MED_SHORT_75_LOWHIGH, aes(x=time, y= Pripc, group = as.factor(B_V))) + 
         geom_line(aes(color = as.factor(B_V)))+ theme_classic() + geom_vline(xintercept = c(30,100))


###
MED_SHORT_75_LOWHIGH_Split <- split(MED_SHORT_75_LOWHIGH, MED_SHORT_75_LOWHIGH$B_V)
for (k in seq(1,3)){
  MED_SHORT_75_LOWHIGH_Split [[k]]$PriPC <- PrI_PC(MED_SHORT_75_LOWHIGH_Split [[k]]$G)
  
  MED_SHORT_75_LOWHIGH_Split [[k]]$cumtrans <- cumsum(  MED_SHORT_75_LOWHIGH_Split [[k]]$PriPC )
  MED_SHORT_75_LOWHIGH_Split [[k]]$Propcumtrans <- MED_SHORT_75_LOWHIGH_Split [[k]]$cumtrans/max(MED_SHORT_75_LOWHIGH_Split [[k]]$cumtrans)
  
  full_30 <-  sum(MED_SHORT_75_LOWHIGH_Split[[k]][MED_SHORT_75_LOWHIGH_Split[[k]]$time <31,]$PriPC)
  full_5 <-  sum(MED_SHORT_75_LOWHIGH_Split[[k]][MED_SHORT_75_LOWHIGH_Split[[k]]$time <6,]$PriPC)
  
  MED_SHORT_75_LOWHIGH_Split [[k]]$Prop30cumtrans <- MED_SHORT_75_LOWHIGH_Split [[k]]$cumtrans/full_30
  MED_SHORT_75_LOWHIGH_Split [[k]]$Prop5cumtrans <- MED_SHORT_75_LOWHIGH_Split [[k]]$cumtrans/full_5
  
  
   }

plot(
     MED_SHORT_75_LOWHIGH_Split[[1]]$time,
     MED_SHORT_75_LOWHIGH_Split[[1]]$PriPC,type = 'l', xlab = 'time', ylab= 'trans. prob',col = 'blue',
     ylim=c(0,1))
lines(
  MED_SHORT_75_LOWHIGH_Split[[2]]$time,
  MED_SHORT_75_LOWHIGH_Split[[2]]$PriPC,type = 'l', col = 'black')
lines(
  MED_SHORT_75_LOWHIGH_Split[[3]]$time,
  MED_SHORT_75_LOWHIGH_Split[[3]]$PriPC,type = 'l', col = 'red')



plot(
  MED_SHORT_75_LOWHIGH_Split[[1]]$time,
  MED_SHORT_75_LOWHIGH_Split[[1]]$Propcumtrans,type = 'l', xlab = 'DPI',
  ylab = 'Cumulative transmission')

lines(
    MED_SHORT_75_LOWHIGH_Split[[2]]$time,
    MED_SHORT_75_LOWHIGH_Split[[2]]$Propcumtrans,type = 'l', col = 'red')
  
)


lines(
  MED_SHORT_75_LOWHIGH_Split[[1]]$time,
  log10(MED_SHORT_75_LOWHIGH_Split[[1]]$G+1),type = 'l', xlab = 'DPI',
  ylab = 'Cumulative transmission', col = 'blue')

lines(
  MED_SHORT_75_LOWHIGH_Split[[1]]$time,
  log10(MED_SHORT_75_LOWHIGH_Split[[2]]$G+1),type = 'l', xlab = 'DPI',
  ylab = 'Cumulative transmission')

plot(
  MED_SHORT_75_LOWHIGH_Split[[3]]$time,
  log10(MED_SHORT_75_LOWHIGH_Split[[3]]$G+1),type = 'l', col = 'red',
  xlab = 'DPI', ylab = 'Gams (Log10)')

)





plot(
  MED_SHORT_75_LOWHIGH_Split[[1]]$time,
  MED_SHORT_75_LOWHIGH_Split[[1]]$Prop5cumtrans,type = 'l', xlab = 'time',
  ylab = 'Proportion of Cumulative Fitness of 5 day infection', xlim = c(0,7),ylim = c(0,1))

lines(
  MED_SHORT_75_LOWHIGH_Split[[2]]$time,
  MED_SHORT_75_LOWHIGH_Split[[2]]$Prop5cumtrans,type = 'l', col = 'red')

)

plot( MED_SHORT_75_LOWHIGH_Split[[3]]$time,  log10(MED_SHORT_75_LOWHIGH_Split[[3]]$G+1),type='l',
      xlab = 'time', ylab = 'R',col = 'orange')
lines( MED_SHORT_75_LOWHIGH_Split[[2]]$time,  log10(MED_SHORT_75_LOWHIGH_Split[[2]]$G+1),type='l',
      xlab = 'time', ylab = 'R',col = 'red')
lines( MED_SHORT_75_LOWHIGH_Split[[1]]$time,  log10(MED_SHORT_75_LOWHIGH_Split[[1]]$G+1),type='l',
      xlab = 'time', ylab = 'R')

###


FIT75 <- subset(Fitness_MODEL_PC_FULL_MED,
                round(Fitness_MODEL_PC_FULL_MED$C_V,4) == 0.89 &
                Fitness_MODEL_PC_FULL_MED$status == 'success')

FIT75$proportionofacutephase30 <- FIT75$endtime/30
FIT75$proportionofacutephase100<- FIT75$endtime/100

ggplot(FIT75,aes(x = B_V, y= proportionofacutephase30))+geom_line() + xlab("Burst size") + 
  ylab("Proportion of the acute phase in a 30 day infection") 

ggplot(FIT75,aes(x = B_V, y= proportionofacutephase100))+geom_line() + xlab("Burst size") + 
  ylab("Proportion of the acute phase in a  100 day infection") 
