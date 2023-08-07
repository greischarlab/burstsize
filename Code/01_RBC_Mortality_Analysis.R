# HEADER --------------------------------------------
#
# Author: Damie Pak
# Copyright (c) Damie Pak, 2023
# Email:  dp495@cornell.edu
#
# Date: July 12, 2023
#
# Script Name: 01_RBC_Mortality_Analysis.R
#
# Script Description: This is the script that analyzes the 
# mouse RBC time series. We use it to determine the mortality
#threshold. This mortality threshold is crucial for all 
#future analysis as it 
#We then use use the figure for a schematic in Figure 2
#and we vary the mortality threshold for a supp figure. 
#
# Notes: 
#

#######################################
###Analysis of the RBC plots of mice###
#######################################
library(here)

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))

###Main data for analysis 
main_dat <- read.csv(here("Data","RBC_plots_MICE.csv"))

###Split based on mouse and treatment
splitted_mt <- split(main_dat, list(main_dat$Mouse,main_dat$Treatment),drop=TRUE)

###Mice that died have NA in their time series and therefore, I figure out which 
###indices belong to either the died/survived groups.

###Mice that died - 15 treatments
which_mort <- which(do.call(rbind,lapply(splitted_mt, function(x) any(is.na(x) == TRUE))) == TRUE)
###Mice that survived - 25 treatments
which_surv <- which(do.call(rbind,lapply(splitted_mt, function(x) any(is.na(x) == TRUE))) == FALSE)

mice_mortality <- splitted_mt[which_mort]
mice_survival <- splitted_mt[which_surv]

###Of the mice that died, I'm looking at the minimum red blood cell
Mice_MORT_MIN <-do.call(rbind,
                        lapply(mice_mortality, function(x) 
                        cbind.data.frame(
                              Min = min(x$RBC,na.rm = TRUE),
                              Mouse = unique(x$Mouse),
                              Treatment = unique(x$Treatment),
                              Mort = "Mortality")))

###Of the mice that survived, I'm looking at the minimum red blood cell
Mice_SURV_MIN <-do.call(rbind,
                        lapply(mice_survival, function(x) 
                        cbind.data.frame(
                                 Min = min(x$RBC,na.rm = TRUE),
                                 Mouse = unique(x$Mouse),
                                 Treatment = unique(x$Treatment),
                                 Mort = 'Survival')))

min(Mice_SURV_MIN$Min)

Mice_ALL <- rbind(Mice_MORT_MIN,Mice_SURV_MIN)

###Setting the mortality threshold by density

MICE_MORTALITY_THRESHOLD <- 
  ggplot(Mice_ALL,(aes(x = Mort, y = Min)))+
  geom_boxplot(width = 0.5)+
  geom_beeswarm(size = 3,cex = 4)+
  geom_hline(yintercept = 0.65, color='#e6154c', size=0.9)+
  xlab("")+
  ylab("Minimum RBC (x 10^6)")+
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 14,color='black'),
        axis.title = element_text(size = 15, color = 'black' ))

###So the threshold of death is  ~6.7 *10^5, rounded to 6.5 * 10^5

ggsave(here("Figures", "01F_MouseMortalityData.pdf"), width = 4, height =5,
       units ='in')


###An aside: seeing what the variability is in the initial RBCs

###Before the infection is t = 0,

R_Initial <- do.call(rbind,lapply(splitted_mt , function(x) x[x$Day.PI == -5, c("RBC","Treatment")]))


R_Initial_NoPHZ <- subset(R_Initial, 
                          R_Initial$Treatment %in% c("G","E","C","A"))


min(R_Initial_NoPHZ $RBC) # 5.413 * 10^6
max(R_Initial_NoPHZ $RBC) # 8.573 * 10^6
mean(R_Initial_NoPHZ$RBC) #6.85 * 10^6

summary(R_Initial_NoPHZ$RBC)


ggplot(R_Initial_NoPHZ, 
       aes(x = 1 , y = RBC)) + 
  geom_boxplot(width = 0.5)+
  geom_beeswarm(size = 4,cex = 4) + 
  xlab("")+
  ylab("Initial RBC (x 10^6)")+
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 14,color='black'),
        axis.title = element_text(size = 15, color = 'black'))

ggsave(here("Figures", "Supp_MouseInitialRBCData.pdf"), width = 4, height =5,
       units ='in')


