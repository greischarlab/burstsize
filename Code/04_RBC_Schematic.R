#######################
###Schematic Figure ###
#######################
library(here)

FULL_MODEL_PC_DT <- as.data.frame(fread(here("Output","Full_Model", 
                                  "FULL_MODEL_PC_DT.csv")))

###This is the schematic figure on the three ways that
###infection can turn out (1) kills the host OR (2) does not
###lead to unestablished infection and (3)

###These are of interest
R_C_V_Examples <- subset(FULL_MODEL_PC_DT,
                         FULL_MODEL_PC_DT$B_V %in% c(13,40) &
                         round(FULL_MODEL_PC_DT$C_V, 2) ==  0.35 & 
                           FULL_MODEL_PC_DT$time < 11)

R_C_V_Examples$Status <- ifelse(R_C_V_Examples$B_V == 13, "success", 'mort')
R_C_V_Examples$I_ALL <- rowSums(R_C_V_Examples [,4:103])

R_C_V_Full <- R_C_V_Examples[,c("B_V","time","R", "I_ALL", "Status")]

R_C_V_Full <- melt(R_C_V_Full, id.vars= c("B_V","Status",'time'))


###Successful Infection AND Killing infection

ggplot(R_C_V_Full, aes(x = time, y = log10(value), 
                       color = variable, group =variable))+
       geom_line() + 
       geom_hline(yintercept = log10(6.5 * 10^5))+
       facet_wrap(~Status, ncol = 1) +
       theme_void()

ggsave(file = here("Figures","Raw",
                   "Schematic_infections.pdf"), 
       width = 4, height=7, units='in')
