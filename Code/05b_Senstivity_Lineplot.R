# HEADER --------------------------------------------
#
# Author: Damie Pak
# Copyright (c) Damie Pak, 2023
# Email:  dp495@cornell.edu
#
# Date: 2023-04-21
#
# Script Name: 01_RBC_Mortality_Analysis.R
#
# Script Description: This is the script that is for 
#Figure 6 which shows how changes in the mortality rate of
#merozoite (uM), the lambda rate, and the initial RBCS.
#
# Notes: 
#


library(here)
### Packages to load
source(here("Code", "Helper_Function_Code", "FUNC_00_Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "FUNC_00_best_long_strategyfinder.R"))
### We need the data here and give everyone the grouping name

###Original model
FITNESS_MODEL <- as.data.frame(fread(here(
  "Output", "Fitness_Model",
  "Fitness_MODEL_PC_med.csv"
)))[,-1]
FITNESS_MODEL$change <- 1
FITNESS_MODEL$grouping <- "Original"

###Blood replenishment model
FITNESS_L <- as.data.frame(fread(here(
  "Output", "Fitness_Model",
  "Fitness_L.csv"
)))
FITNESS_L$grouping <- "lambda"

###Initial RBC model
FITNESS_RI <- as.data.frame(fread(here(
  "Output", "Fitness_Model",
  "Fitness_RI.csv"
)))
FITNESS_RI$grouping <- "RI"

###Merozoite mortality model
FITNESS_UM <- as.data.frame(fread(here(
  "Output", "Fitness_Model",
  "Fitness_UM.csv"
)))
FITNESS_UM$grouping <- "UM"

##################
### Day Labeler###
##################
Parameter_Names <- c(
  "lambda" = "A. RBC replenishment",
  "RI" = "B. Initial RBC",
  "UM" = "C. Merozoite mortality rate"
)

FITNESS_ALL <- rbind(FITNESS_MODEL, FITNESS_L[,-1], FITNESS_RI[,-1], FITNESS_UM[,-1])
FITNESS_ALL_SPLIT <- split(FITNESS_ALL, FITNESS_ALL$change)
FITNESS_MAX <- do.call(rbind, by(FITNESS_ALL, FITNESS_ALL$change, function(x) x[which.max(x$end_fitness), ],
  simplify = FALSE
))

########################################
### CONSTANT TRANSMISSION INVESTMENT ####
########################################

Optimal_Transmission_Investment = 0.76


FITNESS_OPT<- subset(FITNESS_ALL, FITNESS_ALL$C_V == Optimal_Transmission_Investment &
  FITNESS_ALL$change != 1 &
  FITNESS_ALL$status == "success")

na.omit(by(FITNESS_OPT , list(FITNESS_OPT$change,
                             FITNESS_OPT $grouping),
   Best_Strategy_Finder))



FITNESS_OPT_O <- subset(FITNESS_ALL, FITNESS_ALL$C_V == 0.76 &
  FITNESS_ALL$change == 1 &
  FITNESS_ALL$status == "success")[, c("B_V", "end_fitness")]



###Have to do to a silly convoluted way because facet_wrap and
###would not let me have x-axis on both facets without letting me 
###changing through scale_x_continuous


Fitness_CosntantC_Lineplot_Merozoite_Mortality <- 
  ggplot(
  subset(FITNESS_OPT,
         FITNESS_OPT$grouping == 'UM'),
  aes(
    x = B_V,
    y = end_fitness,
    group = as.factor(change),
    color = as.factor(change)))+
  geom_line(size = 1.0) 
  scale_x_continuous(breaks = seq(12,32,2), limits=c(12,32))+
  annotate("line",
    x = FITNESS_OPT_O$B_V,
    y = FITNESS_OPT_O$end_fitness,
    size = 1.2
  ) +
  scale_color_manual(values = c(
    "#ef2366", "#1199ee"
  )) +
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Cumulative \ntransmission potential") +
  ggtitle("A. Merozoite mortality, uM")+
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "off",
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 15))


Fitness_CosntantC_Lineplot_Initial_RI <- 
  ggplot(
    subset(FITNESS_OPT,
           FITNESS_OPT$grouping == 'RI'),
    aes(
      x = B_V,
      y = end_fitness,
      group = as.factor(change),
      color = as.factor(change)))+
  geom_line(size = 1.0) +
  scale_x_continuous(breaks = seq(12,32,2), limits=c(12,32))+
  
  annotate("line",
           x = FITNESS_OPT_O$B_V,
           y = FITNESS_OPT_O$end_fitness,
           size = 1.2
  ) +
  scale_color_manual(values = c(
     "#1199ee","#ef2366"
  )) +
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Cumulative \ntransmission potential") +
  ggtitle("B. Initial red blood cells, R(0)")+
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "off",
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 15))

Fitness_CosntantC_Lineplot_Merozoite_Mortality/Fitness_CosntantC_Lineplot_Initial_RI


ggsave(here("Figures", "Raw", "Supp_Optimal_75_SA.pdf"),
  height = 6, width = 5, units = "in"
)

###Supplementary:

Fitness_CosntantC_Lineplot_Lambda <- 
  ggplot(
    subset(FITNESS_OPT,
           FITNESS_OPT$grouping == "lambda"),
    aes(
      x = B_V,
      y = end_fitness,
      group = as.factor(change),
      color = as.factor(change)))+
  geom_line(size = 1.0) +
  scale_x_continuous(breaks = seq(12,32,2),limits=c(12,32))+
  
  annotate("line",
           x = FITNESS_OPT_O$B_V,
           y = FITNESS_OPT_O$end_fitness,
           size = 1.2
  ) +
  scale_color_manual(values = c(
    "#ef2366", "#1199ee"
  )) +
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Cumulative \ntransmission potential (f)") +
  ggtitle("RBC replenishment, uM")+
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "off",
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 15))

ggsave(here("Figures", "Raw", "Supp_Line_Lambda.pdf"),
       height = 3, width = 5, units = "in"
)
