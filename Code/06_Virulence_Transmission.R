# HEADER --------------------------------------------
#
# Author: Damie Pak
# Copyright (c) Damie Pak, 2023
# Email:  dp495@cornell.edu
#
# Date: 2023-04-21
#
# Script Name: 06_Virulence_Transmission .R
#
# Script Description: This is the script that looks at the fitness
# output and investigate how different burst size affects the 
# virulence, the total cumulative gametocyctes, and the acute phase
# duration
#
# Notes: 
#

########################################################
### Understanding how the viable burst size combinations#
### influence virulence, transmission, and such.     ####
########################################################

library(here)

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "02_Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "Virulence_Transmission_BurstSize.R"))
source(here("Code","Simulator_Code", "Simulator_Main_Gflux.R"))
###

FULL_MODEL_PC_DT <- as.data.frame(fread(here(
  "Output", "Full_Model",
  "FULL_MODEL_PC_DT.csv"
)))
Fitness_MODEL_PC_FULL <- as.data.frame(fread(here(
  "Output", "Fitness_Model",
  "Fitness_MODEL_PC_FULL.csv"
)))

### burst size for a transmission investment of 54%
BV_54 <- subset(Fitness_MODEL_PC_FULL,
                Fitness_MODEL_PC_FULL$C_V == 0.54 &
                Fitness_MODEL_PC_FULL$status != 'Fail')

BV_54_Burstsize <- BV_54$B_V
BV_54_CV <- BV_54$C_V
BV_54_Endtime <- BV_54$endtime

### We need to rerun the Fitness_Model for the Gflux
Fitness_MODEL_PC_GFlux <- mcmapply(Simulator_MalariaPC_DDE_BC_GFLUX,
                             c(BV_54_Burstsize),
                             c(0.54),
                             c(BV_54_Endtime),
                             mc.cores = 1,
                             SIMPLIFY = FALSE)

for (strain in seq(1,nrow(BV_54))){
  Fitness_MODEL_PC_GFlux[[strain]]$status <- BV_54$status[[strain]]
}

Gflux_DF <- do.call(rbind,
        lapply(Fitness_MODEL_PC_GFlux, Virulence_Gam_Finder))

Gflux_DF$optimal <- ifelse(Gflux_DF$B_V == 8, 1,NA)
### Proportion of red blood cells compared too burst size
rbcmin_to_burst_size_GG <-
  ggplot(
    Gflux_DF,
    aes(
      x = B_V,
      y = RMin,
      group = status,
    )
  ) +
  geom_line(aes(color = status),
            size = 1.2
  ) +
  scale_color_manual(values = c(
    "success" = "black",
    "mort" = "#83d4a5"
  )) +
  new_scale_color() +
  geom_point(
    aes(
      color = as.factor(optimal),
      shape = as.factor(optimal),
    ),
    size = 3.5,
    stroke = 1
  ) +
  scale_y_log10(limits = c(10^4.5,10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_color_manual(values = c(
    `NA` = NA,
    `1` = "#FF52A3"
  ))  +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Minimum RBC density (log10)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black")
  )


### Maximum gametocyte density to burst size
gcum_to_burst_size_GG  <-
  ggplot(
    Gflux_DF,
    aes(
      x = B_V,
      y = GMax,
      group = status,
    )
  ) +
  geom_line(aes(color = status),
            size = 1.2
  ) +
  scale_color_manual(values = c(
    "success" = "black",
    "mort" = "#83d4a5"
  )) +
  new_scale_color() +
  geom_point(
    aes(
      color = as.factor(optimal),
      shape = as.factor(optimal),
    ),
    size = 3.5,
    stroke = 1
  ) +
  scale_color_manual(values = c(
    `NA` = NA,
    `1` = "#FF52A3"
  )) +
  scale_y_log10(limits = c(10^4.5,10^7),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Cumulative gametocyte density (log10)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black")
  )

##################################################
### Length of the acute phase to the burst size ##
##################################################
acute_phase_to_burst_size_GG  <-
  ggplot(
    Gflux_DF,
    aes(
      x = B_V,
      y = endtime,
      group = status,
    )
  ) +
  geom_line(aes(color = status),
    size = 1.2
  ) +
  scale_color_manual(values = c(
    "success" = "black",
    "mort" = "#83d4a5"
  )) +
  new_scale_color() +
  geom_point(
    aes(
      color = as.factor(optimal),
      shape = as.factor(optimal),
    ),
    size = 3.5,
    stroke = 1
  ) +
  scale_color_manual(values = c(
    `NA` = NA,
    `1` = "#FF52A3"
  )) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(expression(paste("Burst size", "( ", beta, ")"))) +
  ylab("Acute phase (days)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black")
  )


rbcmin_to_burst_size_GG / gcum_to_burst_size_GG / acute_phase_to_burst_size_GG +
    plot_annotation(tag_levels = "A")

ggsave(here("Figures", "Raw", "Virulence_Transmission.pdf"), height = 8, width = 4,
      units = 'in')

  