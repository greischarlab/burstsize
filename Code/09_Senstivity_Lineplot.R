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
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))

### We need the data here and give everyone the grouping name

###Original model
FITNESS_MODEL <- as.data.frame(fread(here(
  "Output", "Fitness_Model",
  "Fitness_MODEL_PC_med.csv"
)))
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

FITNESS_ALL <- rbind(FITNESS_MODEL, FITNESS_L, FITNESS_RI, FITNESS_UM)
FITNESS_ALL_SPLIT <- split(FITNESS_ALL, FITNESS_ALL$change)
FITNESS_MAX <- do.call(rbind, by(FITNESS_ALL, FITNESS_ALL$change, function(x) x[which.max(x$end_fitness), ],
  simplify = FALSE
))
FITNESS_MAX$RC <- FITNESS_MAX$B_V * (1 - FITNESS_MAX$C_V)
FITNESS_MAX$increase <- ifelse(FITNESS_MAX$RC > 3.625, "increase", "decrease")
FITNESS_MAX$cum_trans <- ifelse(FITNESS_MAX$end_fitness > 28.04553,
  "Fitter", "Not_Fitter"
)

########################################
### OPTIMAL TRANSMISSION INVESTMENT #####
########################################

FITNESS_ALL_OPTIMAL_CV <- NULL
for (k in seq(1, nrow(FITNESS_MAX))) {
  tmp <- FITNESS_MAX[k, ]
  tmp_fitness <- FITNESS_ALL_SPLIT[[k]]

  optimal_cv <- subset(tmp_fitness, tmp_fitness$C_V == tmp$C_V &
    tmp_fitness$status == "success")

  FITNESS_ALL_OPTIMAL_CV[[k]] <- cbind.data.frame(
    B_V = optimal_cv$B_V,
    C_V = optimal_cv$C_V,
    change = optimal_cv$change,
    grouping = optimal_cv$grouping,
    endfitness = optimal_cv$end_fitness
  )
}

FITNESS_ALL_OPTIMAL_CV_F <- do.call(rbind, FITNESS_ALL_OPTIMAL_CV)

FITNESS_ALL_OPTIMAL_CV_1 <- subset(
  FITNESS_ALL_OPTIMAL_CV_F,
  FITNESS_ALL_OPTIMAL_CV_F$change != 1
)


FITNESS_ALL_OPTIMAL_CV_2 <- subset(FITNESS_ALL_OPTIMAL_CV_F, FITNESS_ALL_OPTIMAL_CV_F$C_V == 0.75 &
  FITNESS_ALL_OPTIMAL_CV_F$change == "Original")[, c("B_V", "endfitness")]


OPT_CV_GG <- ggplot(
  FITNESS_ALL_OPTIMAL_CV_1,
  aes(
    x = B_V,
    y = endfitness,
    group =as.factor(change),
    color = as.factor(change)
  )
) +
  geom_line(size = 1.2) +
  facet_wrap(~grouping,
    ncol = 1,
    labeller = as_labeller(Parameter_Names)
  ) +
  annotate("line",
    x = FITNESS_ALL_OPTIMAL_CV_2$B_V,
    y = FITNESS_ALL_OPTIMAL_CV_2$endfitness,
    size = 1.2
  ) +
  scale_color_manual(values = c(
    "#ef233c", "#118ab2",
    "#118ab2", "#ef233c",
    "#ef233c", "#118ab2"
  )) +
  xlab("Burst size") +
  ylab("Cumulative transmission potential") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "off",
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 15)
  )

########################################
### Optimal burst size              ####
########################################
FITNESS_ALL_OPTIMAL_BV <- NULL
for (k in seq(1, nrow(FITNESS_MAX))) {
  tmp <- FITNESS_MAX[k, ]
  tmp_fitness <- FITNESS_ALL_SPLIT[[k]]

  optimal_bv <- subset(tmp_fitness, tmp_fitness$B_V == tmp$B_V &
    tmp_fitness$status == "success")



  FITNESS_ALL_OPTIMAL_BV[[k]] <- cbind.data.frame(
    B_V = optimal_bv$B_V,
    C_V = optimal_bv$C_V,
    change = optimal_bv$change,
    grouping = optimal_bv$grouping,
    endfitness = optimal_bv$end_fitness
  )
}

FITNESS_ALL_OPTIMAL_BV_F <- do.call(rbind, FITNESS_ALL_OPTIMAL_BV)

FITNESS_ALL_OPTIMAL_BV_1 <- subset(
  FITNESS_ALL_OPTIMAL_BV_F,
  FITNESS_ALL_OPTIMAL_BV_F$change != 1
)


FITNESS_ALL_OPTIMAL_BV_2 <- subset(
  FITNESS_ALL_OPTIMAL_BV_F,
  FITNESS_ALL_OPTIMAL_BV_F$B_V == 8 &
    FITNESS_ALL_OPTIMAL_BV_F$change == "Original"
)[, c("C_V", "endfitness")]




OPT_BV_GG <- ggplot(
  FITNESS_ALL_OPTIMAL_BV_1,
  aes(
    x = C_V,
    y = endfitness,
    group = as.factor(change),
    color = as.factor(change)
  )
) +
  geom_line(size = 1.2) +
  facet_wrap(~grouping,
    ncol = 1,
    labeller = as_labeller(Parameter_Names)
  ) +
  annotate("line",
    x = FITNESS_ALL_OPTIMAL_BV_2$C_V,
    y = FITNESS_ALL_OPTIMAL_BV_2$endfitness,
    size = 1.2
  ) +
  scale_color_manual(values = c(
    "#ef233c", "#118ab2",
    "#118ab2", "#ef233c",
    "#ef233c", "#118ab2"
  )) +
  xlab("Transmission investment") +
  ylab("Cumulative transmission potential") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "off",
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 15)
  )


OPT_CV_GG + OPT_BV_GG


ggsave(here("Figures", "Raw", "Supp_Optimal_BV_CV_SA.pdf"),
  height = 8, width = 9, units = "in"
)

########################################
### CONSTANT TRANSMISSION INVESTMENT ####
########################################
FITNESS_75 <- subset(FITNESS_ALL, FITNESS_ALL$C_V == 0.75 &
  FITNESS_ALL$change != 1 &
  FITNESS_ALL$status == "success")

FITNESS_75_O <- subset(FITNESS_ALL, FITNESS_ALL$C_V == 0.75 &
  FITNESS_ALL$change == 1 &
  FITNESS_ALL$status == "success")[, c("B_V", "end_fitness")]


Fitness_CosntantC_Lineplot <- ggplot(
  FITNESS_75,
  aes(
    x = B_V,
    y = end_fitness,
    group = as.factor(change),
    color = as.factor(change)))+
  geom_line(size = 1.2) +
  facet_wrap(~grouping,
    ncol = 1,
    labeller = as_labeller(Parameter_Names)
  ) +
  annotate("line",
    x = FITNESS_75_O$B_V,
    y = FITNESS_75_O$end_fitness,
    size = 1.2
  ) +
  scale_color_manual(values = c(
    "#ef233c", "#118ab2",
    "#118ab2", "#ef233c",
    "#ef233c", "#118ab2"
  )) +
  xlab("Burst size") +
  ylab("Cumulative transmission potential") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "off",
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 15)
  )

ggsave(here("Figures", "Raw", "Supp_Optimal_75_SA.pdf"),
  height = 8, width = 5, units = "in"
)
