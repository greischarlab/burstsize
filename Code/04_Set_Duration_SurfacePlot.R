# HEADER --------------------------------------------
#
# Author: Damie Pak
# Copyright (c) Damie Pak, 2023
# Email:  dp495@cornell.edu
#
# Date: 2023-04-21
#
# Script Name: 04_Set_Duration_SurfacePlot.R
#
# Script Description: This is a script that looks at the cumulative
# transmission potential at specific end points.
#
# Notes:
#
#
#
###################
### Set durations##
###################
library(here)

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "Fitness_SD_Function.R"))
source(here("Code", "Helper_Function_Code", "Fitness_Functions.R"))

FULL_MODEL_PC_DT <- as.data.frame(fread(here(
  "Output", "Full_Model",
  "FULL_MODEL_PC_DT.csv"
)))

FULL_MODEL_PC <- FULL_MODEL_PC_DT [,-1]

FULL_MODEL_PC  <- split(FULL_MODEL_PC , list(FULL_MODEL_PC$B_V,
                                             FULL_MODEL_PC$C_V), 
                                             drop = TRUE)

remove(FULL_MODEL_PC_DT)
gc()



### Burst Size Versus Transmission Investment ###
B_V <- seq(1, 50, 0.5) # Burst size
C_V <- seq(.01, 1, 0.01) # Transmission investment
B_V_C_V <- expand.grid(B_V = B_V, C_V = C_V) # Different combinations

###############################
### The time point of interests#
###############################

time_point_of_interest <- c(10, 20, 30, 40, 50, 60,100)
lifespan_interest <- 100


Set_Duration_DF<- 
  Calculate_Set_Duration_CP(FULL_MODEL_PC,
                            B_V_C_V,
                            time_point_of_interest,
                            "no")
  

Set_Duration_DF_lxmx <- 
  Calculate_Set_Duration_CP(FULL_MODEL_PC,
                            B_V_C_V,
                            lifespan_interest, 
                            "yes"
)

Grapher_Mortality_DF <- NULL
for (k in seq(1, length(Set_Duration_DF))) {
  Grapher_Mortality_DF[[k]] <-
    cbind.data.frame(grapher_mortality_boundary(Set_Duration_DF[[k]]),
      end_time = time_point_of_interest[[k]]
    )
}



Grapher_Mortality_DF <- do.call(rbind, Grapher_Mortality_DF)
Set_Duration_DF <- do.call(rbind, Set_Duration_DF)


### Finding the maximum fitness points for each of the time point of interest

Max_Point_DF <- do.call(rbind, by(Set_Duration_DF, Set_Duration_DF$end_time, function(x) {
  x[which.max(x$fitness), ]
}, simplify = FALSE))

### Day Labeler
Day_Names <- c(
  `10` = "A. Day 10",
  `20` = "B. Day 20",
  `30` = "C. Day 30",
  `40` = "D. Day 40",
  `50` = "E. Day 50",
  `60` = "F. Day 60"
)

GG_Fitness_Set_Duration <-
  ggplot(
    tmp2,
    aes(x = round(B_V,3), y = round(C_V,3))
  ) +
  geom_raster(aes(fill = fitness)) +
  geom_raster(
    data = subset(
      tmp2,
      tmp2$status == "Fail"
    ),
    aes(x = B_V, y = C_V), fill = "#d1dbe8"
  ) +
  geom_segment(
    data =
      Grapher_Mortality_DF,
    aes(
      x = x, xend = xend,
      y = y, yend = yend, color = id
    ),
    size = 1.2, lineend = "round"
  ) +
  geom_point(
    data = Max_Point_DF, aes(x = B_V, y = C_V),
    size = 2, shape = 21, col = "#ffb5df",
    stroke = 1
  ) +
  facet_wrap(~end_time, labeller = as_labeller(Day_Names)) +
  scale_fill_viridis(option = "inferno", name = "Cumulative \ntransmisison potential") +
  scale_color_manual(
    values = c(
      "fail" = "black",
      "mort" = "#b8fcd5"
    ),
    guide = "none"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression("Burst size" ~ beta)) +
  ylab("Transmission investment (c)") +
  annotate("text",
    x = 12.5, y = 0.92, label = "Unestablished \ninfection",
    size = 5
  ) +
  annotate("text",
    x = 43, y = 0.1, label = "Host \nMortality",
    size = 5, color = "#b8fcd5"
  ) +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    legend.position = "top"
  )


GG_Fitness_Set_Duration

ggsave(
  here(
    "Figures", "Raw",
    "Fitness_Surface_Plot_SD.pdf"
  ),
  width = 12, height = 9
)


####

GG_Fitness_Set_Duration <-
  ggplot(
   
GG_Fitness_Set_Duration <-
  ggplot(Set_Duration_DF_lxmx[[1]],
    aes(x = round(B_V,3), y = round(C_V,3))
  ) +
  geom_raster(aes(fill = fitness)) +
  geom_raster(
    data = subset(
      tmp2,
      tmp2$status == "Fail"
    ),
    aes(x = B_V, y = C_V), fill = "#d1dbe8"
  ) +
  geom_segment(
    data =
      Grapher_Mortality_DF,
    aes(
      x = x, xend = xend,
      y = y, yend = yend, color = id
    ),
    size = 1.2, lineend = "round"
  ) +
  geom_point(
    data = Max_Point_DF, aes(x = B_V, y = C_V),
    size = 2, shape = 21, col = "#ffb5df",
    stroke = 1
  ) +
  facet_wrap(~end_time, labeller = as_labeller(Day_Names)) +
  scale_fill_viridis(option = "inferno", name = "Cumulative \ntransmisison potential") +
  scale_color_manual(
    values = c(
      "fail" = "black",
      "mort" = "#b8fcd5"
    ),
    guide = "none"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression("Burst size" ~ beta)) +
  ylab("Transmission investment (c)") +
  annotate("text",
    x = 12.5, y = 0.92, label = "Unestablished \ninfection",
    size = 5
  ) +
  annotate("text",
    x = 43, y = 0.1, label = "Host \nMortality",
    size = 5, color = "#b8fcd5"
  ) +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    legend.position = "top"
  )

,
    aes(x = round(B_V,3), y = round(C_V,3))
  ) +
  geom_raster(aes(fill = fitness)) +
  geom_raster(
    data = subset(
      tmp2,
      tmp2$status == "Fail"
    ),
    aes(x = B_V, y = C_V), fill = "#d1dbe8"
  ) +
  geom_segment(
    data =
      Grapher_Mortality_DF,
    aes(
      x = x, xend = xend,
      y = y, yend = yend, color = id
    ),
    size = 1.2, lineend = "round"
  ) +
  geom_point(
    data = Max_Point_DF, aes(x = B_V, y = C_V),
    size = 2, shape = 21, col = "#ffb5df",
    stroke = 1
  ) +
  facet_wrap(~end_time, labeller = as_labeller(Day_Names)) +
  scale_fill_viridis(option = "inferno", name = "Cumulative \ntransmisison potential") +
  scale_color_manual(
    values = c(
      "fail" = "black",
      "mort" = "#b8fcd5"
    ),
    guide = "none"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression("Burst size" ~ beta)) +
  ylab("Transmission investment (c)") +
  annotate("text",
           x = 12.5, y = 0.92, label = "Unestablished \ninfection",
           size = 5
  ) +
  annotate("text",
           x = 43, y = 0.1, label = "Host \nMortality",
           size = 5, color = "#b8fcd5"
  ) +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    legend.position = "top"
  )


