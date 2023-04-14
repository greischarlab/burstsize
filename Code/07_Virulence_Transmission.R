########################################################
### Understanding how the viable burst size combinations#
### influence virulence, transmission, and such.     ####
########################################################

library(here)

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "Virulence_Transmission_BurstSize.R"))

###

FULL_MODEL_PC_DT <- as.data.frame(fread(here(
  "Output", "Full_Model",
  "FULL_MODEL_PC_DT.csv"
)))
Fitness_MODEL_PC_FULL <- as.data.frame(fread(here(
  "Output", "Fitness_Model",
  "Fitness_MODEL_PC_FULL.csv"
)))


Main_DF_VT <- Virulence_Transmission_DF_Calculator(
  FULL_MODEL_PC_DT,
  Fitness_MODEL_PC_FULL,
  0.54
)
Main_DF_VT$optimal <- ifelse(Main_DF_VT$CTP == max(Main_DF_VT$CTP), 1, NA)

######################################################################
### Cumulative transmission potential versus peak gametocyte density###
######################################################################

trans_to_gams <-
  ggplot(
    Main_DF_VT,
    aes(
      x = GMax,
      y = CTP,
      group = Status,
      linetype = Status
    )
  ) +
  geom_line(size =1.2) + 
  geom_point(
    aes(
      color = as.factor(optimal),
      shape = as.factor(optimal),
    ),
    size = 4,
    stroke = 1
  ) +
  scale_color_manual(values = c(
    `NA` = NA,
    `1` = "#FF52A3"
  )) +
  scale_linetype_manual(values = c(
    'success' = 1,
    'mort' = 3
  ))+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Peak gametocyte density") +
  ylab("Cumulative transmission potential (f)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black")
  ) 
  



### Proportion of red blood cells compared too burst size
rbcmin_to_burst_size <-
  ggplot(
    Main_DF_VT,
    aes(x = B_V, 
        y = RMin,
        group = Status,
        linetype = Status
    )
  ) +
  geom_line(size =1.2) + 
  geom_point(
    aes(
      color = as.factor(optimal),
      shape = as.factor(optimal),
    ),
    size = 4,
    stroke = 1
  ) +
  scale_color_manual(values = c(
    `NA` = NA,
    `1` = "#FF52A3"
  )) +
  scale_linetype_manual(values = c(
    'success' = 1,
    'mort' = 3
  ))+ 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(expression(paste("Burst size", "( ",beta,")"))) +
  ylab("Minimum RBC density") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black")
  )


### Maximum gametocyte density to burst size 
gammax_to_burst_size <-
  ggplot(
    Main_DF_VT,
    aes(
      x = B_V,
      y = GMax,
      group = Status,
      linetype = Status
    )
  ) +
  geom_line(size =1.2) + 
  geom_point(
    aes(
      color = as.factor(optimal),
      shape = as.factor(optimal),
    ),
    size = 4,
    stroke = 1
  ) +
  scale_color_manual(values = c(
    `NA` = NA,
    `1` = "#FF52A3"
  )) +
  scale_linetype_manual(values = c(
    'success' = 1,
    'mort' = 3
  ))+ 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(expression(paste("Burst size", "( ",beta,")"))) +
  ylab("Peak gametocyte density") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black")
  )

### Length of the acute phase to the burst size 
acute_phase_to_burst_size <-
  ggplot(
    Main_DF_VT,
    aes(
      x = B_V,
      y = End_Time,
      group = Status,
      linetype = Status
    )
  ) +
  geom_line(size =1.2) + 
  geom_point(
    aes(
      color = as.factor(optimal),
      shape = as.factor(optimal),
    ),
    size = 4,
    stroke = 1
  ) +
  scale_color_manual(values = c(
    `NA` = NA,
    `1` = "#FF52A3"
  )) +
  scale_linetype_manual(values = c(
    'success' = 1,
    'mort' = 3
  ))+ 
  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(expression(paste("Burst size", "( ",beta,")"))) +
  ylab("Acute phase") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black")
  )
  xlab(expression(paste("Burst size", "( ",beta,")"))) 


(trans_to_gams) / (rbcmin_to_burst_size + gammax_to_burst_size + acute_phase_to_burst_size) +
    plot_annotation(tag_levels = "A")

  ggsave(here("Figures", "Raw", "Virulence_Transmission.pdf"), height =8, width =10,
         units = 'in')

  