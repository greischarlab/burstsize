library(here)

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "02_Fitness_Functions.R"))
source(here("Code", "Helper_Function_Code", "04_RM_Function.R"))

###

FULL_MODEL_PC_DT <- as.data.frame(fread(here(
  "Output", "Full_Model",
  "FULL_MODEL_PC_DT.csv"
)))
Fitness_MODEL_PC_FULL <- as.data.frame(fread(here(
  "Output", "Fitness_Model",
  "Fitness_MODEL_PC_FULL.csv"
)))


# R_M#

rate_PMR_dataframe1 <- rate_PMR_data(
  FULL_MODEL_PC_DT,
  Fitness_MODEL_PC_FULL,
  c(33,34,35),
  c(0.89),
  48
)

rate_PMR_dataframe2 <- rate_PMR_data(
  FULL_MODEL_PC_DT,
  Fitness_MODEL_PC_FULL,
  c(32 ),
  c(0.88),
  48
)

rate_PMR_dataframe <- rbind(rate_PMR_dataframe1[[1]])
                            #rate_PMR_dataframe2[[1]]

end_points1 <- rate_PMR_dataframe1[[2]]
end_points2 <- rate_PMR_dataframe2[[2]]

end_points <- rbind(end_points1)#end_points2)

### Expected merozoites
mer_GG <-
  ggplot(
    rate_PMR_dataframe,
    aes(
      x = time,
      y = rate,
      fill = as.factor(B_V),
      color = as.factor(B_V),
      group = as.factor(B_V)
    )
  ) +

  geom_line(aes(alpha = which.time), size = 1.2) +
  geom_point(
    data = end_points,
    aes(x = endtime, y = PMR_END),
    color = "black",
    size = 3.2,
    shape = 21
  ) +
  geom_hline(
    yintercept = 1,
    alpha = 1,
    color = "black"
  ) +
  geom_segment(
    data = end_points,
    aes(
      x = PMR_root,
      xend = PMR_root,
      y = 0.7,
      yend = 1,
      color = as.factor(B_V)
    ),
    linetype = 3, size = 1.1
  ) +
  scale_color_manual(
    values = c(
      "black",
      "#3DD8EE",
      "#FD7E6C"
    ),
    name = "Burst size"
  ) +
  scale_fill_manual(
    values = c(
      "black",
      "#3DD8EE",
      "#FD7E6C"
    )
  ) +
  xlab("Days post-infection") +
  ylab("Expected merozoites invading (uncommited) (R_M)") +
  scale_x_continuous(
    breaks = seq(0, 100, 10),
    limits = c(0, 100)
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black")
  )


dailytp_GG <-
  ggplot(
    rate_PMR_dataframe,
    aes(
      x = time, y = Daily_Trans_Prob,
      fill = as.factor(B_V),
      color = as.factor(B_V),
      group = as.factor(B_V)
    )
  ) +
  geom_line(aes(alpha = which.time), size = 1.2) +
  geom_segment(
    data = end_points,
    aes(
      x = Daily_Trans_Root,
      xend = Daily_Trans_Root,
      y = 0,
      yend = 1,
      color = as.factor(B_V)
    ),
    linetype = 3, size = 1.1
  ) +
  geom_point(
    data = end_points,
    aes(
      x = endtime,
      y = Daily_Trans_End
    ),
    color = "black",
    size = 3.2,
    shape = 21
  ) +
  scale_color_manual(
    values =
      c("#3DD8EE", "black", "#FD7E6C"),
    name = "Burst size"
  ) +
  scale_fill_manual(values = c(
    "#3DD8EE",
    "black",
    "#FD7E6C"
  )) +
  ylab("Daily transmission \n probability (P_c)") +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(0, 100, 10),
    limits = c(0, 100)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black")
  )


### CUMULATIVE TRANSMISSION POTENTIAL
cumcp_GG <-
  ggplot(
    rate_PMR_dataframe,
    aes(
      x = time, y = Cum_Trans_Potential,
      fill = as.factor(B_V),
      color = as.factor(B_V),
      group = as.factor(B_V)
    )
  ) +
  geom_line(aes(alpha = which.time), size = 1.2) +
  geom_segment(
    data = end_points,
    aes(
      x = Cum_Trans_Root,
      xend = Cum_Trans_Root,
      y = 0,
      yend = 91,
      color = as.factor(B_V)
    ),
    linetype = 3, size = 1.1
  ) +
  geom_point(
    data = end_points,
    aes(x = endtime, y = Cum_Trans_End),
    color = "black",
    size = 3.2,
    shape = 21
  ) +
  scale_color_manual(
    values = c("#3DD8EE", "black", "#FD7E6C"),
    name = "Burst size"
  ) +
  scale_fill_manual(values = c(
    "#3DD8EE",
    "black",
    "#FD7E6C"
  )) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(0, 100, 10),
    limits = c(0, 100)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  xlab("Days post-infection") +
  ylab("Cumulative \ntransmission potential (f)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black")
  )


mer_GG / dailytp_GG / cumcp_GG + plot_annotation(tag_levels = "A") 




ggsave(file = here("Figures", "Raw", "RM_3_Plot.pdf"), width = 8, height = 9, units = "in")
