


MIN_RBC_GG_Plotter <- function(x){
### Proportion of red blood cells compared too burst size
rbcmin_to_burst_size_GG <-
  ggplot(
   x,
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
  ylab("Minimum RBC \ndensity") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black")
  )
}

### Maximum gametocyte density to burst size
GCUM_GG_Plotter  <- function(x){
  ggplot(
   x,
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
  ylab("Cumulative \ngametocyte density") +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 16, color = "black")
    )
}
##################################################
### Length of the acute phase to the burst size ##
##################################################
ACUTE_PHASE_GG_Plotter <- function(x){
  ggplot(
    x,
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
  ylab("Acute phase \n(days)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black")
  )
}

Triplot_Plotter <- function(x,name){
  
  rbcmin_to_burst_size_GG <- MIN_RBC_GG_Plotter(x)
  gcum_to_burst_size_GG <- GCUM_GG_Plotter(x)
  acute_phase_to_burst_size_GG <- ACUTE_PHASE_GG_Plotter(x)


  rbcmin_to_burst_size_GG / gcum_to_burst_size_GG / acute_phase_to_burst_size_GG +
  plot_annotation(tag_levels = "A")


  ggsave(here("Figures","Raw", paste("Virulence_Transmission",name,".pdf",sep = "")),
       height = 8, width = 4, unit = "in")

  print(rbcmin_to_burst_size_GG / gcum_to_burst_size_GG / acute_phase_to_burst_size_GG +
          plot_annotation(tag_levels = "A"))
}
