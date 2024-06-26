RM_Plotter <- function(x,y){
  
  ### Expected merozoites
  mer_GG <-
    ggplot(
      x,
      aes(
        x = time,
        y = rate,
        fill = as.factor(B_V),
        color = as.factor(B_V),
        group = as.factor(B_V)
      )
    ) +
    
    geom_line(aes(alpha = which.time), size = 1.1) +
    geom_point(
      data = y ,
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
      data = y ,
      aes(
        x = PMR_root,
        xend = PMR_root,
        y = 0.5,
        yend = 1,
        color = as.factor(B_V)
      ),
      linetype = 3, size = 1.1
    ) +
    scale_color_manual(
      values = c(
        "#acdbe8",
        "#FF116B",
        "black"
      ),
      name = "Burst size"
    ) +
    scale_fill_manual(
      values = c(
        "#acdbe8",
        "#FF116B",
        "black"
      ),
    ) +
    xlab("Days post-infection") +
    ylab("Effective \n merozoite number (R_M)") +
    scale_x_continuous(
      expand=c(0,0),
      breaks = seq(0, 50, 5),
      limits = c(0, 50)
    ) +
    scale_y_continuous(
      expand=c(0,0),
      breaks= seq(0.5, 2.25, 0.25)
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
  
  
  dailyG_GG <-
    ggplot(
      x,
      aes(
        x = time, y = log10(Gam_Time+1),
        fill = as.factor(B_V),
        color = as.factor(B_V),
        group = as.factor(B_V)
      )
    ) +
    geom_line(aes(alpha = which.time), size = 1.1) +
    geom_segment(
      data = y,
      aes(
        x = Daily_Trans_Root,
        xend = Daily_Trans_Root,
        y = 0,
        yend = 8,
        color = as.factor(B_V)
      ),
      linetype = 3, size = 1.1
    ) +
    geom_point(
      data = y,
      aes(
        x = endtime,
        y = log10(Gam_End+1)
      ),
      color = "black",
      size = 3.2,
      shape = 21
    ) +
    scale_color_manual(
      values = c(
        "#acdbe8",
        "#FF116B",
        "black"
      ),
      name = "Burst size"
    ) +
    scale_fill_manual(
      values = c(
        "#acdbe8",
        "#FF116B",
        "black"
      ),
    ) +
    ylab("Gametocyte density") +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = seq(0, 50, 5),
      limits = c(0, 50)
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
  
  
  
  dailytp_GG <-
    ggplot(
      x,
      aes(
        x = time, y = Daily_Trans_Prob,
        fill = as.factor(B_V),
        color = as.factor(B_V),
        group = as.factor(B_V)
      )
    ) +
    geom_line(aes(alpha = which.time), size = 1.1) +
    geom_segment(
      data = y,
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
      data = y,
      aes(
        x = endtime,
        y = Gam_End
      ),
      color = "black",
      size = 3.2,
      shape = 21
    ) +
    scale_color_manual(
      values = c(
        "#acdbe8",
        "#FF116B",
        "black"
      ),
      name = "Burst size"
    ) +
    scale_fill_manual(
      values = c(
        "#acdbe8",
        "#FF116B",
        "black"
      ),
    )  +
    ylab("Daily transmission \n probability (P_c)") +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = seq(0, 50, 5),
      limits = c(0, 50)
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.1)) +
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
      x,
      aes(
        x = time, y = Daily_Trans_Prob,
        fill = as.factor(B_V),
        color = as.factor(B_V),
        group = as.factor(B_V)
      )
    ) +
    geom_line(aes(alpha = which.time), size = 1.1) +
    geom_segment(
      data = y,
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
      data = y,
      aes(
        x = endtime,
        y = Daily_Trans_End
      ),
      color = "black",
      size = 3.2,
      shape = 21
    ) +
    scale_color_manual(
      values = c(
        "#acdbe8",
        "#FF116B",
        "black"
      ),
      name = "Burst size"
    ) +
    scale_fill_manual(
      values = c(
        "#acdbe8",
        "#FF116B",
        "black"
      ),
    )  +
    xlab("Days post-infection")+
    ylab("Daily transmission \n probability (P_c)") +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = seq(0, 50, 5),
      limits = c(0, 50)
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.1)) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 13, color = "black")
    )
  
  
  
  mer_GG /dailyG_GG/ dailytp_GG + plot_annotation(tag_levels = "A") 
  
  
}
