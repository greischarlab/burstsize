library(here)

###Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))
source(here("Code", "Helper_Function_Code", "RM_Functions.R"))

###

FULL_MODEL_PC_DT <- as.data.frame(fread(here("Output","Full_Model", 
                                             "FULL_MODEL_PC_DT.csv")))
Fitness_MODEL_PC_FULL <- as.data.frame(fread(here("Output","Fitness_Full",
                                                  "Fitness_MODEL_PC_FULL.csv")))


#R_M#

rate_PMR_dataframe <- rate_PMR_data(FULL_MODEL_PC_DT, 
                 c(5.5,6,10), 0.35, 48)


tmp_time_interest <- subset(rate_PMR_dataframe[[1]], 
                            rate_PMR_dataframe[[1]]$time <= time_of_interest)
  
end_points <- rate_PMR_dataframe[[2]]
  
  mer_GG <- ggplot(tmp_time_interest,
                   aes(x = time,
                       y = rate,
                       fill = as.factor(B_V),
                       color = as.factor(B_V),
                       group = as.factor(B_V))) +
    annotate('rect', xmin = 0, 
             xmax = time_of_interest, 
             ymin = 0,
             ymax = 1, 
             alpha = 0.3, 
             fill='grey')+
    geom_line(size=1.2) +
    geom_point(data = end_points,
               aes(x = endtime, y = PMR_END),
               color = 'black',
               size = 2.2, 
               shape = 21)+
    geom_hline(yintercept = 1,alpha = 0.2)+
    geom_segment(aes(x= 6.5, xend = 6.5,
                     y = 0, yend = 1),
                 color = "#FD7E6C",
                 size = 1.0,
                 linetype = 3)+
    geom_segment(aes(x = 14.1, xend = 14.1,
                     y = 0, yend = 1),
                 color = "black",
                 size = 1.0, 
                 linetype = 3)+
    geom_segment(aes(x= 15.8, 
                     xend = 15.8,
                     y = 0, 
                     yend = 1),
                 color = "#3DD8EE", 
                 size = 1.0, 
                 linetype = 3)+
    scale_color_manual(
      values = c("#3DD8EE",
                 "black",
                 "#FD7E6C"),
      name='Burst size')+
    scale_fill_manual(
      values = c("#3DD8EE",
                 "black",
                 "#FD7E6C"))+
    theme_bw()+
    xlab("Days post-infection")+
    ylab("Expected merozoites invading (uncommited) (R_M)")+
    scale_x_continuous(expand = c(0,0),
                       breaks = seq(0,35,5))+
    scale_y_continuous(expand = c(0,0))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          # axis.title.x = element_blank(),
          # axis.text.x = element_blank(),
          #  axis.ticks.x = element_blank(),
          legend.position = "none", 
          text = element_text(size = 14, color = 'black'))
  
  
  dailytp_GG <- 
    ggplot(tmp_time_interest, 
           aes(x = time, y = Daily_Trans_Prob, 
               fill = as.factor(B_V),
               color = as.factor(B_V),
               group = as.factor(B_V))) +
    geom_line(size = 1.2) +
    geom_segment(aes(x= 6.5+2, xend = 6.5+2,
                     y =0, yend = 0.9970716),
                 color="#FD7E6C", 
                 size=1.0,linetype = 3)+
    geom_segment(aes(x= 14.1 + 2, xend = 14.1 + 2,
                     y = 0, yend = 0.9312130),
                 color="black", 
                 size=1.0, linetype = 3)+
    geom_segment(aes(x = 15.8+2, xend = 15.8+2,
                     y = 0, yend = 0.7816508), 
                 color = "#3DD8EE", size = 1.0,
                 linetype = 3)+
    geom_point(data = end_points,
               aes(x = endtime, y = Daily_Trans_End),
               color = 'black',
               size = 2.2, 
               shape = 21)+
    scale_color_manual(values = 
                         c("#3DD8EE","black","#FD7E6C"),
                       name='Burst size')+
    scale_fill_manual(values = c("#3DD8EE",
                                 "black",
                                 "#FD7E6C"))+
    theme_bw()+
    scale_x_continuous(expand=c(0,0),breaks=seq(0,35,5))+
    scale_y_continuous(expand=c(0,0),limits=c(0,1.05))+
    xlab("Days post-infection")+
    ylab("Daily \ntransmission potential")+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          #axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          legend.position = "none",
          text = element_text(size = 14, color = 'black'))
  
  
  
  ###CUMULATIVE TRANSMISSION POTENTIAL
  cumcp_GG <-
    ggplot(tmp_time_interest, 
           aes(x = time, y = Cum_Trans_Potential, 
               fill = as.factor(B_V),
               color = as.factor(B_V),
               group = as.factor(B_V))) +
    geom_line(size=1.2)+
    geom_segment(aes(x= 6.5+2, xend = 6.5+2,
                     y =0, yend =36),color="#FD7E6C",size=1.0,linetype=3)+
    
    geom_segment(aes(x= 14.1+2, xend = 14.1+2,
                     y =0, yend =36),color="black",size=1.0,linetype=3)+
    
    geom_segment(aes(x= 15.8+2, xend = 15.8+2,
                     y =0, yend =36),color="#3DD8EE",size=1.0,linetype=3)+
    geom_point(data = end_points,
               aes(x = endtime, y = Cum_Trans_End),
               color = 'black',
               size = 2.2, 
               shape = 21)+
    scale_color_manual(
      values=c("#3DD8EE","black","#FD7E6C"),
      name='Burst size')+
    scale_fill_manual(values = c("#3DD8EE",
                                 "black",
                                 "#FD7E6C"))+
    theme_bw()+
    scale_x_continuous(expand=c(0,0),breaks=seq(0,35,5))+
    scale_y_continuous(expand=c(0,0))+
    theme_bw()+
    xlab("Days post-infection")+
    ylab("Cumulative \ntransmission potential")+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none",
          text = element_text(size = 14, color = 'black'))
  
  mer/dailytp/cumcp + plot_layout(guides='collect')
  



ggsave(file = here("Figures","mer_dal_cumcp_3.pdf"),width = 8, height=9,units='in')

