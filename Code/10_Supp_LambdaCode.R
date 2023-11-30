###############################
###Supplemental lambda plot###
##############################
library(here)

source(here("Code","Helper_Function_Code",
            "07_sensitivity_analysis_functions.R"))
       
B_V = c(15.5)
C_V = 0.76
change = c(277500, 370000, 462500)
B_V_C_V <- expand.grid(B_V = B_V, C_V = C_V, change = change)

Lambda_Supp_DF <-mcmapply(Simulator_Malaria_BC_Lambda,
         c(B_V_C_V$B_V),
         c(B_V_C_V$C_V),
         c(4385.96491),
         c(B_V_C_V$change),
         mc.cores = 5,
         SIMPLIFY = FALSE
)

for (k in seq(1, length(Lambda_Supp_DF ))){
  Lambda_Supp_DF[[k]]$change <- B_V_C_V [k,]$change
}




Lambda_Supp_DF  <- do.call(rbind,Lambda_Supp_DF )


ggplot(subset(Lambda_Supp_DF,
              Lambda_Supp_DF$time < 50),
       aes(x = time, y = R, color = as.factor(change), 
           group = as.factor(change)))+
  geom_line(linewidth= 0.89)+
  xlab("Days post-infection")+
  ylab("RBC abundance (log10)") + 
  scale_color_manual(name = "RBC replenishment",
                       values = c('#1199ee','black','#ef2366'))+
  scale_y_log10(
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 15),
        axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 13)) 


ggsave(here("Figures","Raw", "lambda_RBC_Supp.pdf"),
       units = 'in', height = 4, width = 5)  
  
