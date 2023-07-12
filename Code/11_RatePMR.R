
library(here)

# This is the main simulation assuming nothing has changed
# about the initial red blood cell density or the replenishment rate.
# This is one of the longer code to source

### Packages to load
source(here("Code", "Helper_Function_Code", "Packages_Loader.R"))

### Main modeling code
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_CUT.cpp"))
sourceCpp(here("Code", "RCPP_Code", "rcpp_malaria_dynamics_UNCUT.cpp"))
source(here("Code", "Simulator_Code", "Simulator_Main.R"))



HighB_HighC <- Simulator_Malaria_BC_FULL(49,0.82) #Optimal for 50 days
HighB_HighC$probG <- PrI_PC(HighB_HighC$G)

HighB_HighC_SAMERC <- Simulator_Malaria_BC(10,0.118) #Optimal for 50 days
HighB_HighC_SAMERC$probG <- PrI_PC(HighB_HighC_SAMERC $G)




LowB_MedC <- Simulator_Malaria_BC_FULL(14.5,0.72) #optimal for 25 days
LowB_MedC $probG <- PrI_PC(LowB_MedC $G)

a<- ggplot(HighB_HighC, aes(x = M, y=R
                        ,color = time)) + 
  
  geom_path(size = 2,lineend = 'round')+
  scale_color_viridis(option = 'turbo',name = 'Time') + 
  geom_path(data = LowB_MedC, aes(x= M, y= R),size =2,lineend='round')+


  xlab("Merozoite abundance")+
  ylab("RBC abundance")+
  theme_bw()+
  theme(panel.grid = element_blank())


b<-ggplot(HighB_HighC, aes(x = G, y=R
                        ,color = time)) + 
  
  geom_path(size = 2,lineend = 'round')+
  scale_color_viridis(option = 'turbo',name = 'Time') + 
  geom_path(data = LowB_MedC, aes(x= G, y= R),size =2,lineend='round')+
  
  
  ylab("RBC abundance")+
  xlab("Gam abundance")+
  theme_bw()+
  theme(panel.grid = element_blank())


a+b+plot_layout(guides = 'collect')







GAM_TS <- ggplot(HighB_HighC, aes(x = time, y= log10(G+1))) + 
  geom_line()+
  geom_line(data = LowB_MedC, aes(x= time, y= log10(G+1)),col ='red')+
  geom_hline(yintercept = log10(1270))+
  xlab("Days post-infection")+
  ylab("Gametocyte abundance (log10)")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = c(5,25,50), linetype = 2)+
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14))



GAM_Prob_TS <- ggplot(HighB_HighC, aes(x = time, y=probG )) + 
  geom_line()+
  geom_line(data = LowB_MedC, aes(x= time, y= probG),col ='red')+
  geom_line(data =HighB_HighC_SAMERC, aes(x= time, y= probG),col ='blue')+
  
  xlab("Days post-infection")+
  ylab("Transmission probability")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = c(5,25,50), linetype = 2)+
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14))


RBC_TS <- ggplot(HighB_HighC, aes(x = time, y= log10(R+1))) + 
  geom_line()+
  geom_line(data = LowB_MedC, aes(x= time, y= log10(R+1)),col ='red')+
  geom_line(data = HighB_HighC_SAMERC, aes(x= time, y= log10(R+1)),col ='blue')+
  
  xlab("Days post-infection")+
  ylab("RBC abundance (log10)")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = c(5,25,50), linetype = 2)+
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14))



p =  2.5e-6

HighB_HighC_PMR <- data.frame(rate_PMR(
  HighB_HighC ,49,0.82,48))
       
LowB_MedC_PMR <- data.frame(rate_PMR(
  LowB_MedC  , 14.5,0.72,48))

HighB_HighC_SAMERC_PMR <- data.frame(rate_PMR(
  HighB_HighC_SAMERC ,10,0.118,48))



PMR<- ggplot(
  HighB_HighC_PMR, aes(x = time, y= rate)) + 
  geom_line(col = 'green')+
  geom_line(data= LowB_MedC_PMR , 
            aes(x= time, y= rate),col='black')+
  geom_line(data = 
              HighB_HighC_SAMERC_PMR,
            aes(x = time, y= rate),col = 'blue')+
  xlab("Days post-infection")+
  ylab("Total merozoites invading successfully")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = c(5,25,50), linetype = 2)+
  geom_hline(yintercept =1,linetype = 3, color = 'red')+
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14))

plot(HighB_HighC_PMR $rate, HighB_HighC$G)
points(LowB_MedC_PMR$rate, LowB_MedC$G, col='red')



ggsave(here("Figures", "Raw", "RBC_GAM_TS.pdf"), height = 4, width = 4,
       units = 'in')