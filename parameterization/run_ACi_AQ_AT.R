library(reshape2)
library(ggplot2)
library(BioCro)
rm(list=ls())
source("../helper_functions/colorblind_palette.R")
source("my_functions/biocro_FvCB.R")
use_default_Q10s = TRUE 

my_palette = cb8()
keyword = "ld11"
Tgrowth = 24 
new_gamma_star = FALSE 
new_Vcmax_norm = FALSE
new_Jmax_norm  = FALSE 

if(!use_default_Q10s){
  Q10s = read.csv(paste0('outputs/ePhotosynthesis_optQ10_ACi_',keyword,'.csv'))
  Q10s = Q10s[-1]
  print(Q10s)
}

prefix = c("ACi","AQ","AT")
PAR = 2000 
executable_path <- "./myephoto.exe" 
##get alpha1 & alpha2
para_output = read.csv(paste0('outputs/ePhotosynthesis_optimal_alpha1_alpha2_',keyword,'.csv'))
alpha1 = para_output$alpha1
alpha2 = para_output$alpha2
total_P_sf = para_output$total_P_sf
##get FvCB parameters
Vcmax25      <- para_output$vcmax25 
Jmax25       <- para_output$jmax25 
Rd25         <- para_output$Rd25 
TPU25        <- para_output$TPU25
print(para_output)
if("alpha_old" %in% colnames(para_output)){
  alpha_old  <- para_output$alpha_old
}else{
  alpha_old  <- 0 
}

curve_options = 1:3 #1: A-Ci;2: A-Q; 3: A-T 
for (curve_option in curve_options){ 
  #call ephoto c++
  if(use_default_Q10s){
    args = c(alpha1,alpha2,total_P_sf,PAR,curve_option)
    output_figure_name = paste0("figs_25C/",prefix[curve_option],"_Q",PAR,"_",keyword,"_CTL.pdf")
  }else{
    args = c(alpha1,alpha2,total_P_sf,PAR,curve_option,Q10s)
    output_figure_name = paste0("figs_25C/",prefix[curve_option],"_Q",PAR,"_",keyword,"_optQ10.pdf")
  }
  args = as.numeric(args)
  print(length(args))
  system2(executable_path, args = args)
  #read in ephoto results
  ephoto = read.csv("output.data",header=FALSE)
  colnames(ephoto) = c("PAR","Tleaf","Ci","An")
  
  #run Farquhar
  An_farquhar = NA*(1:dim(ephoto)[1])
  
  for (i in 1:dim(ephoto)[1]){
    #For Farquhar, we need the total Q, as there's a calculation of
    #absorption inside the Farquhar function
    output_farquhar  = BioCro_FvCB(ephoto$PAR[i],ephoto$Tleaf[i], ephoto$Ci[i], Vcmax25, Jmax25, Rd25, TPU25,alpha_old,Tgrowth,new_gamma_star,new_Vcmax_norm,new_Jmax_norm)
    An_farquhar[i]   = output_farquhar$An
    #get ephoto's An
    Rd = Rd25 * arrhenius_exponential(18.72, 46.39e3, ephoto$Tleaf[i]+273.15)
    ephoto$An[i] = ephoto$An[i] - Rd
  }
  
  df_for_plot = cbind(ephoto,An_FvCB=An_farquhar)
  colnames(df_for_plot)[colnames(df_for_plot)=="An"] = "An_ePhoto"
  df_for_plot = melt(df_for_plot, id.vars=c("PAR","Tleaf","Ci"))
  
  if(curve_option==1){
    myplot<-ggplot(data=df_for_plot,aes(x=Ci, y=value,colour=variable)) +
         geom_line()  +  geom_point()+
      labs(x = bquote(Ci~(ppm)),
           y = bquote(An~(mu*mol ~ m^-2 ~ s^-1))
      ) +
      theme_bw() +
      theme(text = element_text(size = 16),
            legend.position = c(0.8, 0.2),
            legend.title = element_blank(),
            plot.margin = margin(0.5, 1, 0.5, 0.5, "cm") #(top, right, bottom, left)  
           )
  }else if(curve_option==2){
    myplot<-ggplot(data=df_for_plot,aes(x=PAR, y=value,colour=variable)) +
         geom_line()  +  geom_point()+
      labs(x = bquote(Q~(mu*mol ~ m^-2 ~ s^-1)),
           y = bquote(An~(mu*mol ~ m^-2 ~ s^-1))
      ) +
      theme_bw() +
      theme(text = element_text(size = 16),
            legend.position = c(0.8, 0.2),
            legend.title = element_blank(),
            plot.margin = margin(0.5, 1, 0.5, 0.5, "cm") #(top, right, bottom, left)  
           )
  }else{
    myplot<-ggplot(data=df_for_plot,aes(x=Tleaf, y=value,colour=variable)) +
         geom_line()  +  geom_point()+
      labs(x = bquote(T~(degree*C)),
           y = bquote(An~(mu*mol ~ m^-2 ~ s^-1))
      ) +
      scale_color_manual(values=my_palette[c(3,2)])+
      theme_bw() +
      theme(text = element_text(size = 16),
            legend.position = c(0.8, 0.2),
            legend.title = element_blank(),
            plot.margin = margin(0.5, 1, 0.5, 0.5, "cm") #(top, right, bottom, left)  
           )
  }
  # Save ggplot to PDF
  ggsave(output_figure_name, plot = myplot, width = 6, height = 6)
}
