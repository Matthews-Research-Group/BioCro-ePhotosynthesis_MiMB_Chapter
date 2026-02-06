library(reshape2)
library(ggplot2)
library(BioCro)
library(PhotoGEA)
rm(list=ls())
source("my_functions/biocro_FvCB.R")
plot_type = 1 #1:ACi; 2: AQ
prefix = c("ACi","AQ")

use_default_Q10s = FALSE 
Tgrowth = 24 
new_gamma_star =  FALSE 
new_Vcmax_norm =  FALSE
new_Jmax_norm  =  FALSE 

if(use_default_Q10s){
  exp_type = "CTL"
}else{
  exp_type = "Q10opt"
}
keyword = "ld11"
#output figure name. Make sure to change before running
output_figure_name = paste0("figs/",prefix[plot_type],"_with_OBS_",exp_type,"_",keyword,".pdf")

executable_path <- "./myephoto_single.exe" 

alpha1_alpha2 = read.csv(paste0('outputs/ePhotosynthesis_optimal_alpha1_alpha2_',keyword,'.csv'))
alpha1 = alpha1_alpha2[2]
alpha2 = alpha1_alpha2[3]

if(!use_default_Q10s){
  Q10s = read.csv(paste0('outputs/ePhotosynthesis_optQ10_ACi_',keyword,'.csv'))
  Q10s = Q10s[-1]
}

Vcmax25      <- alpha1_alpha2$vcmax25
Jmax25       <- alpha1_alpha2$jmax25
Rd25         <- alpha1_alpha2$Rd25
TPU25        <- alpha1_alpha2$TPU25
if('alpha_old' %in% alpha1_alpha2$parameter){
  alpha_old   <- alpha1_alpha2$alpha_old
}else{
  alpha_old   <- 0 
}
print(alpha1_alpha2)

if(plot_type==1){
  licor_data <- read.csv.exdf(paste0("my_functions/",keyword,"_aci.csv")) 
  
  Tleaf_all  <- licor_data[, 'TleafCnd'] 
  Qin_all    <- licor_data[, 'Qin']
  Ci_all     <- licor_data[, 'Ci']
  
  date = substr(licor_data[, 'date'], 1, 8)
  curve_identifier = paste(date,licor_data[, 'instrument'],licor_data[, 'plot'],sep='-')
  An_obs = licor_data[,'A']
  
  An_FvCB     = NA* 1:length(Tleaf_all)
  An_ePhoto   = NA* 1:length(Tleaf_all)
   
  for (i in 1:length(Tleaf_all)){ 
    PAR   = Qin_all[i]
    Tleaf = Tleaf_all[i]
    Ci    = Ci_all[i]
    #For Farquhar, we need the total Q, as there's a calculation of
    #absorption inside the Farquhar function
    output_farquhar  = BioCro_FvCB(PAR,Tleaf, Ci, Vcmax25, Jmax25, Rd25, TPU25,alpha_old,Tgrowth,new_gamma_star,new_Vcmax_norm,new_Jmax_norm)
    An_FvCB[i]   = output_farquhar$An
  
    #call ephoto c++
    if(use_default_Q10s){
      args = c(alpha1,alpha2,PAR,Tleaf, Ci)
    }else{
      args = c(alpha1,alpha2,PAR,Tleaf, Ci,Q10s)
    }
    args = as.numeric(args)
    output  <-system2(executable_path, args = args, stdout = TRUE, stderr = TRUE)
    if((length(output))>1) output = output[3]
    numbers <- as.numeric(unlist(strsplit(output, ",")))
    ephoto =  numbers[1]
    Rd = Rd25 * arrhenius_exponential(18.72, 46.39e3, Tleaf+273.15)
    An_ePhoto[i] = ephoto - Rd
  }
  An_ePhoto = as.numeric(An_ePhoto)
  df = data.frame(Tleaf=Tleaf_all,PAR=Qin_all,Ci=Ci_all,An_ePhoto,An_FvCB,An_obs,curve_identifier)

  #plot average line
  unique_identifier = unique(df$curve_identifier)
  xx = 0
  An_ePhoto  = c()
  An_FvCB    = c()
  An_obs     = c()
  total_id = 0
  for (id in unique_identifier){
    tmp = df[df$curve_identifier == id,1:6]
    xx = xx + tmp
    An_ePhoto = cbind(An_ePhoto,tmp$An_ePhoto)
    An_FvCB   = cbind(An_FvCB  ,tmp$An_FvCB)
    An_obs    = cbind(An_obs   ,tmp$An_obs)
    total_id = total_id + 1
  }
  sd_ePhoto   = apply(An_ePhoto,1,sd)
  sd_FvCB     = apply(An_FvCB,1,sd)
  sd_obs      = apply(An_obs,1,sd)
  df_avg      = xx/total_id
  df_avg2     = df_avg[,c('Ci','An_ePhoto','An_FvCB','An_obs')]
  df_avg2_sub     = df_avg2[df_avg2$Ci>100,]

  df_melt = reshape2::melt(df_avg2,id.vars = c("Ci"),variable.name = "variable", value.name = "value")
  df_melt = cbind(df_melt,sd = c(sd_ePhoto,sd_FvCB,sd_obs))
  myplot<-ggplot(df_melt, aes(x = Ci, y = value, color = variable)) +
    geom_line() + geom_point(size=3) +
    scale_color_manual(values=c("green","blue", "black"))+
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(0.05)) +
    labs(x = bquote(Ci~(ppm)),
         y = bquote(An~(mu*mol ~ m^-2 ~ s^-1))
    ) +
    theme_bw() +
    theme(text = element_text(size = 20),legend.position = c(0.8, 0.2))
# Save ggplot to PDF
  ggsave(output_figure_name, plot = myplot, width = 6, height = 6)

}else if(plot_type==2){
  licor_data <- read.csv.exdf("for_calibration_with_OBS/my_scripts/ld11_bb.csv") #this comes from "get_ld11_bb.R"

  Tleaf_all  <- licor_data[, 'TleafCnd'] 
  Qin_all    <- licor_data[, 'Qin']
  Ci_all     <- licor_data[, 'Ci']

  date = substr(licor_data[, 'date'], 1, 8)
  
  curve_identifier = paste(date,licor_data[, 'instrument'],licor_data[, 'plot'],licor_data[, 'replicate'],sep='-')
  
  An_obs = licor_data[,'A']

  An_FvCB     = NA* 1:length(Tleaf_all)
  An_ePhoto   = NA* 1:length(Tleaf_all)
   
  for (i in 1:length(Tleaf_all)){ 
    PAR = Qin_all[i]
    Tleaf = Tleaf_all[i]
    Ci  = Ci_all[i]
    #For Farquhar, we need the total Q, as there's a calculation of
    #absorption inside the Farquhar function
    output_farquhar  = BioCro_FvCB(PAR,Tleaf, Ci, Vcmax25, Jmax25, Rd25, TPU25)
    An_FvCB[i]   = output_farquhar$An
  
    #call ephoto c++
    args = c(alpha1,alpha2,PAR,Tleaf, Ci,Q10s)
    args = as.numeric(args)
    ephoto<-system2(executable_path, args = args, stdout = TRUE, stderr = TRUE)
    if(length(ephoto)!=1) {
      print(ephoto)
      ephoto = tail(ephoto,1)
    }
    ephoto =  as.numeric(ephoto)
    Rd = Rd25 * arrhenius_exponential(18.72, 46.39e3, Tleaf+273.15)
    An_ePhoto[i] = ephoto - Rd
  }
  An_ePhoto = as.numeric(An_ePhoto)

  df = data.frame(Qin=Qin_all,An_ePhoto,An_FvCB,An_obs,curve_identifier)

  #plot average line
  unique_identifier = unique(df$curve_identifier)
  xx = 0
  An_ePhoto  = c()
  An_FvCB    = c()
  An_obs     = c()
  for (id in unique_identifier){
    tmp = df[df$curve_identifier == id,1:4]
    xx = xx + tmp
    An_ePhoto = cbind(An_ePhoto,tmp$An_ePhoto)
    An_FvCB   = cbind(An_FvCB  ,tmp$An_FvCB)
    An_obs    = cbind(An_obs   ,tmp$An_obs)
  }
  sd_ePhoto   = apply(An_ePhoto,1,sd)
  sd_FvCB     = apply(An_FvCB,1,sd)
  sd_obs      = apply(An_obs,1,sd)
  df_avg = xx/length(unique_identifier)
  df_melt = reshape2::melt(df_avg,id.vars = c("Qin"),variable.name = "variable", value.name = "value")
  df_melt = cbind(df_melt,sd = c(sd_ePhoto,sd_FvCB,sd_obs))
  myplot<-ggplot(df_melt, aes(x = Qin, y = value, color = variable)) +
    geom_line() + geom_point(size=3) +
    scale_color_manual(values=c("green","blue", "black"))+
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(0.05)) +
    labs(x = bquote(Q~(mu*mol ~ m^-2 ~ s^-1)),
         y = bquote(An~(mu*mol ~ m^-2 ~ s^-1))
    ) +
    theme_bw() +
    theme(text = element_text(size = 20),legend.position = c(0.8, 0.2))
# Save ggplot to PDF
  ggsave(output_figure_name, plot = myplot, width = 6, height = 6)
}
