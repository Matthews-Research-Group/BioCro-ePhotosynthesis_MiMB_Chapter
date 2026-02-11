library(ggplot2)
library(gridExtra)
rm(list=ls())
source('merge_rds_results.R')
out_keys   = c("ePhoto","FvCB")
use_ePhoto = FALSE
years      = 2014
folder_path = paste0("Results_Bernacchi/year",years,"site1")

if(use_ePhoto){
  exp_id  = 1 #this will use the specific enzyme opt experiment's vmax scaling factors
}else{
  exp_id  = 0 #FvCB does not need this. Always set to 0
}

result_list = list()
for (i in 1:length(out_keys)){
  out_key = out_keys[i]
  site_id      = 1
  result_list[[i]] <- merge_results(site_id,years,exp_id,out_key,folder_path)
}

time_stamp = result_list[[1]][[1]]$doy + result_list[[1]][[1]]$hour/24
ephoto = result_list[[1]][[1]]
fvcb   = result_list[[2]][[1]]
ephoto$time = time_stamp
fvcb$time   = time_stamp
ephoto$cum_sunlit_A = cumsum(ephoto$sunlit_Assim_layer_0)
ephoto$cum_shaded_A = cumsum(ephoto$shaded_Assim_layer_0)
ephoto$cum_canopy_A = cumsum(ephoto$canopy_assimilation_rate)
fvcb$cum_sunlit_A = cumsum(fvcb$sunlit_Assim_layer_0)
fvcb$cum_shaded_A = cumsum(fvcb$shaded_Assim_layer_0)
fvcb$cum_canopy_A = cumsum(fvcb$canopy_assimilation_rate)
vars2plot = c('sunlit_Assim_layer_0','shaded_Assim_layer_0','cum_sunlit_A','cum_shaded_A',
              'lai','Grain','canopy_assimilation_rate','cum_canopy_A')
source('plot_functions.R')
timeseries(ephoto,fvcb,vars2plot)

df2plot = data.frame(year=years,models=out_keys,TotalDays=NA,Seed=NA,CanopyA=NA,Tmean=NA,Qmean=NA,WSmean=NA)

for (i in 1:length(out_keys)){
  tmp = result_list[[i]][[1]]
  cumA = cumsum(tmp$canopy_assimilation_rate)
  df2plot$Seed[i]      = tail(tmp$Grain,1)
  df2plot$CanopyA[i]   = tail(cumA,1)
  df2plot$TotalDays[i] = dim(tmp)[1]/24
  df2plot$Tmean[i]     = mean(tmp$temp)
  df2plot$Qmean[i]     = mean(tmp$par_incident_direct + tmp$par_incident_diffuse)
  df2plot$WSmean[i]     = mean(tmp$StomataWS)
} 


if(TRUE){
  # ggplot(df2plot, aes(x = year, y = Seed, color = models)) +
  ggplot(df2plot, aes(x = models, y = Seed)) +
    geom_bar(stat = "identity", fill = "skyblue") +      # Use lines for each group
    labs(x = "models", y = "Seed yield (t/ha)") +
    theme_minimal(base_size = 20)
  ggsave(paste0("figs/compare_models_","_newWater.png"),plot = last_plot(),
         width = 12, height = 8, units = "in", dpi = 300)
  # Arrange the plots in a grid and save to a PNG
  # png(paste0("barplot_site",site_id,"_",out_key,".png"), width = 12, height = 8, units = "in", res = 300) # Open PNG device
  # grid.arrange(p1, p2, p3, p4,ncol = 2) # Arrange the plots in a 2-column grid
  dev.off() # Close PNG device
}
