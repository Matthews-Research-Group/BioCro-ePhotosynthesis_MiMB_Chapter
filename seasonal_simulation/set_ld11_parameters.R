set_ld11_parameters<- function(year,soybean_parameters0){
#get FvCB parameters
  soybean_fvcb_parameters = read.csv('../parameterization/FvCB_parameters//ld11_aci_fit_parameters_avg.csv')
  # From Ball-Berry curves
  soybean_parameters0$b0 = 0.085           # estimated from average of 2021 and 2022 data 
  soybean_parameters0$b1 = 5.375           # estimated from average of 2021 and 2022 data 

  # From A-Ci curves
  soybean_parameters0$Vcmax_at_25 = soybean_fvcb_parameters$mean[soybean_fvcb_parameters$parameter=="Vcmax_at_25"]
  soybean_parameters0$Jmax_at_25  = soybean_fvcb_parameters$mean[soybean_fvcb_parameters$parameter=="Jmax_at_25"]
  soybean_parameters0$RL_at_25    = soybean_fvcb_parameters$mean[soybean_fvcb_parameters$parameter=="RL_at_25"]
  soybean_parameters0$Tp_at_25   = soybean_fvcb_parameters$mean[soybean_fvcb_parameters$parameter=="TPU_at_25"] 
  #iSp
  soybean_parameters0$iSp = 3

#get atm CO2
  CO2_series = get_CO2()
  #atmospheric CO2
  soybean_parameters0$Catm = CO2_series$co2[CO2_series$year==year] 

  return(soybean_parameters0)
}

get_CO2<-function(){
#co2 ppm data of SSP5-8.5 from 1970-2060. Yufeng
years = 1970:2060
co2_level =
        c(325.58,326.23,328.21,330.81,331.67,331.74,332.50,334.34,336.02,337.66,
          339.75,341.08,341.65,343.24,344.83,346.31,347.84,349.58,352.26,354.00,
          355.14,356.49,357.22,357.89,359.44,361.50,363.41,364.58,367.28,369.33,
          370.70,372.26,374.41,377.16,378.64,380.61,382.83,384.31,386.50,388.02,
          390.78,393.04,395.04,397.71,399.59,401.68,405.11,407.85,410.85,413.92,
          417.06,420.29,423.62,427.07,430.63,434.30,438.08,441.97,445.97,450.08,
          454.29,458.63,463.10,467.70,472.43,477.29,482.28,487.41,492.67,498.06,
          503.60,509.27,515.09,521.07,527.19,533.47,539.91,546.50,553.25,560.16,
          567.24,574.48,581.92,589.55,597.38,605.40,613.63,622.05,630.68,639.51,648.54)
co2_out = data.frame(year = years,co2=co2_level)
return(co2_out)
}
