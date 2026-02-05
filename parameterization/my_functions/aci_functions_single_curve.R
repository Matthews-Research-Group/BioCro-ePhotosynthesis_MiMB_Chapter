#This converts a simple dataframe of A-Ci data to
#the exdf format that is needed to run PhotoGEA
get_vcmax_jmax <- function(A_Ci_df,my_fit_option,Tgrowth,cj_crossover_max,new_gamma_star,new_Vcmax_norm,new_Jmax_norm)
{
  aci_data      = A_Ci_df
  aci_data_exdf = exdf(aci_data)

  aci_data_exdf$units$Tleaf = "degrees C"
  aci_data_exdf$units$Ci = "micromol mol^(-1)"
  aci_data_exdf$units$A = "micromol m^(-2) s^(-1)"
  
  if("Ca" %in% names(A_Ci_df)){
    aci_data_exdf$units$Ca = "micromol mol^(-1)"
  }else{
    # #add a bogus constant Ca
    aci_data_exdf <- set_variable(
      aci_data_exdf,
      'Ca',
      value = 1,
      units = 'micromol mol^(-1)'
    )
  }
  
  #assume 1 atmos pressure
  aci_data_exdf <- set_variable(
    aci_data_exdf,
    'total_pressure',
    value = 1,
    units = 'bar'
  )
  aci_data_exdf <- set_variable(
    aci_data_exdf,
    'Oxygen',
    value = 21,
    units = 'percent'
  )
  #calculate Kc, J_norm, Ko,Vcmax
  my_c3_arrehenius = c3_temperature_param_bernacchi
  if(new_gamma_star){
    my_c3_arrehenius$Gamma_star_norm$c  = 17.9  #Tobacco:https://doi.org/10.1111/pbi.13750
    my_c3_arrehenius$Gamma_star_norm$Ea = 34.8
  }
  if(new_Jmax_norm){
    my_c3_arrehenius$J_norm$c  = 20.99  #Tobacco:https://doi.org/10.1111/pbi.13750
    my_c3_arrehenius$J_norm$Ea = 52
  }
  aci_data_exdf <- calculate_temperature_response(aci_data_exdf, my_c3_arrehenius,tleaf_column_name = 'Tleaf')

  #overwrite arrehenius with a new vcmax temperature response
  if(new_Vcmax_norm){
    aci_data_exdf[, 'Vcmax_norm'] <- Vcmax_multiplier(aci_data_exdf[, 'Tleaf']+273.15,Tgrowth)
  }
  
  # Calculate Cc
  aci_data_exdf <- set_variable(
    aci_data_exdf,
    'gmc_at_25',
    'mol m^(-2) s^(-1) bar^(-1)',
    'process_soybean_aci',
    Inf
  )
  aci_data_exdf <- apply_gm(aci_data_exdf)
  
#  # We can fit just one curve from the data set, although it is rare to do this

  c3_aci_results <- fit_c3_aci(
    aci_data_exdf,
    a_column_name = 'A',
    Ca_atmospheric = 420,
    fit_options = my_fit_option, 
    cj_crossover_max = cj_crossover_max,
    Wj_coef_C = ELECTRONS_PER_CARBOXYLATION,
    Wj_coef_Gamma_star = ELECTRONS_PER_OXYGENATION * 2
  )
  
  #here, no average is needed
  c3_aci_averages <- c3_aci_results$parameters
  
  # Compile table of parameter values to use with BioCro
  fvcb_parameters <- list(
    cultivar = 'ld11',
    electrons_per_carboxylation = ELECTRONS_PER_CARBOXYLATION,
    electrons_per_oxygenation = ELECTRONS_PER_OXYGENATION,
    Vcmax = c3_aci_averages[, 'Vcmax_at_25'],
    J     = c3_aci_averages[, 'J_at_25'],
    Rd    = c3_aci_averages[, 'RL_at_25'],
    TPU   = c3_aci_averages[, 'Tp_at_25'],
    alpha_old = c3_aci_averages[, 'alpha_old']
  )
  
  # One complication is that PhotoGEA returns a value for J, but not Jmax. In
  # BioCro, values of PhiPSII, Q, and Jmax are used to determine J. Here, we will
  # solve those equations for Jmax.
  leaf_reflectance   = 0.1
  leaf_transmittance = 0.05
  fvcb_parameters$Jmax <- get_jmax25(
    soybean$parameters$theta_0,            # dimensionless
    soybean$parameters$beta_PSII,          # dimensionless
    fvcb_parameters$J,                     # micromol / m^2 / s
    leaf_reflectance,                      # dimensionless
    leaf_transmittance,                    # dimensionless
    mean(aci_data[, 'Tleaf']),             # degrees C
    mean(aci_data_exdf[, 'PAR'])           # micromol / m^2 / s
  )


  fvcb_parameters$TPU25 <- c3_aci_averages[,'Tp_at_25']
  #rescale the TPU to 25 C
  # average_Tleaf = mean(aci_data[, 'Tleaf']) 
  # fvcb_parameters$TPU25 <- get_TPU_at_25(
  #        c3_aci_averages[,'Tp'],
  #        average_Tleaf
  #        )

  return(c(fvcb_parameters$Vcmax,fvcb_parameters$Jmax,fvcb_parameters$TPU25,fvcb_parameters$Rd,fvcb_parameters$alpha_old))
}

# dev.new()
# print(xyplot(
#   A + Ac + Aj + A_fit ~ Cc,
#   data = c3_aci_results$fits$main_data,
#   type = 'b',
#   pch = 16,
#   auto.key = list(space = 'right'),
#   grid = TRUE,
#   xlab = paste('Chloroplast CO2 concentration [', c3_aci_results$fits$units$Ci, ']'),
#   ylab = paste('Net CO2 assimilation rate [', c3_aci_results$fits$units$A, ']'),
#   par.settings = list(
#     superpose.line = list(col = multi_curve_colors()),
#     superpose.symbol = list(col = multi_curve_colors(), pch = 16)
#   )
# ))
