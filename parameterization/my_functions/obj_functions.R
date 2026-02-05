# Define the RMSE function
rmse <- function(observed, predicted) {
  # Check if the lengths of the observed and predicted vectors are the same
  if (length(observed) != length(predicted)) {
    stop("The lengths of observed and predicted values must be the same.")
  }
  
  # Calculate the squared differences
  squared_diff <- (observed - predicted)^2
  
  # Calculate the mean of the squared differences
  mean_squared_diff <- mean(squared_diff)
  
  # Calculate the square root of the mean squared differences (RMSE)
  rmse_value <- sqrt(mean_squared_diff)
  
  return(rmse_value)
}

#Vcmax25,Jmax25,TPU25 are not used for ePhoto's run
obj_Q10<-function(Q10s,curve_type,alpha1,alpha2,Vcmax25,
                  Jmax25,Rd25,TPU25,
                  Tleaf_all,Qin_all,Ci_all,curve_identifier,An_obs,targets_all,targets)
#Because each replicate has slightly different Ci (depsite being close),
#here we run ePhoto for each replicate's condition first.
#then the RMSE will be calculated by using their averages
#(Why not calculate RMSE for each replicate? I think it might lead to over-fitting)
#often time, each replicate's Ci sequence is close to each other
#Since RMSE is applied to averages, why not averaging the OBS and then run ePhoto?
#See obj_Q10_avg for this implementation
{
  avg_rmse =  TRUE#whether to calcualte RMSE based on avg of replicates
  if(curve_type==1){
    An_ePhoto   = NA* 1:length(Tleaf_all)
    for (i in 1:length(Tleaf_all)){ 
      PAR = Qin_all[i]
      Tleaf = Tleaf_all[i]
      Ci  = Ci_all[i]
      
      Q10s_new = targets_all
      Q10s_new[names(targets_all) %in% names(targets)] = Q10s
      # Arguments for the executable
      args <- c(as.numeric(alpha1), as.numeric(alpha2), PAR, Tleaf,Ci,Q10s_new)
      
      # Run the executable and capture the output
      output <- system2(executable_path, args = args, stdout = TRUE, stderr = TRUE)
      numbers <- as.numeric(unlist(strsplit(output, ",")))
      ephoto =  numbers[1]

      Rd = Rd25 * arrhenius_exponential(18.72, 46.39e3, Tleaf+273.15)
      An_ePhoto[i] = ephoto - Rd
    }
    An_ePhoto = as.numeric(An_ePhoto)
    df = data.frame(Ci=Ci_all,An_ePhoto,An_obs,curve_identifier)
    #calculate RMSE by averages by curve_id
    if(avg_rmse){
        unique_identifier = unique(df$curve_identifier)
        xx = 0
        An_ePhoto  = c()
        An_obs     = c()
        for (id in unique_identifier){
          tmp = df[df$curve_identifier == id,c("Ci","An_ePhoto","An_obs")]
          xx = xx + tmp
          An_ePhoto = cbind(An_ePhoto,tmp$An_ePhoto)
          An_obs    = cbind(An_obs   ,tmp$An_obs)
        }
        df_avg = xx/length(unique_identifier)
        #select only points with Ci larger than some threshold
        df_avg = df_avg[df_avg$Ci>100,] 
        err = rmse(df_avg$An_obs,df_avg$An_ePhoto) 
    }else{
        #check the main script, the Ci may have been sub setted there
        err = rmse(df$An_obs,df$An_ePhoto) 
    }
    return(err)

  }else if(curve_type==2){
    stop('under construction...')
  }
}

obj_Q10_avg<-function(Q10s,curve_type,alpha1,alpha2,Vcmax25,
                  Jmax25,Rd25,TPU25,
                  Tleaf_all,Qin_all,Ci_all,curve_identifier,An_obs,targets_all,targets)
{
  if(curve_type==1){
#first, we average the OBS data
    df = data.frame(Ci=Ci_all,Tleaf=Tleaf_all, Qin=Qin_all,An_obs,curve_identifier)
    #average by curve_id
    unique_identifier = unique(df$curve_identifier)
    xx = 0
    total_id = 1
    for (id in unique_identifier){
      if(id=="20210805-ripe3-4") next #this replicate's Ci seq is strange for both soy and tobacco 
      tmp = df[df$curve_identifier == id,c("Ci","Tleaf","Qin","An_obs")]
      xx = xx + tmp
      total_id = total_id + 1
    }
    df_avg = xx/total_id
    #select only points with Ci larger than some threshold
    df_avg = df_avg[df_avg$Ci>100,] 
#second, we run ePhoto for the avarges
    An_ePhoto   = NA* 1:dim(df_avg)[1]
    for (i in 1:dim(df_avg)[1]){ 
      PAR   = df_avg$Qin[i]
      Tleaf = df_avg$Tleaf[i]
      Ci    = df_avg$Ci[i]
      
      Q10s_new = targets_all
      Q10s_new[names(targets_all) %in% names(targets)] = Q10s
      # Arguments for the executable
      args <- c(as.numeric(alpha1), as.numeric(alpha2), PAR, Tleaf,Ci,Q10s_new)
      
      # Run the executable and capture the output
      ephoto <- system2(executable_path, args = args, stdout = TRUE, stderr = TRUE)
      if(length(ephoto)!=1) {
        # print(ephoto)
        ephoto = tail(ephoto,1)
      }
      ephoto = as.numeric(ephoto)
      Rd = Rd25 * arrhenius_exponential(18.72, 46.39e3, Tleaf+273.15)
      An_ePhoto[i] = ephoto - Rd
    }
    An_ePhoto = as.numeric(An_ePhoto)
    err = rmse(df_avg$An_obs,An_ePhoto) 
    return(err)

  }else if(curve_type==2){
    stop('under construction...')
  }
}

An_Q10<-function(Q10s,curve_type,alpha1,alpha2,Rd25,Tleaf,PAR,Ci)
{
  if(curve_type==1){
      # Arguments for the executable
      args <- c(as.numeric(alpha1), as.numeric(alpha2), PAR, Tleaf,Ci,Q10s)
      
      # Run the executable and capture the output
      ephoto <- system2(executable_path, args = args, stdout = TRUE, stderr = TRUE)
      if(length(ephoto)!=1) {
        # print(ephoto)
        ephoto = tail(ephoto,1)
      }
      ephoto = as.numeric(ephoto)
      Rd = Rd25 * arrhenius_exponential(18.72, 46.39e3, Tleaf+273.15)
      An_ePhoto = ephoto - Rd
      An_ePhoto = as.numeric(An_ePhoto)
      return(An_ePhoto)
  }else if(curve_type==2){
    stop('under construction...')
  }
}

# Define the LASSO regularized objective function with target values
lasso_obj_function <- function(params, constant_list, lambda) {
  # Calculate the original objective
  original_value <- obj_Q10(params,
                            constant_list$curve_type,
                            constant_list$alpha1,
                            constant_list$alpha2,
                            constant_list$Vcmax25,
                            constant_list$Jmax25,
                            constant_list$Rd25,
                            constant_list$TPU25,
                            constant_list$Tleaf_all,
                            constant_list$Qin_all,
                            constant_list$Ci_all,
                            constant_list$curve_identifier,
                            constant_list$An_obs,
                            constant_list$targets_all,
                            constant_list$targets
                            )
  
  # Calculate the modified LASSO penalty
  lasso_penalty <- lambda * sum(abs(params - constant_list$targets))
  
  print(params)
  
  # Return the combined objective function value
  original_value + lasso_penalty
}

# Define the LASSO regularized objective function with target values
l2_obj_function <- function(params, constant_list, lambda) {
  # Calculate the original objective
  original_value <- obj_Q10(params,
                            constant_list$curve_type,
                            constant_list$alpha1,
                            constant_list$alpha2,
                            constant_list$Vcmax25,
                            constant_list$Jmax25,
                            constant_list$Rd25,
                            constant_list$TPU25,
                            constant_list$Tleaf_all,
                            constant_list$Qin_all,
                            constant_list$Ci_all,
                            constant_list$curve_identifier,
                            constant_list$An_obs,
                            constant_list$targets_all,
                            constant_list$targets
                           )
  
  # Calculate the modified L2 penalty
  l2_penalty <- lambda * sum((params - constant_list$targets)^2)
  
  # Return the combined objective function value
  original_value + l2_penalty
}

An_obj_function <- function(params, constant_list) {
  # Calculate the original objective
  model_An <- An_Q10(params,
                            constant_list$curve_type,
                            constant_list$alpha1,
                            constant_list$alpha2,
                            constant_list$Rd25,
                            constant_list$Tleaf,
                            constant_list$PAR,
                            constant_list$Ci
  )
  
  # Return the combined objective function value
  model_An
}
