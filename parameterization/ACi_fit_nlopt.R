#This is for fitting the A-Ci curves of ePhotosynthesis
#to match target Vcmax and Jmax
# Clear the workspace
rm(list=ls())
# Load required packages
# library(lattice)
library(nloptr)
library(PhotoGEA)
library(BioCro)
source('my_functions/aci_defaults.R') #get_jmax25 function
source('my_functions/biocro_FvCB.R')  #Vcmax_multiplier function
source('my_functions/aci_functions_single_curve.R')

keyword = "ld11"
FvCB_parameter_file = paste0("FvCB_parameters/",keyword,"_aci_fit_parameters_avg.csv")
aci_fit_results <- read.csv(FvCB_parameter_file)

Tgrowth = 24 
cj_crossover_max = NA #set to NA otherwise
new_gamma_star = FALSE #if using Amanda's Gamma_star
new_Vcmax_norm = FALSE
new_Jmax_norm  = FALSE 

#for PhotoGEA's outputs, we do a bit re-arrange of the results
if(TRUE){
  x <-  read.csv(FvCB_parameter_file)
  aci_fit_results =  matrix(x$mean,nrow=1)
  colnames(aci_fit_results) = x$parameter
  aci_fit_results = as.data.frame(aci_fit_results)
  if(length(colnames(aci_fit_results))==6){
    colnames(aci_fit_results) = c('vcmax25','j25','jmax25','rd25','tpu25','alpha_old')
  }else if(length(colnames(aci_fit_results))==5){
    colnames(aci_fit_results) = c('vcmax25','j25','jmax25','rd25','tpu25')
  }
}
vcmax_target <- aci_fit_results$vcmax25
jmax_target  <- aci_fit_results$jmax25
Rd_at_25     <- aci_fit_results$rd25 
tpu_target   <- aci_fit_results$tpu25 
print(aci_fit_results)

if("alpha_old" %in% colnames(aci_fit_results)){
  alpha_old  <- aci_fit_results$alpha_old
}else{
  alpha_old  <- 0
}

#If you used PhotoGEA to estimate the FvCB parameters, 
#make sure this option here is consistent with that
my_fit_option = list(RL_at_25  = Rd_at_25,
#Since ePhoto cannot account for a decline in A, i always set alpha_old to 0 
                     alpha_old = 0
#                     Tp   = tpu_target 
                     )

obj_func<-function(x){
  # Specify FvCB parameter values
  # To-do: it's better to pass these as function parameters instead of 
  # reading it in as something out of scope. Still works in R, but not good
  vcmax_target <- aci_fit_results$vcmax25
  jmax_target  <- aci_fit_results$jmax25
  Rd_at_25     <- aci_fit_results$rd25 
  tpu_target   <- aci_fit_results$tpu25 
#set high light to derive Vcmax and Jmax
  PAR = 1800
  system(paste("./myephoto.exe",x[1],x[2],x[3],PAR,1))
  ePhoto_result = read.table("output.data",header=FALSE,sep = ",")
  A_Ci_df = as.data.frame(ePhoto_result)
  colnames(A_Ci_df) = c("PAR","Tleaf","Ci","A")
  #ephoto's assimilation is the gross,which needs to substract Rd to get An
  Rd = Rd_at_25 * arrhenius_exponential(18.72, 46.39e3, A_Ci_df$Tleaf+273.15) 
  A_Ci_df$A = A_Ci_df$A - Rd 
  output = get_vcmax_jmax(A_Ci_df,my_fit_option,Tgrowth,cj_crossover_max,new_gamma_star,new_Vcmax_norm,new_Jmax_norm)
  vcmax25 = output[1]
  jmax25  = output[2]
  tpu25   = output[3]
  rd25    = output[4]
 
  #minimize vcmax and jmax distances to the targets 
  error = (vcmax25 - vcmax_target)^2 + (jmax25 - jmax_target)^2 #+ (tpu25 - tpu_target)^2
  if(is.na(error)) error=999
  return(error)
}

# These bounds are subject to changes. My strategy is to compare the 
# optimized alpha1&alpha2 values with the bounds.
# If either is very close to bounds, it suggests that there might be   
# a better solution if extending the bounds.
lbs = c(0.5,0.6,1.0) #alpha1, alpha2
ubs = c(1.0,2.5,1.0)
init_guess = c(0.7,1.0,1.0) 

opt_results = nloptr(x0 = init_guess, eval_f = obj_func, 
                     lb = lbs,
                     ub = ubs,
                     opts = list("algorithm" = "NLOPT_LN_SBPLX",
                                 "xtol_rel"  = 1.0e-4,
                                 "maxeval"   = 200,
				 "print_level"=3)
 		     )
alpha1_alpha2 = as.data.frame(opt_results$solution)
rownames(alpha1_alpha2) = c('alpha1','alpha2','total_P_sf')
alpha1_alpha2 = t(alpha1_alpha2)
output_data = data.frame(alpha1_alpha2,vcmax25=vcmax_target,jmax25=jmax_target,Rd25=Rd_at_25,TPU25=tpu_target,alpha_old=alpha_old)
write.csv(output_data,paste0("outputs/ePhotosynthesis_optimal_alpha1_alpha2_",keyword,".csv"))
