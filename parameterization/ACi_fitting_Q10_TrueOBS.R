library(BioCro)
library(PhotoGEA)
library(nloptr)
source("my_functions/biocro_FvCB.R") # for Arrehnius equation
source("my_functions/obj_functions.R")
curve_type = 1 #1:ACi; 2: AQ
prefix = c("ACi","AQ")
outfile_prefix = prefix[curve_type]

keyword = "ld11" #check if the species changed Pi value in ePhoto's exe

all_enzymes = c( 'Q10_1',  'Q10_2', 'Q10_3',
                 'Q10_5',  'Q10_6', 'Q10_7',
                 'Q10_8',  'Q10_9', 'Q10_10',
                 'Q10_13', 'Q10_23')
#the default values of these enzymes
targets_all = c(1.93,rep(2.0,length(all_enzymes)-1))
names(targets_all) = all_enzymes
 
enzymes_to_optimize = c( 'Q10_1',  'Q10_2', 'Q10_3',
                         'Q10_5',  'Q10_6', 'Q10_7',
                         'Q10_8',  'Q10_9', 'Q10_10',
                         'Q10_13', 'Q10_23')
#enzymes_to_optimize = c( 'Q10_1','Q10_5','Q10_9','Q10_23') 
number_of_enzymes = length(enzymes_to_optimize) #total number of Q10 to be optimized
#the default values of these enzymes
targets = targets_all[all_enzymes%in%enzymes_to_optimize]
names(targets) = enzymes_to_optimize

executable_path <- "./myephoto_single.exe"

alpha1_alpha2 = read.csv(paste0('outputs/ePhotosynthesis_optimal_alpha1_alpha2_',keyword,'.csv'))
print(alpha1_alpha2)
alpha1 = alpha1_alpha2[2]
alpha2 = alpha1_alpha2[3]

Vcmax25      <- alpha1_alpha2$vcmax25
Jmax25       <- alpha1_alpha2$jmax25
Rd25         <- alpha1_alpha2$Rd25
TPU25        <- alpha1_alpha2$TPU25

if(curve_type==1){
    licor_data <- read.csv.exdf(paste0("my_functions/",keyword,"_aci.csv")) 
    
    Tleaf_all  <- licor_data[, 'TleafCnd'] 
    Qin_all    <- licor_data[, 'Qin']
    Ci_all     <- licor_data[, 'Ci']
    
    date = substr(licor_data[, 'date'], 1, 8)
    curve_identifier = paste(date,licor_data[, 'instrument'],licor_data[, 'plot'],sep='-')
    An_obs = licor_data[,'A']

    #remove invalid replicate with Jmax of NA
    all_replicate_Vcmax_Jmax = read.csv.exdf(paste0("my_functions/",keyword,"_aci_fit_parameters.csv"))
    all_curve_id  = paste(all_replicate_Vcmax_Jmax[,'day'],all_replicate_Vcmax_Jmax[,'instrument'],
                          all_replicate_Vcmax_Jmax[,'plot'],sep='-')
    invalid_curve = all_curve_id[is.na(all_replicate_Vcmax_Jmax[,'Jmax_at_25'])]

    index_to_remove <- which(curve_identifier==invalid_curve)
    if (length(index_to_remove) > 0) {
      Tleaf_all = Tleaf_all[-index_to_remove]
      Qin_all   = Qin_all[-index_to_remove]
      Ci_all    = Ci_all[-index_to_remove]
      curve_identifier    = curve_identifier[-index_to_remove]
      An_obs    = An_obs[-index_to_remove]
    }
}else if(curve_type==2){
  licor_data <- read.csv.exdf("for_calibration_with_OBS/my_scripts/ld11_bb.csv") #this comes from "get_ld11_bb.R"

  Tleaf_all  <- licor_data[, 'TleafCnd'] 
  Qin_all    <- licor_data[, 'Qin']
  Ci_all     <- licor_data[, 'Ci']

  date = substr(licor_data[, 'date'], 1, 8)
  
  curve_identifier = paste(date,licor_data[, 'instrument'],licor_data[, 'plot'],licor_data[, 'replicate'],sep='-')
  
  An_obs = licor_data[,'A']
}


lbs = rep(1.5,number_of_enzymes) 
ubs = rep(2.5,number_of_enzymes)
init_guess = rep(2.0,length(lbs))

#You can set lambda to 0 to NOT use LASSO
lambda = 0.5

if(keyword == "tobacco_Amanda"){
  subset_condition = which(Ci_all>100 & Ci_all<1000)
  split_obs <- split(obs_num[subset_condition], curve_identifier[subset_condition])
  print(split_obs)
  #check if all replicates have the same obs number sequences
  #by comparing to the first vector
  all_equal <- all(sapply(split_obs, function(x) setequal(x, split_obs[[1]])))
  if(!all_equal) stop("Yufeng: subset condition error!")
}else{
  #use only data with Ci larger than 100ppm
  subset_condition = which(Ci_all>100)
}

constants_list = list(curve_type    = curve_type,
                 alpha1        = alpha1,
                 alpha2        = alpha2,
                 Vcmax25       = Vcmax25,
                 Jmax25        = Jmax25,
                 Rd25          = Rd25,
                 TPU25         = TPU25,
                 Tleaf_all     = Tleaf_all[subset_condition],
                 Qin_all       = Qin_all[subset_condition],
                 Ci_all        = Ci_all[subset_condition],
                 curve_identifier = curve_identifier[subset_condition],
                 An_obs        = An_obs[subset_condition],
                 targets_all   = targets_all,
                 targets       = targets
                 )
# Define the objective function for nloptr, including constant parameters
objective_function <- function(params) {
  l2_obj_function(params, constants_list, lambda)
}

opt_results = nloptr(x0 = init_guess, eval_f = objective_function, 
                     lb = lbs,
                     ub = ubs,
                     #NLOPT_LN_SBPLX, NLOPT_GN_ISRES
                     opts = list("algorithm" = "NLOPT_LN_SBPLX","xtol_rel"=1.0e-4,
                                 "maxeval"   = 200,
                                 "print_level"= 1)
                     )

Q10s  = as.data.frame(opt_results$solution)
rownames(Q10s) = names(targets)
Q10s  = t(Q10s)
write.csv(Q10s,paste0("outputs/ePhotosynthesis_optQ10_",outfile_prefix,"_",keyword,".csv"))

