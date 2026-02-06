rm(list=ls())

use_default_Q10s = TRUE 

if(use_default_Q10s){
  exp_type = "CTL"
}else{
  exp_type = "Q10opt"
}
keyword = "ld11"

executable_path <- "./myephoto_single.exe" 

alpha1_alpha2 = read.csv(paste0('../parameterization/outputs/ePhotosynthesis_optimal_alpha1_alpha2_',keyword,'.csv'))
alpha1 = alpha1_alpha2[2]
alpha2 = alpha1_alpha2[3]

if(!use_default_Q10s){
  Q10s = read.csv(paste0('../parameterization/outputs/ePhotosynthesis_optQ10_ACi_',keyword,'.csv'))
  Q10s = Q10s[-1]
}

print(alpha1_alpha2)
PAR   = 2000
Tleaf = 25
Ci    = 400
#call ephoto c++
if(use_default_Q10s){
  args = c(alpha1,alpha2,PAR,Tleaf, Ci)
}else{
  args = c(alpha1,alpha2,PAR,Tleaf, Ci,Q10s)
}
args = as.numeric(args)
system(sprintf(
  "%s %s > log 2>&1",
  shQuote("./myephoto_single.exe"),
  paste(args, collapse = " ")
))


