merge_results<-function(site_id,years,exp_id,out_key,folder_path){
  result_list = list()
  for (i in 1:length(years)){
    year = years[i]
    output_folder = paste0("../",folder_path,"/results_site",site_id,"_exp",exp_id,"_",year,"_",out_key)  #the folder to save daily outputs
    print(output_folder)
    all_files    = dir(path=output_folder,pattern="*.rds",full.names = TRUE)
    result_oneyear = c()
    for (fn in all_files){
      result = readRDS(fn)
      #remove the last record since it's the extra hour saved
      result = result[1:24,]
      result_oneyear = rbind(result_oneyear,result)
    }
    result_list[[i]] = result_oneyear
  }
  return(result_list)
}