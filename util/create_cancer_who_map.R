create_cancer_who_map<- function (target_cancer_map,who_ped_words){
  target_cancer_map_compressed <- target_cancer_map %>% distinct(diseases)
  
  for (iter in 1: dim(target_cancer_map_compressed)[1]){
    print(iter)
    cancer_name <- target_cancer_map_compressed$diseases[iter]
    if (length(which(grepl(cancer_name,who_ped_words$who_paediatric_cancers,fixed = TRUE)))>0){
      target_cancer_map_compressed$is_pedcancer_who[iter] <- "Yes"
    } 
    else{
      if(!is.na(cancer_name)){
        fuzzy_values <- pedcan_fuzzy_match(cancer_name,who_ped_words)
        target_cancer_map_compressed$is_pedcancer_who[iter]<-fuzzy_values[[1]]
        target_cancer_map_compressed$fuzzy_distance_ped[iter]<-fuzzy_values[[2]]
        target_cancer_map_compressed$matched_who_name_ped[iter] <- fuzzy_values[[3]]
      }
    }
  }
  return(target_cancer_map_compressed)
}