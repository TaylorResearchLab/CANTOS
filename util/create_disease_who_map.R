create_disease_who_map<- function (target_disease_map,who_general_words){
  target_disease_map_compressed <- target_disease_map %>% distinct(diseases)
  
  for (iter in 1: dim(target_disease_map_compressed)[1]){
    print(iter)
    disease_name <- target_disease_map_compressed$diseases[iter]
    #if (length(which(grepl(disease_name,who_general_words$who_cancers,ignore.case = TRUE)))>0){
    if (length(which(grepl(disease_name,who_general_words$who_cancers,fixed=TRUE)))>0){
      target_disease_map_compressed$is_cancer_who[iter] <- "Yes"
    } 
    else{
      if(!is.na(disease_name)){
        fuzzy_values <- cancer_fuzzy_match(disease_name,who_general_words)
        target_disease_map_compressed$is_cancer_who[iter]<-fuzzy_values[[1]]
        target_disease_map_compressed$fuzzy_distance[iter]<-fuzzy_values[[2]]
        target_disease_map_compressed$matched_who_name[iter] <- fuzzy_values[[3]]
      }
    }
  }
  return(target_disease_map_compressed)
}