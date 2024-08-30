is_pedcan_who <- function(target_cancer_map,who_ped_words){
  target_cancer_map$is_pedcancer_who <-NA
  target_cancer_map$fuzzy_distance_ped<-NA
  target_cancer_map$matched_who_name_ped <- NA
  target_cancer_map_compressed<-create_cancer_who_map(target_cancer_map,who_ped_words)
  
  for (iter in 1: dim(target_cancer_map)[1]){
    cancer_name <- target_cancer_map$diseases[iter]
      if(!is.na(cancer_name)){
        index_cancer <- which(target_cancer_map_compressed==cancer_name)
        target_cancer_map$is_pedcancer_who[iter]<-target_cancer_map_compressed$is_pedcancer_who[index_cancer]
        target_cancer_map$fuzzy_distance_ped[iter]<-target_cancer_map_compressed$fuzzy_distance_ped[index_cancer]
        target_cancer_map$matched_who_name_ped[iter]<- target_cancer_map_compressed$matched_who_name_ped[index_cancer]
      }
    }
  
  return(target_cancer_map)
}

