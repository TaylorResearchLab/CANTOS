is_cancer_who <- function(target_disease_map,who_general_words){
  target_disease_map$is_cancer_who <-NA_character_
  target_disease_map$fuzzy_distance<-NA_character_
  target_disease_map$matched_who_name <- NA_character_
  
  target_disease_map_compressed<-create_disease_who_map(target_disease_map,who_general_words)
  
  
  for (iter in 1: dim(target_disease_map)[1]){
    print(iter)
    disease_name <- target_disease_map$diseases[iter]
    if(!is.na(disease_name)){
    index_disease <- which(target_disease_map_compressed==disease_name)
    target_disease_map$is_cancer_who[iter]<-target_disease_map_compressed$is_cancer_who[index_disease]
    target_disease_map$fuzzy_distance[iter]<-target_disease_map_compressed$fuzzy_distance[index_disease]
    target_disease_map$matched_who_name[iter]<- target_disease_map_compressed$matched_who_name[index_disease]
  }
  }
  return(target_disease_map)
}