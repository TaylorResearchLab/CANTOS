pedcan_fuzzy_match <- function(cancer_name,who_ped_words){
  
  Match_Flag <- "No"
  fuzzy_distance <- NA
  cancer_names<-NA
  for (distance in c(0.1,0.2)){
    if(length(which(agrepl(cancer_name,who_ped_words$who_paediatric_cancers,ignore.case = TRUE,max.distance = distance)))>0){
      Match_Flag <- "Yes"
      fuzzy_distance =distance
      cancer_names = paste(who_ped_words$who_paediatric_cancers[which(agrepl(cancer_name,who_ped_words$who_paediatric_cancers,
                                                                             ignore.case = TRUE,max.distance = distance))],collapse = ",")
      break
    }
  }
  
  return (list(Match_Flag,fuzzy_distance,cancer_names)) 
}