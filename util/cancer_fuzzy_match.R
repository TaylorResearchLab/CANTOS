cancer_fuzzy_match <- function(disease_name,who_general_words){
  
  Match_Flag <- "No"
  fuzzy_distance <- NA
  cancer_names<-NA
  for (distance in c(0.1,0.2)){
    if(length(which(agrepl(disease_name,who_general_words$who_cancers,ignore.case = TRUE,max.distance = distance)))>0){
      Match_Flag <- "Yes"
      fuzzy_distance =distance
      cancer_names = paste(who_general_words$who_cancers[which(agrepl(disease_name,who_general_words$who_cancers,ignore.case = TRUE,max.distance = distance))],collapse = ",")
      break
    }
  }
  
 return (list(Match_Flag,fuzzy_distance,cancer_names)) 
}
  
