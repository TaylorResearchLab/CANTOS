# This script is used for splitting conjoined drug names from clinical trials into individual drug names.

split_drugs <- function(drug_list_refined){
  
  # Remove explanation in ()
  drug_list_refined$Intervention_CT<-gsub("\\s*\\([^\\)]+\\)","",as.character(drug_list_refined$Intervention_CT))
  
  
  #Remove eplanation in []
  drug_list_refined$Intervention_CT<-gsub("\\[[^][]*]", "",drug_list_refined$Intervention_CT)
  
  # Remove copyright symbol
  
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "®", "")
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "™", "")
  
 
  # Remove the + sign 
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "\\+", ",") # Replace + with ,
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "\\bplus\\b", ",") # Replace + with ,
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "\\bPLUS\\b", ",") # Replace + with ,
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "\\bPlus\\b", ",") # Replace + with ,
  
  
  drug_list_refined <- separate_rows(drug_list_refined, Intervention_CT, sep = ",")
  
  
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "\\band\\b", ",")# Replace ands with ,
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "\\bAnd\\b", ",")# Replace ands with ,
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "\\bAND\\b", ",")# Replace ands with ,
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "\\&", ",") # Replace + with ,
  
  drug_list_refined <- separate_rows(drug_list_refined, Intervention_CT, sep = ",")
  
  
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "\\bwith\\b", ",")# Replace ands with ,
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "\\bWith\\b", ",")# Replace ands with ,
  drug_list_refined$Intervention_CT<-str_replace_all(drug_list_refined$Intervention_CT, "\\bWITH\\b", ",")# Replace ands with ,
  
  drug_list_refined <- separate_rows(drug_list_refined, Intervention_CT, sep = ",")
  
  # Split based on ;
  drug_list_refined <- separate_rows(drug_list_refined, Intervention_CT, sep = ";")
  
  
  #Remove preceeding white space
  drug_list_refined$Intervention_CT<-trimws(drug_list_refined$Intervention_CT)
  
  # Remove blank and NA 
  drug_list_refined<-drug_list_refined[!(is.na(drug_list_refined$Intervention_CT) | drug_list_refined$Intervention_CT==""), ]
  
  # make lower case
  drug_list_refined$Intervention_CT<-tolower(drug_list_refined$Intervention_CT)
  
  # Uncomment if you are using CT data directly from CT website and not ACCT
  # corticosteroid:
  # drug_list_refined$Intervention_CT <- gsub(".*corticosteroid:", "",drug_list_refined$Intervention_CT)  
  # drug_list_refined$Intervention_CT<-trimws(drug_list_refined$Intervention_CT) # remove white space before corticosteroid
  
  # Remove (numeric) mg/ml 
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) mg/ml.*", "", drug_list_refined$Intervention_CT)
  # 1200 mg/kg
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) mg/kg.*", "", drug_list_refined$Intervention_CT)
  #375 mg/m2
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) mg/m2.*", "", drug_list_refined$Intervention_CT)
  #25 mg/1 ml
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) mg/([0-9]+) ml.*", "", drug_list_refined$Intervention_CT)
  
  #0.5 mcg/kg
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) mcg/kg.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)mcg/kg.*", "", drug_list_refined$Intervention_CT)
  
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) mcg/ml.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)mcg/ml.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) mcg.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)mcg.*", "", drug_list_refined$Intervention_CT)
  
  #5 µg / 0.5 ml
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) \U03BCg/([0-9]+) ml.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) \U03BCg/ ([0-9]+) ml.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)\U03BCg/ ([0-9]+) ml.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)\U03BCg/([0-9]+) ml.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)\U03BCg/([0-9]+)ml.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)\U03BCg.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) \U03BCg.*", "", drug_list_refined$Intervention_CT)
  
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)µg.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT <-  gsub("([0-9]+) µg.*", "", drug_list_refined$Intervention_CT)
  
  
  # [0-9] mg descriptions 
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) mg.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)mg.*", "", drug_list_refined$Intervention_CT)
  
  
  # 10000u/ml
  
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)u/ml.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) u/ml.*", "", drug_list_refined$Intervention_CT)
  
  # gm 
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) gm.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)gm.*", "", drug_list_refined$Intervention_CT)
  
  # Injection
  drug_list_refined$Intervention_CT<- gsub("Injection.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("injection.*", "", drug_list_refined$Intervention_CT)
  
  
  
  #1.0ug/kg
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)ug/kg.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) ug/kg.*", "", drug_list_refined$Intervention_CT)
  
  # ml/kg
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) ml/kg.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)ml/kg.*", "", drug_list_refined$Intervention_CT)
  
  #ml
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) ml.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)ml.*", "", drug_list_refined$Intervention_CT)
  
  
  # 100u / 0.5ml
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)u / ([0-9]+)ml.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)u /([0-9]+)ml.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)u/([0-9]+)ml.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)u/ ([0-9]+)ml.*", "", drug_list_refined$Intervention_CT)
  
  
  
  
  # Remove %
  drug_list_refined$Intervention_CT<- gsub("([0-9]+) %.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("([0-9]+)%.*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("%([0-9]+).*", "", drug_list_refined$Intervention_CT)
  drug_list_refined$Intervention_CT<- gsub("% ([0-9]+).*", "", drug_list_refined$Intervention_CT)
  
  
  return(drug_list_refined)
}