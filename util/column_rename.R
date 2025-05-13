column_rename <- function(col_names){
  
for(iter in 1:length(col_names)){
  if(col_names[iter]=="valid_euclidean_dist_v3" ){
    col_names[iter]="LTE-3+Euclidean Distance"
  }else if(col_names[iter]=="valid_af_v3"){
    col_names[iter]="LTE-3+AP"
    
  }else if(col_names[iter]=="valid_euclidean_dist_ada2"){
    col_names[iter]="ADA-002+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_af_ad2"){
    col_names[iter]="ADA-002+AP"
    
  }else if(col_names[iter]=="valid_kmeans_v3"){
    col_names[iter]="LTE-3+KMeans"
    
  }else if(col_names[iter]=="valid_kmeans_ada2"){
    col_names[iter]="ADA-002+KMeans"
    
  }else if(col_names[iter]=="valid_lv_match"){
    col_names[iter]="Levenshtein"
    
  }else if(col_names[iter]=="valid_af_lv"){
    col_names[iter]="Levenshtein+AP"
    
  }else if(col_names[iter]=="valid_jw_match"){
    col_names[iter]="Jaro Winkler"
  }else if(col_names[iter]=="valid_af_jw"){
    col_names[iter]="Jaro Winkler+AP"
    
  }else if(col_names[iter]=="valid_cosine_match"){
    col_names[iter]="Cosine"
    
  }else if(col_names[iter]=="valid_af_cosine"){
    col_names[iter]="Cosine+AP"
    
  }else if(col_names[iter]=="valid_euclidean_dist_llama"){
    col_names[iter]="LLama3.0+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_biobert"){
    col_names[iter]="BioBERT+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_medllama"){
    col_names[iter]="MedLlama_13B+Euclidean Dist"
    
  }else if(col_names[iter]=="valid_euclidean_dist_pubmedbert"){
    col_names[iter]="PubMedBERT+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_modernbert"){
    col_names[iter]="ModernBERT+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_medllama_7b"){
    col_names[iter]="MedLlama2+Euclidean Distance"
 
  }else if(col_names[iter]=="valid_euclidean_dist_llama_32_3b"){
    col_names[iter]="LLama3.2_3B+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_phi4"){
    col_names[iter]="Phi-4+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_llama_33_70b"){
    col_names[iter]="Llama3.3_70B+Euclidean Distance"
  }else if(col_names[iter]=="valid_euclidean_dist_MiniLM_L6_v2"){
    col_names[iter]="all-MiniLM-L6-v2+Euclidean Distance"
  }else if(col_names[iter]=="valid_euclidean_mpnet_base"){
    col_names[iter]="all-mpnet-base-v2+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_e5_large"){
    col_names[iter]="e5-large+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_gtr_t5_large"){
    col_names[iter]="gtr-t5-large+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_roberta"){
    col_names[iter]="all-roberta-large-v1+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_MiniLM_L12_v2"){
    col_names[iter]="all-MiniLM-L12-v2+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_labse"){
    col_names[iter]="LaBSE+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_sciBERT"){
    col_names[iter]="sciBERT+Euclidean Distance"
    
  }
  else if(col_names[iter]=="valid_euclidean_dist_sapBERT"){
    col_names[iter]="sapBERT+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_cohere"){
    col_names[iter]="embed-english-v2.0+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_deepseek"){
    col_names[iter]="DeepSeek_8B+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_BioGPT"){
    col_names[iter]="BioGPT+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_clincalBERT"){
    col_names[iter]="ClinicalBERT+Euclidean Distance"
    
  }else if(col_names[iter]=="valid_euclidean_dist_e5large_v2"){
    col_names[iter]="e5large_v2+Euclidean Dist"
    
  }else if(col_names[iter]=="valid_euclidean_dist_nomic"){
    col_names[iter]="nomic+Euclidean Distance"
  }
}
  return(col_names)
}