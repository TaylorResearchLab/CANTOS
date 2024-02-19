distance_clusters<- function (vect1,cutoff=0.0005,tumor_names){

  lower_bound <- quantile(vect1, cutoff)
  outlier_ind <- which((vect1< lower_bound))
  clustered_tumors <- tumor_names[outlier_ind]
  clustered_tumors<-paste(clustered_tumors,collapse = ";")
  return(clustered_tumors)
}
