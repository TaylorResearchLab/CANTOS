string_dissimilarity<- function (S1,S2,meth){
  library(stringdist)
  d<-stringdist(S1, S2, method = meth)
  return(d)
}
  