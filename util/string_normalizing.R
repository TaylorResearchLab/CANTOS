string_normalzing<- function (S1,S2){
  library(tidyverse)
  s1_len= str_length(S1)
  s2_len=str_length(S2)
  return(max(s1_len,s2_len))
}
