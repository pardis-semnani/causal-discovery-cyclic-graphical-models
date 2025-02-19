source("/Users/pardissemnani/Documents/UBC/Msc projects/Causality/Cyclic Causal Inference/paths.R")
library(combinat)




is.ancestor <- function(n,E,s,t,st.paths=list()) #Checks if s is an ancestor of t in G=([n],E)
{
  if(length(st.paths)==0) {st.paths <- st_paths(n,E,s,t)}
  if (length(st.paths)==0) {return(FALSE)}
  for (P in st.paths){
    flag <- TRUE
    for (i in seq(2,length(P)-1,by=2)){
      if (P[[i]][1]!=P[[i-1]]) {
        flag <- FALSE
        break
      }
    }
    if (flag) {return(TRUE)}
  }
  return(FALSE)
}

is.ancestor.set <- function(n,E,sT.paths) #Checks if s is an ancestor of T in G=([n],E) #sT.paths: list of all the paths from s to T in G
{
  if (length(sT.paths)==0) {return(FALSE)}
  for (P in sT.paths){
    flag <- TRUE
    for (i in seq(2,length(P)-1,by=2)){
      if (P[[i]][1]!=P[[i-1]]) {
        flag <- FALSE
        break
      }
    }
    if (flag) {return(TRUE)}
  }
  return(FALSE)
}


seq_int <- function(from,to,by=1){
  if (from > to) {
    return(c())
    }
  else {return(seq(from,to,by))}
}



d_connections <- function(n,E) #Outputs the list of d-connections in G=([n],E)
{
  L <- list()
  for (s in seq_int(1,n-1)) {
    for (t in seq_int(s+1,n)) {  
      
      PP <- st_paths(n,E,s,t)
      my_vec <- c(seq_int(1,s-1),seq_int(s+1,t-1),seq_int(t+1,n))
      subsets <- if(length(my_vec)<=1) {list(c(),my_vec)} else {c(list(c()),unlist(lapply(1:length(my_vec),combinat::combn,x = my_vec,simplify = FALSE),recursive = FALSE))}
      for (S in subsets){
        
        for (P in PP) {
          flag <- TRUE
          for (i in seq_int(3,length(P)-2,by=2)) {
            if (P[[i]] %in% S & P[[i-1]][2]!=P[[i+1]][2]) { #if a vertex of S exists on P, but is not a collider on P
              flag <- FALSE
              break
            }
            else if (P[[i-1]][2]==P[[i+1]][2] & !(P[[i]] %in% S) & !(TRUE %in% sapply(S,is.ancestor,n=n,E=E,s=P[[i]]))) {
              flag <- FALSE
              break
            }
          }
          if (flag) {L[[length(L)+1]] <- list(s,t,as.double(S)); break}
        }
      }
    }
  }
  return(L)
}





