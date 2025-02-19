paths <- function(n,E,m,m.minus.1_paths) #n: number of vertices, E: set of edges, m: length of paths, m.minus.1_paths: list of paths in G=([n],E) of length m-1
{
  if (m==1) {
    L <- list()
    for (e in E) {
      L[[length(L)+1]] <- list(e[1],e,e[2])
      L[[length(L)+1]] <- list(e[2],e,e[1])
    }
    return(L)
  }
  else {
    L <- list()
    for (P in m.minus.1_paths) {
      end <- P[[length(P)]]
      for (v in seq_int(1,n)){
        if (!(v %in% P) & (list(c(end,v)) %in% E)) {
          P_copy <- P
          P_copy[[length(P_copy)+1]] <- c(end,v)
          P_copy[[length(P_copy)+1]] <- v
          L[[length(L)+1]] <- P_copy
        }
        if (!(v %in% P) & (list(c(v,end)) %in% E)) {
          P_copy <- P
          P_copy[[length(P_copy)+1]] <- c(v,end)
          P_copy[[length(P_copy)+1]] <- v
          L[[length(L)+1]] <- P_copy
        }
      }
    }
    return(L)
  }
}
  
  
  
  st_paths <- function(n,E,s,t) #Outputs a list of all st-paths in G=([n],E) 
  {
    L <- list()
    m.minus.1_paths <- list()
    for (m in 1:n){
      m.minus.1_paths <- paths(n,E,m,m.minus.1_paths)
      for (P in m.minus.1_paths) {
        if (P[[1]]==s & P[[length(P)]]==t){
          L[[length(L)+1]] <- P
        }
      }
    }
    return(L)
  }
  
  
  
  all_paths <- function(n,E) #Outputs a list of all paths in G=([n],E) 
  {
    L <- list()
    m.minus.1_paths <- list()
    for (m in 1:length(E)){
      m.minus.1_paths <- paths(n,E,m,m.minus.1_paths)
      L <- c(L,m.minus.1_paths)
    }
    return(L)
  }
  
  all_paths_separated <- function(n,E) #Outputs the list of all st-paths for all s,t in [n]
  {
    M = matrix(rep(list(),n^2),nrow=n,ncol=n)
    PATHS <- all_paths(n,E)
    for (P in PATHS) {
      T= append(M[P[[1]],P[[2]]],P)
      M[[P[[1]],P[[length(P)]]]][[length(M[[P[[1]],P[[length(P)]]]])+1]] <- P
    }
    return(M)
  }