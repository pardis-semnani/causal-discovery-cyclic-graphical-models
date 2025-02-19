source("d_connections.R")
source("five_properties.R")

library(igraph)
library(R.utils)



interruptor <- function(FUN,args, time.limit, ALTFUN){
  
  results <- 
    tryCatch({
      withTimeout({FUN(args)}, timeout=time.limit)
    }, error = function(e){
      if(grepl("reached elapsed time limit",e$message))
        ALTFUN(args) else
          paste(e$message,"EXTRACTERROR")
    })
  
  if(grepl("EXTRACTERROR",results)){
    print(gsub("EXTRACTERROR","",results))
    results <- NULL
  } 
  
  return(results)
}   


is.ancestor_cycle <- function(n,E,Reachable,C1,C2) {
  for (u in C1) {
    for (v in C2) {
      if (v %in% Reachable[[u]]) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}


check_algorithm <- function(n,p,N){#n: number of vertices, #p: how likely an edge is to appear in the random graph, #N:maximum number of iterations allowed for each cycle
  #random graph
  g <- erdos.renyi.game(n, p, type = "gnp",directed=TRUE)
  #defining the edge set E
  E_mat <- get.edgelist(g)
  E <- list()
  for (i in 1:dim(E_mat)[1]) {E[[i]] <- c(E_mat[i,])}
  print(c("our edge set is",E))
  # E=list(c(3,1),c(6,1),c(4,2),c(7,2),c(6,3),c(9,3),c(3,5),c(10,5),c(14,5),c(16,5),c(2,6),c(8,6),c(17,6),c(3,7),c(5,7),c(6,7),
  #        c(3,8),c(7,9),c(15,10),c(18,11),c(13,12),c(2,13),c(7,13),c(12,13),c(9,14),c(6,15),c(9,15),c(12,15),c(18,15),c(16,20),
  #        c(2,17),c(4,17),c(18,20),c(20,18),c(6,19),c(7,19))
  # E = list(c(8,1),c(3,4),c(5,4),c(9,4),c(2,6),c(8,6),c(3,7),c(3,8),c(4,8),c(8,10),c(2,9),c(3,9),c(5,9),c(8,9))
  # E = list(c(2,1),c(7,1),c(9,1),c(5,2),c(10,2),c(2,3),c(8,4),c(10,5),c(3,7),c(6,7),c(3,9),c(6,9))
  #defining the vertex set V
  V <- list()
  for (v in 1:n) {V[[v]] <- v}
  #finding the corresponding partition
  Reachable <- list()
  for (s in 1:n) {
    Reachable[[s]] <- c(s)
  }
  for (e in E) {
    for (s in 1:n) {
      if (e[1] %in% Reachable[[s]]) {
        Reachable[[s]] <- unique(c(Reachable[[s]],Reachable[[e[2]]]))
      }
    }
  }
  
  P <- list()
  for (s in 1:n) {P[[s]] <- c(s)}
  for (s in seq_int(1,n-1)) {
    for (t in seq_int(s+1,n)) {
      if (s %in% Reachable[[t]] & t %in% Reachable[[s]]) {
        i <- section(s,P)
        j <- section(t,P)
        if (i!=j) {P[[i]] <- unique(c(P[[i]],P[[j]]));P[[j]] <- NULL}
      }
    }
  }
  #determining the order of the partition
  for (i in seq_int(1,length(P)-1)) {
    for (j in seq_int(i+1,length(P))) {
      if (is.ancestor_cycle(n,E,Reachable,P[[j]],P[[i]])) {
        a <- P[[j]]
        P[[j]] <- P[[i]]
        P[[i]] <- a
      }
    }
  }
  print(c("length of Partition is",length(P)))
  print(P)
  #determining the descendants of common children in for each pair of vertices
  Common_Children <- list()
  for (s in seq_int(1,n-1)) {
    for (t in seq_int(s+1,n)) {
      index <- seq_int(max(section(s,P),section(t,P))+1,length(P))
      for (i in index) {
        for (v in P[[i]]) {
          if (list(c(s,v)) %in% E & list(c(t,v)) %in% E) {
            Common_Children[[as.character(i)]] <- unique(append(Common_Children[[as.character(i)]],list(c(s,t))))
            for (j in seq_int(i+1,length(P))) {
              if (is.ancestor_cycle(n,E,Reachable,P[[i]],P[[j]])) {
                Common_Children[[as.character(j)]] <- unique(append(Common_Children[[as.character(j)]],list(c(s,t))))
                # index <- index[-j]
              }
            }
            break
          }
        }
      }
    }
  }
  print(c("common children ",Common_Children))
  D <- list()
  #start recovering cycle by cycle
  for(i in 1:length(P)) {
    #defining ComCh for cycle P[[i]]
    ComCh <- list()
    for (p in Common_Children[[as.character(i)]]) {
      flag <- TRUE
      for (j in seq_int(1,i-1)) {
        if (is.ancestor_cycle(n,E,Reachable,P[[j]],P[[i]]) & list(p) %in% Common_Children[[as.character(j)]]) {
          flag <- FALSE
          break
        }
      }
      if (flag) {ComCh[[length(ComCh)+1]] <- p}
    }
    ComCh <- unique(ComCh)
    print(c("ComChi is",ComCh))
    #defining NoComCh for cycle P[[i]]
    NoComCh <- list()
    for (s in seq_int(1,n-1)) {
      for (t in seq_int(s+1,n)) {
        if (!(list(c(s,t)) %in% Common_Children[[as.character(i)]]) & section(s,P)<i & section(t,P)<i) {
          NoComCh <- append(NoComCh,list(c(s,t)))
        }
      }
    }
    print(c("NoComChi is",NoComCh))
    #defining V1
    V1 <- list()
    for (v in P[[i]]) {
      V1 <- append(V1, list(v))
    }
    #defining E1 (inner p-adjacencies) and H1 (outer-p-adjacencies) for each cycle
    E1 <- list()
    E1_directed <- list()
    H1 <- list()
    for (e in E) {
      if (e[1] %in% V1 & e[2] %in% V1) {
        E1 <- append(E1,list(sort(e)))
        E1_directed <- append(E1_directed,list(e))
      }
      else if(e[2] %in% V1) {
        H1 <- append(H1,list(e))
      }
    }
    for (e1 in H1) {
      for (e2 in E1_directed) {
        if (e1[2]==e2[2]) {
          H1 <- append(H1,list(c(e1[1],e2[1])))
        }
      }
    }
    for (e1 in E1_directed) {
      for (e2 in E1_directed) {
        if (e1[1]!=e2[1] & e1[2]==e2[2]) {
          E1 <- append(E1,list(sort(c(e1[1],e2[1]))))
        }
      }
    }
    E1 <- unique(E1)
    H1 <- unique(H1)
    print(c("E1",E1))
    print(c("H1",H1))
    print(c("We are recovering cycle no.",i))
    
    I <<- 0
    A <- recover(V=V1,N=N,E=E1,H=H1,ComCh=ComCh,NoComCh=NoComCh)
    if (is.character(A)) {return("failed")}
    D <- append(D,A)
  }
  
  # plot(graph(c(unlist(unique(D)))))
  return(c(unique(D),length(P)))
}

# recover <- function(V,N,E,H,ComCh,NoComCh,D=list(),D_copy=list(),Cause=list(),avoid=0,a=0)