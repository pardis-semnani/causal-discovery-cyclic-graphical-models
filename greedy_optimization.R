source("five_properties.R")
source("d_connections.R")


g <- function(n,d_connections,P) {#Assume n>=3
   b = c(0 , n^3 , n^(3:n) , n^6 , n^3*sum(n^(3:n)))
   #m=n+2
   s <- E6(n,d_connections,P) + b[n+2]*E4(n,d_connections,P)
   for (k in seq_int(1,n-2)) {
     s <- s - b[n+2]*prod(b[seq_int(k+3,n+1)]+1)*length(mutually_exclusive(n,d_connections,P,k))
   }
   s <- s +  b[n+2]*prod(b[3:n+1]+1)*length(perfect_nonconductors(n,d_connections,P)) + b[n+2]*prod(b[2:n+1]+1)*length(p_adjacencies(n,d_connections,P))
   return(s)
}



greedy_optimization <- function(n,d_connections,P0) {
  Track <- list()
  L <- list(P0)
  i <- 1
  base <- g(n,d_connections,P0)
  K <- list(P0)
  while (length(L)>=i) {
    flag <- FALSE #true if a non-increasing neighbor is found
    P <- L[[i]]
    print(P)
    print(base)
    Track[[length(Track)+1]] <- P
    for (j in seq_int(1,length(P))) {
      if (j < length(P)) {
        #swaps P[j] and P[j+1]
        P_new <- P
        place_holder <- P[[j+1]]
        P_new[[j+1]] <- P[[j]]
        P_new[[j]] <- place_holder
        switch <- TRUE #true if P[j] and P[j+1] are d-separated given what comes before so far
        for (s in P[[j]]) {
          if (switch) {
            for (t in P[[j+1]]) {
              S <- as.double(seq(1,n)[lapply(seq(1,n),section,P=P) <= j+1 & seq(1,n)!=s & seq(1,n)!=t])
              if (list(list(min(s,t),max(s,t),S)) %in% d_connections) {
                switch <- FALSE
                if (!(list(P_new) %in% K)) {
                  base_new <- g(n,d_connections,P_new)
                  K[[length(K)+1]] <- P_new
                  if (base_new <= base) {
                    L <- list(P_new)
                    i <- 1
                    base <- base_new
                    flag <- TRUE
                  }
                }
                break
              }
            }
          }
        }
        if (flag) {break}
        if (switch & !(list(P_new) %in% L)) {L[[length(L)+1]] <- P_new}
        #adds elements of P[j] to P[j+1]
        for (k in 1:length(P[[j]])) {
          P_new <- P
          P_new[[j]] <- P[[j]][-k]
          P_new[[j+1]] <- sort(c(P[[j+1]],P[[j]][k]))
          if (length(P_new[[j]])==0) {P_new[[j]] <- NULL}
          if (!(list(P_new) %in% K)) {
            base_new <- g(n,d_connections,P_new)
            K[[length(K)+1]] <- P_new
            if (base_new <= base) {
              L <- list(P_new)
              i <- 1
              base <- base_new
              flag <- TRUE
              break
            }
          }
        }
        if (flag) {break}
        #adds elements of P[j+1] to P[j]
        for (k in 1:length(P[[j+1]])) {
          P_new <- P
          P_new[[j+1]] <- P[[j+1]][-k]
          P_new[[j]] <- sort(c(P[[j]],P[[j+1]][k]))
          if (length(P_new[[j+1]])==0) {P_new[[j+1]] <- NULL}
          if (!(list(P_new) %in% K)) {
            base_new <- g(n,d_connections,P_new)
            K[[length(K)+1]] <- P_new
            if (base_new <= base) {
              L <- list(P_new)
              i <- 1
              base <- base_new
              flag <- TRUE
              break
            }
          }
        }
        if (flag) {break}
      }
      #adds elements of P[j] to an empty set before it 
      for (k in 1:length(P[[j]])) {
        P_new <- c(P[seq_int(1,j-1)],list(c(P[[j]][k])),list(P[[j]][-k]),P[seq_int(j+1,length(P))])
        if (length(P_new[[j+1]])==0) {P_new[[j+1]] <- NULL}
        if (!(list(P_new) %in% K)) {
          base_new <- g(n,d_connections,P_new)
          K[[length(K)+1]] <- P_new
          if (base_new <= base) {
            L <- list(P_new)
            i <- 1
            base <- base_new
            flag <- TRUE
            break
          }
        }
      }
      if (flag) {break}
      #adds elements of P[j] to an empty set after it 
      for (k in 1:length(P[[j]])) {
        P_new <- c(P[seq_int(1,j-1)],list(c(P[[j]][-k])),list(P[[j]][k]),P[seq_int(j+1,length(P))])
        if (length(P_new[[j]])==0) {P_new[[j]] <- NULL}
        if (!(list(P_new) %in% K)) {
          base_new <- g(n,d_connections,P_new)
          K[[length(K)+1]] <- P_new
          if (base_new <= base) {
            L <- list(P_new)
            i <- 1
            base <- base_new
            flag <- TRUE
            break
          }
        }
      }
      if (flag) {break}
    }
    if (!flag) {i <- i+1}
  }
  return(Track)
}


connected <- function(V,E) {#Checks if the graph (V,E) is connected. V: list of vertices, E: list of vectors(edges)
  i <- 1
  L <- V[1]
  while (length(L) >=i) {
    v <- L[[i]]
    for (w in V) {
      if (list(c(min(v,w),max(v,w))) %in% E & !(list(w) %in% L)) {
        L[[length(L)+1]] <- w
      }
    }
    i <- i+1
  }
  if (length(L) == length(V)) {
    return (TRUE)
  } else {
    return (FALSE)
  }
}


two_edge_connected <- function(V,E) {#Checks if the graph (V,E) is 2-edge-connected. V: list of vertices, E: list of vectors(edges)
  for (k in seq_int(1,length(E))) {
    if (list(E[[k]][1]) %in% V & list(E[[k]][2]) %in% V) {
      if (!connected(V,E[-k])) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

check_connectivity <- function(n,d_connectons,P0) {
  Track <- greedy_optimization(n,d_connections,P0)
  for (i in length(Track):1) {
    P <- Track[[i]]
    E <- p_adjacencies(n,d_connections,P)
    flag <- TRUE
    for (section in P) {
      V <- list()
      for (v in section) {
        V[[length(V)+1]] <- v
      }
      if (!two_edge_connected(V,E)) {
        flag <- FALSE
        break
      }
    }
    if (flag) {
      return(P)
    }
  }
  return("failed!")
}

#You should be able to add elements to empty set, too - check
#You shouldn't get back to the neighbors already checked - check
#Each P has to be a list of vectors. - check
#The necessary condition for an ordered partition to work should be used. - check


#December 12:
#Can the algorithm get more efficient by analyzing g part by part?
#Can the search be over the space of graphs rather than ordered partitions? In this case, how should we define neighborhoods?