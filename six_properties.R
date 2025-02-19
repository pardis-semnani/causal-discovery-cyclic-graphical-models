source("d_connections.R")
library(cIRT)



check_d_connections <- function(s,t,S,d_connections=list(),PATHS=list(),A="a") { #A = (I-C)^(-1) (I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM
  if (length(PATHS)==0 & is.character(A)) {
    if (list(list(s,t,S)) %in% d_connections) {return(TRUE)}
    else {return(FALSE)}
  }
  else if (is.character(A)) {
    for (P in PATHS[[s,t]]) {
      flag <- TRUE
      for (i in seq_int(3,length(P)-2,by=2)) {
        if (P[[i]] %in% S & P[[i-1]][2]!=P[[i+1]][2]) { #if a vertex of S exists on P, but is not a collider on P
          flag <- FALSE
          break
        }
        else if (P[[i-1]][2]==P[[i+1]][2] & !(P[[i]] %in% S)){
          S.PATHS <- list() 
          for (s in S) {
            S.PATHS <- append(S.PATHS,PATHS[[P[[i]],s]])
          }
          if (!is.ancestor.set(n,E,S.PATHS)) {
            flag <- FALSE
            break
          }
        }
      }
      if (flag) {return(TRUE)}
    }
    return(FALSE)
  }
  else {
    alp = 0
    Sig <- A
    
    if ((length(S)==0 && abs(Sig[s,t])<=alp) || (length(S)>1 && abs(det(Sig[c(unlist(S),s),c(unlist(S),t)]))<=alp)) {
      return(FALSE)
    }
    else {return(TRUE)}
  }
}

section <- function(s,P) #Outputs i such that P[[i]] is the section of the ordered partition P containing s
{
  return(as.double(which(lapply(P,match,x=c(s),nomatch=0)!=0)))
}

##############################################

p_adjacencies <- function(n,d_connections=list(),P,order,PATHS=list(),A="a") #Outputs the list of p-adjacencies associated to the partially ordered partition (P,order), i.e. E^(1)_(P,order)
{#PATHS: all st-separated paths in the original graph in the form of a matrix - d_connections: all d_connections in the original graph / Exactly one of them should be given to the graph
 #A = (I-C)^(-1) (I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM
  L <- list()
  for (s in seq_int(1,n-1)) {
    for (t in seq_int(s+1,n)) {
      S <- c()
      for (z in seq_int(1,n)) {
        if ((list(c(section(z,P),section(s,P))) %in% order | list(c(section(z,P),section(t,P))) %in% order) & z!=s & z!=t) {
          S <- unique(append(S,z))
        }
      }
      S <- as.double(S)
      if (check_d_connections(s,t,S,d_connections,PATHS,A)) {
        L[[length(L)+1]] <- as.double(c(s,t))
      }
    }
  }
  return(L)
}

##############################################

p_adjacent <- function(n,d_connections=list(),P,order,s,t,PATHS=list(),A="a") #Checks whether s and t are p-adjacent with respect to (P,order)
{#A = (I-C)^(-1) (I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM
  S <- c()
  for (z in seq_int(1,n)) {
    if ((list(c(section(z,P),section(s,P))) %in% order | list(c(section(z,P),section(t,P))) %in% order) & z!=s & z!=t) {
      S <- unique(append(S,z))
    }
  }
  S <- as.double(S)
  if (s!=t & check_d_connections(min(s,t),max(s,t),S,d_connections,PATHS,A)) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

##############################################

itineraries <- function(n,d_connections=list(),P,order,S,k,kminus1.itineraries="a",PATHS=list(),A="a") #Outputs all the uncovered itineraries in S of length k with respect to the p-adjacencies associated with the partially ordered partition (P,order)
  #S is a vector or list
  #A = (I-C)^(-1) (I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM
{
  L <- list()
  if (k==1) {
    S <- unlist(S)
    for (a in S) {
      L[[length(L)+1]] <- list(a)
    }
    return(L)
  } 
  else {
    S <- unlist(S)
    if (is.character(kminus1.itineraries)) {kminus1.itineraries <- itineraries(n=n,d_connections=d_connections,P=P,order=order,S=S,k=k-1,PATHS=PATHS,A=A)}
    for (I in kminus1.itineraries) {
      for (a in S) {
        if ((p_adjacent(n,d_connections,P,order,a,I[[1]],PATHS,A)) & !(TRUE %in% lapply(I[seq_int(2,k-1)],p_adjacent,n=n,d_connections=d_connections,P=P,order=order,s=a,PATHS=PATHS,A=A)) & !(a %in% I)) {
          L[[length(L)+1]] <- c(list(a),I)
        }
      }
    }
    return(L)
  }
}

##############################################

mutually_exclusive <- function(n,d_connections=list(),P,order,k,kminus1.itineraries="a",PATHS=list(),A="a") #Outputs the list of mutually exclusive conductors <x_0,...,x_{k+1}> associated to the partially ordered partition (P,order)
{#A = (I-C)^(-1) (I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM
  L <- list()
  IT <- list()
  for(i in 1:length(P)) {
    if (is.character(kminus1.itineraries)) {it <- itineraries(n,d_connections,P,order,P[[i]],k,kminus1.itineraries,PATHS,A)}
    else {it <- itineraries(n,d_connections,P,order,P[[i]],k,kminus1.itineraries[[i]],PATHS,A)}
    IT [[length(IT) + 1]] <- it
    for (I in it) {
      for (a in seq_int(1,n)) {
        if (!(list(c(i,section(a,P))) %in% order) & p_adjacent(n,d_connections,P,order,a,I[[1]],PATHS,A) & !(TRUE %in% lapply(I[seq_int(2,k)],p_adjacent,n=n,d_connections=d_connections,P=P,order=order,s=a,PATHS=PATHS,A=A))) {
          for (b in seq_int(1,n)) {
            if (!(list(c(i,section(b,P))) %in% order) & p_adjacent(n,d_connections,P,order,b,I[[k]],PATHS,A) & !(TRUE %in% lapply(I[seq_int(1,k-1)],p_adjacent,n=n,d_connections=d_connections,P=P,order=order,s=b,PATHS,A=A)) & a!=b & !p_adjacent(n,d_connections,P,order,a,b,PATHS,A=A)) {
              L [[length(L)+1]] <- c(list(a),I,list(b)) 
            }
          }
        }
      }
    }
  }
  return(list(L,IT))
}

mutually_exclusive_efficient <- function(n,d_connections=list(),P,order,k,kminus1.itineraries="a",PATHS=list(),p_adjacencies,A="a") #Outputs the list of mutually exclusive conductors <x_0,...,x_{k+1}> associated to the partially ordered partition (P,order)
{
  L <- list()
  IT <- list()
  is.itineraries.checked <- list()
  for (i in seq_int(1,length(P))) {
    IT[[i]] <- list()
    is.itineraries.checked[[i]] <- FALSE
  }
  for (a in seq_int(1,n-1)) {
    for (b in seq_int(a+1,n)) {
      if (!(list(c(a,b)) %in% p_adjacencies)) {
        for (i in 1:length(P)) {
          if (!(list(c(i,section(a,P))) %in% order) & !(list(c(i,section(b,P))) %in% order)) {
            if (!is.itineraries.checked[[i]]) {
              if (is.character(kminus1.itineraries)) {it <- itineraries(n,d_connections,P,order,P[[i]],k,kminus1.itineraries,PATHS,A)}
              else {it <- itineraries(n,d_connections,P,order,P[[i]],k,kminus1.itineraries[[i]],PATHS,A)}
              IT [[i]] <- it
              is.itineraries.checked[[i]] <- TRUE
            }
            for (I in IT[[i]]) {
              if (list(c(min(a,I[[1]]),max(a,I[[1]]))) %in% p_adjacencies & list(c(min(b,I[[k]]),max(b,I[[k]]))) %in% p_adjacencies) {
                flag <- TRUE
                for (j in seq_int(2,k)) {
                  if (list(c(min(a,I[[j]]),max(a,I[[j]]))) %in% p_adjacencies) {
                    flag <- FALSE
                    break
                  }
                }
                if (flag) {
                  for (j in seq_int(1,k-1)) {
                    if (list(c(min(b,I[[j]]),max(b,I[[j]]))) %in% p_adjacencies) {
                      flag <- FALSE
                      break
                    }
                  }
                }
                if (flag) {L [[length(L)+1]] <- c(list(a),I,list(b)) }
              }
            }
          }
        }
      }
    }
  }
  return(list(L,IT))
}

##############################################

perfect_nonconductors <- function(n,d_connections=list(),P,order,PATHS=list(),A="a") {#Outputs all the unshielded perfect nonconductors in the graph associated with (P,order) - A: coefficient matrix of the linear Gaussian SEM
  L <- list()
  for (a in seq_int(1,n)) {
    for (c in seq_int(1,n)) {
      if (a!=c & !p_adjacent(n,d_connections,P,order,a,c,PATHS,A)) {
        for (b in seq_int(1,n)) {
          if (!(list(c(section(b,P),section(a,P))) %in% order) & !(list(c(section(b,P),section(c,P))) %in% order) & p_adjacent(n,d_connections,P,order,a,b,PATHS,A) & p_adjacent(n,d_connections,P,order,b,c,PATHS,A)) {
            S <- c()
            for (z in seq_int(1,n)) {
              if ((list(c(section(z,P),section(a,P))) %in% order | list(c(section(z,P),section(b,P))) %in% order | list(c(section(z,P),section(c,P))) %in% order) & z!=a & z!=c) {
                S <- append(S,z)
              }
            }
            S <- as.double(S)
            if (check_d_connections(min(a,c),max(a,c),S,d_connections,PATHS,A)) {
              L[[length(L)+1]] <- list(a,b,c) 
            }
          }
        }
      }
    }
  }
  return(L)
}

perfect_nonconductors_efficient <- function(n,d_connections=list(),P,order,PATHS=list(),non_conductors,A="a") {#Outputs all the unshielded perfect nonconductors in the graph associated with (P,order) - #A = (I-C)^(-1) (I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM
  L <- list()
  for (triple in non_conductors) {
    a <- triple[[1]]
    b <- triple [[2]]
    c <- triple [[3]]
    S <- c()
    for (z in seq_int(1,n)) {
      if ((list(c(section(z,P),section(a,P))) %in% order | list(c(section(z,P),section(b,P))) %in% order | list(c(section(z,P),section(c,P))) %in% order) & z!=a & z!=c) {
        S <- unique(append(S,z))
      }
    }
    S <- as.double(S)
    if (check_d_connections(min(a,c),max(a,c),S,d_connections,PATHS,A)) {
      L[[length(L)+1]] <- list(a,b,c) 
    }
  }
  return(L)
}
##############################################

E4 <- function(n,d_connections=list(),P,order,PATHS=list(),A="a") #Property 4 of Richardson's Theorem
{
  counter <- 0
  for (a in seq_int(1,n)) {
    for (c in seq_int(1,n)) {
      if (a!=c & !p_adjacent(n,d_connections,P,order,a,c,PATHS,A)) {
        for (b1 in seq_int(1,n)) {
          S <- c()
          for (z in seq_int(1,n)) {
            if ((list(c(section(z,P),section(a,P))) %in% order | list(c(section(z,P),section(b1,P))) %in% order | list(c(section(z,P),section(c,P))) %in% order) & z!=a & z!=c) {
              S <- append(S,z)
            }
          }
          S <- as.double(S)
          if (b1!=a & b1!=c & p_adjacent(n,d_connections,P,order,a,b1,PATHS,A) & p_adjacent(n,d_connections,P,order,c,b1,PATHS,A)
              & !(list(c(section(b1,P),section(a,P))) %in% order) & !(list(c(section(b1,P),section(c,P))) %in% order) & !check_d_connections(min(a,c),max(a,c),S,d_connections,PATHS,A)) {
            for (b2 in seq_int(1,n)) {
              S <- c()
              for (z in seq_int(1,n)) {
                if ((list(c(section(z,P),section(a,P))) %in% order | list(c(section(z,P),section(b2,P))) %in% order | list(c(section(z,P),section(c,P))) %in% order) & z!=a & z!=c) {
                  S <- append(S,z)
                }
              }
              S <- as.double(S)
              if (b2!=a & b2!=c & b1!=b2 & p_adjacent(n,d_connections,P,order,a,b2,PATHS,A) & p_adjacent(n,d_connections,P,order,c,b2,PATHS,A)
                  & !(list(c(section(b2,P),section(a,P))) %in% order) & !(list(c(section(b2,P),section(c,P))) %in% order) & !check_d_connections(min(a,c),max(a,c),S,d_connections,PATHS,A)
                  & list(c(section(b1,P),section(b2,P))) %in% order) {
                counter <- counter + 1
              }
            }
          }
        }
      }
    }
  }
  return(counter)
}


E4_efficient <- function(n,d_connections=list(),P,order,PATHS=list(),nonconductors,perfect_nonconductors,A="a") {#A = (I-C)^(-1) (I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM
  counter <- 0
  for (triple1 in nonconductors) {
    for (triple2 in nonconductors) {
      if (triple1[[1]]==triple2[[1]] & triple1[[3]]==triple2[[3]] & triple1[[2]]!=triple2[[2]] & list(c(section(triple1[[2]],P),section(triple2[[2]],P))) %in% order & !(list(triple1) %in% perfect_nonconductors) & !(list(triple2) %in% perfect_nonconductors)) {
        counter <- counter+1
      }
    }
  }
  return(counter)
}

##############################################

E6 <- function(n,d_connections=list(),P,order,m_exclusive=list(),perfect_nonconductors=list(),PATHS=list(),A="a") {#Property 6 of Richardson's Theorem - #A = (I-C)^(-1) (I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM
  counter <- 0
  if (length(m_exclusive)==0) {
    for (k in seq_int(1,n-2)) {
      for (ME in mutually_exclusive(n=n,d_connections=d_connections,P=P,order=order,k=k,PATHS=PATHS,A)[[1]]) {
        a <- ME[[1]]
        c <- ME[[k+2]]
        for (b in seq_int(1,n)) {
          S <- c()
          for (z in seq_int(1,n)) {
            if ((list(c(section(z,P),section(a,P))) %in% order | list(c(section(z,P),section(b,P))) %in% order | list(c(section(z,P),section(c,P))) %in% order) & z!=a & z!=c) {
              S <- append(S,z)
            }
          }
          S <- as.double(S)
          if (b!=a & b!=c & p_adjacent(n,d_connections,P,order,a,b,PATHS,A) & p_adjacent(n,d_connections,P,order,b,c,PATHS,A) 
              & !(list(c(section(b,P),section(a,P))) %in% order) & !(list(c(section(b,P),section(c,P))) %in% order)
              & !check_d_connections(min(a,c),max(a,c),S,d_connections,PATHS,A)
              & list(c(section(ME[[2]],P),section(b,P))) %in% order) {
            counter <- counter+1
          }
        }
      }
    }
  }
  else {
    for (k in seq_int(1,n-2)) {
      for (ME in m_exclusive[[k]]) {
        a <- ME[[1]]
        c <- ME[[k+2]]
        for (triple in m_exclusive[[1]]) {
          if (!(list(triple) %in% perfect_nonconductors) & triple[[1]]==a & triple[[3]]==c
              & list(c(section(ME[[2]],P),section(triple[[2]],P))) %in% order) {
            counter <- counter+1
          }
        }
      }
    }
  }
  return(counter)
}



good_with_p_adjacencies <- function(n,d_connections=list(),P,order,PATHS=list(),p_adjacencies=list(),A="a") { #check if every consecutive pair of sections are connected with a p-adjacency, and every p-adjacency between two sections makes them comparable - #A = (I-C)^(-1) (I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM
  if (length(p_adjacencies)==0) {
    p_adjacencies <- p_adjacencies(n,d_connections,P,order,PATHS,A)
  }
  for (e in p_adjacencies) {
    if (!(list(c(section(e[1],P),section(e[2],P))) %in% order) & !(list(c(section(e[2],P),section(e[1],P))) %in% order)) {
      return (FALSE)
    }
  }
  for (e in order) {
    order_new <- order[-which(order %in% list(e))]
    if (is.partial.order(length(P),order_new)) {
      flag <- FALSE
      for (a in P[[e[1]]]) {
        for (b in P[[e[2]]]) {
          if (list(c(min(a,b),max(a,b))) %in% p_adjacencies) {
            flag <- TRUE
            break
          }
        }
        if (flag) {break}
      }
      if (!flag) {return(FALSE)}
    }
  }
  return(TRUE)
}

#############################################
unshielded_conductors <- function(n,d_connections=list(),P,order,PATHS=list(),A="a") {#Outputs all the unshielded conductors in the graph associated with (P,order) - A: coefficient matrix of the linear Gaussian SEM
  L <- list()
  for (a in seq_int(1,n)) {
    for (c in seq_int(a,n)) {
      if (a!=c & !p_adjacent(n,d_connections,P,order,a,c,PATHS,A)) {
        for (b in seq_int(1,n)) {
          if (((list(c(section(b,P),section(a,P))) %in% order) | (list(c(section(b,P),section(c,P))) %in% order)) & p_adjacent(n,d_connections,P,order,a,b,PATHS,A) & p_adjacent(n,d_connections,P,order,b,c,PATHS,A)) {
            L[[length(L)+1]] <- list(a,b,c) 
          }
        }
      }
    }
  }
  return(L)
}


##################################################
# Partially ordered partitions used in the running example of Section 3 of the paper:
P1 = list(c(1),c(2),c(3,4,5))
order1 = list(c(1,1),c(2,2),c(3,3),c(1,3),c(2,3))

P2 = list(c(1,2),c(3,4,5))
order2 = list(c(1,1),c(2,2),c(1,2))

P3 = list(c(1,2),c(3,4),c(5))
order3 = list(c(1,1),c(2,2),c(3,3),c(1,2))

P4 = list(c(1,2),c(3,4,5))
order4 = list(c(1,1),c(2,2),c(2,1))







