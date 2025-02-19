# source("d_connections.R")
# source("paths.R")
# source("test.R")
# source("greedy_optimization_partial_orders.R")
library(cIRT)




section <- function(s,P) #Outputs i such that P[[i]] is the section of the ordered partition P containing s
{
  return(as.double(which(lapply(P,match,x=c(s),nomatch=0)!=0)))
}

##############################################

p_adjacencies_graph <- function(n,E,P,order) #Outputs the list of p-adjacencies in the graph (G,[n]) whose associated partially ordered partition is (P,order).
{
  L <- list()
  for (s in seq_int(1,n-1)) {
    for (t in seq_int(s+1,n)) {
      if (list(c(s,t)) %in% E || list(c(t,s)) %in% E) {L[[length(L)+1]] <- as.double(c(s,t))}
      else {
        for (z in seq_int(1,n)) {
          if(z!=s && z!=t && list(c(s,z)) %in% E && list(c(t,z)) %in% E && (section(z,P)==section(s,P) || section(z,P)==section(t,P))) {
            L[[length(L)+1]] <- as.double(c(s,t))
            break
          }
        }
      }
    }
  }
  return(L)
}

##############################################


itineraries_graph <- function(n,E,P,order,S,k,kminus1.itineraries="a",p_adjacencies) #Outputs all the uncovered itineraries in S of length k in ([n],E), whose associated partially ordered partition is (P,order).
  #S is a vector or list
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
    if (is.character(kminus1.itineraries)) {kminus1.itineraries <- itineraries_graph(n=n,E=E,P=P,order=order,S=S,k=k-1,p_adjacencies=p_adjacencies)}
    for (I in kminus1.itineraries) {
      for (a in S) {
        if (list(c(min(a,I[[1]]),max(a,I[[1]]))) %in% p_adjacencies & !(a %in% I)) {
          flag <- TRUE
          for (i in seq_int(2,k-1)) {
            if (list(c(min(a,I[[i]]),max(a,I[[i]]))) %in% p_adjacencies) {
              flag <- FALSE
              break
            }
          }
          if (flag) {L[[length(L)+1]] <- c(list(a),I)}
        }
      }
    }
    return(L)
  }
}

##############################################

mutually_exclusive_efficient_graph <- function(n,E,P,order,k,kminus1.itineraries="a",p_adjacencies) #Outputs the list of mutually exclusive conductors <x_0,...,x_{k+1}> in ([n],E) with associated partially ordered partition (P,order)
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
              if (is.character(kminus1.itineraries)) {it <- itineraries_graph(n,E,P,order,P[[i]],k,kminus1.itineraries,p_adjacencies)}
              else {it <- itineraries_graph(n,E,P,order,P[[i]],k,kminus1.itineraries[[i]],p_adjacencies)}
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

perfect_nonconductors_efficient_graph <- function(n,E,P,order,non_conductors) {#Outputs all the unshielded perfect nonconductors in the graph ([n],E) with associated partially ordered partiton (P,order) 
  L <- list()
  for (triple in non_conductors) {
    a <- triple[[1]]
    b <- triple [[2]]
    cc <- triple [[3]]
    for (z in seq_int(1,n)) {
      if (list(c(a,z)) %in% E && list(c(cc,z)) %in% E && list(c(section(z,P), section(b,P))) %in% order) {
        L[[length(L)+1]] <- list(a,b,cc) 
        break
      }
    }
  }
  return(L)
}
##############################################


E4_efficient_graph <- function(n,E,P,order,nonconductors,perfect_nonconductors) {
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

E6_graph <- function(n,E,P,order,m_exclusive=list(),perfect_nonconductors) {#Property 6 of Richardson's Theorem 
  counter <- 0
  if (length(m_exclusive)==0) {
    for (k in seq_int(1,n-2)) {
      for (ME in mutually_exclusive(n=n,d_connections=d_connections,P=P,order=order,k=k,PATHS=PATHS,A)) {
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


##############################################

true_score_graph <- function(n,E,P,order) { # Outputs the optimal score of the graph ([n],E) whose associated partially ordered partition is (P,order).
  base_new <- c()
  P_ADJ <- p_adjacencies_graph(n,E,P,order)
  base_new[1] <- length(P_ADJ)
  M <- mutually_exclusive_efficient_graph(n=n,E=E,P=P,order=order,k=1,p_adjacencies=P_ADJ)
  m_exclusive <- list()
  m_exclusive[[1]] <- M[[1]]
  it <- M[[2]]
  base_new[2] <- length(m_exclusive[[1]])
  PN <- perfect_nonconductors_efficient_graph(n,E,P,order,m_exclusive[[1]])
  base_new[3] <- length(PN)
  base_new[4] <- E4_efficient_graph(n,E,P,order,m_exclusive[[1]],PN) 
  for (t in seq_int(2,n-2)) {
    M <-  mutually_exclusive_efficient_graph(n,E,P,order,t,it,P_ADJ)
    m_exclusive[[t]] <- M[[1]]
    it <- M[[2]]
    base_new[t+3] <- length(m_exclusive[[t]])
  }
  base_new[n+2] <- E6_graph(n,E,P,order,m_exclusive,PN)
  return(base_new)
}




######################### test
n = 7
p = 0.6
#defining the edge set E
g <- erdos.renyi.game(n, p, type = "gnp",directed=TRUE)
E_mat <- get.edgelist(g)
E <- list()
for (i in seq_int(1,dim(E_mat)[1])) {E[[i]] <- c(E_mat[i,])}


#finding the right partially ordered partition 
PPP <- find_partially_ordered_partition(n,E)
P <- PPP[[1]]
order <- PPP[[2]]

#finding all the paths in the graph
M <- all_paths_separated(n,E)


#finding the optimal score
base <- c(Inf,-1,Inf,Inf,rep(-1,n-3),Inf)
compare_posets(n,list(),P,order,base,M,"a")

true_score_graph(n,E,P,order)
# Result: This code works!!