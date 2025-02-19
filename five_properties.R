source("d_connections.R")

section <- function(s,P) #Outputs i such that P[[i]] is the section of the ordered partition P containing s
{
  return(which(lapply(P,match,x=c(s),nomatch=0)!=0))
}

##############################################

p_adjacencies <- function(n,d_connections,P) #Outputs the list of p-adjacencies associated to the ordered partition P, i.e. E^(1)_P
{
  L <- list()
  for (s in seq(1,n-1)) {
    for (t in seq(s+1,n)) {
      S <- seq(1,n)[lapply(seq(1,n),section,P=P) <= max(section(s,P),section(t,P)) & seq(1,n)!=s & seq(1,n)!=t]
      S <- as.double(S)
      if (list(list(s,t,S)) %in% d_connections) {
        L[[length(L)+1]] <- as.double(c(s,t))
      }
    }
  }
  return(L)
}

##############################################

p_adjacent <- function(n,d_connections,P,s,t) #Checks whether s and t are p-adjacent with respect to P
{
  if (list(as.double(c(min(s,t),max(s,t)))) %in% p_adjacencies(n,d_connections,P)) {return(TRUE)} else {return(FALSE)}
}

##############################################

itineraries <- function(n,d_connections,P,S,k) #Outputs all the uncovered itineraries in S of length k with respect to the p-adjacencies associated with the ordered partition P
#S is a vector or list
{
  p_adjacencies <- p_adjacencies(n,d_connections,P)
  L <- list()
  if (k==1) {
    S <- unlist(S)
    for (a in S) {
      L[[length(L)+1]] <- list(a)
    }
    return(L)
  } 
  else {
    for (I in itineraries(n,d_connections,P,S,k-1)) {
      for (a in S) {
        if ((p_adjacent(n,d_connections,P,a,I[[1]])) & !(TRUE %in% lapply(I[seq_int(2,k-1)],p_adjacent,n=n,d_connections=d_connections,P=P,s=a)) & !(a %in% I)) {
          L[[length(L)+1]] <- c(list(a),I)
        }
      }
    }
    return(L)
  }
}

##############################################

mutually_exclusive <- function(n,d_connections,P,k) #Outputs the list of mutually exclusive conductors <x_0,...,x_{k+1}> associated to the ordered partition P 
{
  L <- list()
  for(i in 1:length(P)) {
    for (I in itineraries(n,d_connections,P,P[[i]],k)) {
      for (a in 1:n) {
        if (section(a,P)<i & p_adjacent(n,d_connections,P,a,I[[1]]) & !(TRUE %in% lapply(I[seq_int(2,k)],p_adjacent,n=n,d_connections=d_connections,P=P,s=a))) {
          for (b in 1:n) {
            if (section(b,P)<i & p_adjacent(n,d_connections,P,b,I[[k]]) & !(TRUE %in% lapply(I[seq_int(1,k-1)],p_adjacent,n=n,d_connections=d_connections,P=P,s=b)) & a!=b & !p_adjacent(n,d_connections,P,a,b)) {
              L [[length(L)+1]] <- c(list(a),I,list(b))
            }
          }
        }
      }
    }
  }
  return(L)
}

##############################################

perfect_nonconductors <- function(n,d_connections,P) {#Outputs all the unshielded perfect nonconductors in the graph associated with P
  L <- list()
  for (a in seq(1,n)) {
    for (c in seq(1,n)) {
      if(a!=c & !p_adjacent(n,d_connections,P,a,c)) {
        for (b in seq(1,n)) {
          S <- seq(1,n)[lapply(seq(1,n),section,P=P) <= max(section(a,P),section(b,P),section(c,P)) & seq(1,n)!=a & seq(1,n)!=c]
          S <- as.double(S)
          if (p_adjacent(n,d_connections,P,a,b) & p_adjacent(n,d_connections,P,b,c) & list(list(min(a,c),max(a,c),S)) %in% d_connections) {
            L[[length(L)+1]] <- list(a,b,c) 
          }
        }
      }
    }
  }
  return(L)
}

##############################################

E4 <- function(n,d_connections,P) #Property 4 of Richardson's Theorem
{
  counter <- 0
  for (a in seq(1,n)) {
    for (c in seq(1,n)) {
      if (a!=c & !p_adjacent(n,d_connections,P,a,c)) {
        for (b1 in seq(1,n)) {
          S <- seq(1,n)[lapply(seq(1,n),section,P=P) <= max(section(a,P),section(b1,P),section(c,P)) & seq(1,n)!=a & seq(1,n)!=c]
          S <- as.double(S)
          if (p_adjacent(n,d_connections,P,a,b1) & p_adjacent(n,d_connections,P,c,b1)
              & max(section(a,P),section(c,P)) < section(b1,P) & !(list(list(min(a,c),max(a,c),S)) %in% d_connections)) {
            for (b2 in seq(1,n)) {
              S <- seq(1,n)[lapply(seq(1,n),section,P=P) <= max(section(a,P),section(b2,P),section(c,P)) & seq(1,n)!=a & seq(1,n)!=c]
              S <- as.double(S)
              if (b1!=b2 & p_adjacent(n,d_connections,P,a,b2) & p_adjacent(n,d_connections,P,c,b2)
                  & max(section(a,P),section(c,P)) < section(b2,P) & !(list(list(min(a,c),max(a,c),S)) %in% d_connections)
                  & section(b1,P)<=section(b2,P)) {
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

##############################################

E6 <- function(n,d_connections,P) { #Property 5 of Richardson's Theorem
  counter <- 0
  for (k in seq(1,n)) {
    for (ME in mutually_exclusive(n,d_connections,P,k)) {
      a <- ME[[1]]
      c <- ME[[k+2]]
      for (b in seq(1,n)) {
        S <- seq(1,n)[lapply(seq(1,n),section,P=P) <= max(section(a,P),section(b,P),section(c,P)) & seq(1,n)!=a & seq(1,n)!=c]
        S <- as.double(S)
        if (p_adjacent(n,d_connections,P,a,b) & p_adjacent(n,d_connections,P,b,c) 
            & max(section(a,P),section(c,P))<section(b,P)
            & !(list(list(min(a,c),max(a,c),S)) %in% d_connections & a!=c)
            & section(ME[[2]],P)<=section(b,P)) {
          counter <- counter+1
        }
      }
    }
  }
  return(counter)
}

