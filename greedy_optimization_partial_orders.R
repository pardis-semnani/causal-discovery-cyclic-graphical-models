source("six_properties.R")
source("d_connections.R")


compare_posets <- function(n,d_connections=list(),P,order,base,PATHS=list(),A="a") {
  #A = (I-C)^(-1)(I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM 
  #0 when (P,order) is increasing compared to base, 
  #base_new when (P,order) is decreasing compared to base,
  #1 when base_new is the same as base.
  flag <- FALSE #true if (P,order) has a strictly better score than base
  base_new <- c()
  P_ADJ <- p_adjacencies(n,d_connections,P,order,PATHS,A)
  base_new[1] <- length(P_ADJ)
  if (base_new[1] > base[1]) {
    return(0)
  }
  # else if (!good_with_p_adjacencies(n,d_connections,P,order,PATHS,P_ADJ)) {
  #   return(0)
  # }
  else if (base_new[1] < base[1]) {
    flag <- TRUE
  }
  M <- mutually_exclusive_efficient(n=n,d_connections=d_connections,P=P,order=order,k=1,PATHS=PATHS,p_adjacencies=P_ADJ,A=A)
  m_exclusive <- list()
  m_exclusive[[1]] <- M[[1]]
  it <- M[[2]]
  base_new[2] <- length(m_exclusive[[1]])
  if (!flag & base_new[2] < base[2]) {
    return(0)
  }
  else if (!flag & base_new[2] > base[2]) {
    flag <- TRUE
  }
  PN <- perfect_nonconductors_efficient(n,d_connections,P,order,PATHS,m_exclusive[[1]],A)
  base_new[3] <- length(PN)
  if (!flag & base_new[3] > base[3]) {
    return(0)
  }
  else if (!flag & base_new[3] < base[3]) {
    flag <- TRUE
  }
  base_new[4] <- E4_efficient(n,d_connections,P,order,PATHS,m_exclusive[[1]],PN,A) 
  if (!flag & base_new[4] > base[4]) {
    return(0)
  }
  if (!flag & base_new[4] < base[4]) {
    flag <- TRUE
  }
  for (t in seq_int(2,n-2)) {
    M <-  mutually_exclusive_efficient(n,d_connections,P,order,t,it,PATHS,P_ADJ,A)
    m_exclusive[[t]] <- M[[1]]
    it <- M[[2]]
    base_new[t+3] <- length(m_exclusive[[t]])
    if (!flag & base_new[t+3] < base[t+3]) {
      return(0)
    }
    if (!flag & base_new[t+3] > base[t+3]) {
      flag <- TRUE
    }
  }
  base_new[n+2] <- E6(n,d_connections,P,order,m_exclusive,PN,PATHS,A)
  if (!flag & base_new[n+2] > base[n+2]) {
    return(0)
  }
  else if (!flag & base_new[n+2] < base[n+2]) {
    flag <- TRUE
  }
  if (flag) {
    return(base_new)
  }
  return(1)
}


#The following function has not been updated with PATHS and A options.
compare_posets_total_orders <- function(n,d_connections,P,order,base) {
  #0 when (P,order) is increasing compared to base, 
  #base_new when (P,order) is decreasing compared to base,
  #1 when base_new is the same as base.
  flag <- FALSE #true if (P,order) has a strictly better score than base
  base_new <- c()
  base_new[1] <- length(p_adjacencies(n,d_connections,P,order))
  if (base_new[1] > base[1]) {
    return(0)
  }
  else if (base_new[1] < base[1]) {
    flag <- TRUE
  }
  M <- mutually_exclusive(n,d_connections,P,order,1)
  it <- M[[2]]
  base_new[2] <- length(M[[1]])
  if (!flag & base_new[2] < base[2]) {
    return(0)
  }
  else if (!flag & base_new[2] > base[2]) {
    flag <- TRUE
  }
  for (t in seq_int(2,n-2)) {
    M <-  mutually_exclusive(n,d_connections,P,order,t,it)
    it <- M[[2]]
    base_new[t+1] <- length(M[[1]])
    if (!flag & base_new[t+1] < base[t+1]) {
      return(0)
    }
    if (!flag & base_new[t+1] > base[t+1]) {
      flag <- TRUE
    }
  }
  if (flag) {
    return(base_new)
  }
  return(1)
}


is.partial.order <- function(n,order) {
  for (a in seq_int(1,n)) {
    if (!(list(c(a,a)) %in% order)) {
      return(FALSE)
    }
  }
  for (a in seq_int(1,n)) {
    for (b in seq_int(1,n)) {
      for (c in seq_int(1,n)) {
        if (list(c(a,b)) %in% order & list(c(b,c)) %in% order & !(list(c(a,c)) %in% order)) {
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)
}

has.been.checked <- function(n,P,order,K) {#checks whether (P,order) is equivalent to any of the partially ordered partitions in K
  for (i in seq_int(1,length(K))) {
    flag <- TRUE #true if (P1,order1) is equivalent to (P,order)
    P1 <- K[[i]][[1]]
    order1 <- K[[i]][[2]]
    for (a in seq_int(1,n)) {
      for (b in seq_int(1,n)) {
        if (list(c(section(a,P1),section(b,P1))) %in% order1 & !(list(c(section(a,P),section(b,P))) %in% order)) {
          flag <- FALSE
          break
        }
        else if (list(c(section(a,P),section(b,P))) %in% order & !(list(c(section(a,P1),section(b,P1))) %in% order1)) {
          flag <- FALSE
          break
        }
      }
      if (!flag) {break}
    }
    if (flag) {return(TRUE)}
  }
  return(FALSE)
}




################################################
#greedy optimization that finds a path with at most N partially ordered partitions with the same score.
greedy_optimization_partial_orders_depth_first <- function(n,d_connections=list(),P0,order0,N=Inf,L=list(),PATHS=list(),A="a") {#L: list of partially ordered partitions to avoid
  # PATHS: list of all st-separated paths in the original graph in the form of a matrix / #A = (I-C)^(-1)(I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM  / Exactly one of PATHS or d_connections or A should be given to the function
  # Track <- list()
  K <- c(L,list(list(P0,order0))) #list of checked partially ordered partitions
  base <- c(Inf,-1,Inf,Inf,rep(-1,n-3),Inf)
  print("hello1")
  base <- compare_posets(n,d_connections,P0,order0,base,PATHS,A)
  II <- c(length(K)) #a vector consisting of the (index of) elements in K that have equal score to the current (P,order)
  I <- c(length(K)) #a vector consisting of the (index of) elements in K that have equal score to the current (P,order) in the current path
  # L <- list(list(P0,order0)) #list of equal partially ordered partitions
  # i <- 1
  flag <- TRUE
  P <- P0
  order <- order0
  # counter <- 0
  print("hello")
  while (TRUE) {
    while (flag & length(I)<=N) {
      # if (i==1) {print(base)}
      print(length(K))
      print(base)
      print(c("counter",length(I)))
      # print(K[[length(K)]])
      flag <- FALSE #if a non-increasing neighbor is found
      # P <- L[[i]][[1]]
      # order <- L[[i]][[2]]
      cons <- list() #list of consecutive pairs of sections in order
      # Track [[length(Track)+1]] <- L[[i]]
      
      #finding neighbors by changing the partial order
      for (a in seq_int(1,length(P))) {
        for (b in seq_int(1,length(P))) {
          #removing from partial order
          if (a!=b & list(c(a,b)) %in% order) {
            order_new <- order[-which(order %in% list(c(a,b)))]
            if (is.partial.order(length(P),order_new)) {
              cons [[length(cons)+1]] <- c(a,b)
              if(!has.been.checked(n,P,order_new,K)) {
                # if (!(list(list(P,order_new)) %in% K)) {
                base_new <- compare_posets(n,d_connections,P,order_new,base,PATHS,A)
                K [[length(K)+1]] <- list(P,order_new)
                if (length(base_new) > 1) {
                  I <- c(length(K))
                  II <- c(length(K))
                  base <- base_new
                  order <- order_new
                  flag <- TRUE
                  break
                }
                else if (base_new==1) {
                  I [length(I)+1] <- length(K)
                  II [length(II)+1] <- length(K)
                  order <- order_new
                  flag <- TRUE
                  break
                }
              }
            }
          }
          #adding to partial order
          else if (a!=b & !(list(c(b,a)) %in% order)) {
            order_new <- append(order,list(c(a,b)))
            if (is.partial.order(length(P),order_new) & !has.been.checked(n,P,order_new,K)) {
              base_new <- compare_posets(n,d_connections,P,order_new,base,PATHS,A)
              K [[length(K)+1]] <- list(P,order_new)
              if (length(base_new) > 1) {
                I <- c(length(K))
                II <- c(length(K))
                base <- base_new
                order <- order_new
                flag <- TRUE
                break
              }
              else if (base_new==1) {
                I[length(I)+1] <- length(K)
                II[length(II)+1] <- length(K)
                order <- order_new
                flag <- TRUE
                break
              }
            }
          }
        }
        if (flag) {break}
      }
      if (flag) {next}
      
      #########################################################
      #add one element of one section to a consecutive section
      for (pair in cons) {
        a <- pair[1]; b <- pair[2]
        #adding an element of P[[b]] to P[[a]]
        for (k in seq_int(1,length(P[[b]]))) {
          P_new <- P
          P_new[[b]] <- P[[b]][-k]
          P_new[[a]] <- sort(c(P[[a]],P[[b]][k]))
          order_new <- order
          order_new1 <- order_new
          if (length(P_new[[b]])==0) {
            for (pair in order_new1) {
              if (pair[1]==b | pair[2]==b) {
                order_new <- order_new[-which(order_new %in% list(pair))]
              }
            }
            for (j in seq_int(1,length(order_new))) {
              pair <- order_new[[j]]
              if (pair[1]>b & pair[2]>b) {
                order_new[[j]] <- c(pair[1]-1,pair[2]-1)
              }
              else if (pair[1]>b) {
                order_new[[j]] <- c(pair[1]-1,pair[2])
              }
              else if (pair[2]>b) {
                order_new[[j]] <- c(pair[1],pair[2]-1)
              }
            }
            P_new[[b]] <- NULL
          }
          if(!has.been.checked(n,P_new,order_new,K)) {
            # if (!(list(list(P_new,order_new)) %in% K)) {
            base_new <- compare_posets(n,d_connections,P_new,order_new,base,PATHS,A)
            K [[length(K)+1]] <- list(P_new,order_new)
            if (length(base_new) > 1) {
              I <- c(length(K))
              II <- c(length(K))
              base <- base_new
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
            else if (base_new==1) {
              I[length(I)+1] <- length(K)
              II[length(II)+1] <- length(K)
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
          }
        }
        if (flag) {break}
        
        #adding an element of P[[a]] to P[[b]]
        for (k in seq_int(1,length(P[[a]]))) {
          P_new <- P
          P_new[[a]] <- P[[a]][-k]
          P_new[[b]] <- sort(c(P[[b]],P[[a]][k]))
          order_new <- order
          order_new1 <- order_new
          if (length(P_new[[a]])==0) {
            for (pair in order_new1) {
              if (pair[1]==a | pair[2]==a) {
                order_new <- order_new[-which(order_new %in% list(pair))]
              }
            }
            for (j in seq_int(1,length(order_new))) {
              pair <- order_new[[j]]
              if (pair[1]>a & pair[2]>a) {
                order_new[[j]] <- c(pair[1]-1,pair[2]-1)
              }
              else if (pair[1]>a) {
                order_new[[j]] <- c(pair[1]-1,pair[2])
              }
              else if (pair[2]>a) {
                order_new[[j]] <- c(pair[1],pair[2]-1)
              }
            }
            P_new[[a]] <- NULL
          }
          if(!has.been.checked(n,P_new,order_new,K)) {
            # if (!(list(list(P_new,order_new)) %in% K)) {
            base_new <- compare_posets(n,d_connections,P_new,order_new,base,PATHS,A)
            K [[length(K)+1]] <- list(P_new,order_new)
            if (length(base_new) > 1) {
              I <- c(length(K))
              II <- c(length(K))
              base <- base_new
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
            else if (base_new==1) {
              I [length(I)+1] <- length(K)
              II [length(II)+1] <- length(K)
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
          }
        }
        if (flag) {break}
      }
      if (flag) {next}
      ################################################
      #add one element of one section to an empty set 
      for (a in seq_int(1,length(P))) {
        #add the element as a section right after the current section
        for (k in seq_int(1,length(P[[a]]))) {
          P_new <- c(P[seq_int(1,a-1)],list(P[[a]][-k]),P[seq_int(a+1,length(P))],list(c(P[[a]][k])))
          order_new <- append(order,list(c(as.double(length(P_new)),as.double(length(P_new)))))
          for (b in seq_int(1,length(P))) {
            if (list(c(b,a)) %in% order) {order_new <- append(order_new,list(c(as.double(b),as.double(length(P_new)))))}
            if (list(c(a,b)) %in% order & b!=a) {order_new <- append(order_new,list(c(as.double(length(P_new)),as.double(b))))}
          }
          order_new1 <- order_new
          if (length(P_new[[a]])==0) {
            # if(i==92 & length(L)==123) {print("order_new1 is");print(order_new1)}
            for (pair in order_new1) {
              # if(i==92 & length(L)==123) {print("hello");print(pair)}
              if (pair[1]==a | pair[2]==a) {
                order_new <- order_new[-which(order_new %in% list(pair))]
              }
            }
            for (j in seq_int(1,length(order_new))) {
              pair <- order_new[[j]]
              if (pair[1]>a & pair[2]>a) {
                order_new[[j]] <- c(pair[1]-1,pair[2]-1)
              }
              else if (pair[1]>a) {
                order_new[[j]] <- c(pair[1]-1,pair[2])
              }
              else if (pair[2]>a) {
                order_new[[j]] <- c(pair[1],pair[2]-1)
              }
              if(i==92 & length(L)==123) {print("bye");print(order_new)}
            }
            P_new[[a]] <- NULL
          }
          if(!has.been.checked(n,P_new,order_new,K)) {
            # if (!(list(list(P_new,order_new)) %in% K)) {
            base_new <- compare_posets(n,d_connections,P_new,order_new,base,PATHS,A)
            K [[length(K)+1]] <- list(P_new,order_new)
            if (length(base_new) > 1) {
              I <- c(length(K))
              II <- c(length(K))
              base <- base_new
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
            else if (base_new==1) {
              I [length(I)+1] <- length(K)
              II [length(II)+1] <- length(K)
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
          }
        }
        if (flag) {break}
        
        #add the element as a section right before the current section
        for (k in seq_int(1,length(P[[a]]))) {
          P_new <- c(P[seq_int(1,a-1)],list(P[[a]][-k]),P[seq_int(a+1,length(P))],list(c(P[[a]][k])))
          order_new <- append(order,list(c(as.double(length(P_new)),as.double(length(P_new)))))
          for (b in seq_int(1,length(P))) {
            if (list(c(b,a)) %in% order & b!=a) {order_new <- append(order_new,list(c(as.double(b),as.double(length(P_new)))))}
            if (list(c(a,b)) %in% order) {order_new <- append(order_new,list(c(as.double(length(P_new)),as.double(b))))}
          }
          order_new1 <- order_new
          if (length(P_new[[a]])==0) {
            for (pair in order_new1) {
              if (pair[1]==a | pair[2]==a) {
                order_new <- order_new[-which(order_new %in% list(pair))]
              }
            }
            for (j in seq_int(1,length(order_new))) {
              pair <- order_new[[j]]
              if (pair[1]>a & pair[2]>a) {
                order_new[[j]] <- c(pair[1]-1,pair[2]-1)
              }
              else if (pair[1]>a) {
                order_new[[j]] <- c(pair[1]-1,pair[2])
              }
              else if (pair[2]>a) {
                order_new[[j]] <- c(pair[1],pair[2]-1)
              }
            }
            P_new[[a]] <- NULL
          }
          if(!has.been.checked(n,P_new,order_new,K)) {
            # if (!(list(list(P_new,order_new)) %in% K)) {
            base_new <- compare_posets(n,d_connections,P_new,order_new,base,PATHS,A)
            K [[length(K)+1]] <- list(P_new,order_new)
            if (length(base_new) > 1) {
              I <- c(length(K))
              II <- c(length(K))
              base <- base_new
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
            else if (base_new==1) {
              I [length(I)+1] <- length(K)
              II [length(II)+1] <- length(K)
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
          }
        }
        if (flag) {break}
      }
      if (flag) {next}
    }
    if (length(I) > N) {return(list(K,II))}
    I <- I[-length(I)]
    if (length(I)==0) {break}
    P <- K[[I[length(I)]]][[1]]
    order <- K[[I[length(I)]]][[2]]
    flag <- TRUE
  }
  return(list(K,II))
}

########################################################
greedy_optimization_partial_orders_depth_first_modified <- function(n,d_connections=list(),P0,order0,N=Inf,L=list(),PATHS=list(),A="a") {#L: list of partially ordered partitions to avoid
  # PATHS: list of all st-separated paths in the original graph in the form of a matrix / #A = (I-C)^(-1)(I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM  / Exactly one of PATHS or d_connections or A should be given to the function
  # Track <- list()
  K <- c(L,list(list(P0,order0))) #list of checked partially ordered partitions
  base <- c(Inf,-1,Inf,Inf,rep(-1,n-3),Inf)
  print("hello1")
  base <- compare_posets(n,d_connections,P0,order0,base,PATHS,A)
  II <- c(length(K)) #a vector consisting of the (index of) elements in K that have equal score to the current (P,order)
  I <- c(length(K)) #a vector consisting of the (index of) elements in K that have equal score to the current (P,order) in the current path
  # L <- list(list(P0,order0)) #list of equal partially ordered partitions
  # i <- 1
  flag <- TRUE
  P <- P0
  order <- order0
  # counter <- 0
  print("hello")
  while (TRUE) {
    while (flag & length(I)<=N) {
      # if (i==1) {print(base)}
      print(length(K))
      print(base)
      print(c("counter",length(I)))
      # print(K[[length(K)]])
      flag <- FALSE #if a non-increasing neighbor is found
      # P <- L[[i]][[1]]
      # order <- L[[i]][[2]]
      cons <- list() #list of consecutive pairs of sections in order
      # Track [[length(Track)+1]] <- L[[i]]
      
      #finding neighbors by changing the partial order
      for (a in seq_int(1,length(P))) {
        for (b in seq_int(1,length(P))) {
          #removing from partial order
          if (a!=b & list(c(a,b)) %in% order) {
            order_new <- order[-which(order %in% list(c(a,b)))]
            if (is.partial.order(length(P),order_new)) {
              cons [[length(cons)+1]] <- c(a,b)
              if(!has.been.checked(n,P,order_new,K)) {
                # if (!(list(list(P,order_new)) %in% K)) {
                base_new <- compare_posets(n,d_connections,P,order_new,base,PATHS,A)
                K [[length(K)+1]] <- list(P,order_new)
                if (length(base_new) > 1) {
                  I <- c(length(K))
                  II <- c(length(K))
                  base <- base_new
                  order <- order_new
                  flag <- TRUE
                  break
                }
                else if (base_new==1) {
                  I [length(I)+1] <- length(K)
                  II [length(II)+1] <- length(K)
                  order <- order_new
                  flag <- TRUE
                  break
                }
              }
            }
          }
          #adding to partial order
          else if (a!=b & !(list(c(b,a)) %in% order)) {
            order_new <- append(order,list(c(a,b)))
            if (is.partial.order(length(P),order_new) & !has.been.checked(n,P,order_new,K)) {
              base_new <- compare_posets(n,d_connections,P,order_new,base,PATHS,A)
              K [[length(K)+1]] <- list(P,order_new)
              if (length(base_new) > 1) {
                I <- c(length(K))
                II <- c(length(K))
                base <- base_new
                order <- order_new
                flag <- TRUE
                break
              }
              else if (base_new==1) {
                I[length(I)+1] <- length(K)
                II[length(II)+1] <- length(K)
                order <- order_new
                flag <- TRUE
                break
              }
            }
          }
        }
        if (flag) {break}
      }
      if (flag) {next}
      
      #########################################################
      #add one element of one section to a consecutive section
      for (pair in cons) {
        a <- pair[1]; b <- pair[2]
        #adding an element of P[[b]] to P[[a]]
        for (k in seq_int(1,length(P[[b]]))) {
          P_new <- P
          P_new[[b]] <- P[[b]][-k]
          P_new[[a]] <- sort(c(P[[a]],P[[b]][k]))
          order_new <- order
          order_new1 <- order_new
          if (length(P_new[[b]])==0) {
            for (pair in order_new1) {
              if (pair[1]==b | pair[2]==b) {
                order_new <- order_new[-which(order_new %in% list(pair))]
              }
            }
            for (j in seq_int(1,length(order_new))) {
              pair <- order_new[[j]]
              if (pair[1]>b & pair[2]>b) {
                order_new[[j]] <- c(pair[1]-1,pair[2]-1)
              }
              else if (pair[1]>b) {
                order_new[[j]] <- c(pair[1]-1,pair[2])
              }
              else if (pair[2]>b) {
                order_new[[j]] <- c(pair[1],pair[2]-1)
              }
            }
            P_new[[b]] <- NULL
          }
          if(!has.been.checked(n,P_new,order_new,K)) {
            # if (!(list(list(P_new,order_new)) %in% K)) {
            base_new <- compare_posets(n,d_connections,P_new,order_new,base,PATHS,A)
            K [[length(K)+1]] <- list(P_new,order_new)
            if (length(base_new) > 1) {
              I <- c(length(K))
              II <- c(length(K))
              base <- base_new
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
            else if (base_new==1) {
              I[length(I)+1] <- length(K)
              II[length(II)+1] <- length(K)
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
          }
        }
        if (flag) {break}
        
        #adding an element of P[[a]] to P[[b]]
        for (k in seq_int(1,length(P[[a]]))) {
          P_new <- P
          P_new[[a]] <- P[[a]][-k]
          P_new[[b]] <- sort(c(P[[b]],P[[a]][k]))
          order_new <- order
          order_new1 <- order_new
          if (length(P_new[[a]])==0) {
            for (pair in order_new1) {
              if (pair[1]==a | pair[2]==a) {
                order_new <- order_new[-which(order_new %in% list(pair))]
              }
            }
            for (j in seq_int(1,length(order_new))) {
              pair <- order_new[[j]]
              if (pair[1]>a & pair[2]>a) {
                order_new[[j]] <- c(pair[1]-1,pair[2]-1)
              }
              else if (pair[1]>a) {
                order_new[[j]] <- c(pair[1]-1,pair[2])
              }
              else if (pair[2]>a) {
                order_new[[j]] <- c(pair[1],pair[2]-1)
              }
            }
            P_new[[a]] <- NULL
          }
          if(!has.been.checked(n,P_new,order_new,K)) {
            # if (!(list(list(P_new,order_new)) %in% K)) {
            base_new <- compare_posets(n,d_connections,P_new,order_new,base,PATHS,A)
            K [[length(K)+1]] <- list(P_new,order_new)
            if (length(base_new) > 1) {
              I <- c(length(K))
              II <- c(length(K))
              base <- base_new
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
            else if (base_new==1) {
              I [length(I)+1] <- length(K)
              II [length(II)+1] <- length(K)
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
          }
        }
        if (flag) {break}
      }
      if (flag) {next}
      ################################################
      #add one element of one section to an empty set 
      for (a in seq_int(1,length(P))) {
        #add the element as a section right after the current section
        for (k in seq_int(1,length(P[[a]]))) {
          P_new <- c(P[seq_int(1,a-1)],list(P[[a]][-k]),P[seq_int(a+1,length(P))],list(c(P[[a]][k])))
          order_new <- append(order,list(c(as.double(length(P_new)),as.double(length(P_new)))))
          for (b in seq_int(1,length(P))) {
            if (list(c(b,a)) %in% order) {order_new <- append(order_new,list(c(as.double(b),as.double(length(P_new)))))}
            if (list(c(a,b)) %in% order & b!=a) {order_new <- append(order_new,list(c(as.double(length(P_new)),as.double(b))))}
          }
          order_new1 <- order_new
          if (length(P_new[[a]])==0) {
            # if(i==92 & length(L)==123) {print("order_new1 is");print(order_new1)}
            for (pair in order_new1) {
              # if(i==92 & length(L)==123) {print("hello");print(pair)}
              if (pair[1]==a | pair[2]==a) {
                order_new <- order_new[-which(order_new %in% list(pair))]
              }
            }
            for (j in seq_int(1,length(order_new))) {
              pair <- order_new[[j]]
              if (pair[1]>a & pair[2]>a) {
                order_new[[j]] <- c(pair[1]-1,pair[2]-1)
              }
              else if (pair[1]>a) {
                order_new[[j]] <- c(pair[1]-1,pair[2])
              }
              else if (pair[2]>a) {
                order_new[[j]] <- c(pair[1],pair[2]-1)
              }
              if(i==92 & length(L)==123) {print("bye");print(order_new)}
            }
            P_new[[a]] <- NULL
          }
          if(!has.been.checked(n,P_new,order_new,K)) {
            # if (!(list(list(P_new,order_new)) %in% K)) {
            base_new <- compare_posets(n,d_connections,P_new,order_new,base,PATHS,A)
            K [[length(K)+1]] <- list(P_new,order_new)
            if (length(base_new) > 1) {
              I <- c(length(K))
              II <- c(length(K))
              base <- base_new
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
            else if (base_new==1) {
              I [length(I)+1] <- length(K)
              II [length(II)+1] <- length(K)
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
          }
        }
        if (flag) {break}
        
        #add the element as a section right before the current section
        for (k in seq_int(1,length(P[[a]]))) {
          P_new <- c(P[seq_int(1,a-1)],list(P[[a]][-k]),P[seq_int(a+1,length(P))],list(c(P[[a]][k])))
          order_new <- append(order,list(c(as.double(length(P_new)),as.double(length(P_new)))))
          for (b in seq_int(1,length(P))) {
            if (list(c(b,a)) %in% order & b!=a) {order_new <- append(order_new,list(c(as.double(b),as.double(length(P_new)))))}
            if (list(c(a,b)) %in% order) {order_new <- append(order_new,list(c(as.double(length(P_new)),as.double(b))))}
          }
          order_new1 <- order_new
          if (length(P_new[[a]])==0) {
            for (pair in order_new1) {
              if (pair[1]==a | pair[2]==a) {
                order_new <- order_new[-which(order_new %in% list(pair))]
              }
            }
            for (j in seq_int(1,length(order_new))) {
              pair <- order_new[[j]]
              if (pair[1]>a & pair[2]>a) {
                order_new[[j]] <- c(pair[1]-1,pair[2]-1)
              }
              else if (pair[1]>a) {
                order_new[[j]] <- c(pair[1]-1,pair[2])
              }
              else if (pair[2]>a) {
                order_new[[j]] <- c(pair[1],pair[2]-1)
              }
            }
            P_new[[a]] <- NULL
          }
          if(!has.been.checked(n,P_new,order_new,K)) {
            # if (!(list(list(P_new,order_new)) %in% K)) {
            base_new <- compare_posets(n,d_connections,P_new,order_new,base,PATHS,A)
            K [[length(K)+1]] <- list(P_new,order_new)
            if (length(base_new) > 1) {
              I <- c(length(K))
              II <- c(length(K))
              base <- base_new
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
            else if (base_new==1) {
              I [length(I)+1] <- length(K)
              II [length(II)+1] <- length(K)
              P <- P_new
              order <- order_new
              flag <- TRUE
              break
            }
          }
        }
        if (flag) {break}
      }
      if (flag) {next}
    }
    I <- I[-length(I)]
    if (length(I)==0) {break}
    P <- K[[I[length(I)]]][[1]]
    order <- K[[I[length(I)]]][[2]]
    flag <- TRUE
  }
  return(list(K,II))
}











#Not updated with PATHS and A options
greedy_optimization_partial_orders <- function(n,d_connections,P0,order0) {
  # Track <- list()
  K <- list(list(P0,order0)) #list of checked partially ordered partitions
  base <- c(Inf,-1,Inf,Inf,rep(-1,n-3),Inf)
  base <- compare_posets(n,d_connections,P0,order0,base)
  L <- list(list(P0,order0)) #list of equal partially ordered partitions
  i <- 1
  print("hello")
  while (length(L) >= i) {
    if (i==1) {print(base)}
    print(length(L))
    flag <- FALSE #if a decreasing neighbor is found
    P <- L[[i]][[1]]
    order <- L[[i]][[2]]
    cons <- list() #list of consecutive pairs of sections in order
    # Track [[length(Track)+1]] <- L[[i]]
    
    #finding neighbors by changing the partial order
    for (a in seq_int(1,length(P))) {
      for (b in seq_int(1,length(P))) {
        b <- 
          #removing from partial order
          if (a!=b & list(c(a,b)) %in% order) {
            order_new <- order[-which(order %in% list(c(a,b)))]
            if (is.partial.order(length(P),order_new)) {
              cons [[length(cons)+1]] <- c(a,b)
              if(!has.been.checked(n,P,order_new,K)) {
                # if (!(list(list(P,order_new)) %in% K)) {
                base_new <- compare_posets(n,d_connections,P,order_new,base)
                K [[length(K)+1]] <- list(P,order_new)
                if (length(base_new) > 1) {
                  base <- base_new
                  L <- list(list(P,order_new))
                  i <- 1
                  flag <- TRUE
                  break
                }
                else if (base_new==1) {
                  L[[length(L)+1]] <- list(P,order_new)
                }
              }
            }
          }
        #adding to partial order
        else if (a!=b & !(list(c(b,a)) %in% order)) {
          order_new <- append(order,list(c(a,b)))
          if (is.partial.order(length(P),order_new) & !has.been.checked(n,P,order_new,K)) {
            base_new <- compare_posets(n,d_connections,P,order_new,base)
            K [[length(K)+1]] <- list(P,order_new)
            if (length(base_new) > 1) {
              base <- base_new
              L <- list(list(P,order_new))
              i <- 1
              flag <- TRUE
              break
            }
            else if (base_new==1) {
              L[[length(L)+1]] <- list(P,order_new)
            }
          }
        }
      }
      if (flag) {break}
    }
    if (flag) {next}
    
    #########################################################
    #add one element of one section to a consecutive section
    for (pair in cons) {
      a <- pair[1]; b <- pair[2]
      #adding an element of P[[b]] to P[[a]]
      for (k in seq_int(1,length(P[[b]]))) {
        P_new <- P
        P_new[[b]] <- P[[b]][-k]
        P_new[[a]] <- sort(c(P[[a]],P[[b]][k]))
        order_new <- order
        order_new1 <- order_new
        if (length(P_new[[b]])==0) {
          for (pair in order_new1) {
            if (pair[1]==b | pair[2]==b) {
              order_new <- order_new[-which(order_new %in% list(pair))]
            }
          }
          for (j in seq_int(1,length(order_new))) {
            pair <- order_new[[j]]
            if (pair[1]>b & pair[2]>b) {
              order_new[[j]] <- c(pair[1]-1,pair[2]-1)
            }
            else if (pair[1]>b) {
              order_new[[j]] <- c(pair[1]-1,pair[2])
            }
            else if (pair[2]>b) {
              order_new[[j]] <- c(pair[1],pair[2]-1)
            }
          }
          P_new[[b]] <- NULL
        }
        if(!has.been.checked(n,P_new,order_new,K)) {
          # if (!(list(list(P_new,order_new)) %in% K)) {
          base_new <- compare_posets(n,d_connections,P_new,order_new,base)
          K [[length(K)+1]] <- list(P_new,order_new)
          if (length(base_new) > 1) {
            base <- base_new
            L <- list(list(P_new,order_new))
            i <- 1
            flag <- TRUE
            break
          }
          else if (base_new==1) {
            L[[length(L)+1]] <- list(P_new,order_new)
          }
        }
      }
      if (flag) {break}
      
      #adding an element of P[[a]] to P[[b]]
      for (k in seq_int(1,length(P[[a]]))) {
        P_new <- P
        P_new[[a]] <- P[[a]][-k]
        P_new[[b]] <- sort(c(P[[b]],P[[a]][k]))
        order_new <- order
        order_new1 <- order_new
        if (length(P_new[[a]])==0) {
          for (pair in order_new1) {
            if (pair[1]==a | pair[2]==a) {
              order_new <- order_new[-which(order_new %in% list(pair))]
            }
          }
          for (j in seq_int(1,length(order_new))) {
            pair <- order_new[[j]]
            if (pair[1]>a & pair[2]>a) {
              order_new[[j]] <- c(pair[1]-1,pair[2]-1)
            }
            else if (pair[1]>a) {
              order_new[[j]] <- c(pair[1]-1,pair[2])
            }
            else if (pair[2]>a) {
              order_new[[j]] <- c(pair[1],pair[2]-1)
            }
          }
          P_new[[a]] <- NULL
        }
        if(!has.been.checked(n,P_new,order_new,K)) {
          # if (!(list(list(P_new,order_new)) %in% K)) {
          base_new <- compare_posets(n,d_connections,P_new,order_new,base)
          K [[length(K)+1]] <- list(P_new,order_new)
          if (length(base_new) > 1) {
            base <- base_new
            L <- list(list(P_new,order_new))
            i <- 1
            flag <- TRUE
            break
          }
          else if (base_new==1) {
            L[[length(L)+1]] <- list(P_new,order_new)
          }
        }
      }
      if (flag) {break}
    }
    if (flag) {next}
    ################################################
    #add one element of one section to an empty set 
    for (a in seq_int(1,length(P))) {
      #add the element as a section right after the current section
      for (k in seq_int(1,length(P[[a]]))) {
        P_new <- c(P[seq_int(1,a-1)],list(P[[a]][-k]),P[seq_int(a+1,length(P))],list(c(P[[a]][k])))
        order_new <- append(order,list(c(as.double(length(P_new)),as.double(length(P_new)))))
        for (b in seq_int(1,length(P))) {
          if (list(c(b,a)) %in% order) {order_new <- append(order_new,list(c(as.double(b),as.double(length(P_new)))))}
          if (list(c(a,b)) %in% order & b!=a) {order_new <- append(order_new,list(c(as.double(length(P_new)),as.double(b))))}
        }
        order_new1 <- order_new
        if (length(P_new[[a]])==0) {
          if(i==92 & length(L)==123) {print("order_new1 is");print(order_new1)}
          for (pair in order_new1) {
            if(i==92 & length(L)==123) {print("hello");print(pair)}
            if (pair[1]==a | pair[2]==a) {
              order_new <- order_new[-which(order_new %in% list(pair))]
            }
          }
          for (j in seq_int(1,length(order_new))) {
            pair <- order_new[[j]]
            if (pair[1]>a & pair[2]>a) {
              order_new[[j]] <- c(pair[1]-1,pair[2]-1)
            }
            else if (pair[1]>a) {
              order_new[[j]] <- c(pair[1]-1,pair[2])
            }
            else if (pair[2]>a) {
              order_new[[j]] <- c(pair[1],pair[2]-1)
            }
            # if(i==92 & length(L)==123) {print("bye");print(order_new)}
          }
          P_new[[a]] <- NULL
        }
        if(!has.been.checked(n,P_new,order_new,K)) {
          # if (!(list(list(P_new,order_new)) %in% K)) {
          base_new <- compare_posets(n,d_connections,P_new,order_new,base)
          K [[length(K)+1]] <- list(P_new,order_new)
          if (length(base_new) > 1) {
            base <- base_new
            L <- list(list(P_new,order_new))
            i <- 1
            flag <- TRUE
            break
          }
          else if (base_new==1) {
            L[[length(L)+1]] <- list(P_new,order_new)
          }
        }
      }
      if (flag) {break}
      
      #add the element as a section right before the current section
      for (k in seq_int(1,length(P[[a]]))) {
        P_new <- c(P[seq_int(1,a-1)],list(P[[a]][-k]),P[seq_int(a+1,length(P))],list(c(P[[a]][k])))
        order_new <- append(order,list(c(as.double(length(P_new)),as.double(length(P_new)))))
        for (b in seq_int(1,length(P))) {
          if (list(c(b,a)) %in% order & b!=a) {order_new <- append(order_new,list(c(as.double(b),as.double(length(P_new)))))}
          if (list(c(a,b)) %in% order) {order_new <- append(order_new,list(c(as.double(length(P_new)),as.double(b))))}
        }
        order_new1 <- order_new
        if (length(P_new[[a]])==0) {
          for (pair in order_new1) {
            if (pair[1]==a | pair[2]==a) {
              order_new <- order_new[-which(order_new %in% list(pair))]
            }
          }
          for (j in seq_int(1,length(order_new))) {
            pair <- order_new[[j]]
            if (pair[1]>a & pair[2]>a) {
              order_new[[j]] <- c(pair[1]-1,pair[2]-1)
            }
            else if (pair[1]>a) {
              order_new[[j]] <- c(pair[1]-1,pair[2])
            }
            else if (pair[2]>a) {
              order_new[[j]] <- c(pair[1],pair[2]-1)
            }
          }
          P_new[[a]] <- NULL
        }
        if(!has.been.checked(n,P_new,order_new,K)) {
          # if (!(list(list(P_new,order_new)) %in% K)) {
          base_new <- compare_posets(n,d_connections,P_new,order_new,base)
          K [[length(K)+1]] <- list(P_new,order_new)
          if (length(base_new) > 1) {
            base <- base_new
            L <- list(list(P_new,order_new))
            i <- 1
            flag <- TRUE
            break
          }
          else if (base_new==1) {
            L[[length(L)+1]] <- list(P_new,order_new)
          }
        }
      }
      if (flag) {break}
    }
    if (flag) {next}
    i <- i+1 
  }
  return(L)
}









#Not updated with PATHS and A options
greedy_optimization_total_orders_depth_first <- function(n,d_connections,P0,order0) {
  # Track <- list()
  K <- list(list(P0,order0)) #list of checked partially ordered partitions
  base <- c(Inf,rep(-1,n-2))
  print("hello1")
  base <- compare_posets_total_orders(n,d_connections,P0,order0,base)
  L <- list() #list of strictly increasing neighbors that were checked
  # i <- 1
  flag <- TRUE
  P <- P0
  order <- order0
  print("hello")
  while (flag) {
    # if (i==1) {print(base)}
    print(length(K))
    print(base)
    flag <- FALSE #if a non-increasing neighbor is found
    # P <- L[[i]][[1]]
    # order <- L[[i]][[2]]
    cons <- list() #list of consecutive pairs of sections in order
    # Track [[length(Track)+1]] <- L[[i]]
    
    #finding neighbors by swapping two consecutive sections of the partition
    for (a in seq_int(1,length(P))) {
      for (b in seq_int(1,length(P))) {
        if (a!=b & list(c(a,b)) %in% order) {
          order_new <- order
          order_new[[which(order %in% list(c(a,b)))]] <- c(b,a)
          if (is.partial.order(length(P),order_new)) {
            cons [[length(cons)+1]] <- c(a,b)
            if(!has.been.checked(n,P,order_new,K)) {
              # if (!(list(list(P,order_new)) %in% K)) {
              base_new <- compare_posets_total_orders(n,d_connections,P,order_new,base)
              K [[length(K)+1]] <- list(P,order_new)
              if (length(base_new) > 1) {
                base <- base_new
                order <- order_new
                flag <- TRUE
                break
              }
              else if (base_new==1) {
                order <- order_new
                flag <- TRUE
                break
              }
              else {
                L [[length(L)+1]] <- list(P,order_new)
              }
            }
          }
        }
      }
      if (flag) {break}
    }
    if (flag) {next}
    
    #########################################################
    #add one element of one section to a consecutive section
    for (pair in cons) {
      a <- pair[1]; b <- pair[2]
      #adding an element of P[[b]] to P[[a]]
      for (k in seq_int(1,length(P[[b]]))) {
        P_new <- P
        P_new[[b]] <- P[[b]][-k]
        P_new[[a]] <- sort(c(P[[a]],P[[b]][k]))
        order_new <- order
        order_new1 <- order_new
        if (length(P_new[[b]])==0) {
          for (pair in order_new1) {
            if (pair[1]==b | pair[2]==b) {
              order_new <- order_new[-which(order_new %in% list(pair))]
            }
          }
          for (j in seq_int(1,length(order_new))) {
            pair <- order_new[[j]]
            if (pair[1]>b & pair[2]>b) {
              order_new[[j]] <- c(pair[1]-1,pair[2]-1)
            }
            else if (pair[1]>b) {
              order_new[[j]] <- c(pair[1]-1,pair[2])
            }
            else if (pair[2]>b) {
              order_new[[j]] <- c(pair[1],pair[2]-1)
            }
          }
          P_new[[b]] <- NULL
        }
        if(!has.been.checked(n,P_new,order_new,K)) {
          # if (!(list(list(P_new,order_new)) %in% K)) {
          base_new <- compare_posets_total_orders(n,d_connections,P_new,order_new,base)
          K [[length(K)+1]] <- list(P_new,order_new)
          if (length(base_new) > 1) {
            base <- base_new
            P <- P_new
            order <- order_new
            flag <- TRUE
            break
          }
          else if (base_new==1) {
            P <- P_new
            order <- order_new
            flag <- TRUE
            break
          }
          else {
            L [[length(L)+1]] <- list(P_new,order_new)
          }
        }
      }
      if (flag) {break}
      
      #adding an element of P[[a]] to P[[b]]
      for (k in seq_int(1,length(P[[a]]))) {
        P_new <- P
        P_new[[a]] <- P[[a]][-k]
        P_new[[b]] <- sort(c(P[[b]],P[[a]][k]))
        order_new <- order
        order_new1 <- order_new
        if (length(P_new[[a]])==0) {
          for (pair in order_new1) {
            if (pair[1]==a | pair[2]==a) {
              order_new <- order_new[-which(order_new %in% list(pair))]
            }
          }
          for (j in seq_int(1,length(order_new))) {
            pair <- order_new[[j]]
            if (pair[1]>a & pair[2]>a) {
              order_new[[j]] <- c(pair[1]-1,pair[2]-1)
            }
            else if (pair[1]>a) {
              order_new[[j]] <- c(pair[1]-1,pair[2])
            }
            else if (pair[2]>a) {
              order_new[[j]] <- c(pair[1],pair[2]-1)
            }
          }
          P_new[[a]] <- NULL
        }
        if(!has.been.checked(n,P_new,order_new,K)) {
          # if (!(list(list(P_new,order_new)) %in% K)) {
          base_new <- compare_posets_total_orders(n,d_connections,P_new,order_new,base)
          K [[length(K)+1]] <- list(P_new,order_new)
          if (length(base_new) > 1) {
            base <- base_new
            P <- P_new
            order <- order_new
            flag <- TRUE
            break
          }
          else if (base_new==1) {
            P <- P_new
            order <- order_new
            flag <- TRUE
            break
          }
          else {
            L [[length(L)+1]] <- list(P_new,order_new)
          }
        }
      }
      if (flag) {break}
    }
    if (flag) {next}
    ################################################
    #add one element of one section to an empty set 
    for (a in seq_int(1,length(P))) {
      #add the element as a section right after the current section
      for (k in seq_int(1,length(P[[a]]))) {
        P_new <- c(P[seq_int(1,a-1)],list(P[[a]][-k]),P[seq_int(a+1,length(P))],list(c(P[[a]][k])))
        order_new <- append(order,list(c(as.double(length(P_new)),as.double(length(P_new)))))
        for (b in seq_int(1,length(P))) {
          if (list(c(b,a)) %in% order) {order_new <- append(order_new,list(c(as.double(b),as.double(length(P_new)))))}
          if (list(c(a,b)) %in% order & b!=a) {order_new <- append(order_new,list(c(as.double(length(P_new)),as.double(b))))}
        }
        order_new1 <- order_new
        if (length(P_new[[a]])==0) {
          if(i==92 & length(L)==123) {print("order_new1 is");print(order_new1)}
          for (pair in order_new1) {
            if(i==92 & length(L)==123) {print("hello");print(pair)}
            if (pair[1]==a | pair[2]==a) {
              order_new <- order_new[-which(order_new %in% list(pair))]
            }
          }
          for (j in seq_int(1,length(order_new))) {
            pair <- order_new[[j]]
            if (pair[1]>a & pair[2]>a) {
              order_new[[j]] <- c(pair[1]-1,pair[2]-1)
            }
            else if (pair[1]>a) {
              order_new[[j]] <- c(pair[1]-1,pair[2])
            }
            else if (pair[2]>a) {
              order_new[[j]] <- c(pair[1],pair[2]-1)
            }
            if(i==92 & length(L)==123) {print("bye");print(order_new)}
          }
          P_new[[a]] <- NULL
        }
        if(!has.been.checked(n,P_new,order_new,K)) {
          # if (!(list(list(P_new,order_new)) %in% K)) {
          base_new <- compare_posets_total_orders(n,d_connections,P_new,order_new,base)
          K [[length(K)+1]] <- list(P_new,order_new)
          if (length(base_new) > 1) {
            base <- base_new
            P <- P_new
            order <- order_new
            flag <- TRUE
            break
          }
          else if (base_new==1) {
            P <- P_new
            order <- order_new
            flag <- TRUE
            break
          }
          else {
            L [[length(L)+1]] <- list(P_new,order_new)
          }
        }
      }
      if (flag) {break}
      
      #add the element as a section right before the current section
      for (k in seq_int(1,length(P[[a]]))) {
        P_new <- c(P[seq_int(1,a-1)],list(P[[a]][-k]),P[seq_int(a+1,length(P))],list(c(P[[a]][k])))
        order_new <- append(order,list(c(as.double(length(P_new)),as.double(length(P_new)))))
        for (b in seq_int(1,length(P))) {
          if (list(c(b,a)) %in% order & b!=a) {order_new <- append(order_new,list(c(as.double(b),as.double(length(P_new)))))}
          if (list(c(a,b)) %in% order) {order_new <- append(order_new,list(c(as.double(length(P_new)),as.double(b))))}
        }
        order_new1 <- order_new
        if (length(P_new[[a]])==0) {
          for (pair in order_new1) {
            if (pair[1]==a | pair[2]==a) {
              order_new <- order_new[-which(order_new %in% list(pair))]
            }
          }
          for (j in seq_int(1,length(order_new))) {
            pair <- order_new[[j]]
            if (pair[1]>a & pair[2]>a) {
              order_new[[j]] <- c(pair[1]-1,pair[2]-1)
            }
            else if (pair[1]>a) {
              order_new[[j]] <- c(pair[1]-1,pair[2])
            }
            else if (pair[2]>a) {
              order_new[[j]] <- c(pair[1],pair[2]-1)
            }
          }
          P_new[[a]] <- NULL
        }
        if(!has.been.checked(n,P_new,order_new,K)) {
          # if (!(list(list(P_new,order_new)) %in% K)) {
          base_new <- compare_posets_total_orders(n,d_connections,P_new,order_new,base)
          K [[length(K)+1]] <- list(P_new,order_new)
          if (length(base_new) > 1) {
            base <- base_new
            P <- P_new
            order <- order_new
            flag <- TRUE
            break
          }
          else if (base_new==1) {
            P <- P_new
            order <- order_new
            flag <- TRUE
            break
          }
          else {
            L [[length(L)+1]] <- list(P_new,order_new)
          }
        }
      }
      if (flag) {break}
    }
    if (flag) {next}
  }
  return(list(P,order,L))
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

#Not updated with PATHS and A
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