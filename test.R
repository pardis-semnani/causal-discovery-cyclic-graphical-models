source("check_algorithm.R")



find_partially_ordered_partition <- function(n,E) {
  Reachable <- list()
  for (s in seq_int(1,n)) {
    Reachable[[s]] <- c(s)
  }
  for (e in E) {
    for (s in seq_int(1,n)) {
      if (e[1] %in% Reachable[[s]]) {
        Reachable[[s]] <- unique(c(Reachable[[s]],Reachable[[e[2]]]))
      }
    }
  }
  
  P <- list()
  for (s in seq_int(1,n)) {P[[s]] <- c(s)}
  for (s in seq_int(1,n-1)) {
    for (t in seq_int(s+1,n)) {
      if (s %in% Reachable[[t]] & t %in% Reachable[[s]]) {
        i <- section(s,P)
        j <- section(t,P)
        if (i!=j) {P[[i]] <- unique(c(P[[i]],P[[j]]));P[[j]] <- NULL}
      }
    }
  }
  order <- list()
  #determining the order of the partition
  for (i in seq_int(1,length(P))) {
    for (j in seq_int(1,length(P))) {
      if (is.ancestor_cycle(n,E,Reachable,P[[i]],P[[j]])) {
        order <- unique(append(order,list(c(i,j))))
      }
    }
  }
  return(list(P,order))
}


structures <- function(n,E) {
  PPP <- find_partially_ordered_partition(n,E)
  P <- PPP[[1]]
  order <- PPP[2]
  
  #finding p-adjacencies
  P_ADJ <- list()
  for (a in seq_int(1,n-1)) {
    for (b in seq_int(a+1,n)) {
      if (list(c(a,b)) %in% E | list(c(b,a)) %in% E) {
        P_ADJ[[length(P_ADJ)+1]] <- c(a,b)
      }
      else {
        for (c in seq_int(1,n)) {
          if (c!=a & c!=b & list(c(a,c)) %in% E & list(c(b,c)) %in% E & (section(a,P)==section(c,P) | section(b,P)==section(c,P))) {
            P_ADJ[[length(P_ADJ)+1]] <- c(a,b)
            break
          }
        }
      }
    }
  }
  
  Non_conductors <- list()
  for (a in seq_int(1,n-1)) {
    for (c in seq_int(a+1,n)) {
      if (!(list(c(a,c)) %in% P_ADJ)) {
        for (b in seq_int(1,n)) {
          if (b!=a & b!=c & list(c(min(a,b),max(a,b))) %in% P_ADJ & list(c(min(c,b),max(c,b))) %in% P_ADJ & list(c(section(a,P),section(b,P))) %in% order & list(c(section(c,P),section(b,P))) %in% order) {
            Non_conductors[[length(Non_conductors)+1]] <- list(a,b,c)
          }
        }
      }
    }
  }
  
  PN <- list()
  for (triple in Non_conductors) {
    a <- triple[[1]]
    b <- triple[[2]]
    c <- triple[[3]]
    if (list(c(a,b)) %in% E & list(c(c,b)) %in% E) {
      PN[[length(PN)+1]] <- list(a,b,c)
      for (triple2 in Non_conductors) {
        if (triple2[[1]]==a & triple2[[3]]==c & triple2[[2]]!=b & list(c(section(b,P),section(triple2[[2]],P))) %in% order) {
          PN[[length(PN)+1]] <- triple2
        }
      }
    }
  }
  PN <- unique(PN)
  return(list(P_ADJ,Non_conductors,PN))
}


graph_discovery_with_original_POP <- function(n,E_original,N_SCCR=Inf,N_counter=Inf) { #N_SCCR: how many reruns should be done for each component of each partially ordered partition
  
  #graph discovery
  counter <- 0
  while (TRUE) {
    counter <- counter+1
    if (counter>N_counter) {
      return(list(FALSE,counter))
    }
    PPP <- find_partially_ordered_partition(n,E_original) 
    P_hat <- PPP[[1]]
    order_hat <- PPP[[2]]
    
    print(c("counter is",counter))
    
    markov_equivalent <- TRUE #TRUE if a Markov equivalent graph is constructed
    stru <- structures(n,E_original)
    P_ADJ <- stru[[1]]
    if (!good_with_p_adjacencies(n,list(),P_hat,order_hat,list(),P_ADJ)) {
      print("not good with p_adjacencies")
      markov_equivalent <- FALSE
      next
    }
    Non_conductors <- stru[[2]]
    PN <- stru[[3]]
    
    
    
    #constructing the graph component by component
    E <- list()
    for (i in seq_int(1,length(P_hat))) {
      #defining A,B
      A <- list()
      B <- list()
      for (e in P_ADJ) {
        if (section(e[1],P_hat)==i & section(e[2],P_hat)==i) {
          A[[length(A)+1]] <- e
        }
        else if (section(e[1],P_hat)==i & list(c(section(e[2],P_hat),i)) %in% order_hat) {
          B[[length(B)+1]] <- c(e[2],e[1])
        }
        else if (section(e[2],P_hat)==i & list(c(section(e[1],P_hat),i)) %in% order_hat) {
          B[[length(B)+1]] <- e
        }
      }
      #defining ComCh,NoComCh
      ComCh <- list()
      NoComCh <- list()
      for (triple in Non_conductors) {
        a <- triple[[1]]
        b <- triple[[2]]
        c <- triple[[3]]
        if (b %in% P_hat[[i]]) {
          if (!(list(triple) %in% PN)) {
            NoComCh[[length(NoComCh)+1]] <- c(min(a,c),max(a,c))
          }
          else {
            flag <- TRUE
            for (triple2 in PN) {
              if (triple2[[1]]==a & triple2[[3]]==c & section(triple[[2]],P_hat)!=i & list(c(section(triple[[2]],P_hat),i)) %in% order_hat) {
                flag <-FALSE
                break
              }
            }
            if (flag) {
              ComCh[[length(ComCh)+1]] <- c(min(a,c),max(a,c))
            }
          }
        }
      }
      ComCh <- unique(ComCh)
      NoComCh <- unique(NoComCh)
      
      #running SCCR algorithm
      C <- list()
      for (v in P_hat[[i]]) {C[[length(C)+1]] <- v}
      I <<- 0
      R <- recover(V=C,N=N_SCCR,E=A,H=B,ComCh=ComCh,NoComCh=NoComCh)
      if (is.character(R)) {
        print("failed in algorithm")
        markov_equivalent <- FALSE
        break
        
      }
      else { #check if C is a strongly connected component.
        Reachable <- list()
        for (s in C) {
          Reachable[[as.character(s)]] <- c(s) 
        }
        for (e in R) {
          for (s in C) {
            if (e[1] %in% Reachable[[as.character(s)]]) {
              Reachable[[as.character(s)]] <- sort(unique(c(Reachable[[as.character(s)]],Reachable[[as.character(e[2])]])))
            }
          }
        }
        for (s in C) {
          for (t in C) {
            if (!t %in% Reachable[[as.character(s)]]) {
              print("no cycle")
              markov_equivalent <- FALSE
              break
              
            }
          }
          if (!markov_equivalent) {break}
        }
        if (!markov_equivalent) {break}
      }#else
      E <- append(E,R)
    }
    if (markov_equivalent) {return(list(TRUE,P_hat,order_hat,E,counter))}
  }
}


test_graph_discovery_with_original_POP<- function(n,N_SCCR=Inf,N_counter=Inf,N_test,m_min,m_max) {
  number_of_edges <- c()
  time <- c()
  accuracy <- c()
  reason_of_failure <- c()
  how_many_tries <- c()
  
  #determining m
  m = runif(N_test,min=m_min,max=m_max)
  
  for (i in 1:N_test) {
    #defining the edge set E
    g <- erdos.renyi.game(n, m[i], type = "gnp",directed=TRUE)
    E_mat <- get.edgelist(g)
    E <- list()
    for (i in seq_int(1,dim(E_mat)[1])) {E[[i]] <- c(E_mat[i,])}
    number_of_edges = c(number_of_edges,length(E))
    print(E)
    
    # #finding the right partially ordered partition 
    # PPP <- find_partially_ordered_partition(n,E)
    # P <- PPP[[1]]
    # order <- PPP[[2]]
    # if (check_the_right_pop) {
    #   initial <- list(PPP)
    #   N_optimization <- 1
    #   N_counter <- 1
    # }
    
    # #finding all the paths in the graph
    # M = all_paths_separated(n,E)
    # 
    # #finding the optimal score
    # base <- c(Inf,-1,Inf,Inf,rep(-1,n-3),Inf)
    # base_right <- compare_posets(n,list(),P,order,base,M)
    
    #running the algorithm
    start_time <- Sys.time()
    E_recovered <- graph_discovery_with_original_POP(n,E,N_SCCR,N_counter)
    end_time <- Sys.time()
    time = c(time,difftime(end_time,start_time,units="secs"))
    
    #checking the accuracy
    if (!E_recovered[[1]]) {
      accuracy <- c(accuracy,0)
      reason_of_failure <- c(reason_of_failure,"graph failed")
    }
    # else if (compare_posets(n,list(),E_recovered[[2]],E_recovered[[3]],base_right,M)==0) {
    #   accuracy <- c(accuracy,0)
    #   reason_of_failure <- c(reason_of_failure,"class failed")
    # }
    else {
      accuracy <- c(accuracy,1)
      reason_of_failure <- c(reason_of_failure,"success")
    }
    how_many_tries <- c(how_many_tries,E_recovered[[length(E_recovered)]])
    print(data.frame(number_of_edges,time,accuracy,reason_of_failure,how_many_tries))
  }
  return(data.frame(number_of_edges,time,accuracy,reason_of_failure,how_many_tries))
}
