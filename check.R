library(cIRT)


# Generates n samples from the uniform distribution on (xMin,gapMin)U(gapMax,xMax)
runifGap <- function(n, xMin, xMax, gapMin, gapMax) {
  dGap <- gapMax- gapMin
  x <- runif(n, xMin, xMax - dGap)
  i <- x > gapMin
  x[i] <- x[i] + dGap
  x
}


compare_greedy_optimization <- function (n,p_min,p_max,m,N1,N2,linear_Gaussian=FALSE) { #n: number of vertices, p: prob. of an edge appearing, m: number of random graphs, N1: depth of greedy_optimization_partial_order_depth_first, N2: depth of ..._modified
  #If linear_Gaussian=TRUE, then d_connections are read off of the covariance matrix
  time1 = c()
  time2 = c()
  prob = c()
  number_of_edges = c()
  accuracy1 = c()
  accuracy2 = c()
  for (j in seq_int(1,m)) {
    print(j)
    #determining p
    p = runif(1,min=p_min,max=p_max)
    prob = c(prob,p)
    
    #defining the edge set E
    g <- erdos.renyi.game(n, p, type = "gnp",directed=TRUE)
    E_mat <- get.edgelist(g)
    E <- list()
    for (i in seq_int(1,dim(E_mat)[1])) {E[[i]] <- c(E_mat[i,])}
    number_of_edges = c(number_of_edges,length(E))
    
    
    #defining the coefficient matrix corresponding to the edge set E (for a linear SEM)
    if (linear_Gaussian) {
      coeff <- runifGap(length(E),-1,1,-0.25,0.25)
      A <- matrix(0,n,n)
      counter <- 0
      for (e in E) {
        counter <- counter + 1
        A[e[[2]],e[[1]]] <- coeff[counter]
      }
      A <- solve(diag(n)-A)
      A <- A %*% t(A)
    }
    else {A <- "a"}
    
    
    
    #finding the right partially ordered partition 
    PPP <- find_partially_ordered_partition(n,E)
    P <- PPP[[1]]
    order <- PPP[[2]]
    
    #finding all the paths in the graph
    if (linear_Gaussian) {M <- list()}
    else {M <- all_paths_separated(n,E)}
   
    
    
    #finding the optimal score
    base <- c(Inf,-1,Inf,Inf,rep(-1,n-3),Inf)
    # base_right <- compare_posets(n,list(),P,order,base,M,A)
    base_right <- true_score_graph(n,E,P,order)
    
    
    
    
    
    #defining the initial partially ordered partition
    P0 = list(seq_int(1,n))
    order0 = list(c(1,1))
    
    P1 = list()
    for (t in seq_int(1,n)) {P1[[length(P1)+1]] <- c(t)}
    order1 = list()
    for (t in seq_int(1,n)) {order1[[length(order1)+1]] <- c(t,t)}
    
    P2 = list(seq_int(1,floor(n/2)),seq_int(floor(n/2)+1,n))
    order2 = list(c(1,1),c(2,2),c(1,2))
    
    #executing greedy_optimization_partial_orders_depth_first
    start_time <- Sys.time()
    L0 =greedy_optimization_partial_orders_depth_first(n=n,P=P0,order=order0,N=N1,PATHS=M,A=A)
    L1 =greedy_optimization_partial_orders_depth_first(n=n,P=P1,order=order1,N=N1,PATHS=M,A=A)
    L2 =greedy_optimization_partial_orders_depth_first(n=n,P=P2,order=order2,N=N1,PATHS=M,A=A)
    end_time <- Sys.time()
    time1 = c(time1,difftime(end_time,start_time,units="secs"))
    K0 <- L0[[1]]
    II0 <- L0[[2]]
    P_result0 <- K0[[II0[length(II0)]]][[1]]
    order_result0 <- K0[[II0[length(II0)]]][[2]]
    
    K1 <- L1[[1]]
    II1 <- L1[[2]]
    P_result1 <- K1[[II1[length(II1)]]][[1]]
    order_result1 <- K1[[II1[length(II1)]]][[2]]
    
    K2 <- L2[[1]]
    II2 <- L2[[2]]
    P_result2 <- K2[[II2[length(II2)]]][[1]]
    order_result2 <- K2[[II2[length(II2)]]][[2]]
    if(length(compare_posets(n,list(),P_result0,order_result0,base_right,M,A))>1 || length(compare_posets(n,list(),P_result1,order_result1,base_right,M,A))>1 || length(compare_posets(n,list(),P_result2,order_result2,base_right,M,A))>1) {accuracy1 = c(accuracy1,2)} # Code 2 denotes that we have found a partially ordered partition with a score better that the score of the original partially ordered partition, which would not have been the case if our d_connection tests were accurate.}
    else if (compare_posets(n,list(),P_result0,order_result0,base_right,M,A)==1 | compare_posets(n,list(),P_result1,order_result1,base_right,M,A)==1 | compare_posets(n,list(),P_result2,order_result2,base_right,M,A)==1){
      accuracy1 = c(accuracy1,1)
    } 
    else {
      accuracy1 = c(accuracy1,0)
    }
    
    #executing greedy_optimization_partial_orders_depth_first_modified
    # start_time <- Sys.time()
    # L2 = greedy_optimization_partial_orders_depth_first_modified(n=n,P=P0,order=order0,N=N2,PATHS=M)
    # end_time <- Sys.time()
    # time2 = c(time2,difftime(end_time,start_time,units="secs"))
    # K <- L2[[1]]
    # II <- L2[[2]]
    # P_result <- K[[II[length(II)]]][[1]]
    # order_result <- K[[II[length(II)]]][[1]]
    # accuracy2 = c(accuracy2,compare_posets(n,list(),P_result,order_result,base_right,M))
    print(data.frame(time1,accuracy1,prob,number_of_edges))
  }
  # df <- data.frame(time1,accuracy1,time2,accuracy2)
  df <- data.frame(time1,accuracy1,prob,number_of_edges)
  return(df)
}



markov_equivalent_graph <- function(n,d_connections=list(),PATHS=list(),N_optimization=Inf,N_SCCR=Inf,N_counter=Inf,initial=list(),A="a") { #N_counter: how many partially ordered partitions should be checked before we stop. N_SCCR: how many reruns should be done for each component of each partially ordered partition
  #A = (I-C)^(-1)(I-C)^(-T), where C is the coefficient matrix of the linear Gaussian SEM 
  
  
  #defining the initial partially ordered partitions
  if (length(initial)==0) {
    initial <- list(list(list(seq_int(1,n)),list(c(1,1))))
    
    P = list()
    for (t in seq_int(1,n)) {P[[length(P)+1]] <- c(t)}
    order = list()
    for (t in seq_int(1,n)) {order[[length(order)+1]] <- c(t,t)}
    initial[[2]] <- list(P,order)
    
    initial[[3]] <- list(list(seq_int(1,floor(n/2)),seq_int(floor(n/2)+1,n)),list(c(1,1),c(2,2),c(1,2)))
  }
  
  #optimizing the graphical score over the set of partially ordered partitions
  base <- c(Inf,-1,Inf,Inf,rep(-1,n-3),Inf)
  for (i in seq_int(1,length(initial))) {
    L <- greedy_optimization_partial_orders_depth_first(n=n,d_connections=d_connections,P=initial[[i]][[1]],order=initial[[i]][[2]],N=N_optimization,PATHS=PATHS,A=A)
    base_new <- compare_posets(n,d_connections,L[[1]][[L[[2]][1]]][[1]],L[[1]][[L[[2]][1]]][[2]],base,PATHS,A)
    if (length(base_new) > 1) {
      base <- base_new
      K <- L[[1]]
      II <- L[[2]]
    }
    else if (base_new==1) {
      II <- c(II,L[[2]]+length(K))
      K <- c(K,L[[1]])
    }
  }
  
  #graph discovery
  counter <- 0
  while (TRUE) {
    counter <- counter+1
    if (counter > N_counter) {
      return(FALSE)
    }
    if (length(II)>0) {
      P_hat <- K[[II[1]]][[1]]
      order_hat <- K[[II[1]]][[2]]
    }
    else {
      L <- greedy_optimization_partial_orders_depth_first(n,d_connections,P_hat,order_hat,N=1,L=K,PATHS=PATHS,A=A)
      if (length(L[[2]])<2) {return(FALSE)} #maybe improve this
      P_hat <- L[[1]][[L[[2]][2]]][[1]]
      order_hat <- L[[1]][[L[[2]][2]]][[2]]
      K <- L[[1]]
    }
    
    print(c("counter",counter))
    print(c("the partition is ",P_hat))
    print(c("the order is", order_hat))
    
    markov_equivalent <- TRUE #TRUE if a Markov equivalent graph is constructed
    P_ADJ <- p_adjacencies(n,d_connections,P_hat,order_hat,PATHS,A)
    if (!good_with_p_adjacencies(n,d_connections,P_hat,order_hat,PATHS,P_ADJ,A)) {
      print("not good with p_adjacencies")
      II <- II[-1]
      markov_equivalent <- FALSE
      next
    }
    Non_conductors <- mutually_exclusive_efficient(n=n,d_connections=d_connections,P=P_hat,order=order_hat,k=1,PATHS=PATHS,p_adjacencies=P_ADJ,A=A)[[1]]
    PN <- perfect_nonconductors_efficient(n,d_connections,P_hat,order_hat,PATHS,Non_conductors,A)
    
    
    # #finding a linear extension of (P_hat,order_hat)
    # for (i in seq_int(1,length(P_hat)-1)) {
    #   for (j in seq_int(i+1,length(P_hat))) {
    #     if (list(c(j,i)) %in% order_hat) {
    #       a <- P_hat[[j]]
    #       P_hat[[j]] <- P[[i]]
    #       P_hat[[i]] <- a
    #       order_hat[[which(order_hat %in% list(c(j,i)))]] <- c(i,j)
    #     }
    #   }
    # }
    
    
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
        II <- II[-1]
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
              II <- II[-1]
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
    if (markov_equivalent) {return(list(TRUE,P_hat,order_hat,E))}
  }
}


test_markov_equivalent_graph <- function(n,N_optimization=Inf,N_SCCR=Inf,N_counter=Inf,N_test,initial=list(),m_min,m_max,check_the_right_pop=FALSE,linear_Gaussian=FALSE) {
  #If linear_Gaussian=TRUE, then d_connections are read off of the covariance matrix
  
  number_of_edges <- c()
  time <- c()
  accuracy <- c()
  reason_of_failure <- c()
  POP <- c()
  
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
    # E = list(c(7,1),c(8,1),c(1,2),c(2,8),c(6,2),c(4,3),c(6,3),c(8,3),c(3,6),c(7,6),c(5,7),c(7,8))
    
    
    
    #defining the coefficient matrix corresponding to the edge set E (for a linear SEM)
    if (linear_Gaussian) {
      coeff <- runifGap(length(E),-1,1,-0.25,0.25)
      A <- matrix(0,n,n)
      counter <- 0
      for (e in E) {
        counter <- counter + 1
        A[e[[2]],e[[1]]] <- coeff[counter]
      }
      A <- solve(diag(n)-A)
      A <- A %*% t(A)
    }
    else {A <- "a"}
    
    #finding the right partially ordered partition 
    PPP <- find_partially_ordered_partition(n,E)
    P <- PPP[[1]]
    order <- PPP[[2]]
    if (check_the_right_pop) {
      initial <- list(PPP)
      N_optimization <- 1
      N_counter <- 1
    }
    
    
    #finding all the paths in the graph
    if (linear_Gaussian) {M <- list()}
    else {M <- all_paths_separated(n,E)}
    
    #finding the optimal score
    base <- c(Inf,-1,Inf,Inf,rep(-1,n-3),Inf)
    base_right <- compare_posets(n,list(),P,order,base,M,A)
    
    #running the algorithm
    start_time <- Sys.time()
    E_recovered <- markov_equivalent_graph(n,list(),M,N_optimization,N_SCCR,N_counter,initial,A)
    end_time <- Sys.time()
    time = c(time,difftime(end_time,start_time,units="secs"))
    
    #does SCCR work on the right ordered partition
    E_POP <- graph_discovery_with_original_POP(n,E,N_SCCR,1)
    if (!E_POP[[1]]) {
      POP <- c(POP,0)
    }
    else {
      POP <- c(POP,1)
    }
    
    #checking the accuracy
    if (!E_recovered[[1]]) {
      accuracy <- c(accuracy,0)
      reason_of_failure <- c(reason_of_failure,"graph failed")
    }
    else if (compare_posets(n,list(),E_recovered[[2]],E_recovered[[3]],base_right,M,A)==0) {
      accuracy <- c(accuracy,0)
      reason_of_failure <- c(reason_of_failure,"class failed")
    }
    else {
      accuracy <- c(accuracy,1)
      reason_of_failure <- c(reason_of_failure,"success")
    }
    print(data.frame(number_of_edges,time,accuracy,reason_of_failure,POP))
  }
  return(data.frame(number_of_edges,time,accuracy,reason_of_failure,POP))
}
