source("d_connections.R")



mod <- function(a,b) {
  if (a %% b==0) {return(b)}
  else {return(a %% b)}
}



check_parents <- function(e,Parents,V,E,H,NoComCh) {
  flag <- TRUE
  for (v in Parents[[as.character(e[2])]]) {
    if(v!=e[1] & list(sort(c(v,e[1]))) %in% NoComCh) {return(FALSE)}
    else if (v!=e[1] & (v %in% V | e[1] %in% V) & !(list(sort(c(v,e[1]))) %in% E) & !(list(c(v,e[1])) %in% H) & !(list(c(e[1],v)) %in% H)) {return(FALSE)}
  }
  return(TRUE)
}

recover <- function(V,N,E,H,ComCh,NoComCh,D=list(),D_copy=list(),Cause=list(),avoid=0,a=0){ #In ComCh and NoComCh, the two elements in each pair are sorted.
  D_initial <- D; D_copy_initial <- D_copy; Cause_initial <- Cause; avoid_initial <- avoid; a_initial <- a
  print(c("Our rerun no. is",I))
  if (I==N) {return("falied")}
  I <<- I+1
  plot(graph(c(unlist(unique(D[!(D_copy %in% list(0))])))))
  first_edge <- TRUE
  #defining Checked and Parent
  #Checked contains sorted edges
  Checked <- list()
  Parents <- list()
  for (i in seq_int(1,length(D))) {
    e <- D[[i]]
    if (!(list(sort(e)) %in% Checked)) {
      Checked <- append(Checked,list(sort(e)))
    }
    if (!(D_copy[i] %in% list(0))) {
      Parents[[as.character(e[2])]] <- append(Parents[[as.character(e[2])]],list(e[1]))
    }
  }
  
  #choosing the common children
  if (length(which(D_copy %in% list(2))) < 2*length(ComCh)) {
    #L: list of pairs in ComCh which are already assigned a common child in D
    L <- list()
    for (i in which(D_copy %in% list(2))) {
      L[[length(L)+1]] <- sort(c(D[[i]][1],Cause[[i]]))
    }
    for (p in ComCh) {
      if (!(list(p) %in% L)) {
        flag <- FALSE
        for (v in V) {
          if (list(c(p[1],v)) %in% H & list(c(p[2],v)) %in% H & !(list(avoid) %in% list(c(p[1],v))) & !(list(avoid) %in% list(c(p[2],v))) & check_parents(c(p[1],v),Parents,V,E,H,NoComCh) & check_parents(c(p[2],v),Parents,V,E,H,NoComCh)) {
            D <- append(D,list(c(p[1],v),c(p[2],v)))
            D_copy <- append(D_copy,list(2,2))
            Cause <- append(Cause,list(p[2],p[1]))
            Parents[[as.character(v)]] <- append(Parents[[as.character(v)]],list(p[1],p[2]))
            Checked <- unique(append(Checked,list(sort(c(p[1],v)),sort(c(p[2],v)))))
            flag <- TRUE
            plot(graph(c(unlist(unique(D[!(D_copy %in% list(0))])))))
            break
          }
        }
        if (!flag) {return("failure")} #maybe improve this later
      }
    }
  }
  
  #defining a
  if (a==0) {
    for (i in 1:length(V)) {
      s <- 0
      for (e in E) {
        if (e[1]==V[[i]] | e[2]==V[[i]]) {s <- s+1}
      }
      for (e in H) {
        if (e[2]==V[[i]]) {s <- s+1}
      }
      if(i==1 || s<M[2]) {M <- c(V[[i]],s)}
    }
    a <- M[1]
  }
  Rep <- 0
  Rep_compare <- length(D)
  Record <- c(rep(list(0),length(V)))
  #adding the next edge
  while (length(Checked)!=length(E)+length(H) || (length(V)==2 && length(E)>0 && (!(list(c(V[[1]],V[[2]])) %in% D[!(D_copy %in% list(0))]) | !(list(c(V[[2]],V[[1]])) %in% D[!(D_copy %in% list(0))])))) {
    if (Rep > Rep_compare) {print("too many steps");return(recover(V,N,sample(E),sample(H),ComCh,NoComCh,D_initial,D_copy_initial,Cause_initial,avoid_initial,a_initial))}
    outsider <- FALSE
    flag <- FALSE
    reserved <- list()
    for (e in c(E,H)) {
      b <- 0
      if (e[2]==a & list(e) %in% E & (!first_edge | !(list(avoid) %in% list(c(e[2],e[1]))))) {b <- e[1];reserved <- append(reserved,e[1])}
      else if (e[2]==a & list(e) %in% H & (!first_edge | !(list(avoid) %in% list(e)))) {b <- e[1]}
      else if (e[1]==a & (!first_edge | !(list(avoid) %in% list(e)))) {b <- e[2];reserved <- append(reserved,e[2])}
      if (b!=0 & !(list(sort(c(a,b))) %in% Checked)) {
        flag <- TRUE
        break
      }
    }
    if (!flag) {
      if(length(reserved)==0) {return("failed")}#???
      Record[[which(V %in% list(a))]] <- Record[[which(V %in% list(a))]] + 1
      b <- reserved[[mod(Record[[which(V %in% list(a))]],length(reserved))]]
    }
    cause_num <- length(D)+1
    if (b %in% V) {
      D <- append(D,list(c(a,b)))
      D_copy <- append(D_copy,list(c(a,b)))
      Cause <- append(Cause,list(cause_num))
      Parents[[as.character(b)]] <- append(Parents[[as.character(b)]],list(a))
      if (list(sort(c(a,b))) %in% Checked) {Rep <- Rep+1}
      else {Rep <- 0;Rep_compare <- length(D)}
      Checked <- unique(append(Checked,list(sort(c(a,b)))))
      
      plot(graph(c(unlist(unique(D[!(D_copy %in% list(0))])))))
    }
    else {
    outsider <- TRUE
    D <- append(D,list(c(b,a)))
    D_copy <- append(D_copy,list(c(b,a)))
    Cause <- append(Cause,list(cause_num))
    Parents[[as.character(a)]] <- append(Parents[[as.character(a)]],list(b))
    if (list(sort(c(a,b))) %in% Checked) {Rep <- Rep+1}
    else {Rep <- 0;Rep_compare <- length(D)}
    Checked <- unique(append(Checked,list(sort(c(b,a)))))
    plot(graph(c(unlist(unique(D[!(D_copy %in% list(0))])))))
    }
    if (b %in% V & check_parents(c(a,b),Parents,V,E,H,NoComCh)) {a <- b}
    else if (!(b %in% V) & check_parents(c(b,a),Parents,V,E,H,NoComCh)) {a <- a}
    #the added edge is problematic
    else {
      i=length(D_copy)
      Potential_Problems <- list()
      flag <- FALSE
      vertex_repeat <- list(0,Inf)
      while (TRUE) {
        if (vertex_repeat[[1]]>vertex_repeat[[2]]) {
          print("OHH!
                ")
          return("failed")
        }
        e <- D[[i]]
        if  (e[1] %in% V) {
          #can we delete the problematic edge?
          for (v in V) {
            if (list(sort(c(v,e[1]))) %in% E & list(sort(c(v,e[2]))) %in% E & check_parents(c(e[1],v),Parents,V,E,H,NoComCh) & check_parents(c(e[2],v),Parents,V,E,H,NoComCh) & (!first_edge | !(list(avoid) %in% list(c(e[1],v),c(e[2],v))))) {#avoid??? 
              #in D, we put the reverse of a deleted edge
              D[[i]] <- c(e[2],e[1])
              D <- append(D,list(c(e[1],v),c(e[2],v)))
              D_copy[[i]] <- 0
              D_copy <- append(D_copy,list(1,1))
              Cause[[i]] <- cause_num #maybe refine this by determining which addition of (a,b) caused the problem.
              Cause[[length(Cause)+1]] <- cause_num
              Cause[[length(Cause)+1]] <- cause_num
              Parents[[as.character(v)]] <- append(Parents[[as.character(v)]],list(e[1],e[2]))
              Parents[[as.character(e[2])]][max(which(Parents[[as.character(e[2])]]==e[1]))] <- NULL
              if (list(sort(c(e[1],v))) %in% Checked & list(sort(c(e[2],v))) %in% Checked) {Rep <- Rep}
              else {Rep <- 0}
              Checked <- unique(append(Checked,list(sort(c(e[1],v)),sort(c(e[2],v)))))
              plot(graph(c(unlist(unique(D[!(D_copy %in% list(0))])))))
              a <- v
              flag <- TRUE
              break
            }
          }
          if (flag) {break}
          #couldn't delete the edge, so let's reverse it.
          D[[i]] <- c(e[2],e[1])
          D_copy[[i]] <- -1
          Cause[[i]] <- cause_num
          Parents[[as.character(e[1])]] <- append(Parents[[as.character(e[1])]],list(e[2]))
          Parents[[as.character(e[2])]][max(which(Parents[[as.character(e[2])]]==e[1]))] <- NULL
          plot(graph(c(unlist(unique(D[!(D_copy %in% list(0))])))))
          if (check_parents(c(e[2],e[1]),Parents,V,E,H,NoComCh)) {
            a <- e[1]
            flag <- TRUE
            break
          }
          else {
            #Reversing doesn't work. Let's see why.
            #L is the least of all problematic edges.
            L <- list()
            for (v in Parents[[as.character(e[1])]]) {
              if (v!=e[2] & !(list(sort(c(v,e[2]))) %in% E) & !(list(c(v,e[2])) %in% H)) {
                L[[length(L)+1]] <- c(v,e[1])
              }
            }
            #if there is a non-touchable problematic edge, then rerun the algorithm.
            if (TRUE %in% (D %in% L & D_copy %in% list(-1,0,1))) {
              j <- min(which(D %in% D[max(which(D %in% L & D_copy %in% list(-1,0,1)))])) #first appearance of the most recent untouchable problem
              # t <- which(D %in% L & D_copy %in% list(-1,0,1))
              # j <- t[ceiling(runif(1,min=0.0001,length(t)))]
              avoid <- D[[j]]
              if (D_copy[j] %in% list(2)) {
                if (j %% 2==1) {t <- c(seq_int(1,j-1),seq_int(j+2,2*length(ComCh)))}
                else {t <- c(seq_int(1,j-2),seq_int(j+1,2*length(ComCh)))}
                D_new <- D[t]
                D_copy_new <- D_copy[t]
                Cause_new <- Cause[t]
                a <- 0
                E <- sample(E)
                H <- sample(H)
                return(recover(V,N,E,H,ComCh,NoComCh,D_new,D_copy_new,Cause_new,avoid,a))
              }
              else {
                j <- Cause[[j]]
                D_new <- D[seq_int(1,j-1)]
                D_copy_new <- D_copy[seq_int(1,j-1)]
                Cause_new <- Cause[seq_int(1,j-1)]
                for (k in seq_int(1,j-1)) {
                  if (D_copy_new[k] %in% list(-1,0) & Cause[[k]]>=j) {
                    D_new[[k]] <- c(D_new[[k]][2],D_new[[k]][1])
                    D_copy_new[[k]] <- D_new[[k]]
                    Cause[[k]] <- k
                  }
                }
                if (D_copy[j] %in% list(-1,0) & list(sort(D[[j]])) %in% E) {a <- D[[j]][2]}
                else if(D_copy[j] %in% list(-1,0)) {a <- D[[j]][1]}
                else if(list(sort(D[[j]])) %in% E) {a <- D[[j]][1]}
                else {a <- D[[j]][2]}
                E <- sample(E)
                H <- sample(H)
                return(recover(V,N,E,H,ComCh,NoComCh,D_new,D_copy_new,Cause_new,avoid,a))
              }
            }
            #if all the problematic edges are touchable, then start with the most recent one. At the end, check if all the problems are solved.
            else {
              Potential_Problems <- append(Potential_Problems,list(c(e[2],e[1])))
              i <- max(which(D %in% L & !(D_copy %in% list(0))))
              vertex_repeat[[1]] <- 0
              next
            }
          }#if(reversing doesn't work)
        }#if(e[1] %in% V)
        else {
          #delete the edge that's coming from outside of the cycle.
          for (v in V) {
            
            if (list(c(e[1],v)) %in% H & list(sort(c(v,e[2]))) %in% E & check_parents(c(e[1],v),Parents,V,E,H,NoComCh) & check_parents(c(e[2],v),Parents,V,E,H,NoComCh) & (!first_edge | !(list(avoid) %in% list(c(e[1],v),c(e[2],v))))) { #avoid changed??
              if (!(D_copy[i] %in% list(2))) {
                D[[i]] <- c(e[2],e[1])
                D_copy[[i]] <- 0
                Cause[[i]] <- cause_num #maybe refine this by determining which addition of (a,b) caused the problem.
              }
              D <- append(D,list(c(e[1],v),c(e[2],v)))
              D_copy <- append(D_copy,list(1,1))
              Cause[[length(Cause)+1]] <- cause_num
              Cause[[length(Cause)+1]] <- cause_num
              Parents[[as.character(v)]] <- append(Parents[[as.character(v)]],list(e[1],e[2]))
              Parents[[as.character(e[2])]][max(which(Parents[[as.character(e[2])]]==e[1]))] <- NULL
              if (list(sort(c(e[1],v))) %in% Checked & list(sort(c(e[2],v))) %in% Checked) {Rep <- Rep}
              else {Rep <- 0;Rep_compare <- length(D)}
              Checked <- unique(append(Checked,list(sort(c(e[1],v)),sort(c(e[2],v)))))
              plot(graph(c(unlist(unique(D[!(D_copy %in% list(0))])))))
              aa <- v
              flag <- TRUE
              break
            }
          }
          #check if the deleted edge also holds a common child
          L <- which(D %in% list(e) & D_copy %in% list(2))
          if (flag & length(L)>0) {
            for (j in L) {
              flag <- FALSE
              for (v in V) {
                if (list(c(e[1],v)) %in% H & list(c(Cause[[j]],v)) %in% H & check_parents(c(e[1],v),Parents,V,E,H,NoComCh) & check_parents(c(Cause[[j]],v),Parents,V,E,H,NoComCh)) {
                 D[[j]] <- c(e[1],v)
                 D[[ifelse(c(j %% 2==1),c(j+1),c(j-1))]] <- c(Cause[[j]],v)
                 Parents[[as.character(v)]] <- append(Parents[[as.character(v)]],list(e[1],Cause[[j]]))
                 Parents[[as.character(e[2])]][min(which(Parents[[as.character(e[2])]]==e[1]))] <- NULL
                 Parents[[as.character(e[2])]][min(which(Parents[[as.character(e[2])]]==Cause[[j]]))] <- NULL
                 if (list(sort(c(e[1],v))) %in% Checked & list(sort(Cause[[j]],v)) %in% Checked) {Rep <- Rep}
                 else {Rep <- 0;Rep_compare <- length(D)}
                 Checked <- unique(append(Checked,list(sort(c(e[1],v)),sort(Cause[[j]],v))))
                 plot(graph(c(unlist(unique(D[!(D_copy %in% list(0))])))))
                 flag <- TRUE 
                 break
                }
              }
              if (!flag) {break}
            }
          }
          if (flag) {
            a <- aa
            break
          }
          #couldn't delete the edge. so, let's see if we can reverse something.
          #create a list of all problematic edges
          L <- list()
          for (v in Parents[[as.character(e[2])]]) {
            if ((v %in% V & !(list(c(e[1],v)) %in% H)) | (!(v %in% V) & v!=e[1] & list(sort(c(e[1],v))) %in% NoComCh)) {
              L[[length(L)+1]] <- c(v,e[2])
            }
          }
          #if there is a non-touchable problematic edge, then rerun the algorithm.
          if (TRUE %in% (D %in% L & D_copy %in% list(-1,0,1))) {
            #j <- min(which(D %in% D[max(which(D %in% L & D_copy %in% list(-1,0,1)))]))#first appearance of the most recent untouchable problem
            t <- which(D %in% L & D_copy %in% list(-1,0,1))
            j <- t[ceiling(runif(1,min=0.0001,length(t)))]
            avoid <- D[[j]]
            if (D_copy[j] %in% list(2)) {
              if (j %% 2==1) {t <- c(seq_int(1,j-1),seq_int(j+2,2*length(ComCh)))}
              else {t <- c(seq_int(1,j-2),seq_int(j+1,2*length(ComCh)))}
              D_new <- D[t]
              D_copy_new <- D_copy[t]
              Cause_new <- Cause[t]
              a <- 0
              E <- sample(E)
              H <- sample(H)
              return(recover(V,N,E,H,ComCh,NoComCh,D_new,D_copy_new,Cause_new,avoid,a))
            }
            else {
              j <- Cause[[j]]
              D_new <- D[seq_int(1,j-1)]
              D_copy_new <- D_copy[seq_int(1,j-1)]
              Cause_new <- Cause[seq_int(1,j-1)]
              for (k in seq_int(1,j-1)) {
                if (D_copy_new[k] %in% list(-1,0) & Cause[[k]]>=j) {
                  D_new[[k]] <- c(D_new[[k]][2],D_new[[k]][1])
                  D_copy_new[[k]] <- D_new[[k]]
                  Cause[[k]] <- k
                }
              }
              if (D_copy[j] %in% list(-1,0) & list(sort(D[[j]])) %in% E) {a <- D[[j]][2]}
              else if(D_copy[j] %in% list(-1,0)) {a <- D[[j]][1]}
              else if(list(sort(D[[j]])) %in% E) {a <- D[[j]][1]}
              else {a <- D[[j]][2]}
              E <- sample(E)
              H <- sample(H)
              return(recover(V,N,E,H,ComCh,NoComCh,D_new,D_copy_new,Cause_new,avoid,a))
            }
          }
          #if all the problematic edges are touchable, then start with the most recent one. At the end, check if all the problems are solved.
          else {
            Potential_Problems <- append(Potential_Problems,list(c(e[1],e[2])))
            i <- max(which(D %in% L & !(D_copy %in% list(0))))
            if (vertex_repeat[[1]]==0) {
              vertex_repeat[[2]] <- 0
              for (edge in D) {
                if(edge[2]==e[2]) {vertex_repeat[[2]] <- vertex_repeat[[2]]+1}
              }
            }
            vertex_repeat[[1]] <- vertex_repeat[[1]]+1
            next
          }
          # if (outsider) {A <- a;a <- b;b <- A}
          # problem <- list()
          # for (v in Parents[[as.character(b)]]) {
          #   if (v!=a & (((v %in% V | a %in% V) & !(list(sort(c(v,a))) %in% E) & !(list(c(v,a)) %in% H) & !(list(c(a,v)) %in% H)) | list(sort(c(a,v))) %in% NoComCh)) {
          #     problem[[length(problem)+1]] <- v
          #   }
          # }
          # print(c("length of problem is",length(problem),"outsider"))
          # v <- problem[[ceiling(runif(1,min=0.0001,max=length(problem)))]]
          # j <- min(which(D %in% list(c(v,b))))
          # avoid <- D[[j]]
          # if (D_copy[j] %in% list(2)) {
          #   if (j %% 2==1) {t <- c(seq_int(1,j-1),seq_int(j+2,2*length(ComCh)))}
          #   else {t <- c(seq_int(1,j-2),seq_int(j+1,2*length(ComCh)))}
          #   D_new <- D[t]
          #   D_copy_new <- D_copy[t]
          #   Cause_new <- Cause[t]
          #   a <- 0
          #   E <- sample(E)
          #   H <- sample(H)
          #   return(recover(V,N,E,H,ComCh,NoComCh,D_new,D_copy_new,Cause_new,avoid,a))
          # }
          # else {
          #   j <- Cause[[j]]
          #   D_new <- D[seq_int(1,j-1)]
          #   D_copy_new <- D_copy[seq_int(1,j-1)]
          #   Cause_new <- Cause[seq_int(1,j-1)]
          #   for (k in seq_int(1,j-1)) {
          #     if (D_copy_new[k] %in% list(-1,0) & Cause[[k]]>=j) {
          #       D_new[[k]] <- c(D_new[[k]][2],D_new[[k]][1])
          #       D_copy_new[[k]] <- D_new[[k]]
          #       Cause[[k]] <- k
          #     }
          #   }
          #   if (D_copy[j] %in% list(-1,0) & list(sort(D[[j]])) %in% E) {a <- D[[j]][2]}
          #   else if(D_copy[j] %in% list(-1,0)) {a <- D[[j]][1]}
          #   else if(list(sort(D[[j]])) %in% E) {a <- D[[j]][1]}
          #   else {a <- Cause[[j]][2]}
          #   E <- sample(E)
          #   H <- sample(H)
          #   return(recover(V,N,E,H,ComCh,NoComCh,D_new,D_copy_new,Cause_new,avoid,a))
          # }
        }#if (!(e[1] %in% V))
      }#while(TRUE)
      for (e in Potential_Problems) {
        if (!check_parents(e,Parents,V,E,H,NoComCh)) {
          #if (outsider) {A <- a;a <- b;b <- A}
          problem <- list()
          for (v in Parents[[as.character(e[2])]]) {
            if (v!=e[1] & (((v %in% V | e[1] %in% V) & !(list(sort(c(v,e[1]))) %in% E) & !(list(c(v,e[1])) %in% H) & !(list(c(e[1],v)) %in% H)) | list(sort(c(e[1],v))) %in% NoComCh)) {
              problem[[length(problem)+1]] <- v
            }
          }
          print(c("length of problem is",length(problem),"potential problems"))
          v <- problem[[ceiling(runif(1,min=0.0001,max=length(problem)))]]
          j <- min(which(D %in% list(c(v,e[2]))))
          avoid <- D[[j]]
          if (D_copy[j] %in% list(2)) {
            if (j %% 2==1) {t <- c(seq_int(1,j-1),seq_int(j+2,2*length(ComCh)))}
            else {t <- c(seq_int(1,j-2),seq_int(j+1,2*length(ComCh)))}
            D_new <- D[t]
            D_copy_new <- D_copy[t]
            Cause_new <- Cause[t]
            a <- 0
            E <- sample(E)
            H <- sample(H)
            return(recover(V,N,E,H,ComCh,NoComCh,D_new,D_copy_new,Cause_new,avoid,a))
          }
          else {
            j <- Cause[[j]]
            D_new <- D[seq_int(1,j-1)]
            D_copy_new <- D_copy[seq_int(1,j-1)]
            Cause_new <- Cause[seq_int(1,j-1)]
            for (k in seq_int(1,j-1)) {
              if (D_copy_new[k] %in% list(-1,0) & Cause[[k]]>=j) {
                D_new[[k]] <- c(D_new[[k]][2],D_new[[k]][1])
                D_copy_new[[k]] <- D_new[[k]]
                Cause[[k]] <- k
              }
            }
            if (D_copy[j] %in% list(-1,0) & list(sort(D[[j]])) %in% E) {a <- D[[j]][2]}
            else if(D_copy[j] %in% list(-1,0)) {a <- D[[j]][1]}
            else if(list(sort(D[[j]])) %in% E) {a <- D[[j]][1]}
            else {a <- D[[j]][2]}
            E <- sample(E)
            H <- sample(H)
            return(recover(V,N,E,H,ComCh,NoComCh,D_new,D_copy_new,Cause_new,avoid,a))
          }
        }
      }#for(e in Potential_Problems)
    }#if(the added edge is problematic)
    first_edge <- FALSE
  }#while(not all the edges have been checked)
  return(D[!(D_copy %in% list(0))])
}
