ordered_partitions <- function(n) {
  if (n==1) {
    return(list(list(c(1))))
  }
  else {
    K <- list()
    for (P in ordered_partitions(n-1)) {
      for (i in 1:length(P)) {
        #Add n+1 to the ith section of the ordered partition of P
        P_copy = P
        P_copy[[i]] <- c(P_copy[[i]],n)
        K[[length(K)+1]] <- P_copy
        #Add (n+1) as a separate section between the (i-1)th and ith section of P
        P_copy = P
        for (j in seq(length(P)+1,i+1,by=-1)) {
          P_copy[[j]] <- P_copy[[j-1]]
        }
        P_copy[[i]] <- c(n)
        K[[length(K)+1]] <- P_copy
      }
      #Add (n+1) as the last section to the ordered partition P
      P[[length(P)+1]] <- c(n) 
      K[[length(K)+1]] <- P 
    }
    return(K)
  }
}