source("five_properties.R")
source("ordered_partitions.R")

optimize <- function(n,d_connections) {
  P <- ordered_partitions(n)
  # optimizing the p-adjacencies
  A <- lapply(lapply(P,p_adjacencies,n=n,d_connections=d_connections),length)
  P <- P[A == min(unlist(A))]
  # optimizing perfect unshielded nonconductors
  A <- lapply(lapply(P,perfect_nonconductors,n=n,d_connections=d_connections),length)
  P <- P[A == min(unlist(A))]
  # optimizing mutually exclusive conductors
  for (k in 1:n) {
    A <- lapply(lapply(P,mutually_exclusive,n=n,d_connections=d_connections,k=k),length)
    P <- P[A == max(unlist(A))]
  }
  # optimizing condition 4 of Richardson's Theorem
  A <- lapply(P,E4,n=n,d_connections=d_connections)
  P <- P[A == min(unlist(A))]
  # optimizing condition 5 of Richardson's Theorem
  A <- lapply(P,E6,n=n,d_connections=d_connections)
  P <- P[A == min(unlist(A))]
  return(P)
}
