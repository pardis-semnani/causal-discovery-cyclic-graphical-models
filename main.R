test_markov_equivalent_graph(n=8,N_optimization=1,N_SCCR=100,N_counter=300,N_test=1,initial=list(list(list(c(1,2,3,6,8),c(4,5,7)),list(c(1,1),c(2,2),c(2,1)))),m_min=0.2,m_max=0.2,check_the_right_pop=FALSE)

result71 <- test_markov_equivalent_graph(n=7,N_optimization=30,N_SCCR=100,N_counter=300,N_test=30,m_min=0.2,m_max=0.2)
test =markov_equivalent_graph(7,d_connections=list(),PATHS=M,N_optimization=30,N_SCCR=100,N_counter=300,initial=list()) 
n=7
E=list(c(1,3),c(2,1),c(2,6),c(3,1),c(3,7),c(3,4),c(3,5),c(3,6),c(4,1),c(4,2),c(4,7),c(4,6),c(6,3),c(7,6))


  # P1 = list(c(10),c(1,2,3,4,5,7,8,9),c(6))
# order1 = list(c(1,1),c(2,2),c(3,3),c(2,1))
# compare_posets(n,list(),P1,order1,base,M)
source("d_connections.R")
source("five_properties.R")
source("optimization.R")
source("greedy_optimization.R")
source("recovering_cycles.R")
source("check_algorithm.R")
source("greedy_optimization_partial_orders.R")
source("test.R")
source("six_properties.R")
source("check.R")
source("six_properties_from_original_graph.R")
# # d_connections(5,list(c(1,2),c(2,3),c(3,1),c(2,4),c(3,5)))
# E = list(c(1,2),c(2,3),c(3,1),c(4,1),c(5,2),c(4,6),c(6,5),c(5,4),c(7,1),c(8,1))
# E = list(c(1,2),c(2,3),c(3,1),c(4,1),c(5,2),c(4,6),c(6,5),c(5,4))
# E = list(c(1,2),c(1,3))
# E = list(c(4,1),c(5,2),c(3,1),c(1,2),c(2,3))
# E = list(c(1,3),c(2,3),c(3,4),c(6,3),c(4,5),c(5,6),c(5,7),c(1,7))
test_markov_equivalent_graph(n=8,N_optimization=1,N_SCCR=100,N_counter=300,N_test=1,initial=list(list(list(c(1,2,3,6,8),c(4,5,7)),list(c(1,1),c(2,2),c(2,1)))),m_min=0.2,m_max=0.2,check_the_right_pop=FALSE)

# E = list(c(1,2),c(2,3),c(3,4),c(5,4),c(2,5),c(4,2),c(6,7),c(7,6),c(6,8),c(3,8))
# E6(n,d_connections1,list(c(7),c(8),c(1,2,3),c(4,5,6)),list(c(1,1),c(2,2),c(3,3),c(4,4),c(1,3),c(2,3),c(4,3)),0)
n = 10
p = 0.2
# random graph
g <- erdos.renyi.game(n, p, type = "gnp",directed=TRUE)
#defining the edge set E
E_mat <- get.edgelist(g)
E <- list()
for (i in 1:dim(E_mat)[1]) {E[[i]] <- c(E_mat[i,])}
print(E)
####
PPP <- find_partially_ordered_partition(n,E)
P <- PPP[[1]]
order <- PPP[[2]]
print(PPP)
# d_connections = d_connections(n,E)
M = all_paths_separated(n,E)
# check_d_connections(s,t,c(),PATHS=M)
# for (s in seq_int(1,n-1)) {
#   for (t in seq_int(s+1,n)) {
#     my_vec <- c(seq_int(1,s-1),seq_int(s+1,t-1),seq_int(t+1,n))
#     subsets <- if(length(my_vec)<=1) {list(c(),my_vec)} else {c(list(c()),unlist(lapply(1:length(my_vec),combinat::combn,x = my_vec,simplify = FALSE),recursive = FALSE))}
#     for (S in subsets) { 
#       A = list(list(s,t,as.double(S))) %in% d_connections
#       B = check_d_connections(s,t,S,PATHS=M)
#       if (A!=B) {
#         print(c(s,t,length(S)))
#         print(c(A,B))
#       }
#     }
#   }
# }
# d_connections3 = d_connections(n,E)
# d_connections1 = d_connections(n,E)
print("hello")
# P = list(c(1,2,3),c(4,5,6),c(7),c(8))
# order = list(c(3,3),c(4,4),c(1,1),c(2,2),c(2,1),c(4,1),c(3,1))
# p_adjacencies(n,d_connections1,P,order)
# p_adjacent(n,d_connections1,P,order,1,6) 
# itineraries(n,d_connections1,P,order,P[1],2)
# mutually_exclusive(n,d_connections1,P,order,2)
# perfect_nonconductors(n,d_connections1,P,order)
# E4(n,d_connections1,P,order)
# E6 (n,d_connections1,P,order)
# P0 = list(c(1,2,3,4,5,6))
# # Track = greedy_optimization(n,d_connections,P0)
# # length(p_adjacencies(n,d_connections,P)) == length(p_adjacencies(n,d_connections,L[[1]]))
# # length(mutually_exclusive(n,d_connections,P,1)) == length(mutually_exclusive(n,d_connections,L[[1]],1))
# # length(mutually_exclusive(n,d_connections,P,2)) == length(mutually_exclusive(n,d_connections,L[[1]],2))
# # length(mutually_exclusive(n,d_connections,P,3)) == length(mutually_exclusive(n,d_connections,L[[1]],3))
# # length(perfect_nonconductors(n,d_connections,P)) == length(perfect_nonconductors(n,d_connections,L[[1]]))
# # E4(n,d_connections,P) == E4(n,d_connections,L[[1]])
# # E6(n,d_connections,P) == E6(n,d_connections,L[[1]])
# # K = optimize(n,d_connections)
# # g(n,d_connections,P)
# # Track
# # two_edge_connected(list(1,2,3),list(c(1,2),c(2,5),c(5,4),c(2,3),c(1,3)))
# P = check_connectivity(n,d_connections,P0)
# I = 0
# V = list(1,2,3,4,5,6)
# N=100
# E=list(c(1,3),c(3,5),c(1,5),c(2,4),c(2,6),c(4,6),c(1,4),c(2,5))
# H=list()
# # ComCh=list()
# # recover(V=V,N=N,E=E,H=H,ComCh=list(),NoComCh=list())
# H1=list(c(7,1),c(8,1),c(9,6),c(7,4),c(8,4),c(9,2))
# ComCh=list(c(7,8))
# NoComCh=list(c(7,9))
# I = 0
# recover(V=V,N=N,E=E,H=H1,ComCh=ComCh,NoComCh=NoComCh)

# E=list(c(1,2),c(1,4),c(3,4),c(3,5),c(2,5),c(5,6),c(1,6),c(3,6),c(4,6),c(2,3),c(2,6),c(1,3),c(4,5),c(1,5))
# recover(V=V,N=N,E=E,H=H,ComCh=list(),NoComCh=list())
# g <- erdos.renyi.game(6, 0.3, type = "gnp",directed=TRUE)
# plot(g)
# I=0

# E1 = list(c(2,7),c(3,6),c(3,9),c(3,5),c(5,10),c(5,14),c(2,6),c(6,8),c(6,17),c(3,7),c(5,7),c(6,7),c(3,8),c(7,9),c(10,15),c(12,13),c(2,13),c(7,13),c(9,14),c(6,15),c(9,15),c(12,15),c(2,17),c(6,9),c(3,10),c(3,14),c(10,14),c(2,8),c(8,17),c(5,6),c(2,12),c(7,12),c(6,12),c(9,12))
# H1=list(c(4,2),c(16,5),c(18,15),c(4,17),c(4,7),c(16,3),c(16,10),c(16,14),c(18,6),c(18,9),c(18,12))
# NoComCh=list(c(4,11),c(4,16),c(4,18),c(4,20),c(11,16),c(11,18),c(11,20),c(16,18),c(16,20),c(18,20))
# ComCh=list()
# V1=list(2,3,5,6,7,8,9,10,12,13,14,15,17)
# N=100
# recover(V1,N,E1,H1,ComCh,NoComCh,D=list(),D_copy=list(),Cause=list(),avoid=0,a=0)

# check_algorithm(30,0.3,1000)
# P=list(c(1),c(2),c(3))
# order=list(c(1,2),c(1,3))
# base=c(3,0,0,0,0)
# compare_posets(n,d_connections1,P,order,base)
base <- c(Inf,-1,Inf,Inf,rep(-1,n-3),Inf)
compare_posets(n,list(),P,order,base,M)
#####
P0 = list(seq_int(1,n))
order0 = list(c(1,1))
# P0 = list(c(1,2),c(3,4,5,6,7))
# order0=list(c(1,1),c(2,2),c(1,2))
# P0 = list(c(1),c(3),c(2))
# order0 = list(c(1,1),c(2,2),c(3,3),c(1,3),c(2,3),c(1,2))
# P0 = list(c(1),c(2),c(3),c(4),c(5))
# order0=list(c(1,1),c(2,2),c(3,3),c(4,4),c(5,5))
# L=greedy_optimization_partial_orders(n,d_connections1,P0,order0)
# L3=greedy_optimization_total_orders_depth_first(n,d_connections1,P0,order0)
# L1 = greedy_optimization_partial_orders_depth_first(n=n,P=P0,order=order0,N=1000,PATHS=M)
# L2 = greedy_optimization_partial_orders_depth_first_modified(n=n,P=P0,order=order0,N=6,PATHS=M)
# list(list(P,order)) %in% L


# E6(n,d_connections1,list(c(1,2,3),c(4),c(5)),list(c(1,1),c(2,2),c(3,3),c(2,1),c(3,1)))
# PP = list(c(2,3),c(5),c(1),c(4))
# oo = list(c(1,1),c(2,2),c(2,1),c(3,3),c(3,1),c(2,3),c(4,4))
# p_adjacencies(n,d_connections1,PP,oo)
# df1 <- compare_greedy_optimization(10,0.15,100,100,3)
# df2 <- compare_greedy_optimization(10,0.2,100,100,3)
# df3 <- compare_greedy_optimization(10,0.25,100,100,3)
# df4 <- compare_greedy_optimization(10,0.15,100,50,3)
# 
# df5 <- compare_greedy_optimization(10,0,0.25,100,50,3)
E_check <- markov_equivalent_graph(n,list(),M,N_optimization=1,N_SCCR=Inf,initial=list(PPP))
print(E_check)
E_recovered <- markov_equivalent_graph(n,list(),M,N_optimization=30,N_SCCR=100)



test_markov_equivalent_graph(n=10,N_optimization=50,N_SCCR=100,N_counter=100,N_test=20,m_min=1,m_max=40,check_the_right_pop =TRUE)
result7 <- test_markov_equivalent_graph(n=7,N_optimization=30,N_SCCR=100,N_counter=300,N_test=20,m_min=1,m_max=42)
result8 <-test_markov_equivalent_graph(n=8,N_optimization=30,N_SCCR=100,N_counter=500,N_test=30,m_min=1,m_max=42)
result9 <-test_markov_equivalent_graph(n=9,N_optimization=40,N_SCCR=100,N_counter=100,N_test=30,m_min=10,m_max=30)
result10 <-test_markov_equivalent_graph(n=10,N_optimization=50,N_SCCR=100,N_counter=100,N_test=10,m_min=10,m_max=30)

SCCR10 <- test_graph_discovery_with_original_POP(n=10,N_SCCR=100,N_counter=100,N_test=30,m_min=1,m_max=90)
SCCR11 <- test_graph_discovery_with_original_POP(n=11,N_SCCR=100,N_counter=100,N_test=30,m_min=1,m_max=110)
SCCR12 <- test_graph_discovery_with_original_POP(n=12,N_SCCR=100,N_counter=100,N_test=30,m_min=1,m_max=132)
SCCR9 <- test_graph_discovery_with_original_POP(n=9,N_SCCR=100,N_counter=700,N_test=30,m_min=1,m_max=72)
SCCR8 <- test_graph_discovery_with_original_POP(n=8,N_SCCR=100,N_counter=100,N_test=30,m_min=1,m_max=56)
SCCR7 <- test_graph_discovery_with_original_POP(n=7,N_SCCR=100,N_counter=100,N_test=30,m_min=1,m_max=42)
SCCR100 <- test_graph_discovery_with_original_POP(n=100,N_SCCR=100,N_counter=100,N_test=5,m_min=1,m_max=9900)

result71 <- test_markov_equivalent_graph(n=7,N_optimization=30,N_SCCR=100,N_counter=300,N_test=16,m_min=0.2,m_max=0.2)
result72 <- test_markov_equivalent_graph(n=7,N_optimization=30,N_SCCR=100,N_counter=300,N_test=30,m_min=0.4,m_max=0.4)
result73 <- test_markov_equivalent_graph(n=7,N_optimization=30,N_SCCR=100,N_counter=300,N_test=30,m_min=0.6,m_max=0.6)
result74 <- test_markov_equivalent_graph(n=7,N_optimization=30,N_SCCR=100,N_counter=300,N_test=30,m_min=0.8,m_max=0.8)
result75 <- test_markov_equivalent_graph(n=7,N_optimization=30,N_SCCR=100,N_counter=300,N_test=30,m_min=0.81,m_max=1)

result82 <- test_markov_equivalent_graph(n=8,N_optimization=30,N_SCCR=100,N_counter=300,N_test=30,m_min=0.2,m_max=0.2) #min and max?
result92 <- test_markov_equivalent_graph(n=9,N_optimization=40,N_SCCR=100,N_counter=300,N_test=30,m_min=0.2,m_max=0.2)
result102 <- test_markov_equivalent_graph(n=10,N_optimization=50,N_SCCR=100,N_counter=300,N_test=30,m_min=0.2,m_max=0.2)
result103 <- test_markov_equivalent_graph(n=10,N_optimization=50,N_SCCR=100,N_counter=300,N_test=30,m_min=0.3,m_max=0.3)
result733 <- test_markov_equivalent_graph(n=7,N_optimization=30,N_SCCR=100,N_counter=300,N_test=30,m_min=0.3,m_max=0.3)
resul833 <- test_markov_equivalent_graph(n=8,N_optimization=30,N_SCCR=100,N_counter=300,N_test=30,m_min=0.3,m_max=0.3)
result933 <- test_markov_equivalent_graph(n=9,N_optimization=40,N_SCCR=100,N_counter=300,N_test=30,m_min=0.3,m_max=0.3)
result1033 <- test_markov_equivalent_graph(n=10,N_optimization=50,N_SCCR=100,N_counter=300,N_test=13,m_min=0.3,m_max=0.3)
result92 <- test_markov_equivalent_graph(n=9,N_optimization=40,N_SCCR=100,N_counter=300,N_test=30,m_min=0.2,m_max=0.2)


opt71 <- compare_greedy_optimization(n=7,p_min=0.2,p_max=0.2,m=30,N1=30,N2=1)
opt72 <- compare_greedy_optimization(n=7,p_min=0.4,p_max=0.4,m=30,N1=30,N2=1)
opt73 <- compare_greedy_optimization(n=7,p_min=0.6,p_max=0.6,m=30,N1=30,N2=1)
opt74 <- compare_greedy_optimization(n=7,p_min=0.8,p_max=0.8,m=30,N1=30,N2=1)

opt82 <- compare_greedy_optimization(n=8,p_min=0.2,p_max=0.2,m=30,N1=30,N2=1) #min and max?
opt92 <- compare_greedy_optimization(n=9,p_min=0.2,p_max=0.2,m=30,N1=40,N2=1)
opt102 <- compare_greedy_optimization(n=10,p_min=0.2,p_max=0.2,m=30,N1=50,N2=1)

opt733 <- compare_greedy_optimization(n=7,p_min=0.3,p_max=0.3,m=30,N1=30,N2=1)
opt83 <- compare_greedy_optimization(n=8,p_min=0.3,p_max=0.3,m=30,N1=30,N2=1) #min and max?
opt93 <- compare_greedy_optimization(n=9,p_min=0.3,p_max=0.3,m=30,N1=40,N2=1)
opt103 <- compare_greedy_optimization(n=10,p_min=0.3,p_max=0.3,m=6,N1=50,N2=1)


##NEW
opt151 <- compare_greedy_optimization(n=15,p_min=0.2,p_max=0.2,m=30,N1=30,N2=1,TRUE) #alp = 0.001
opt152 <- compare_greedy_optimization(n=15,p_min=0.4,p_max=0.4,m=30,N1=30,N2=1,TRUE) #alp = 0.001
opt153 <- compare_greedy_optimization(n=15,p_min=0.6,p_max=0.6,m=30,N1=30,N2=1,TRUE) #alp = 0

opt1512 <- compare_greedy_optimization(n=15,p_min=0.2,p_max=0.2,m=30,N1=30,N2=1,TRUE) #alp = 0
opt1522 <- compare_greedy_optimization(n=15,p_min=0.4,p_max=0.4,m=30,N1=30,N2=1,TRUE) #alp = 0

opt1523 <- compare_greedy_optimization(n=15,p_min=0.4,p_max=0.4,m=30,N1=50,N2=1,TRUE) #alp = 0

opt154 <- compare_greedy_optimization(n=15,p_min=0.8,p_max=0.8,m=30,N1=30,N2=1,TRUE) #alp = 0

opt153 <- compare_greedy_optimization(n=20,p_min=0.4,p_max=0.4,m=30,N1=30,N2=1,TRUE) #alp = 0
opt
###


SCCR71 <- test_graph_discovery_with_original_POP(n=7,N_SCCR=100,N_counter=20,N_test=100,m_min=0.2,m_max=0.2)
SCCR72 <- test_graph_discovery_with_original_POP(n=7,N_SCCR=100,N_counter=20,N_test=100,m_min=0.4,m_max=0.4)
SCCR73 <- test_graph_discovery_with_original_POP(n=7,N_SCCR=100,N_counter=20,N_test=100,m_min=0.6,m_max=0.6)
SCCR74 <- test_graph_discovery_with_original_POP(n=7,N_SCCR=100,N_counter=20,N_test=100,m_min=0.8,m_max=0.8)
SCCR81 <- test_graph_discovery_with_original_POP(n=8,N_SCCR=100,N_counter=20,N_test=100,m_min=0.2,m_max=0.2)
SCCR91 <- test_graph_discovery_with_original_POP(n=9,N_SCCR=100,N_counter=20,N_test=100,m_min=0.2,m_max=0.2)
SCCR101 <- test_graph_discovery_with_original_POP(n=10,N_SCCR=100,N_counter=20,N_test=100,m_min=0.2,m_max=0.2)
SCCR722 <- test_graph_discovery_with_original_POP(n=7,N_SCCR=100,N_counter=20,N_test=100,m_min=0.3,m_max=0.3)
SCCR82 <- test_graph_discovery_with_original_POP(n=8,N_SCCR=100,N_counter=20,N_test=100,m_min=0.3,m_max=0.3)
SCCR92 <- test_graph_discovery_with_original_POP(n=9,N_SCCR=100,N_counter=20,N_test=100,m_min=0.3,m_max=0.3)
SCCR102 <- test_graph_discovery_with_original_POP(n=10,N_SCCR=100,N_counter=20,N_test=100,m_min=0.3,m_max=0.3)

opt102 <- compare_greedy_optimization(n=10,p_min=0.2,p_max=0.2,m=18,N1=60,N2=1)
result82 <- test_markov_equivalent_graph(n=8,N_optimization=30,N_SCCR=100,N_counter=300,N_test=30,m_min=0.2,m_max=0.2) #min and max?



#test greedy optimization on dense graphs 
#test greedy optimization on sparse graphs
#test SCCR on dense graphs /DONE
#test SCCR on sparse graphs /DONE
#test the whole thing on dense graphs /DONE
#test the whole thing on sparse graphs
E = list(c(7,1),c(8,1),c(1,2),c(2,8),c(6,2),c(4,3),c(6,3),c(8,3),c(3,6),c(7,6),c(5,7),c(7,8))

V = list(1,2,3,4,5,6)
N=100
H1=list(c(7,1),c(8,1),c(9,6),c(7,4),c(8,4),c(9,2))
E=list(c(1,3),c(3,5),c(1,5),c(2,4),c(2,6),c(4,6),c(1,4),c(2,5))
ComCh=list(c(7,8))
NoComCh=list(c(7,9))
I = 0
recover(V=V,N=N,E=E,H=H1,ComCh=ComCh,NoComCh=NoComCh)

SCCR20 <- test_graph_discovery_with_original_POP(n=20,N_SCCR=100,N_counter=20,N_test=100,m_min=0,m_max=1)

E1=list(c(1,2),c(3,2),c(2,3))
E2=list(c(1,3),c(1,2),c(3,2))
L1=d_connections(3,E1)
L2=d_connections(3,E2)

n=4
E=list(c(1,2),c(2,3),c(3,2),c(4,3))
d_connections(n,E)

n=5
E=list(c(1,2),c(2,3),c(3,4),c(4,5),c(5,1))
d_connections(n,E)
