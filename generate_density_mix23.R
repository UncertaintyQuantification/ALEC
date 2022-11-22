library(Rcpp)
library(RcppEigen)
sourceCpp(file='src/picard.cpp') 
library(lhs)

#functions of getting external potential
#0-Inf V^ext (no other parameters to change)
# get_V_zero_inf = function(a,L,d){
#   t=a/d
#   k = (a+L)/d+1
#   V_zero_inf = rep(Inf,k)
#   V_zero_inf[(1+t):(k-t)] = 0
#   return(V_zero_inf)
# }

#equation (4.3) binary mixture (one more parameter: epsilon)
get_V_binary = function(a,L,d,epsilon){
  t=a/d
  k = (a+L)/d+1
  x =  a/2+((1:(k-2*t))-1)*d
  V_binary = rep(Inf,k)
  V_binary[(1+t):(k-t)] = -epsilon*((a/(x+a/2))^3+(a/(L+a/2-x))^3)
  return(V_binary)
}

#equation (25) gravitational potential (one more parameter: mg)
get_V_grav = function(a,L,d,mg){
  t=a/d
  k = (a+L)/d+1
  x =  a/2+((1:(k-2*t))-1)*d
  V_grav = rep(Inf,k)
  V_grav[(1+t):(k-t)] = mg*x
  return(V_grav)
}

# #equation (53) power-law potential (three more parameters: u0,x0,a0>0)
# get_V_power = function(a,L,d,u0,x0,a0){
#   t=a/d
#   k = (a+L)/d+1
#   x = -(L-a)/2+((1:(k-2*t))-1)*d
#   V_power = rep(Inf,k)
#   V_power[(1+t):(k-t)] = u0*(abs(x/x0))^a0
#   return(V_power)
# }

#global parameters
a = 1
L = 9
d = 0.01
t=a/d  
k = (a+L)/d+1

alpha = 0.02 #don't change
maxIter = 1000
threshold = 10^(-5) # !!threshold for picard iteration


n5 = 2000 # !!num of densities for mixture of group 3 and group 4


# generate external potential and densities for mixture of group 3 and group 4 
be_V_ext_mix_record = matrix(NA,ncol=k,nrow = (n5))
rho_mix_record = matrix(NA,ncol=k,nrow = (n5))
be_Omega_mix_record = rep(NA,times=(n5))

# initial density
rho0 = rep(0,k) 
rho0[(t+1):(k-t)]=0.2

# mixture of group 2 and group 3 
set.seed(1)
para5 = maximinLHS(n=n5, k=4) # !!use LHD to generate parameter(s) for mixture potential
be_mu5 = para5[,1]*3 #scale for beta*mu
epsilon5 = para5[,2]*2.2 # !!scale for epsilon
mg5 = para5[,3]*3 # !!scale for mg
wt5 = para5[,4] #weight of group 3
time_record_5=system.time(
  for(i in 1:n5){
    print(i)
    V_ext_mix = wt5[i]*get_V_binary(a,L,d,epsilon=epsilon5[i])+(1-wt5[i])*get_V_grav(a,L,d,mg=mg5[i])
    be_V_ext_mix_record[(i),] = V_ext_mix
    exp_n_be_V_ext_mix=exp(-V_ext_mix)
    
    rho_res=picard_rho_cpp(rho=rho0, beta_mu=be_mu5[i], a=a, L=L, exp_n_be_V_ext=exp_n_be_V_ext_mix, alpha=alpha, iteration=maxIter, threshold = threshold)
    rho = rho_res[[2]]
    rho_mix_record[(i),] = rho

    beta_Omega = beta_Omega_cpp(rho=rho,a=a,L=L)
    be_Omega_mix_record[i] = beta_Omega
    }
)
sum(is.na(rho_mix_record))

save.image(file = "generate_density_mix23.RData")
