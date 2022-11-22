library(Rcpp)
library(RcppEigen)
sourceCpp(file='src/picard.cpp') 
library(lhs)

#functions of getting external potential
#0-Inf V^ext (no other parameters to change)
get_V_zero_inf = function(a,L,d){
  t=a/d
  k = (a+L)/d+1
  V_zero_inf = rep(Inf,k)
  V_zero_inf[(1+t):(k-t)] = 0
  return(V_zero_inf)
}

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

#equation (53) power-law potential (three more parameters: u0,x0,a0>0)
get_V_power = function(a,L,d,u0,x0,a0){
  t=a/d
  k = (a+L)/d+1
  x = -(L-a)/2+((1:(k-2*t))-1)*d
  V_power = rep(Inf,k)
  V_power[(1+t):(k-t)] = u0*(abs(x/x0))^a0
  return(V_power)
}

#global parameters
a = 1
L = 9
d = 0.01
t=a/d
k = (a+L)/d+1

alpha = 0.02 #don't change
maxIter = 1000
threshold = 10^(-5) # !!threshold for picard iteration


#generate training set
n1 = 2000 # !!num of densities for zero-infinity potential (no parameters to change)
n2 = 2000 # !!num of densities for binary potential (1 parameter to change: epsilon)
n3 = 2000 # !!num of densities for gravitational potential (1 parameter to change: mg)
n4 = 2000 # !!num of densities for power-law potential (3 parameters to change: u0,x0,a0)


# generate external potential and densities for each group
be_V_ext_record = matrix(NA,ncol=k,nrow = (n1+n2+n3+n4))
rho_record = matrix(NA,ncol=k,nrow = (n1+n2+n3+n4))
be_Omega_record = rep(NA,times=(n1+n2+n3+n4))


# initial density
rho0 = rep(0,k) 
rho0[(t+1):(k-t)]=0.2

# zero_infinity potential
set.seed(1)
para1 = maximinLHS(n=n1, k=1) # !!use LHD to generate parameter(s) for zero_infinity potential
be_mu1 = para1*3 #scale for beta*mu
V_zero_inf = get_V_zero_inf(a,L,d) #same for all beta*mu
for(i in 1:n1){
  print(i)
  be_V_ext_record[i,] = V_zero_inf
  exp_n_be_V_ext=exp(-V_zero_inf)
  
  rho_res=picard_rho_cpp(rho=rho0, beta_mu=be_mu1[i], a=a, L=L, exp_n_be_V_ext=exp_n_be_V_ext, alpha=alpha, iteration=maxIter, threshold = threshold)
  rho = rho_res[[2]]
  rho_record[i,]= rho

  beta_Omega = beta_Omega_cpp(rho=rho,a=a,L=L)
  be_Omega_record[i] = beta_Omega
}

# binary mixture potential
set.seed(1)
para2 = maximinLHS(n=n2, k=2) # !!use LHD to generate parameter(s) for binary mixture potential
be_mu2 = para2[,1]*3 #scale for beta*mu
epsilon2 = para2[,2]*2.2 # !!scale for epsilon
for(i in 1:n2){
  print(i)
  V_ext = get_V_binary(a,L,d,epsilon=epsilon2[i])
  be_V_ext_record[(n1+i),] = V_ext
  exp_n_be_V_ext=exp(-V_ext)
  
  rho_res=picard_rho_cpp(rho=rho0, beta_mu=be_mu2[i], a=a, L=L, exp_n_be_V_ext=exp_n_be_V_ext, alpha=alpha, iteration=maxIter, threshold = threshold)
  rho = rho_res[[2]]
  rho_record[(n1+i),]= rho

  beta_Omega = beta_Omega_cpp(rho=rho,a=a,L=L)
  be_Omega_record[n1+i] = beta_Omega
}

# gravitational potential
set.seed(1)
para3 = maximinLHS(n=n3, k=2) # !!use LHD to generate parameter(s) for gravitational mixture potential
be_mu3 = para3[,1]*3 #scale for beta*mu
mg3 = para3[,2]*3 # !!scale for mg
for(i in 1:n3){
  print(i)
  V_ext = get_V_grav(a,L,d,mg=mg3[i])
  be_V_ext_record[(n1+n2+i),] = V_ext
  exp_n_be_V_ext=exp(-V_ext)
  
  rho_res=picard_rho_cpp(rho=rho0, beta_mu=be_mu3[i], a=a, L=L, exp_n_be_V_ext=exp_n_be_V_ext, alpha=alpha, iteration=maxIter, threshold = threshold)
  rho = rho_res[[2]]
  rho_record[(n1+n2+i),] = rho

  beta_Omega = beta_Omega_cpp(rho=rho,a=a,L=L)
  be_Omega_record[n1+n2+i] = beta_Omega
}

# power-law potential
set.seed(1)
para4 = maximinLHS(n=n4, k=4) # !!use LHD to generate parameter(s) for power-law mixture potential
be_mu4 = para4[,1]*3 #scale for beta*mu
u0_4 = para4[,2]*2+1 # !!scale for m0
x0_4 = para4[,3]*2+1 # !!scale for x0
a0_4 = para4[,4]*3+2 # !!scale for a0
for(i in 1:n4){
  print(i)
  V_ext = get_V_power(a,L,d,u0=u0_4[i],x0=x0_4[i],a0=a0_4[i])
  be_V_ext_record[(n1+n2+n3+i),] = V_ext
  exp_n_be_V_ext=exp(-V_ext)
  
  rho_res=picard_rho_cpp(rho=rho0, beta_mu=be_mu4[i], a=a, L=L, exp_n_be_V_ext=exp_n_be_V_ext, alpha=alpha, iteration=maxIter, threshold = threshold)
  rho = rho_res[[2]]
  rho_record[(n1+n2+n3+i),]=rho

  beta_Omega = beta_Omega_cpp(rho=rho,a=a,L=L)
  be_Omega_record[n1+n2+n3+i] = beta_Omega
}

sum(is.na(rho_record))

# matplot(t(rho_record),type="l") #all densities
# matplot(t(rho_record[(1:n1),]),type="l",ylab="rho",main = "zero-inf") #densities of zero-infinity potential
# matplot(t(rho_record[(1:n2+n1),]),type="l",ylab="rho",main = "binary") #densities of binary potential
# matplot(t(rho_record[(1:n3+n1+n2),]),type="l",ylab="rho",main = "gravitational") #densities of gravitational potential
# matplot(t(rho_record[(1:n4+n1+n2+n3),]),type="l",ylab="rho",main = "power") #densities of power-law potential



