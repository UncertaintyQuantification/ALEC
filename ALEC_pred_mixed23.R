load("ini_4potentials_data.RData")
load("mix_23potentials_data.RData")
source("functions.R")
library(RobustGaSP)
library(Rcpp)
library(RcppEigen)
sourceCpp(file='src/picard.cpp')
sourceCpp(file='src/iterative_alg.cpp')

###############################################################
##################### Learn the new group #####################
###############################################################

##### start from the models built before

N=2000 #N for each
a=1
L=9
k=1001

delete_index=which(rho_record[1,]==0)

sd_threshold=0.01
alphaLevel=0.05

input5 = matrix(be_mu5,ncol = k-length(delete_index),nrow = N)-be_V_ext_mix_record[,-delete_index]
output5 = rho_mix_record[,-delete_index]

###################################################
############## use average criterion ##############
###################################################

#mix group 
n=20 #training
n_test = N-n
train_index = 1:n#as.integer(seq(1,N,length.out = n))

m_rho_ori5=build_GP(rho_record=rho_mix_record,be_V_ext_record=be_V_ext_mix_record,be_mu_record=be_mu5,N=N,n=n,
                    delete_index=delete_index,train_index=train_index,be_Omega_record=be_Omega_mix_record)

m_rho_ori5$RMSE

results5 = AL_GP_mean(testing_input=m_rho_ori5$testing_input, testing_output=m_rho_ori5$testing_output,
                      n_test=n_test, GP_before=m_rho_ori5$m_GP, n_before=n, sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results5$n_aug

results5$RMSE

#group 1
#split training and testing
n=20 #training
n_test = N-n
train_index = 1:n#as.integer(seq(1,N,length.out = n))

m_rho_ori1=orig_rho_GP(rho_record=rho_record,be_V_ext_record=be_V_ext_record,be_mu_record=be_mu_record,
                       group=1,N=N,n=n,delete_index=delete_index,train_index=train_index,
                       be_Omega_record=be_Omega_record)
m_rho_ori1$RMSE
m_rho_ori1$coverage
m_rho_ori1$length95

pred_ori1 = energy_pred_NRMSE(n_test=n_test, pred_rho_di=m_rho_ori1$m_pred$mean, 
                              a=a, L=L, k=k, 
                              be_Omega_test=m_rho_ori1$be_Omega_test, delete_index=delete_index, plot=T)
pred_ori1$be_Omega_RMSE

#augmented GP
results1 = AL_GP_mean(testing_input=m_rho_ori1$testing_input, testing_output=m_rho_ori1$testing_output,
                      n_test=n_test, GP_before=m_rho_ori1$m_GP, n_before=n, sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results1$n_aug

results1$RMSE
results1$coverage
results1$length95

results5with1 = AL_GP_mean(testing_input=input5, testing_output=output5,
                           n_test=N, GP_before=results1$m_GP_aug, n_before=n+results1$n_aug, 
                           sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results5with1$n_aug
results5with1$RMSE
results5with1$coverage
results5with1$length95

energy5with1 = energy_pred_NRMSE(n_test=N-results5with1$n_aug, pred_rho_di=results5with1$pred_record[results5with1$pred_index,], 
                                 a=a, L=L, k=k, 
                                 be_Omega_test=be_Omega_mix_record[results5with1$pred_index], delete_index=delete_index, plot=T)
energy5with1$be_Omega_RMSE


#group 2
#split training and testing
n=20 #training
n_test = N-n
train_index = 1:n#as.integer(seq(1,N,length.out = n))

m_rho_ori2=orig_rho_GP(rho_record=rho_record,be_V_ext_record=be_V_ext_record,be_mu_record=be_mu_record,
                       group=2,N=N,n=n,delete_index=delete_index,train_index=train_index,
                       be_Omega_record=be_Omega_record)
m_rho_ori2$RMSE
m_rho_ori2$coverage
m_rho_ori2$length95

#augmented GP
results2 = AL_GP_mean(testing_input=m_rho_ori2$testing_input, testing_output=m_rho_ori2$testing_output,
                            n_test=n_test, GP_before=m_rho_ori2$m_GP, n_before=n, sd_threshold=sd_threshold,alphaLevel=alphaLevel)

results2$n_aug
results2$RMSE
results2$coverage
results2$length95

results5with2 = AL_GP_mean(testing_input=input5, testing_output=output5,
                                 n_test=N, GP_before=results2$m_GP_aug, n_before=n+results2$n_aug, 
                                 sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results5with2$n_aug
results5with2$RMSE
results5with2$coverage
results5with2$length95

energy5with2 = energy_pred_NRMSE(n_test=N-results5with2$n_aug, pred_rho_di=results5with2$pred_record[results5with2$pred_index,], 
                                 a=a, L=L, k=k, 
                                 be_Omega_test=be_Omega_mix_record[results5with2$pred_index], delete_index=delete_index, plot=T)
energy5with2$be_Omega_RMSE

#group 3
#split training and testing
n=20 #training
n_test = N-n
train_index = 1:n#as.integer(seq(1,N,length.out = n))

m_rho_ori3=orig_rho_GP(rho_record=rho_record,be_V_ext_record=be_V_ext_record,be_mu_record=be_mu_record,
                       group=3,N=N,n=n,delete_index=delete_index,train_index=train_index,
                       be_Omega_record=be_Omega_record)
m_rho_ori3$RMSE
m_rho_ori3$coverage
m_rho_ori3$length95

pred_ori3 = energy_pred_NRMSE(n_test=n_test, pred_rho_di=m_rho_ori3$m_pred$mean, 
                              a=a, L=L, k=k, 
                              be_Omega_test=m_rho_ori3$be_Omega_test, delete_index=delete_index, plot=T)
pred_ori3$be_Omega_RMSE

#augmented GP
results3 = AL_GP_mean(testing_input=m_rho_ori3$testing_input, testing_output=m_rho_ori3$testing_output,
                            n_test=n_test, GP_before=m_rho_ori3$m_GP, n_before=n, sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results3$n_aug

results3$RMSE
results3$coverage
results3$length95

results5with3 = AL_GP_mean(testing_input=input5, testing_output=output5,
                                 n_test=N, GP_before=results3$m_GP_aug, n_before=n+results3$n_aug, 
                                 sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results5with3$n_aug
results5with3$RMSE
results5with3$coverage
results5with3$length95

energy5with3 = energy_pred_NRMSE(n_test=N-results5with3$n_aug, pred_rho_di=results5with3$pred_record[results5with3$pred_index,], 
                                 a=a, L=L, k=k,
                                 be_Omega_test=be_Omega_mix_record[results5with3$pred_index], delete_index=delete_index, plot=T)
energy5with3$be_Omega_RMSE

#group 4: power law
#split training and testing
n=20 #training
n_test = N-n
train_index = 1:n#as.integer(seq(1,N,length.out = n))

m_rho_ori4=orig_rho_GP(rho_record=rho_record,be_V_ext_record=be_V_ext_record,be_mu_record=be_mu_record,
                       group=4,N=N,n=n,delete_index=delete_index,train_index=train_index,
                       be_Omega_record=be_Omega_record)
m_rho_ori4$RMSE
m_rho_ori4$coverage
m_rho_ori4$length95

pred_ori4 = energy_pred_NRMSE(n_test=n_test, pred_rho_di=m_rho_ori4$m_pred$mean, 
                              a=a, L=L, k=k, 
                              be_Omega_test=m_rho_ori4$be_Omega_test, delete_index=delete_index, plot=T)
pred_ori4$be_Omega_RMSE

#augmented GP
results4 = AL_GP_mean(testing_input=m_rho_ori4$testing_input, testing_output=m_rho_ori4$testing_output,
                            n_test=n_test, GP_before=m_rho_ori4$m_GP, n_before=n, sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results4$n_aug
results4$RMSE 
results4$coverage 
results4$length95

results5with4 = AL_GP_mean(testing_input=input5, testing_output=output5,
                                 n_test=N, GP_before=results4$m_GP_aug, n_before=n+results4$n_aug, 
                                 sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results5with4$n_aug
results5with4$RMSE
results5with4$coverage
results5with4$length95

energy5with4 = energy_pred_NRMSE(n_test=N-results5with4$n_aug, pred_rho_di=results5with4$pred_record[results5with4$pred_index,], 
                                 a=a, L=L, k=k, 
                                 be_Omega_test=be_Omega_mix_record[results5with4$pred_index], delete_index=delete_index, plot=T)
energy5with4$be_Omega_RMSE

###all groups combined
n=80
n_test = 4*N-n
set.seed(1)
order_ind = sample(4*N)
train_index = 1:n

m_rho_ori=build_GP(rho_record=rho_record[order_ind,],be_V_ext_record=be_V_ext_record[order_ind,],be_mu_record=be_mu_record[order_ind],N=4*N,n=n,
                   delete_index=delete_index,train_index=train_index,be_Omega_record=be_Omega_record[order_ind])

#augmented GP
resultsall = AL_GP_mean(testing_input=m_rho_ori$testing_input, testing_output=m_rho_ori$testing_output,
                              n_test=n_test, GP_before=m_rho_ori$m_GP, n_before=n, sd_threshold=sd_threshold,alphaLevel=alphaLevel)

resultsall$n_aug

results5withall = AL_GP_mean(testing_input=input5, testing_output=output5,
                           n_test=N, GP_before=resultsall$m_GP_aug, n_before=n+resultsall$n_aug, 
                           sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results5withall$n_aug
results5withall$RMSE
results5withall$coverage
results5withall$length95

energy5withall = energy_pred_NRMSE(n_test=N-results5withall$n_aug, pred_rho_di=results5withall$pred_record[results5withall$pred_index,], 
                                 a=a, L=L, k=k, 
                                 be_Omega_test=be_Omega_mix_record[results5withall$pred_index], delete_index=delete_index, plot=T)
energy5withall$be_Omega_RMSE


###############################################
############## use max criterion ##############
###############################################
results5with1max = AL_GP_max(testing_input=input5, testing_output=output5,
                                        n_test=N, GP_before=results1$m_GP_aug, n_before=n+results1$n_aug, 
                                        sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results5with1max$n_aug

results5with2max = AL_GP_max(testing_input=input5, testing_output=output5,
                                 n_test=N, GP_before=results2$m_GP_aug, n_before=n+results2$n_aug, 
                                 sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results5with2max$n_aug

results5with3max = AL_GP_max(testing_input=input5, testing_output=output5,
                                 n_test=N, GP_before=results3$m_GP_aug, n_before=n+results3$n_aug, 
                                 sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results5with3max$n_aug

results5with4max = AL_GP_max(testing_input=input5, testing_output=output5,
                                 n_test=N, GP_before=results4$m_GP_aug, n_before=n+results4$n_aug, 
                                 sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results5with4max$n_aug

results5withallmax = AL_GP_max(testing_input=input5, testing_output=output5,
                                        n_test=N, GP_before=resultsall$m_GP_aug, n_before=4*n+resultsall$n_aug, 
                                        sd_threshold=sd_threshold,alphaLevel=alphaLevel)
results5withallmax$n_aug

###################################
##### create dataset for plot #####
###################################

#empirical cdf (Fig. 9)
#9 left
error5w1 = abs(results5with1$pred_record[results5with1$pred_index,]-output5[results5with1$pred_index,])
error5w1_sort_all=sort(error5w1,index.return = T)
index_plot5w1=c(seq(1,length(error5w1),200),length(error5w1))

error5w2 = abs(results5with2$pred_record[results5with2$pred_index,]-output5[results5with2$pred_index,])
error5w2_sort_all=sort(error5w2,index.return = T)
index_plot5w2=c(seq(1,length(error5w2),200),length(error5w2))

error5w3 = abs(results5with3$pred_record[results5with3$pred_index,]-output5[results5with3$pred_index,])
error5w3_sort_all=sort(error5w3,index.return = T)
index_plot5w3=c(seq(1,length(error5w3),200),length(error5w3))

error5w4 = abs(results5with4$pred_record[results5with4$pred_index,]-output5[results5with4$pred_index,])
error5w4_sort_all=sort(error5w4,index.return = T)
index_plot5w4=c(seq(1,length(error5w4),200),length(error5w4))

error5wall = abs(results5withall$pred_record[results5withall$pred_index,]-output5[results5withall$pred_index,])
error5wall_sort_all=sort(error5wall,index.return = T)
index_plot5wall=c(seq(1,length(error5wall),200),length(error5wall))

error5w1_plot = cbind(error5w1_sort_all$x[index_plot5w1],(1:length(error5w1))[index_plot5w1]/length(error5w1))
error5w2_plot = cbind(error5w2_sort_all$x[index_plot5w2],(1:length(error5w2))[index_plot5w2]/length(error5w2))
error5w3_plot = cbind(error5w3_sort_all$x[index_plot5w3],(1:length(error5w3))[index_plot5w3]/length(error5w3))
error5w4_plot = cbind(error5w4_sort_all$x[index_plot5w4],(1:length(error5w4))[index_plot5w4]/length(error5w4))
error5wall_plot = cbind(error5wall_sort_all$x[index_plot5wall],(1:length(error5wall))[index_plot5wall]/length(error5wall))

xpoints_mix = c(error5w1_sort_all$x[round(0.95*length(error5w1))],
                error5w2_sort_all$x[round(0.95*length(error5w2))],
                error5w3_sort_all$x[round(0.95*length(error5w3))],
                error5w4_sort_all$x[round(0.95*length(error5w4))],
                error5wall_sort_all$x[round(0.95*length(error5wall))])

#9 right
error5w1max = abs(results5with1max$pred_record[results5with1max$pred_index,]-output5[results5with1max$pred_index,])
error5w1_sort_max=sort(error5w1max,index.return = T)
index_plot5w1max=c(seq(1,length(error5w1max),200),length(error5w1max))

error5w2max = abs(results5with2max$pred_record[results5with2max$pred_index,]-output5[results5with2max$pred_index,])
error5w2_sort_max=sort(error5w2max,index.return = T)
index_plot5w2max=c(seq(1,length(error5w2max),200),length(error5w2max))

error5w3max = abs(results5with3max$pred_record[results5with3max$pred_index,]-output5[results5with3max$pred_index,])
error5w3_sort_max=sort(error5w3max,index.return = T)
index_plot5w3max=c(seq(1,length(error5w3max),200),length(error5w3max))

error5w4max = abs(results5with4max$pred_record[results5with4max$pred_index,]-output5[results5with4max$pred_index,])
error5w4_sort_max=sort(error5w4max,index.return = T)
index_plot5w4max=c(seq(1,length(error5w4max),200),length(error5w4max))

error5wallmax = abs(results5withallmax$pred_record[results5withallmax$pred_index,]-output5[results5withallmax$pred_index,])
error5wall_sort_max=sort(error5wallmax,index.return = T)
index_plot5wallmax=c(seq(1,length(error5wallmax),200),length(error5wallmax))

error5w1max_plot = cbind(error5w1_sort_max$x[index_plot5w1max],(1:length(error5w1max))[index_plot5w1max]/length(error5w1max))
error5w2max_plot = cbind(error5w2_sort_max$x[index_plot5w2max],(1:length(error5w2max))[index_plot5w2max]/length(error5w2max))
error5w3max_plot = cbind(error5w3_sort_max$x[index_plot5w3max],(1:length(error5w3max))[index_plot5w3max]/length(error5w3max))
error5w4max_plot = cbind(error5w4_sort_max$x[index_plot5w4max],(1:length(error5w4max))[index_plot5w4max]/length(error5w4max))
error5wallmax_plot = cbind(error5wall_sort_max$x[index_plot5wallmax],(1:length(error5wallmax))[index_plot5wallmax]/length(error5wallmax))

xpoints_mix_max = c(error5w1_sort_max$x[round(0.95*length(error5w1max))],
                    error5w2_sort_max$x[round(0.95*length(error5w2max))],
                    error5w3_sort_max$x[round(0.95*length(error5w3max))],
                    error5w4_sort_max$x[round(0.95*length(error5w4max))],
                    error5wall_sort_max$x[round(0.95*length(error5wallmax))])


#9 middle
error_col1=apply(error5w1,2,function(x) sum(x>sd_threshold)/dim(error5w1)[1])
error_col2=apply(error5w2,2,function(x) sum(x>sd_threshold)/dim(error5w2)[1])
error_col3=apply(error5w3,2,function(x) sum(x>sd_threshold)/dim(error5w3)[1])
error_col4=apply(error5w4,2,function(x) sum(x>sd_threshold)/dim(error5w4)[1])
error_colall=apply(error5wall,2,function(x) sum(x>sd_threshold)/dim(error5wall)[1])

# error_col1max=apply(error5w1max,2,function(x) sum(x>sd_threshold)/dim(error5w1max)[1])
# error_col2max=apply(error5w2max,2,function(x) sum(x>sd_threshold)/dim(error5w2max)[1])
# error_col3max=apply(error5w3max,2,function(x) sum(x>sd_threshold)/dim(error5w3max)[1])
# error_col4max=apply(error5w4max,2,function(x) sum(x>sd_threshold)/dim(error5w4max)[1])
# error_colallmax=apply(error5wallmax,2,function(x) sum(x>sd_threshold)/dim(error5wallmax)[1])
