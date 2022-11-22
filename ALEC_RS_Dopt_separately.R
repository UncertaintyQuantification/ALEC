load("ini_4potentials_data.RData")
source("functions.R")
library(RobustGaSP)
library(Rcpp)
library(RcppEigen)
sourceCpp(file='src/picard.cpp')
sourceCpp(file='src/iterative_alg.cpp')

#######################################################################
##################### Learn each group separately #####################
#######################################################################

##################### RS1, RS2, ALEC #####################

######## 1) for each group, compare AL and GP########
N=2000 #N for each
a=1
L=9
k=1001

delete_index=which(rho_record[1,]==0)

sd_threshold=0.01
alphaLevel=0.05
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

energy1 = energy_pred_NRMSE(n_test=n_test-results1$n_aug, pred_rho_di=results1$pred_record[results1$pred_index,], 
                            a=a, L=L, k=k, 
                            be_Omega_test=m_rho_ori1$be_Omega_test[results1$pred_index], delete_index=delete_index, plot=T)
energy1$be_Omega_RMSE

#compare with GP with same number of training
train_index_comp = 1:(n+results1$n_aug)#as.integer(seq(1,N,length.out = n+results1$n_aug))

m_rho_comp1=orig_rho_GP(rho_record=rho_record,be_V_ext_record=be_V_ext_record,be_mu_record=be_mu_record,
                        group=1,N=N,n=n+results1$n_aug,delete_index=delete_index,train_index=train_index_comp,
                        be_Omega_record=be_Omega_record)
m_rho_comp1$RMSE
m_rho_comp1$coverage 
m_rho_comp1$length95

energy_comp1 = energy_pred_NRMSE(n_test=n_test-results1$n_aug, pred_rho_di=m_rho_comp1$pred_rho_di, 
                                 a=a, L=L, k=k, 
                                 be_Omega_test=m_rho_comp1$be_Omega_test, delete_index=delete_index, plot=T)
energy_comp1$be_Omega_RMSE


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

pred_ori2 = energy_pred_NRMSE(n_test=n_test, pred_rho_di=m_rho_ori2$m_pred$mean, 
                              a=a, L=L, k=k, 
                              be_Omega_test=m_rho_ori2$be_Omega_test, delete_index=delete_index, plot=T)
pred_ori2$be_Omega_RMSE

#augmented GP
results2 = AL_GP_mean(testing_input=m_rho_ori2$testing_input, testing_output=m_rho_ori2$testing_output,
                            n_test=n_test, GP_before=m_rho_ori2$m_GP, n_before=n, sd_threshold=sd_threshold,alphaLevel=alphaLevel)

results2$n_aug

results2$RMSE
results2$coverage
results2$length95

energy2 = energy_pred_NRMSE(n_test=n_test-results2$n_aug, pred_rho_di=results2$pred_record[results2$pred_index,], 
                            a=a, L=L, k=k, 
                            be_Omega_test=m_rho_ori2$be_Omega_test[results2$pred_index], delete_index=delete_index, plot=T)
energy2$be_Omega_RMSE


#compare with GP with same number of training
train_index_comp = 1:(n+results2$n_aug)#as.integer(seq(1,N,length.out = n+results2$n_aug))

m_rho_comp2=orig_rho_GP(rho_record=rho_record,be_V_ext_record=be_V_ext_record,be_mu_record=be_mu_record,
                        group=2,N=N,n=n+results2$n_aug,delete_index=delete_index,train_index=train_index_comp,
                        be_Omega_record=be_Omega_record)
m_rho_comp2$RMSE
m_rho_comp2$coverage 
m_rho_comp2$length95

energy_comp2 = energy_pred_NRMSE(n_test=n_test-results2$n_aug, pred_rho_di=m_rho_comp2$pred_rho_di, 
                                 a=a, L=L, k=k, 
                                 be_Omega_test=m_rho_comp2$be_Omega_test, delete_index=delete_index, plot=T)
energy_comp2$be_Omega_RMSE


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
#matplot(t(m_rho_ori3$rho_train),type="l")
#matplot(t(m_rho_ori3$rho_test),type="l")

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

energy3 = energy_pred_NRMSE(n_test=n_test-results3$n_aug, pred_rho_di=results3$pred_record[results3$pred_index,], 
                            a=a, L=L, k=k,
                            be_Omega_test=m_rho_ori3$be_Omega_test[results3$pred_index], delete_index=delete_index, plot=T)
energy3$be_Omega_RMSE

#compare with GP with same number of training
train_index_comp = 1:(n+results3$n_aug)#as.integer(seq(1,N,length.out = n+results3$n_aug))

m_rho_comp3=orig_rho_GP(rho_record=rho_record,be_V_ext_record=be_V_ext_record,be_mu_record=be_mu_record,
                        group=3,N=N,n=n+results3$n_aug,delete_index=delete_index,train_index=train_index_comp,
                        be_Omega_record=be_Omega_record)
m_rho_comp3$RMSE
m_rho_comp3$coverage 
m_rho_comp3$length95

energy_comp3 = energy_pred_NRMSE(n_test=n_test-results3$n_aug, pred_rho_di=m_rho_comp3$pred_rho_di, 
                                 a=a, L=L, k=k, 
                                 be_Omega_test=m_rho_comp3$be_Omega_test, delete_index=delete_index, plot=T)
energy_comp3$be_Omega_RMSE


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

energy4 = energy_pred_NRMSE(n_test=n_test-results4$n_aug, pred_rho_di=results4$pred_record[results4$pred_index,], 
                            a=a, L=L, k=k, 
                            be_Omega_test=m_rho_ori4$be_Omega_test[results4$pred_index], delete_index=delete_index, plot=T)
energy4$be_Omega_RMSE


#compare with GP with same number of training
train_index_comp = 1:(n+results4$n_aug)#as.integer(seq(1,N,length.out = n+results4$n_aug))

m_rho_comp4=orig_rho_GP(rho_record=rho_record,be_V_ext_record=be_V_ext_record,be_mu_record=be_mu_record,
                        group=4,N=N,n=n+results4$n_aug,delete_index=delete_index,train_index=train_index_comp,
                        be_Omega_record=be_Omega_record)
m_rho_comp4$RMSE
m_rho_comp4$coverage 
m_rho_comp4$length95

energy_comp4 = energy_pred_NRMSE(n_test=n_test-results4$n_aug, pred_rho_di=m_rho_comp4$pred_rho_di, 
                                 a=a, L=L, k=k,
                                 be_Omega_test=m_rho_comp4$be_Omega_test, delete_index=delete_index, plot=T)
energy_comp4$be_Omega_RMSE

# m_rho_comp4$m_GP@beta_hat 
# m_rho_ori4$m_GP@beta_hat 
# results4$m_GP_aug@beta_hat

##################### D-optimality #####################

results_opt1 = augmented_GP_Doptimal_cpp(testing_input=m_rho_ori1$testing_input, testing_output=m_rho_ori1$testing_output,
                                         n_test=n_test, GP_before=m_rho_ori1$m_GP, n_before=n, D_threshold=30)
results_opt1$n_aug
results_opt1$RMSE 
results_opt1$coverage 
results_opt1$length95
energy_opt1 = energy_pred_NRMSE(n_test=n_test-results_opt1$n_aug, pred_rho_di=results_opt1$pred_record[results_opt1$pred_index,], 
                                a=a, L=L, k=k, 
                                be_Omega_test=m_rho_ori1$be_Omega_test[results_opt1$pred_index], delete_index=delete_index, plot=T)
energy_opt1$be_Omega_RMSE

results_opt2 = augmented_GP_Doptimal_cpp(testing_input=m_rho_ori2$testing_input, testing_output=m_rho_ori2$testing_output,
                                         n_test=n_test, GP_before=m_rho_ori2$m_GP, n_before=n, D_threshold=3.4)
results_opt2$n_aug
results_opt2$RMSE 
results_opt2$coverage 
results_opt2$length95
energy_opt2 = energy_pred_NRMSE(n_test=n_test-results_opt2$n_aug, pred_rho_di=results_opt2$pred_record[results_opt2$pred_index,], 
                                a=a, L=L, k=k, 
                                be_Omega_test=m_rho_ori2$be_Omega_test[results_opt2$pred_index], delete_index=delete_index, plot=T)
energy_opt2$be_Omega_RMSE

results_opt3 = augmented_GP_Doptimal_cpp(testing_input=m_rho_ori3$testing_input, testing_output=m_rho_ori3$testing_output,
                                         n_test=n_test, GP_before=m_rho_ori3$m_GP, n_before=n, D_threshold=2.4)
results_opt3$n_aug
results_opt3$RMSE 
results_opt3$coverage 
results_opt3$length95
energy_opt3 = energy_pred_NRMSE(n_test=n_test-results_opt3$n_aug, pred_rho_di=results_opt3$pred_record[results_opt3$pred_index,], 
                                a=a, L=L, k=k,
                                be_Omega_test=m_rho_ori3$be_Omega_test[results_opt3$pred_index], delete_index=delete_index, plot=T)
energy_opt3$be_Omega_RMSE

results_opt4 = augmented_GP_Doptimal_cpp(testing_input=m_rho_ori4$testing_input, testing_output=m_rho_ori4$testing_output,
                                         n_test=n_test, GP_before=m_rho_ori4$m_GP, n_before=n, D_threshold=1.03)
results_opt4$n_aug
results_opt4$RMSE 
results_opt4$coverage 
results_opt4$length95
energy_opt4 = energy_pred_NRMSE(n_test=n_test-results_opt4$n_aug, pred_rho_di=results_opt4$pred_record[results_opt4$pred_index,], 
                                a=a, L=L, k=k, 
                                be_Omega_test=m_rho_ori4$be_Omega_test[results_opt4$pred_index], delete_index=delete_index, plot=T)
energy_opt4$be_Omega_RMSE

###################################
##### create dataset for plot #####
###################################

#for Fig. 4, using fewer points to plot to make it faster
error1 = abs(results1$pred_record[results1$pred_index,]-m_rho_ori1$testing_output[results1$pred_index,])
error1_sort_all=sort(error1,index.return = T)
index_plot1=c(seq(1,length(error1),200),length(error1))
error2 = abs(results2$pred_record[results2$pred_index,]-m_rho_ori2$testing_output[results2$pred_index,])
error2_sort_all=sort(error2,index.return = T)
index_plot2=c(seq(1,length(error2),200),length(error2))
error3 = abs(results3$pred_record[results3$pred_index,]-m_rho_ori3$testing_output[results3$pred_index,])
error3_sort_all=sort(error3,index.return = T)
index_plot3=c(seq(1,length(error3),200),length(error3))
error4 = abs(results4$pred_record[results4$pred_index,]-m_rho_ori4$testing_output[results4$pred_index,])
error4_sort_all=sort(error4,index.return = T)
index_plot4=c(seq(1,length(error4),200),length(error4))

error4_plot=cbind(error4_sort_all$x[index_plot4],(1:length(error4))[index_plot4]/length(error4))
error3_plot=cbind(error3_sort_all$x[index_plot3],(1:length(error3))[index_plot3]/length(error3))
error2_plot=cbind(error2_sort_all$x[index_plot2],(1:length(error2))[index_plot2]/length(error2))
error1_plot=cbind(error1_sort_all$x[index_plot1],(1:length(error1))[index_plot1]/length(error1))
xpoints = c(error1_sort_all$x[round(0.95*length(error1))],
            error2_sort_all$x[round(0.95*length(error2))],
            error3_sort_all$x[round(0.95*length(error3))],
            error4_sort_all$x[round(0.95*length(error4))])


#for Fig. 5
dat_sep = data.frame(RMSE_rho = c(m_rho_ori1$RMSE,results1$RMSE,m_rho_comp1$RMSE,results_opt1$RMSE,
                                  m_rho_ori2$RMSE,results2$RMSE,m_rho_comp2$RMSE,results_opt2$RMSE,
                                  m_rho_ori3$RMSE,results3$RMSE,m_rho_comp3$RMSE,results_opt3$RMSE,
                                  m_rho_ori4$RMSE,results4$RMSE,m_rho_comp4$RMSE,results_opt4$RMSE),
                     RMSE_Omega = c(pred_ori1$be_Omega_RMSE,energy1$be_Omega_RMSE,energy_comp1$be_Omega_RMSE,energy_opt1$be_Omega_RMSE,
                                    pred_ori2$be_Omega_RMSE,energy2$be_Omega_RMSE,energy_comp2$be_Omega_RMSE,energy_opt2$be_Omega_RMSE,
                                    pred_ori3$be_Omega_RMSE,energy3$be_Omega_RMSE,energy_comp3$be_Omega_RMSE,energy_opt3$be_Omega_RMSE,
                                    pred_ori4$be_Omega_RMSE,energy4$be_Omega_RMSE,energy_comp4$be_Omega_RMSE,energy_opt4$be_Omega_RMSE),
                     n_train = c(as.character(n),paste(n,"+",results1$n_aug,sep=""),as.character(n+results1$n_aug),paste(n,"+",results_opt1$n_aug,sep=""),
                                 as.character(n),paste(n,"+",results2$n_aug,sep=""),as.character(n+results2$n_aug),paste(n,"+",results_opt2$n_aug,sep=""),
                                 as.character(n),paste(n,"+",results3$n_aug,sep=""),as.character(n+results3$n_aug),paste(n,"+",results_opt3$n_aug,sep=""),
                                 as.character(n),paste(n,"+",results4$n_aug,sep=""),as.character(n+results4$n_aug),paste(n,"+",results_opt4$n_aug,sep="")),
                     Group = rep(c("group 1","group 2","group 3","group 4"),each=4),
                     Design = factor(rep(c("RS1","ALEC","RS2","D-opt"),4),levels = c("ALEC","RS1","RS2","D-opt")))


#Group 4 comparison (Fig. 6)
#6(a)
RMSE4 = data.frame(RMSE = c(m_rho_ori4$RMSE_each,m_rho_comp4$RMSE_each,results4$RMSE_each,results_opt4$RMSE_each),
                   design = c(rep("RS1",length(m_rho_ori4$RMSE_each)),rep("RS2",length(m_rho_comp4$RMSE_each)),rep("ALEC",length(results4$RMSE_each)),rep("D-opt",length(results_opt4$RMSE_each))))
RMSE4$design <- factor(RMSE4$design , levels=c("ALEC","RS1", "RS2","D-opt"))

#6(b)
LHD2_ind1=596
LHD1_ind1 = which(pred_ori4$be_Omega_test==energy_comp4$be_Omega_test[LHD2_ind1])
AL_ind1 = which(energy4$be_Omega_test==energy_comp4$be_Omega_test[LHD2_ind1])
Dopt_ind1 = which(energy_opt4$be_Omega_test==energy_comp4$be_Omega_test[LHD2_ind1])
LHD2_ind2=798
LHD1_ind2 = which(pred_ori4$be_Omega_test==energy_comp4$be_Omega_test[LHD2_ind2])
AL_ind2 = which(energy4$be_Omega_test==energy_comp4$be_Omega_test[LHD2_ind2])
Dopt_ind2 = which(energy_opt4$be_Omega_test==energy_comp4$be_Omega_test[LHD2_ind2])

energy4_RS1 = cbind(pred_ori4$be_Omega_test,pred_ori4$pred_be_Omega)
energy4_RS2 = cbind(energy_comp4$be_Omega_test,energy_comp4$pred_be_Omega)
energy4_Dopt= cbind(energy_opt4$be_Omega_test,energy_opt4$pred_be_Omega)
energy4_ALEC= cbind(energy4$be_Omega_test,energy4$pred_be_Omega)
select4 = rbind(c(pred_ori4$be_Omega_test[LHD1_ind1],pred_ori4$pred_be_Omega[LHD1_ind1]),
                c(pred_ori4$be_Omega_test[LHD1_ind2],pred_ori4$pred_be_Omega[LHD1_ind2]))

##6(c)
density1 = rbind(m_rho_ori4$rho_test[LHD1_ind1,-c(1:50,952:1001)],
                 pred_ori4$pred_rho[LHD1_ind1,-c(1:50,952:1001)],
                 energy_comp4$pred_rho[LHD2_ind1,-c(1:50,952:1001)],
                 energy_opt4$pred_rho[Dopt_ind1,-c(1:50,952:1001)],
                 energy4$pred_rho[AL_ind1,-c(1:50,952:1001)])
rownames(density1)=c("Truth","RS1","RS2","Dopt","ALEC")

#6(d)
density2 = rbind(m_rho_ori4$rho_test[LHD1_ind2,-c(1:50,952:1001)],
                 pred_ori4$pred_rho[LHD1_ind2,-c(1:50,952:1001)],
                 energy_comp4$pred_rho[LHD2_ind2,-c(1:50,952:1001)],
                 energy_opt4$pred_rho[Dopt_ind2,-c(1:50,952:1001)],
                 energy4$pred_rho[AL_ind2,-c(1:50,952:1001)])
rownames(density2)=c("Truth","RS1","RS2","Dopt","ALEC")

