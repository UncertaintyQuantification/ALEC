load("ini_4potentials_data.RData")
source("functions.R")
library(RobustGaSP)
library(Rcpp)
library(RcppEigen)
sourceCpp(file='src/picard.cpp')
sourceCpp(file='src/iterative_alg.cpp')

#####################################################################
##################### Learn all groups together #####################
#####################################################################

######## build model together, compare ALEC with RS and D-opt########
N=2000 #N for each
a=1
L=9
k=1001

delete_index=which(rho_record[1,]==0)

sd_threshold=0.01
alphaLevel=0.05

n=80
n_test = 4*N-n
#RS1_ind = c(1:20,1:20+N,1:20+2*N,1:20+3*N)
set.seed(1)
order_ind = sample(4*N)
train_index = 1:n

m_rho_ori=build_GP(rho_record=rho_record[order_ind,],be_V_ext_record=be_V_ext_record[order_ind,],be_mu_record=be_mu_record[order_ind],N=4*N,n=n,
                    delete_index=delete_index,train_index=train_index,be_Omega_record=be_Omega_record[order_ind])
m_rho_ori$RMSE
m_rho_ori$coverage
m_rho_ori$length95

pred_ori = energy_pred_NRMSE(n_test=n_test, pred_rho_di=m_rho_ori$m_pred$mean, 
                              a=a, L=L, k=k,  
                              be_Omega_test=m_rho_ori$be_Omega_test, delete_index=delete_index, plot=T)
pred_ori$be_Omega_RMSE


#augmented GP
resultsall = AL_GP_mean(testing_input=m_rho_ori$testing_input, testing_output=m_rho_ori$testing_output,
                        n_test=n_test, GP_before=m_rho_ori$m_GP, n_before=n, sd_threshold=sd_threshold,alphaLevel=alphaLevel)

resultsall$n_aug
resultsall$RMSE
resultsall$coverage
resultsall$length95

energyall = energy_pred_NRMSE(n_test=n_test-resultsall$n_aug, pred_rho_di=resultsall$pred_record[resultsall$pred_index,], 
                               a=a, L=L, k=k, 
                               be_Omega_test=m_rho_ori$be_Omega_test[resultsall$pred_index], delete_index=delete_index, plot=T)

energyall$be_Omega_RMSE

#compare with GP with same number of training
train_index_comp = 1:(n+resultsall$n_aug)#as.integer(seq(1,4*N,length.out = n+resultsall$n_aug))

m_rho_comp=build_GP(rho_record=rho_record[order_ind,],be_V_ext_record=be_V_ext_record[order_ind,],be_mu_record=be_mu_record[order_ind],
                    N=4*N,n=n+resultsall$n_aug,delete_index=delete_index,train_index=train_index_comp,
                    be_Omega_record=be_Omega_record[order_ind])
m_rho_comp$RMSE
m_rho_comp$coverage 
m_rho_comp$length95

energy_comp = energy_pred_NRMSE(n_test=n_test-resultsall$n_aug, pred_rho_di=m_rho_comp$pred_rho_di, 
                                 a=a, L=L, k=k, 
                                 be_Omega_test=m_rho_comp$be_Omega_test, delete_index=delete_index, plot=T)
energy_comp$be_Omega_RMSE


#D-optimality
resultsall_opt = augmented_GP_Doptimal_cpp(testing_input=m_rho_ori$testing_input, testing_output=m_rho_ori$testing_output,
                                           n_test=n_test, GP_before=m_rho_ori$m_GP, n_before=n, D_threshold = 1.12)#1.2

resultsall_opt$n_aug
resultsall_opt$RMSE
resultsall_opt$coverage
resultsall_opt$length95

energyall_opt = energy_pred_NRMSE(n_test=n_test-resultsall_opt$n_aug, pred_rho_di=resultsall_opt$pred_record[resultsall_opt$pred_index,], 
                              a=a, L=L, k=k, 
                              be_Omega_test=m_rho_ori$be_Omega_test[resultsall_opt$pred_index], delete_index=delete_index, plot=T)

energyall_opt$be_Omega_RMSE



#input decomposition
resultsall_decomp = AL_input_decomp(testing_input=m_rho_ori$testing_input, testing_output=m_rho_ori$testing_output,
                              n_test=n_test, GP_before=m_rho_ori$m_GP, n_before=n, sd_threshold=sd_threshold,alphaLevel=alphaLevel)

resultsall_decomp$n_aug
resultsall_decomp$RMSE
resultsall_decomp$coverage
resultsall_decomp$length95

resultsall_decomp$n_cluster
resultsall_decomp$model_list[[1]]@num_obs
resultsall_decomp$model_list[[2]]@num_obs
resultsall_decomp$model_list[[3]]@num_obs
resultsall_decomp$model_list[[4]]@num_obs


energyall_decomp = energy_pred_NRMSE(n_test=n_test-resultsall_decomp$n_aug, pred_rho_di=resultsall_decomp$pred_record[resultsall_decomp$pred_index,], 
                              a=a, L=L, k=k, 
                              be_Omega_test=m_rho_ori$be_Omega_test[resultsall_decomp$pred_index], delete_index=delete_index, plot=T)

energyall_decomp$be_Omega_RMSE

#####################################################
##################### get table #####################
#####################################################

#check the performance in each group

########RS1 (lower case rs represents RS1)
rs_index=order_ind[train_index]
n1=sum(rs_index<=N)
n2=sum(rs_index>N & rs_index<=2*N)
n3=sum(rs_index>2*N & rs_index<=3*N)
n4=sum(rs_index>3*N & rs_index<=4*N)

pred_ind_rs = order_ind[-train_index]
grp1_pred_ind_rs = which(pred_ind_rs<=N)
rmse1_rs = sqrt(mean((m_rho_ori$pred_rho_di[grp1_pred_ind_rs,]-m_rho_ori$testing_output[grp1_pred_ind_rs,])^2))
cov1_rs = sum(m_rho_ori$testing_output[grp1_pred_ind_rs,]>=m_rho_ori$m_pred$lower95[grp1_pred_ind_rs,] & m_rho_ori$testing_output[grp1_pred_ind_rs,] <= m_rho_ori$m_pred$upper95[grp1_pred_ind_rs,])/length(m_rho_ori$m_pred$upper95[grp1_pred_ind_rs,])
l95_1_rs = mean(m_rho_ori$m_pred$upper95[grp1_pred_ind_rs,]-m_rho_ori$m_pred$lower95[grp1_pred_ind_rs,])
rmseO1_rs = sqrt(mean((pred_ori$pred_be_Omega[grp1_pred_ind_rs]-pred_ori$be_Omega_test[grp1_pred_ind_rs])^2))

grp2_pred_ind_rs = which(pred_ind_rs>N&pred_ind_rs<=2*N)
rmse2_rs = sqrt(mean((m_rho_ori$pred_rho_di[grp2_pred_ind_rs,]-m_rho_ori$testing_output[grp2_pred_ind_rs,])^2))
cov2_rs = sum(m_rho_ori$testing_output[grp2_pred_ind_rs,]>=m_rho_ori$m_pred$lower95[grp2_pred_ind_rs,] & m_rho_ori$testing_output[grp2_pred_ind_rs,] <= m_rho_ori$m_pred$upper95[grp2_pred_ind_rs,])/length(m_rho_ori$m_pred$upper95[grp2_pred_ind_rs,])
l95_2_rs = mean(m_rho_ori$m_pred$upper95[grp2_pred_ind_rs,]-m_rho_ori$m_pred$lower95[grp2_pred_ind_rs,])
rmseO2_rs = sqrt(mean((pred_ori$pred_be_Omega[grp2_pred_ind_rs]-pred_ori$be_Omega_test[grp2_pred_ind_rs])^2))

grp3_pred_ind_rs = which(pred_ind_rs>2*N&pred_ind_rs<=3*N)
rmse3_rs = sqrt(mean((m_rho_ori$pred_rho_di[grp3_pred_ind_rs,]-m_rho_ori$testing_output[grp3_pred_ind_rs,])^2))
cov3_rs = sum(m_rho_ori$testing_output[grp3_pred_ind_rs,]>=m_rho_ori$m_pred$lower95[grp3_pred_ind_rs,] & m_rho_ori$testing_output[grp3_pred_ind_rs,] <= m_rho_ori$m_pred$upper95[grp3_pred_ind_rs,])/length(m_rho_ori$m_pred$upper95[grp3_pred_ind_rs,])
l95_3_rs = mean(m_rho_ori$m_pred$upper95[grp3_pred_ind_rs,]-m_rho_ori$m_pred$lower95[grp3_pred_ind_rs,])
rmseO3_rs = sqrt(mean((pred_ori$pred_be_Omega[grp3_pred_ind_rs]-pred_ori$be_Omega_test[grp3_pred_ind_rs])^2))

grp4_pred_ind_rs =  which(pred_ind_rs>3*N&pred_ind_rs<=4*N)
rmse4_rs = sqrt(mean((m_rho_ori$pred_rho_di[grp4_pred_ind_rs,]-m_rho_ori$testing_output[grp4_pred_ind_rs,])^2))
cov4_rs = sum(m_rho_ori$testing_output[grp4_pred_ind_rs,]>=m_rho_ori$m_pred$lower95[grp4_pred_ind_rs,] & m_rho_ori$testing_output[grp4_pred_ind_rs,] <= m_rho_ori$m_pred$upper95[grp4_pred_ind_rs,])/length(m_rho_ori$m_pred$upper95[grp4_pred_ind_rs,])
l95_4_rs = mean(m_rho_ori$m_pred$upper95[grp4_pred_ind_rs,]-m_rho_ori$m_pred$lower95[grp4_pred_ind_rs,])
rmseO4_rs = sqrt(mean((pred_ori$pred_be_Omega[grp4_pred_ind_rs]-pred_ori$be_Omega_test[grp4_pred_ind_rs])^2))

rs1 = matrix(c(n1,n2,n3,n4,
               rmse1_rs,rmse2_rs,rmse3_rs,rmse4_rs,
               cov1_rs,cov2_rs,cov3_rs,cov4_rs,
               l95_1_rs,l95_2_rs,l95_3_rs,l95_4_rs,
               rmseO1_rs,rmseO2_rs,rmseO3_rs,rmseO4_rs),
             nr=4,dimnames = list(c("Grp1", "Grp2","Grp3","Grp4"),
                                  c("n", "RMSE_rho", "coverage", "length","RMSE_Omega")))
rs1

#######ALEC
aug_index = order_ind[-train_index][-resultsall$pred_index]
n_aug1=sum(aug_index<=N)
n_aug2=sum(aug_index>N & aug_index<=2*N)
n_aug3=sum(aug_index>2*N & aug_index<=3*N)
n_aug4=sum(aug_index>3*N & aug_index<=4*N)

pred_ind_ori = order_ind[-train_index][resultsall$pred_index]
grp1_pred_ind = which(pred_ind_ori<=N)
rmse1 = sqrt(mean((resultsall$pred_record[resultsall$pred_index,][grp1_pred_ind,]-m_rho_ori$testing_output[resultsall$pred_index,][grp1_pred_ind,])^2))
cov1 = sum(m_rho_ori$testing_output[resultsall$pred_index,][grp1_pred_ind,]>=resultsall$lb95_record[resultsall$pred_index,][grp1_pred_ind,] & m_rho_ori$testing_output[resultsall$pred_index,][grp1_pred_ind,]<= resultsall$ub95_record[resultsall$pred_index,][grp1_pred_ind,])/length(resultsall$ub95_record[resultsall$pred_index,][grp1_pred_ind,])
l95_1 = mean(resultsall$ub95_record[resultsall$pred_index,][grp1_pred_ind,]-resultsall$lb95_record[resultsall$pred_index,][grp1_pred_ind,])
rmseO1=sqrt(mean((energyall$pred_be_Omega[grp1_pred_ind]-energyall$be_Omega_test[grp1_pred_ind])^2))

grp2_pred_ind = which(pred_ind_ori>N & pred_ind_ori<=2*N)
rmse2 = sqrt(mean((resultsall$pred_record[resultsall$pred_index,][grp2_pred_ind,]-m_rho_ori$testing_output[resultsall$pred_index,][grp2_pred_ind,])^2))
cov2 = sum(m_rho_ori$testing_output[resultsall$pred_index,][grp2_pred_ind,]>=resultsall$lb95_record[resultsall$pred_index,][grp2_pred_ind,] & m_rho_ori$testing_output[resultsall$pred_index,][grp2_pred_ind,]<= resultsall$ub95_record[resultsall$pred_index,][grp2_pred_ind,])/length(resultsall$ub95_record[resultsall$pred_index,][grp2_pred_ind,])
l95_2 = mean(resultsall$ub95_record[resultsall$pred_index,][grp2_pred_ind,]-resultsall$lb95_record[resultsall$pred_index,][grp2_pred_ind,])
rmseO2=sqrt(mean((energyall$pred_be_Omega[grp2_pred_ind]-energyall$be_Omega_test[grp2_pred_ind])^2))

grp3_pred_ind = which(pred_ind_ori>2*N & pred_ind_ori<=3*N)
rmse3 = sqrt(mean((resultsall$pred_record[resultsall$pred_index,][grp3_pred_ind,]-m_rho_ori$testing_output[resultsall$pred_index,][grp3_pred_ind,])^2))
cov3 = sum(m_rho_ori$testing_output[resultsall$pred_index,][grp3_pred_ind,]>=resultsall$lb95_record[resultsall$pred_index,][grp3_pred_ind,] & m_rho_ori$testing_output[resultsall$pred_index,][grp3_pred_ind,]<= resultsall$ub95_record[resultsall$pred_index,][grp3_pred_ind,])/length(resultsall$ub95_record[resultsall$pred_index,][grp3_pred_ind,])
l95_3 = mean(resultsall$ub95_record[resultsall$pred_index,][grp3_pred_ind,]-resultsall$lb95_record[resultsall$pred_index,][grp3_pred_ind,])
rmseO3=sqrt(mean((energyall$pred_be_Omega[grp3_pred_ind]-energyall$be_Omega_test[grp3_pred_ind])^2))

grp4_pred_ind = which(pred_ind_ori>3*N & pred_ind_ori<=4*N)
rmse4 = sqrt(mean((resultsall$pred_record[resultsall$pred_index,][grp4_pred_ind,]-m_rho_ori$testing_output[resultsall$pred_index,][grp4_pred_ind,])^2))
cov4 = sum(m_rho_ori$testing_output[resultsall$pred_index,][grp4_pred_ind,]>=resultsall$lb95_record[resultsall$pred_index,][grp4_pred_ind,] & m_rho_ori$testing_output[resultsall$pred_index,][grp4_pred_ind,]<= resultsall$ub95_record[resultsall$pred_index,][grp4_pred_ind,])/length(resultsall$ub95_record[resultsall$pred_index,][grp4_pred_ind,])
l95_4 = mean(resultsall$ub95_record[resultsall$pred_index,][grp4_pred_ind,]-resultsall$lb95_record[resultsall$pred_index,][grp4_pred_ind,])
rmseO4=sqrt(mean((energyall$pred_be_Omega[grp4_pred_ind]-energyall$be_Omega_test[grp4_pred_ind])^2))

subgroup = matrix(c(n_aug1,n_aug2,n_aug3,n_aug4,
                    rmse1,rmse2,rmse3,rmse4,
                    cov1,cov2,cov3,cov4,
                    l95_1,l95_2,l95_3,l95_4,
                    rmseO1,rmseO2,rmseO3,rmseO4),
                  nr=4,dimnames = list(c("Grp1", "Grp2","Grp3","Grp4"),
                                       c("n_aug", "RMSE_rho", "coverage", "length","RMSE_Omega")))
subgroup

########RS2 (upper case RS represents RS2)
RS_index=order_ind[train_index_comp]
n1_comp=sum(RS_index<=N)
n2_comp=sum(RS_index>N & RS_index<=2*N)
n3_comp=sum(RS_index>2*N & RS_index<=3*N)
n4_comp=sum(RS_index>3*N & RS_index<=4*N)

pred_ind_RS = order_ind[-train_index_comp]
grp1_pred_ind_RS = which(pred_ind_RS<=N)
rmse1_RS = sqrt(mean((m_rho_comp$pred_rho_di[grp1_pred_ind_RS,]-m_rho_comp$testing_output[grp1_pred_ind_RS,])^2))
cov1_RS = sum(m_rho_comp$testing_output[grp1_pred_ind_RS,]>=m_rho_comp$m_pred$lower95[grp1_pred_ind_RS,] & m_rho_comp$testing_output[grp1_pred_ind_RS,] <= m_rho_comp$m_pred$upper95[grp1_pred_ind_RS,])/length(m_rho_comp$m_pred$upper95[grp1_pred_ind_RS,])
l95_1_RS = mean(m_rho_comp$m_pred$upper95[grp1_pred_ind_RS,]-m_rho_comp$m_pred$lower95[grp1_pred_ind_RS,])
rmseO1_RS = sqrt(mean((energy_comp$pred_be_Omega[grp1_pred_ind_RS]-energy_comp$be_Omega_test[grp1_pred_ind_RS])^2))

grp2_pred_ind_RS = which(pred_ind_RS>N&pred_ind_RS<=2*N)
rmse2_RS = sqrt(mean((m_rho_comp$pred_rho_di[grp2_pred_ind_RS,]-m_rho_comp$testing_output[grp2_pred_ind_RS,])^2))
cov2_RS = sum(m_rho_comp$testing_output[grp2_pred_ind_RS,]>=m_rho_comp$m_pred$lower95[grp2_pred_ind_RS,] & m_rho_comp$testing_output[grp2_pred_ind_RS,] <= m_rho_comp$m_pred$upper95[grp2_pred_ind_RS,])/length(m_rho_comp$m_pred$upper95[grp2_pred_ind_RS,])
l95_2_RS = mean(m_rho_comp$m_pred$upper95[grp2_pred_ind_RS,]-m_rho_comp$m_pred$lower95[grp2_pred_ind_RS,])
rmseO2_RS = sqrt(mean((energy_comp$pred_be_Omega[grp2_pred_ind_RS]-energy_comp$be_Omega_test[grp2_pred_ind_RS])^2))

grp3_pred_ind_RS = which(pred_ind_RS>2*N&pred_ind_RS<=3*N)
rmse3_RS = sqrt(mean((m_rho_comp$pred_rho_di[grp3_pred_ind_RS,]-m_rho_comp$testing_output[grp3_pred_ind_RS,])^2))
cov3_RS = sum(m_rho_comp$testing_output[grp3_pred_ind_RS,]>=m_rho_comp$m_pred$lower95[grp3_pred_ind_RS,] & m_rho_comp$testing_output[grp3_pred_ind_RS,] <= m_rho_comp$m_pred$upper95[grp3_pred_ind_RS,])/length(m_rho_comp$m_pred$upper95[grp3_pred_ind_RS,])
l95_3_RS = mean(m_rho_comp$m_pred$upper95[grp3_pred_ind_RS,]-m_rho_comp$m_pred$lower95[grp3_pred_ind_RS,])
rmseO3_RS = sqrt(mean((energy_comp$pred_be_Omega[grp3_pred_ind_RS]-energy_comp$be_Omega_test[grp3_pred_ind_RS])^2))

grp4_pred_ind_RS =  which(pred_ind_RS>3*N&pred_ind_RS<=4*N)
rmse4_RS = sqrt(mean((m_rho_comp$pred_rho_di[grp4_pred_ind_RS,]-m_rho_comp$testing_output[grp4_pred_ind_RS,])^2))
cov4_RS = sum(m_rho_comp$testing_output[grp4_pred_ind_RS,]>=m_rho_comp$m_pred$lower95[grp4_pred_ind_RS,] & m_rho_comp$testing_output[grp4_pred_ind_RS,] <= m_rho_comp$m_pred$upper95[grp4_pred_ind_RS,])/length(m_rho_comp$m_pred$upper95[grp4_pred_ind_RS,])
l95_4_RS = mean(m_rho_comp$m_pred$upper95[grp4_pred_ind_RS,]-m_rho_comp$m_pred$lower95[grp4_pred_ind_RS,])
rmseO4_RS = sqrt(mean((energy_comp$pred_be_Omega[grp4_pred_ind_RS]-energy_comp$be_Omega_test[grp4_pred_ind_RS])^2))

RS2 = matrix(c(n1_comp,n2_comp,n3_comp,n4_comp,
               rmse1_RS,rmse2_RS,rmse3_RS,rmse4_RS,
               cov1_RS,cov2_RS,cov3_RS,cov4_RS,
               l95_1_RS,l95_2_RS,l95_3_RS,l95_4_RS,
               rmseO1_RS,rmseO2_RS,rmseO3_RS,rmseO4_RS),
             nr=4,dimnames = list(c("Grp1", "Grp2","Grp3","Grp4"),
                                  c("n", "RMSE_rho", "coverage", "length","RMSE_Omega")))
RS2

#######D_opt
aug_index_D = order_ind[-train_index][-resultsall_opt$pred_index]
n_aug_D1=sum(aug_index_D<=N)
n_aug_D2=sum(aug_index_D>N & aug_index_D<=2*N)
n_aug_D3=sum(aug_index_D>2*N & aug_index_D<=3*N)
n_aug_D4=sum(aug_index_D>3*N & aug_index_D<=4*N)

pred_ind_ori_D = order_ind[-train_index][resultsall_opt$pred_index]
grp1_pred_ind_D = which(pred_ind_ori_D<=N)
rmse_D1 = sqrt(mean((resultsall_opt$pred_record[resultsall_opt$pred_index,][grp1_pred_ind_D,]-m_rho_ori$testing_output[resultsall_opt$pred_index,][grp1_pred_ind_D,])^2))
cov_D1 = sum(m_rho_ori$testing_output[resultsall_opt$pred_index,][grp1_pred_ind_D,]>=resultsall_opt$lb95_record[resultsall_opt$pred_index,][grp1_pred_ind_D,] & m_rho_ori$testing_output[resultsall_opt$pred_index,][grp1_pred_ind_D,]<= resultsall_opt$ub95_record[resultsall_opt$pred_index,][grp1_pred_ind_D,])/length(resultsall_opt$ub95_record[resultsall_opt$pred_index,][grp1_pred_ind_D,])
l95_D1 = mean(resultsall_opt$ub95_record[resultsall_opt$pred_index,][grp1_pred_ind_D,]-resultsall_opt$lb95_record[resultsall_opt$pred_index,][grp1_pred_ind_D,])
rmse_DO1=sqrt(mean((energyall_opt$pred_be_Omega[grp1_pred_ind_D]-energyall_opt$be_Omega_test[grp1_pred_ind_D])^2))

grp2_pred_ind_D = which(pred_ind_ori_D>N & pred_ind_ori_D<=2*N)
rmse_D2 = sqrt(mean((resultsall_opt$pred_record[resultsall_opt$pred_index,][grp2_pred_ind_D,]-m_rho_ori$testing_output[resultsall_opt$pred_index,][grp2_pred_ind_D,])^2))
cov_D2 = sum(m_rho_ori$testing_output[resultsall_opt$pred_index,][grp2_pred_ind_D,]>=resultsall_opt$lb95_record[resultsall_opt$pred_index,][grp2_pred_ind_D,] & m_rho_ori$testing_output[resultsall_opt$pred_index,][grp2_pred_ind_D,]<= resultsall_opt$ub95_record[resultsall_opt$pred_index,][grp2_pred_ind_D,])/length(resultsall_opt$ub95_record[resultsall_opt$pred_index,][grp2_pred_ind_D,])
l95_D2 = mean(resultsall_opt$ub95_record[resultsall_opt$pred_index,][grp2_pred_ind_D,]-resultsall_opt$lb95_record[resultsall_opt$pred_index,][grp2_pred_ind_D,])
rmse_DO2=sqrt(mean((energyall_opt$pred_be_Omega[grp2_pred_ind_D]-energyall_opt$be_Omega_test[grp2_pred_ind_D])^2))

grp3_pred_ind_D = which(pred_ind_ori_D>2*N & pred_ind_ori_D<=3*N)
rmse_D3 = sqrt(mean((resultsall_opt$pred_record[resultsall_opt$pred_index,][grp3_pred_ind_D,]-m_rho_ori$testing_output[resultsall_opt$pred_index,][grp3_pred_ind_D,])^2))
cov_D3 = sum(m_rho_ori$testing_output[resultsall_opt$pred_index,][grp3_pred_ind_D,]>=resultsall_opt$lb95_record[resultsall_opt$pred_index,][grp3_pred_ind_D,] & m_rho_ori$testing_output[resultsall_opt$pred_index,][grp3_pred_ind_D,]<= resultsall_opt$ub95_record[resultsall_opt$pred_index,][grp3_pred_ind_D,])/length(resultsall_opt$ub95_record[resultsall_opt$pred_index,][grp3_pred_ind_D,])
l95_D3 = mean(resultsall_opt$ub95_record[resultsall_opt$pred_index,][grp3_pred_ind_D,]-resultsall_opt$lb95_record[resultsall_opt$pred_index,][grp3_pred_ind_D,])
rmse_DO3=sqrt(mean((energyall_opt$pred_be_Omega[grp3_pred_ind_D]-energyall_opt$be_Omega_test[grp3_pred_ind_D])^2))

grp4_pred_ind_D = which(pred_ind_ori_D>3*N & pred_ind_ori_D<=4*N)
rmse_D4 = sqrt(mean((resultsall_opt$pred_record[resultsall_opt$pred_index,][grp4_pred_ind_D,]-m_rho_ori$testing_output[resultsall_opt$pred_index,][grp4_pred_ind_D,])^2))
cov_D4 = sum(m_rho_ori$testing_output[resultsall_opt$pred_index,][grp4_pred_ind_D,]>=resultsall_opt$lb95_record[resultsall_opt$pred_index,][grp4_pred_ind_D,] & m_rho_ori$testing_output[resultsall_opt$pred_index,][grp4_pred_ind_D,]<= resultsall_opt$ub95_record[resultsall_opt$pred_index,][grp4_pred_ind_D,])/length(resultsall_opt$ub95_record[resultsall_opt$pred_index,][grp4_pred_ind_D,])
l95_D4 = mean(resultsall_opt$ub95_record[resultsall_opt$pred_index,][grp4_pred_ind_D,]-resultsall_opt$lb95_record[resultsall_opt$pred_index,][grp4_pred_ind_D,])
rmse_DO4=sqrt(mean((energyall_opt$pred_be_Omega[grp4_pred_ind_D]-energyall_opt$be_Omega_test[grp4_pred_ind_D])^2))

subgroup_D = matrix(c(n_aug_D1,n_aug_D2,n_aug_D3,n_aug_D4,
                    rmse_D1,rmse_D2,rmse_D3,rmse_D4,
                    cov_D1,cov_D2,cov_D3,cov_D4,
                    l95_D1,l95_D2,l95_D3,l95_D4,
                    rmse_DO1,rmse_DO2,rmse_DO3,rmse_DO4),
                  nr=4,dimnames = list(c("Grp1", "Grp2","Grp3","Grp4"),
                                       c("n_aug_D", "RMSE_rho", "cov_Derage", "length","RMSE_Omega")))
subgroup_D



###################################
##### create dataset for plot #####
###################################

dat_all = data.frame(RMSE_rho = c(m_rho_ori$RMSE,rmse1_rs,rmse2_rs,rmse3_rs,rmse4_rs,
                                  resultsall$RMSE,rmse1,rmse2,rmse3,rmse4,
                                  m_rho_comp$RMSE,rmse1_RS,rmse2_RS,rmse3_RS,rmse4_RS,
                                  resultsall_opt$RMSE,rmse_D1,rmse_D2,rmse_D3,rmse_D4),
                     RMSE_Omega = c(pred_ori$be_Omega_RMSE,rmseO1_rs,rmseO2_rs,rmseO3_rs,rmseO4_rs,
                                    energyall$be_Omega_RMSE,rmseO1,rmseO2,rmseO3,rmseO4,
                                    energy_comp$be_Omega_RMSE,rmseO1_RS,rmseO2_RS,rmseO3_RS,rmseO4_RS,
                                    energyall_opt$be_Omega_RMSE,rmse_DO1,rmse_DO2,rmse_DO3,rmse_DO4),
                     Design = rep(c("RS1","ALEC","RS2","D-opt"),each=5),
                     Group = factor(rep(c("all groups","group 1","group 2","group 3","group 4"),4),levels = c("all groups","group 1","group 2","group 3","group 4")))


