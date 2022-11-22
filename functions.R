orig_rho_GP = function(rho_record,be_V_ext_record,be_mu_record,group,N,n,
                       delete_index,train_index,be_Omega_record){
  rho_train = rho_record[(1:N+(group-1)*N),][train_index,]
  be_V_train = be_V_ext_record[(1:N+(group-1)*N),][train_index,]
  be_mu_train = be_mu_record[(1:N+(group-1)*N)][train_index]
  be_Omega_train = be_Omega_record[(1:N+(group-1)*N)][train_index]
  
  rho_test = rho_record[(1:N+(group-1)*N),][-train_index,]
  be_V_test = be_V_ext_record[(1:N+(group-1)*N),][-train_index,]
  be_mu_test = be_mu_record[(1:N+(group-1)*N)][-train_index]
  be_Omega_test = be_Omega_record[(1:N+(group-1)*N)][-train_index]
  
  #build GP model
  input=matrix(be_mu_train,ncol = k-length(delete_index),nrow = n)-be_V_train[,-delete_index]
  testing_input=matrix(be_mu_test,ncol = k-length(delete_index),nrow = N-n)-be_V_test[,-delete_index]
  
  output=rho_train[,-delete_index]
  testing_output=rho_test[,-delete_index]
  
  m_GP=ppgasp(design=input,response=output,nugget.est=F,lower_bound=F,
              isotropic = T,optimization="nelder-mead",num_initial_values = 5)
  
  m_pred=predict(m_GP,testing_input)
  pred_rho_di = m_pred$mean #predicted density with deleted index
  
  #check prediction accuracy of rho
  RMSE = sqrt(mean((pred_rho_di-testing_output)^2))
  NRMSE = sqrt(mean((pred_rho_di-testing_output)^2))/sd(testing_output)
  #coverage
  coverage = sum(testing_output>=m_pred$lower95 & testing_output<= m_pred$upper95)/length(m_pred$lower95)
  #length
  length95 = mean(m_pred$upper95-m_pred$lower95)
  #RMSE_each
  RMSE_each=sqrt(rowMeans((pred_rho_di-testing_output)^2))
  return(list(m_GP=m_GP,pred_rho_di=pred_rho_di,m_pred=m_pred,RMSE=RMSE,NRMSE=NRMSE,coverage=coverage,length95=length95,
              rho_train=rho_train, be_V_train=be_V_train, be_mu_train=be_mu_train,
              be_Omega_train=be_Omega_train,
              rho_test=rho_test,be_V_test=be_V_test,be_mu_test=be_mu_test,
              input=input, testing_input=testing_input,output=output,testing_output=testing_output,
              be_Omega_test=be_Omega_test,
              RMSE_each=RMSE_each))
}

AL_GP_mean = function(testing_input, testing_output, n_test, GP_before, n_before, sd_threshold,alphaLevel=0.05){
  pred_record=matrix(NA, n_test,dim(testing_output)[2] )
  sd_record = matrix(NA, n_test,dim(testing_output)[2] )
  lb95_record=matrix(NA, n_test,dim(testing_output)[2] )
  ub95_record=matrix(NA, n_test,dim(testing_output)[2] )
  
  n_aug=0
  m_GP_aug=GP_before
  m_GP_aug@R0[[1]]=NULL
  pred_index=NULL
  for(i_test in 1: n_test){
    print(i_test)
    testing_here=matrix(testing_input[i_test,],nr=1)
    
    m_pred_individual=predict(m_GP_aug,testing_here)
    
    #if(mean(( abs(m_pred_individual$sd)))<sd_threshold){
    if(sqrt(mean((abs(m_pred_individual$sd)^2)))<sd_threshold/qt(1-alphaLevel/2,df=n_test+n_aug)){
      
      pred_index=c(pred_index,i_test)
      
      if(min(m_pred_individual$mean)<0){
        mean_ind = which(m_pred_individual$mean[1,]<0)
        m_pred_individual$mean[1,mean_ind] = 0
      }
      pred_record[i_test,]=m_pred_individual$mean
      sd_record[i_test,] = m_pred_individual$sd
      if(min(m_pred_individual$lower95)<0){
        lb_ind = which(m_pred_individual$lower95[1,]<0)
        m_pred_individual$lower95[1,lb_ind] = 0
      }
      lb95_record[i_test,]=m_pred_individual$lower95
      if(min(m_pred_individual$upper95)<0){
        ub_ind = which(m_pred_individual$upper95[1,]<0)
        m_pred_individual$upper95[1,ub_ind] = 0
      }
      ub95_record[i_test,]=m_pred_individual$upper95
      
      #sqrt(mean( (pred_record[i_test,]-testing_output[i_test,])^2))/sd(testing_output[i_test,])
    }else{
      n_aug=n_aug+1
      output_add=matrix(testing_output[i_test,],nr=1)
      if(n_aug%%50==0 & m_GP_aug@num_obs<=350){
        aug_input = rbind(m_GP_aug@input,testing_here)
        aug_output = rbind(m_GP_aug@output,output_add)
        m_GP_aug=ppgasp(design=aug_input,response=aug_output,nugget.est=F,lower_bound=F,
                        isotropic = T,optimization="nelder-mead",num_initial_values = 5)
      }else{
        #[[1]]num_obs_new; [[2]]input_new; [[3]]output_new; [[4]]X_new; [[5]]L_new
        #[[6]]LX_new; [[7]]theta_hat_new; [[8]]S_2_all_new/(num_obs_new-q)
        update_res = update_model_cpp(m_GP_aug@num_obs,m_GP_aug@input,m_GP_aug@output,m_GP_aug@X, m_GP_aug@q,
                                  m_GP_aug@theta_hat, m_GP_aug@L,m_GP_aug@alpha,beta=m_GP_aug@beta_hat,
                                  input_add=testing_here, output_add=output_add)
        
        
        m_GP_aug@X=update_res[[4]]  ##only hold for constant mean
        m_GP_aug@output= update_res[[3]]
        #m_GP_aug@R0=R0_aug
        
        m_GP_aug@L=update_res[[5]]
        m_GP_aug@LX=update_res[[6]]
        m_GP_aug@theta_hat=update_res[[7]]
        m_GP_aug@sigma2_hat=update_res[[8]] ##post mode only
        m_GP_aug@input=update_res[[2]]
        m_GP_aug@num_obs=update_res[[1]]
        #break
      }
    }
  }
  #RMSE
  RMSE = sqrt(mean((pred_record[pred_index,]-testing_output[pred_index,])^2))
  #NRMSE
  NRMSE = RMSE/sd(testing_output[pred_index,])
  #coverage
  coverage = sum(testing_output[pred_index,]>=lb95_record[pred_index,] & testing_output[pred_index,]<= ub95_record[pred_index,])/length(ub95_record[pred_index,])
  #length
  length95 = mean(ub95_record[pred_index,]-lb95_record[pred_index,])
  #RMSE_each
  RMSE_each=sqrt(rowMeans((pred_record[pred_index,]-testing_output[pred_index,])^2))
  #delta_each
  delta_each = sqrt(rowMeans(sd_record[pred_index,]^2))
  #difference between RMSE and delta
  diff_rmse=RMSE_each-delta_each
  
  return(list(m_GP_aug=m_GP_aug, n_aug=n_aug, pred_index=pred_index,
              pred_record=pred_record, sd_record=sd_record, lb95_record=lb95_record, ub95_record=ub95_record,
              RMSE=RMSE,NRMSE=NRMSE,coverage=coverage,length95=length95,
              RMSE_each=RMSE_each,delta_each=delta_each,diff_rmse=diff_rmse))
}

AL_GP_max = function(testing_input, testing_output, n_test, GP_before, n_before, sd_threshold,alphaLevel=0.05){
  pred_record=matrix(NA, n_test,dim(testing_input)[2] )
  sd_record = matrix(NA, n_test,dim(testing_input)[2] )
  lb95_record=matrix(NA, n_test,dim(testing_input)[2] )
  ub95_record=matrix(NA, n_test,dim(testing_input)[2] )
  
  n_aug=0
  m_GP_aug=GP_before
  m_GP_aug@R0[[1]]=NULL
  pred_index=NULL
  for(i_test in 1: n_test){
    print(i_test)
    testing_here=matrix(testing_input[i_test,],nr=1)
    
    m_pred_individual=predict(m_GP_aug,testing_here)
    
    #if(mean(( abs(m_pred_individual$sd)))<sd_threshold){
    if(sqrt(max((abs(m_pred_individual$sd)^2)))<sd_threshold/qt(1-alphaLevel/2,df=n_test+n_aug)){
      
      pred_index=c(pred_index,i_test)
      
      if(min(m_pred_individual$mean)<0){
        mean_ind = which(m_pred_individual$mean[1,]<0)
        m_pred_individual$mean[1,mean_ind] = 0
      }
      pred_record[i_test,]=m_pred_individual$mean
      sd_record[i_test,] = m_pred_individual$sd
      if(min(m_pred_individual$lower95)<0){
        lb_ind = which(m_pred_individual$lower95[1,]<0)
        m_pred_individual$lower95[1,lb_ind] = 0
      }
      lb95_record[i_test,]=m_pred_individual$lower95
      if(min(m_pred_individual$upper95)<0){
        ub_ind = which(m_pred_individual$upper95[1,]<0)
        m_pred_individual$upper95[1,ub_ind] = 0
      }
      ub95_record[i_test,]=m_pred_individual$upper95
      
      #sqrt(mean( (pred_record[i_test,]-testing_output[i_test,])^2))/sd(testing_output[i_test,])
    }else{
      n_aug=n_aug+1
      output_add=matrix(testing_output[i_test,],nr=1)
      if(n_aug%%50==0 & m_GP_aug@num_obs<=350){
        aug_input = rbind(m_GP_aug@input,testing_here)
        aug_output = rbind(m_GP_aug@output,output_add)
        m_GP_aug=ppgasp(design=aug_input,response=aug_output,nugget.est=F,lower_bound=F,
                        isotropic = T,optimization="nelder-mead",num_initial_values = 5)
      }else{
        #[[1]]num_obs_new; [[2]]input_new; [[3]]output_new; [[4]]X_new; [[5]]L_new
        #[[6]]LX_new; [[7]]theta_hat_new; [[8]]S_2_all_new/(num_obs_new-q)
        update_res = update_model_cpp(m_GP_aug@num_obs,m_GP_aug@input,m_GP_aug@output,m_GP_aug@X, m_GP_aug@q,
                                      m_GP_aug@theta_hat, m_GP_aug@L,m_GP_aug@alpha,beta=m_GP_aug@beta_hat,
                                      input_add=testing_here, output_add=output_add)

        m_GP_aug@X=update_res[[4]]  ##only hold for constant mean
        m_GP_aug@output= update_res[[3]]
        #m_GP_aug@R0=R0_aug
        
        m_GP_aug@L=update_res[[5]]
        m_GP_aug@LX=update_res[[6]]
        m_GP_aug@theta_hat=update_res[[7]]
        m_GP_aug@sigma2_hat=update_res[[8]] ##post mode only
        m_GP_aug@input=update_res[[2]]
        m_GP_aug@num_obs=update_res[[1]]
        #break
        
        

      }
    }
  }
  #RMSE
  RMSE = sqrt(mean((pred_record[pred_index,]-testing_output[pred_index,])^2))
  #NRMSE
  NRMSE = RMSE/sd(testing_output[pred_index,])
  #coverage
  coverage = sum(testing_output[pred_index,]>=lb95_record[pred_index,] & testing_output[pred_index,]<= ub95_record[pred_index,])/length(ub95_record[pred_index,])
  #length
  length95 = mean(ub95_record[pred_index,]-lb95_record[pred_index,])
  #RMSE_each
  RMSE_each=sqrt(rowMeans((pred_record[pred_index,]-testing_output[pred_index,])^2))
  #delta_each
  delta_each = sqrt(rowMeans(sd_record[pred_index,]^2))
  #difference between RMSE and delta
  diff_rmse=RMSE_each-delta_each
  
  return(list(m_GP_aug=m_GP_aug, n_aug=n_aug, pred_index=pred_index,
              pred_record=pred_record, sd_record=sd_record, lb95_record=lb95_record, ub95_record=ub95_record,
              RMSE=RMSE,NRMSE=NRMSE,coverage=coverage,length95=length95,
              RMSE_each=RMSE_each,delta_each=delta_each,diff_rmse=diff_rmse))
}


AL_input_decomp = function(testing_input, testing_output, n_test, GP_before, n_before, sd_threshold,alphaLevel=0.05){
  pred_record=matrix(NA, n_test,dim(testing_input)[2] )
  sd_record = matrix(NA, n_test,dim(testing_input)[2] )
  lb95_record=matrix(NA, n_test,dim(testing_input)[2] )
  ub95_record=matrix(NA, n_test,dim(testing_input)[2] )
  
  n_aug=0
  n_cluster = 1
  model_list = list(1)
  model_list[[1]]=GP_before
  center_matrix = matrix(apply(model_list[[1]]@input,2,mean),nr=1)
  obs_vec = c(model_list[[1]]@num_obs)
  #m_GP_aug@R0[[1]]=NULL
  pred_index=NULL
  for(i_test in 1: n_test){
    print(i_test)
    testing_here=matrix(testing_input[i_test,],nr=1)
    output_add=matrix(testing_output[i_test,],nr=1)
    
    clus_index = which.min(euclidean_distance(center_matrix, testing_here))
    
    m_pred_individual=predict(model_list[[clus_index]],testing_here)
    
    #if(mean(( abs(m_pred_individual$sd)))<sd_threshold){
    if(sqrt(mean((abs(m_pred_individual$sd)^2)))<sd_threshold/qt(1-alphaLevel/2,df=n_test+n_aug)){
      
      pred_index=c(pred_index,i_test)
      
      if(min(m_pred_individual$mean)<0){
        mean_ind = which(m_pred_individual$mean[1,]<0)
        m_pred_individual$mean[1,mean_ind] = 0
      }
      pred_record[i_test,]=m_pred_individual$mean
      sd_record[i_test,] = m_pred_individual$sd
      if(min(m_pred_individual$lower95)<0){
        lb_ind = which(m_pred_individual$lower95[1,]<0)
        m_pred_individual$lower95[1,lb_ind] = 0
      }
      lb95_record[i_test,]=m_pred_individual$lower95
      if(min(m_pred_individual$upper95)<0){
        ub_ind = which(m_pred_individual$upper95[1,]<0)
        m_pred_individual$upper95[1,ub_ind] = 0
      }
      ub95_record[i_test,]=m_pred_individual$upper95
      
      #sqrt(mean( (pred_record[i_test,]-testing_output[i_test,])^2))/sd(testing_output[i_test,])
    }else{
      n_aug=n_aug+1
      if(model_list[[clus_index]]@num_obs>=400){
        n_cluster = n_cluster+1
        aug_input = rbind(model_list[[clus_index]]@input,testing_here)
        aug_output = rbind(model_list[[clus_index]]@output,output_add)
        set.seed(1)
        decomp = kmeans(aug_input,2)
        group1 = which(decomp$cluster==1)
        m_GP1=ppgasp(design=aug_input[group1,],response=aug_output[group1,],nugget.est=F,lower_bound=F,
                        isotropic = T,optimization="nelder-mead",num_initial_values = 5)
        model_list[[clus_index]]=m_GP1
        center_matrix[clus_index,]=decomp$centers[1,]
        obs_vec[clus_index] = length(group1)
        
        m_GP2=ppgasp(design=aug_input[-group1,],response=aug_output[-group1,],nugget.est=F,lower_bound=F,
                     isotropic = T,optimization="nelder-mead",num_initial_values = 5)
        model_list[[n_cluster]]=m_GP2
        center_matrix=rbind(center_matrix,decomp$centers[2,])
        obs_vec[n_cluster] = dim(aug_input)[1]-length(group1)
      }else{
        center_matrix[clus_index,]=center_matrix[clus_index,]*obs_vec[clus_index]/(obs_vec[clus_index]+1)+testing_here/(obs_vec[clus_index]+1)
        obs_vec[clus_index] = obs_vec[clus_index]+1
        if(model_list[[clus_index]]@num_obs%%50==0&model_list[[clus_index]]@num_obs<=350){
          aug_input = rbind(model_list[[clus_index]]@input,testing_here)
          aug_output = rbind(model_list[[clus_index]]@output,output_add)
          m_GP=ppgasp(design=aug_input,response=aug_output,nugget.est=F,lower_bound=F,
                          isotropic = T,optimization="nelder-mead",num_initial_values = 5)
          model_list[[clus_index]]=m_GP
        }else{
          m_GP_aug=model_list[[clus_index]]
          #[[1]]num_obs_new; [[2]]input_new; [[3]]output_new; [[4]]X_new; [[5]]L_new
          #[[6]]LX_new; [[7]]theta_hat_new; [[8]]S_2_all_new/(num_obs_new-q)
          update_res = update_model_cpp(m_GP_aug@num_obs,m_GP_aug@input,m_GP_aug@output,m_GP_aug@X, m_GP_aug@q,
                                        m_GP_aug@theta_hat, m_GP_aug@L,m_GP_aug@alpha,beta=m_GP_aug@beta_hat,
                                        input_add=testing_here, output_add=output_add)
          
          
          m_GP_aug@X=update_res[[4]]  ##only hold for constant mean
          m_GP_aug@output= update_res[[3]]
          #m_GP_aug@R0=R0_aug
          
          m_GP_aug@L=update_res[[5]]
          m_GP_aug@LX=update_res[[6]]
          m_GP_aug@theta_hat=update_res[[7]]
          m_GP_aug@sigma2_hat=update_res[[8]] ##post mode only
          m_GP_aug@input=update_res[[2]]
          m_GP_aug@num_obs=update_res[[1]]
          
          model_list[[clus_index]]=m_GP_aug
          #break
        }
      }
   
    }
  }
  #RMSE
  RMSE = sqrt(mean((pred_record[pred_index,]-testing_output[pred_index,])^2))
  #NRMSE
  NRMSE = RMSE/sd(testing_output[pred_index,])
  #coverage
  coverage = sum(testing_output[pred_index,]>=lb95_record[pred_index,] & testing_output[pred_index,]<= ub95_record[pred_index,])/length(ub95_record[pred_index,])
  #length
  length95 = mean(ub95_record[pred_index,]-lb95_record[pred_index,])
  #RMSE_each
  RMSE_each=sqrt(rowMeans((pred_record[pred_index,]-testing_output[pred_index,])^2))
  #delta_each
  delta_each = sqrt(rowMeans(sd_record[pred_index,]^2))
  #difference between RMSE and delta
  diff_rmse=RMSE_each-delta_each
  
  return(list(n_cluster=n_cluster, model_list=model_list, obs_vec=obs_vec,
              n_aug=n_aug, pred_index=pred_index,
              pred_record=pred_record, sd_record=sd_record, lb95_record=lb95_record, ub95_record=ub95_record,
              RMSE=RMSE,NRMSE=NRMSE,coverage=coverage,length95=length95,
              RMSE_each=RMSE_each,delta_each=delta_each,diff_rmse=diff_rmse))
}


energy_pred_NRMSE = function(n_test, pred_rho_di, a, L, k, be_Omega_test, delete_index, plot=T){
  pred_rho = matrix(0,ncol=k,nrow=n_test)
  pred_rho[,-delete_index]=pred_rho_di
  pred_be_Omega = rep(NA,n_test)
  for(i in 1:(n_test)){
    pred_be_Omega[i]=beta_Omega_cpp(pred_rho[i,],a,L)
  }
  if(plot == T){
    plot(be_Omega_test,pred_be_Omega)
    abline(0,1,col=2,lwd=1.5)
  }
  #check whether there's NA prediction in beta*Omega
  if(sum(is.na(pred_be_Omega))>0){
    print(paste(sum(is.na(pred_be_Omega)),"NA in prediction of beta*Omega produced"))
    Omega_na_index = which(is.na(pred_be_Omega))
    be_Omega_RMSE=sqrt(mean((pred_be_Omega[-Omega_na_index]-be_Omega_test[-Omega_na_index])^2))
    be_Omega_NRMSE=sqrt(mean((pred_be_Omega[-Omega_na_index]-be_Omega_test[-Omega_na_index])^2))/sd(be_Omega_test[-Omega_na_index])
  }else{
    be_Omega_RMSE=sqrt(mean((pred_be_Omega-be_Omega_test)^2))
    be_Omega_NRMSE=sqrt(mean((pred_be_Omega-be_Omega_test)^2))/sd(be_Omega_test)
  }

  return(list(pred_be_Omega=pred_be_Omega,be_Omega_test=be_Omega_test,
              be_Omega_RMSE=be_Omega_RMSE,be_Omega_NRMSE=be_Omega_NRMSE,pred_rho=pred_rho))
}

pred_GP = function(model, testing_input, testing_output){
  m_pred=predict(model,testing_input)
  m_pred$mean[which(m_pred$mean<0)]=0
  m_pred$upper95[which(m_pred$upper95<0)]=0
  m_pred$lower95[which(m_pred$lower95<0)]=0
  
  pred_rho_di = m_pred$mean #predicted density with deleted index
  
  #check prediction accuracy of rho
  RMSE = sqrt(mean((pred_rho_di-testing_output)^2))
  NRMSE = RMSE/sd(testing_output)
  #coverage
  coverage = sum(testing_output>=m_pred$lower95 & testing_output<= m_pred$upper95)/length(m_pred$lower95)
  #length
  length95 = mean(m_pred$upper95-m_pred$lower95)
  return(list(m_pred=m_pred,pred_rho_di=pred_rho_di,RMSE=RMSE,NRMSE=NRMSE,coverage=coverage,length95=length95))
}

build_GP = function(rho_record,be_V_ext_record,be_mu_record,N,n,
                    delete_index,train_index,be_Omega_record){
  rho_train = rho_record[train_index,]
  be_V_train = be_V_ext_record[train_index,]
  be_mu_train = be_mu_record[train_index]
  be_Omega_train = be_Omega_record[train_index]
  
  rho_test = rho_record[-train_index,]
  be_V_test = be_V_ext_record[-train_index,]
  be_mu_test = be_mu_record[-train_index]
  be_Omega_test = be_Omega_record[-train_index]
  
  #build GP model
  input=matrix(be_mu_train,ncol = k-length(delete_index),nrow = n)-be_V_train[,-delete_index]
  testing_input=matrix(be_mu_test,ncol = k-length(delete_index),nrow = N-n)-be_V_test[,-delete_index]
  
  output=rho_train[,-delete_index]
  testing_output=rho_test[,-delete_index]
  
  m_GP=ppgasp(design=input,response=output,nugget.est=F,lower_bound=F,
              isotropic = T,optimization="nelder-mead",num_initial_values = 5)
  
  m_pred=predict(m_GP,testing_input)
  pred_rho_di = m_pred$mean #predicted density with deleted index
  
  #check prediction accuracy of rho
  RMSE = sqrt(mean((pred_rho_di-testing_output)^2))
  NRMSE = sqrt(mean((pred_rho_di-testing_output)^2))/sd(testing_output)
  #coverage
  coverage = sum(testing_output>=m_pred$lower95 & testing_output<= m_pred$upper95)/length(m_pred$lower95)
  #length
  length95 = mean(m_pred$upper95-m_pred$lower95)
  #RMSE_each
  RMSE_each=sqrt(rowMeans((pred_rho_di-testing_output)^2))
  return(list(m_GP=m_GP,pred_rho_di=pred_rho_di,m_pred=m_pred,RMSE=RMSE,NRMSE=NRMSE,coverage=coverage,length95=length95,RMSE_each=RMSE_each,
              rho_train=rho_train, be_V_train=be_V_train, be_mu_train=be_mu_train,
              be_Omega_train=be_Omega_train,
              rho_test=rho_test,be_V_test=be_V_test,be_mu_test=be_mu_test,
              input=input, testing_input=testing_input,output=output,testing_output=testing_output,
              be_Omega_test=be_Omega_test))
}


separate_group = function(rho_record,be_V_ext_record,be_mu_record,
                          be_Omega_record,
                          delete_index,N){
  input_record = matrix(be_mu_record,ncol = k-length(delete_index),nrow = N*4)-be_V_ext_record[,-delete_index]
  output_record = rho_record[,-delete_index]
  
  input1 = input_record[(1:N),]
  output1 = output_record[(1:N),]
  input2 = input_record[(1:N+(2-1)*N),]
  output2 = output_record[(1:N+(2-1)*N),]
  input3 = input_record[(1:N+(3-1)*N),]
  output3 = output_record[(1:N+(3-1)*N),]
  input4 = input_record[(1:N+(4-1)*N),]
  output4 = output_record[(1:N+(4-1)*N),]

  be_Omega1 = be_Omega_record[(1:N)]
  be_Omega2 = be_Omega_record[(1:N+(2-1)*N)]
  be_Omega3 = be_Omega_record[(1:N+(3-1)*N)]
  be_Omega4 = be_Omega_record[(1:N+(4-1)*N)]
  
  return(list(input_record=input_record,
              output_record=output_record,
       
              input1=input1,output1=output1,
              input2=input2,output2=output2,
              input3=input3,output3=output3,
              input4=input4,output4=output4,

              be_Omega1=be_Omega1,
              be_Omega2=be_Omega2,
              be_Omega3=be_Omega3,
              be_Omega4=be_Omega4))
}

augmented_GP_Doptimal_cpp = function(testing_input, testing_output, n_test, GP_before, n_before, D_threshold){
  pred_record=matrix(NA, n_test,dim(testing_input)[2] )
  sd_record = matrix(NA, n_test,dim(testing_input)[2] )
  lb95_record=matrix(NA, n_test,dim(testing_input)[2] )
  ub95_record=matrix(NA, n_test,dim(testing_input)[2] )
  
  n_aug=0
  m_GP_aug=GP_before
  m_GP_aug@R0[[1]]=NULL
  pred_index=NULL
  for(i_test in 1:n_test){
    print(i_test)
    testing_here=matrix(testing_input[i_test,],nr=1)
    #output_add=matrix(testing_output[i_test,],nr=1)
    
    m_pred_individual=predict(m_GP_aug,testing_here)
    
    #if(mean(( abs(m_pred_individual$sd)))<sd_threshold){
    D_max = D_optimality(input_test=testing_here,L=m_GP_aug@L,input=m_GP_aug@input,beta=m_GP_aug@beta_hat)
    #print(D_max)
    if(abs(D_max)<D_threshold){
      
      pred_index=c(pred_index,i_test)
      
      if(min(m_pred_individual$mean)<0){
        mean_ind = which(m_pred_individual$mean[1,]<0)
        m_pred_individual$mean[1,mean_ind] = 0
      }
      pred_record[i_test,]=m_pred_individual$mean
      sd_record[i_test,] = m_pred_individual$sd
      if(min(m_pred_individual$lower95)<0){
        lb_ind = which(m_pred_individual$lower95[1,]<0)
        m_pred_individual$lower95[1,lb_ind] = 0
      }
      lb95_record[i_test,]=m_pred_individual$lower95
      if(min(m_pred_individual$upper95)<0){
        ub_ind = which(m_pred_individual$upper95[1,]<0)
        m_pred_individual$upper95[1,ub_ind] = 0
      }
      ub95_record[i_test,]=m_pred_individual$upper95
      
      #sqrt(mean( (pred_record[i_test,]-testing_output[i_test,])^2))/sd(testing_output[i_test,])
    }else{
      #print('add')
      n_aug=n_aug+1
      output_add=matrix(testing_output[i_test,],nr=1)
      if(n_aug%%50==0 & m_GP_aug@num_obs<=350){
        aug_input = rbind(m_GP_aug@input,testing_here)
        aug_output = rbind(m_GP_aug@output,output_add)
        m_GP_aug=ppgasp(design=aug_input,response=aug_output,nugget.est=F,lower_bound=F,
                        isotropic = T,optimization="nelder-mead",num_initial_values = 5)
      }else{
        #[[1]]num_obs_new; [[2]]input_new; [[3]]output_new; [[4]]X_new; [[5]]L_new
        #[[6]]LX_new; [[7]]theta_hat_new; [[8]]S_2_all_new/(num_obs_new-q)
        update_res = update_model_cpp(num_obs=m_GP_aug@num_obs,input=m_GP_aug@input,output=m_GP_aug@output,X=m_GP_aug@X, m_GP_aug@q,
                                      m_GP_aug@theta_hat, L=m_GP_aug@L,alpha=m_GP_aug@alpha,beta=m_GP_aug@beta_hat,
                                      input_add=testing_here, output_add=output_add)
        
        if(sum(is.na(update_res[[5]]))==0){
          m_GP_aug@X=update_res[[4]]  ##only hold for constant mean
          m_GP_aug@output= update_res[[3]]
          #m_GP_aug@R0=R0_aug
          
          m_GP_aug@L=update_res[[5]]
          m_GP_aug@LX=update_res[[6]]
          m_GP_aug@theta_hat=update_res[[7]]
          m_GP_aug@sigma2_hat=update_res[[8]] ##post mode only
          m_GP_aug@input=update_res[[2]]
          m_GP_aug@num_obs=update_res[[1]]
          #break
        }else{
          update_res = update_model_w_nug_cpp(m_GP_aug@num_obs,m_GP_aug@input,m_GP_aug@output,m_GP_aug@X, m_GP_aug@q,
                                              m_GP_aug@theta_hat, m_GP_aug@L,m_GP_aug@alpha,beta=m_GP_aug@beta_hat,
                                              input_add=testing_here, output_add=output_add)
          m_GP_aug@X=update_res[[4]]  ##only hold for constant mean
          m_GP_aug@output= update_res[[3]]
          #m_GP_aug@R0=R0_aug
          
          m_GP_aug@L=update_res[[5]]
          m_GP_aug@LX=update_res[[6]]
          m_GP_aug@theta_hat=update_res[[7]]
          m_GP_aug@sigma2_hat=update_res[[8]] ##post mode only
          m_GP_aug@input=update_res[[2]]
          m_GP_aug@num_obs=update_res[[1]]
        }
      }
    }
  }
  #RMSE
  RMSE = sqrt(mean((pred_record[pred_index,]-testing_output[pred_index,])^2))
  #NRMSE
  NRMSE = RMSE/sd(testing_output[pred_index,])
  #coverage
  coverage = sum(testing_output[pred_index,]>=lb95_record[pred_index,] & testing_output[pred_index,]<= ub95_record[pred_index,])/length(ub95_record[pred_index,])
  #length
  length95 = mean(ub95_record[pred_index,]-lb95_record[pred_index,])
  #RMSE_each
  RMSE_each=sqrt(rowMeans((pred_record[pred_index,]-testing_output[pred_index,])^2))
  #delta_each
  delta_each = sqrt(rowMeans(sd_record[pred_index,]^2))
  #difference between RMSE and delta
  diff_rmse=RMSE_each-delta_each
  
  return(list(m_GP_aug=m_GP_aug, n_aug=n_aug, pred_index=pred_index,
              pred_record=pred_record, sd_record=sd_record, lb95_record=lb95_record, ub95_record=ub95_record,
              RMSE=RMSE,NRMSE=NRMSE,coverage=coverage,length95=length95,
              RMSE_each=RMSE_each,delta_each=delta_each,diff_rmse=diff_rmse))
}

