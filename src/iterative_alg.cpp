
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*- 

// we only include RcppEigen.h which pulls Rcpp.h in for us 
#include <RcppEigen.h> 
#include <Rcpp.h> 
#include <cmath> 
#include "ctools.h"

// [[Rcpp::depends(RcppEigen)]] 

using namespace Rcpp;
using namespace std;
using namespace Eigen; 


//#define PI 3.14159265359

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// via the exports attribute we tell Rcpp to make this function
// available from R
//

///////kernel functions
// [[Rcpp::export]]
Eigen::MatrixXd matern_5_2_funct (const Eigen::MatrixXd d, double beta_i){
  //inline static Mat matern_5_2_funct (const Eigen::Map<Eigen::MatrixXd> & d, double beta_i){
  const double cnst = sqrt(5.0);
  Eigen::MatrixXd matOnes = Eigen::MatrixXd::Ones(d.rows(),d.cols());
  Eigen::MatrixXd result = cnst*beta_i*d;
  return ((matOnes + result +
          result.array().pow(2.0).matrix()/3.0).cwiseProduct((-result).array().exp().matrix()));
  
}



// [[Rcpp::export]]
Eigen::MatrixXd euclidean_dist(const Eigen::MatrixXd input1,const Eigen::MatrixXd input2){
  //input are n by p, where p is larger than n
  
  int num_obs1 = input1.rows();
  int num_obs2 = input2.rows();
  
  Eigen::MatrixXd R0=R0.Ones(num_obs1,num_obs2);
  
  for (int i = 0; i < num_obs1; i++){
    
    for (int j = 0; j < num_obs2; j++){
      R0(i,j)=sqrt((input1.row(i)-input2.row(j)).array().pow(2.0).sum());
    }
  }
  return R0;
}

///////update L (cholesky decomposition for R)
// [[Rcpp::export]]
MatrixXd update_L(const int num_obs, const Eigen::MatrixXd L, const Eigen::MatrixXd input_add, const  Eigen::MatrixXd input,
                  double beta, double alpha){
  
  Eigen::MatrixXd L_new = Eigen::MatrixXd::Zero(num_obs+1,num_obs+1); //new L with 1 new input
  L_new.topLeftCorner(num_obs,num_obs) = L; // first num_obs rows are the same as L
  
  Eigen::MatrixXd r0 = Eigen::MatrixXd::Zero(num_obs+1,1);
  r0.topRows(num_obs) = euclidean_dist(input,input_add);
  r0(num_obs,0)=0;
  
  //Eigen::MatrixXd r_new = Eigen::MatrixXd::Zero(num_obs+1,1); //extended row of R
  //r_new(num_obs,0)=1.0; //last entry is 1
  Eigen::MatrixXd r_new = matern_5_2_funct(r0,beta);//isotropic_kernel(r0, beta,kernel_type, alpha);
  
  int k;
  double sum;
  for(int j = 0; j <= num_obs; j++){
    sum = 0.0;
    if(j == num_obs){
      for(k = 0;k<j;k++){
        sum += pow(L_new(j,k), 2);
      }
      // if(sum > 1.0){
      //   sum -= pow(10,-8);
      //   //L_new(num_obs,j)=sum;
      // }//else{
      L_new(num_obs,j) = sqrt(r_new(j,0) - sum);//}
    }else{
      for(k = 0;k<j;k++){
        sum += (L_new(num_obs,k) * L_new(j,k));
      }
      L_new(num_obs,j) = (r_new(j,0) - sum) /L_new(j,j);
    }
  }
  return L_new;
}

// [[Rcpp::export]]
MatrixXd Chol(const MatrixXd &R){
  
  LLT<MatrixXd> lltOfR(R);             // compute the cholesky decomposition of R called lltofR
  MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition
  return L;
}


////Jun28
// [[Rcpp::export]]
List update_model_cpp(int num_obs, MatrixXd input, MatrixXd output, MatrixXd X, int q,
                  VectorXd theta_hat, MatrixXd L, double alpha, double beta,
                  MatrixXd input_add, MatrixXd output_add){
  int p = input.cols();
  int k = output.cols();
  List return_list(8);
  int num_obs_new = num_obs+1;
  return_list[0] = num_obs_new;
  
  MatrixXd input_new = MatrixXd::Zero(num_obs_new,p);
  input_new.topRows(num_obs) = input;
  input_new.bottomRows(1) = input_add;
  return_list[1] = input_new;
  
  MatrixXd output_new = MatrixXd::Zero(num_obs_new,k);
  output_new.topRows(num_obs) = output;
  output_new.bottomRows(1) = output_add;
  return_list[2] = output_new;
  
  MatrixXd X_new = MatrixXd::Zero(num_obs_new,1);
  X_new.topRows(num_obs)=X;
  X_new(num_obs,0)=1.0;
  
  return_list[3]=X_new;
  
  MatrixXd L_new = update_L(num_obs, L, input_add, input, beta, alpha);
  return_list[4]=L_new;
  
  MatrixXd L_inv_X_new;
  MatrixXd R_inv_X_new;
  L_inv_X_new=L_new.triangularView<Lower>().solve(X_new);
  R_inv_X_new=L_new.transpose().triangularView<Upper>().solve(L_inv_X_new); //one forward and one backward to compute R.inv%*%X
  
  //return_list[5]=L_inv_X_new;
  MatrixXd Xt_R_inv_X_new=X_new.transpose()*R_inv_X_new;
  MatrixXd LX_new = sqrt(Xt_R_inv_X_new.array());
  return_list[5]=LX_new;
  
  MatrixXd L_inv_y_new;
  MatrixXd R_inv_y_new;
  L_inv_y_new=L_new.triangularView<Lower>().solve(output_new);
  R_inv_y_new = L_new.transpose().triangularView<Upper>().solve(L_inv_y_new); 
  
  //return_list[7]=L_inv_y_new;
  //MatrixXd yt_R_inv_new= (L_new.transpose().triangularView<Upper>().solve(L_new.triangularView<Lower>().solve(output_new))).transpose(); 
  MatrixXd Xt_R_inv_y_new= X_new.transpose()*R_inv_y_new;
  MatrixXd theta_hat_new=LX_new.transpose().triangularView<Upper>().solve(LX_new.triangularView<Lower>().solve(Xt_R_inv_y_new)); 
  return_list[6]=theta_hat_new;
  
  VectorXd S_2_all_new=VectorXd::Zero(k);
  
  for(int loc_i=0;loc_i<k;loc_i++){
    //S_2_all_new[loc_i]=(yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0);
    S_2_all_new[loc_i]=(output_new.col(loc_i).transpose()*R_inv_y_new.col(loc_i))(0,0)-(Xt_R_inv_y_new.col(loc_i).transpose()*Xt_R_inv_y_new.col(loc_i))(0,0)/Xt_R_inv_X_new(0,0);   
  }
  return_list[7]=S_2_all_new/(num_obs_new-q);
  
  return return_list;
}



// [[Rcpp::export]]
double D_optimality(MatrixXd input_test, MatrixXd L, MatrixXd input, double beta){
  Eigen::MatrixXd r0 = euclidean_dist(input,input_test);
  Eigen::MatrixXd r_star = matern_5_2_funct(r0,beta);
  Eigen::MatrixXd R_inv_r_star = L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(r_star)); 
  return R_inv_r_star.maxCoeff();
}

//For D-opt
// [[Rcpp::export]]
List update_model_w_nug_cpp(int num_obs, MatrixXd input, MatrixXd output, MatrixXd X, int q,
                            VectorXd theta_hat, MatrixXd L, double alpha, double beta,
                            MatrixXd input_add, MatrixXd output_add){
  int p = input.cols();
  int k = output.cols();
  List return_list(8);
  int num_obs_new = num_obs+1;
  return_list[0] = num_obs_new;
  
  MatrixXd input_new = MatrixXd::Zero(num_obs_new,p);
  input_new.topRows(num_obs) = input;
  input_new.bottomRows(1) = input_add;
  return_list[1] = input_new;
  
  MatrixXd output_new = MatrixXd::Zero(num_obs_new,p);
  output_new.topRows(num_obs) = output;
  output_new.bottomRows(1) = output_add;
  return_list[2] = output_new;
  
  MatrixXd X_new = MatrixXd::Zero(num_obs_new,1);
  X_new.topRows(num_obs)=X;
  X_new(num_obs,0)=1.0;
  
  return_list[3]=X_new;
  
  MatrixXd R0_new = euclidean_dist(input_new,input_new);
  MatrixXd R_new = matern_5_2_funct(R0_new,beta);
  MatrixXd R_tilde_new=R_new+MatrixXd::Identity(num_obs_new,num_obs_new)*pow(10,-8);
  //MatrixXd L_new = Chol(R_tilde_new);
  LLT<MatrixXd> lltOfR(R_tilde_new);             // compute the cholesky decomposition of R called lltofR
  MatrixXd L_new = lltOfR.matrixL();   //retrieve factor L  in the decomposition
  return_list[4]=L_new;
  
  MatrixXd L_inv_X_new;
  MatrixXd R_inv_X_new;
  L_inv_X_new=L_new.triangularView<Lower>().solve(X_new);
  R_inv_X_new=L_new.transpose().triangularView<Upper>().solve(L_inv_X_new); //one forward and one backward to compute R.inv%*%X
  
  //return_list[5]=L_inv_X_new;
  MatrixXd Xt_R_inv_X_new=X_new.transpose()*R_inv_X_new;
  MatrixXd LX_new = sqrt(Xt_R_inv_X_new.array());
  return_list[5]=LX_new;
  
  MatrixXd L_inv_y_new;
  MatrixXd R_inv_y_new;
  L_inv_y_new=L_new.triangularView<Lower>().solve(output_new);
  R_inv_y_new = L_new.transpose().triangularView<Upper>().solve(L_inv_y_new); 
  
  //return_list[7]=L_inv_y_new;
  //MatrixXd yt_R_inv_new= (L_new.transpose().triangularView<Upper>().solve(L_new.triangularView<Lower>().solve(output_new))).transpose(); 
  MatrixXd Xt_R_inv_y_new= X_new.transpose()*R_inv_y_new;
  MatrixXd theta_hat_new=LX_new.transpose().triangularView<Upper>().solve(LX_new.triangularView<Lower>().solve(Xt_R_inv_y_new)); 
  return_list[6]=theta_hat_new;
  
  VectorXd S_2_all_new=VectorXd::Zero(k);
  
  for(int loc_i=0;loc_i<k;loc_i++){
    //S_2_all_new[loc_i]=(yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0);
    S_2_all_new[loc_i]=(output_new.col(loc_i).transpose()*R_inv_y_new.col(loc_i))(0,0)-(Xt_R_inv_y_new.col(loc_i).transpose()*Xt_R_inv_y_new.col(loc_i))(0,0)/Xt_R_inv_X_new(0,0);   
  }
  return_list[7]=S_2_all_new/(num_obs_new-q);
  
  return return_list;
}


