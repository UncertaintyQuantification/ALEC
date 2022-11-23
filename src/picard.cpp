
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*- 

// we only include RcppEigen.h which pulls Rcpp.h in for us 
#include <RcppEigen.h> 
#include <Rcpp.h> 
#include <cmath> 
#include "ctools.h"
using namespace Rcpp;
using namespace std;
using namespace Eigen; 
// [[Rcpp::depends(RcppEigen)]] 


// [[Rcpp::export]]
double integral_approx(VectorXd f, int lb, int ub, float d){ //lb and ub - 1 
  VectorXd v =f.segment(lb,ub-lb+1); 
  double res=(v.array().sum()*2.0-v[0]-v[v.size()-1])*d/2.0;
  return res;
}


// [[Rcpp::export]]
List picard_rho_cpp(VectorXd rho, double beta_mu, int a, int L, 
                        VectorXd exp_n_be_V_ext, float alpha, 
                        int iteration, double threshold){
  int k = rho.size();
  float d = (L+a)/(float(k-1));
  int t = a/d;
  double integral1;
  double p1;
  double p2;
  double integral3;
  VectorXd rho_new=VectorXd::Zero(k);
  VectorXd int2_vec=VectorXd::Zero(t+1);
  VectorXd error = VectorXd::Zero(k);
  VectorXd error_max = VectorXd::Zero(iteration);
  List return_list(3);
  int ite = 1;
  while(ite<=iteration && pow((rho_new-rho).array(),2).sum()>threshold ){
    for(int j=t;j<(k-t);j++){
      integral1 = integral_approx(rho,j-t,j,d);
      p1 = log(1.0-integral1);
      
      for (int i=0; i<=t; i++){
        integral3 = integral_approx(rho,j+i-t,j+i,d);
        int2_vec[i] = rho[j+i]/(1.0-integral3);
      }
      p2 = integral_approx(int2_vec,0,int2_vec.size()-1,d);
      
      rho_new[j] = exp(beta_mu+p1-p2)*exp_n_be_V_ext[j];
      error[j]=rho[j]-rho_new[j];
    }
    error_max[ite]=error.array().abs().matrix().maxCoeff();
    rho = (1-alpha)*rho+alpha*rho_new;
    ite++;
  }
  if (ite >= iteration){
    Rcout << "picard did not converge." << endl;
  }
  return_list[0]=ite;
  return_list[1]=rho;
  return_list[2]=error_max.head(ite);
  
  return return_list;
}

//close form from Percus, only related to density
// [[Rcpp::export]]
double beta_Omega_cpp(VectorXd rho, int a, int L){
  int k = rho.size();
  float d = (L+a)/(float(k-1));
  int t = a/d;
  double integral_denominator;
  
  VectorXd integral_vec = VectorXd::Zero(k-t);
  for (int j=t;j<k;j++){
    integral_denominator = integral_approx(rho,j-t,j,d);
    
    integral_vec[j-t]=rho[j-t]/(1.0-integral_denominator);
  }
  
  double res = integral_approx(integral_vec,0,integral_vec.size()-1,d);
  return -res;
}

