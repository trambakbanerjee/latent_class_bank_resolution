//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
data {
  int<lower=2> J;
  int<lower=0> N;
  int<lower=1> p;
  int<lower=1> k;
  int<lower=1,upper=J> y[N,1];
  matrix[N,p] X; // predictor matrix for the outcome
  matrix[N,k] W; // predictor matrix for latent classes
    real alpha0_mean;
  real alpha0_scale;
  real alpha1_mean;
  real alpha1_scale;
}

parameters {
  real alpha0; //intercept in the class membership model
  vector[k] alpha1; //coefficients in the class membership model
  vector[p] beta1; //coefficients in latent class 1
  vector[p] beta2; //coefficients in latent class 2
  ordered[J-1] c; //cut-points
}

model {
alpha0 ~ normal(alpha0_mean, alpha0_scale);
alpha1 ~ normal(alpha1_mean, alpha1_scale);
    vector[J] theta1;
    vector[J] theta2;
    vector[J] theta;
    vector[N] mu1a;
    vector[N] mu1;
    mu1a = alpha0 + (W*alpha1);
    mu1 = Phi_approx(mu1a);    
for (i in 1:N) {
    real eta1; real eta2;
    eta1 = X[i,] * beta1; 
    eta2 = X[i,] * beta2;
    theta1[1] = 1 - Phi_approx(eta1 - c[1]);
    theta2[1] = 1 - Phi_approx(eta2 - c[1]);
    theta[1] = mu1[i]*theta1[1]+(1-mu1[i])*theta2[1]; 
    for (j in 2:(J-1)){
    theta1[j] = Phi_approx(eta1 - c[j-1]) - Phi_approx(eta1 - c[j]);
    theta2[j] = Phi_approx(eta2 - c[j-1]) - Phi_approx(eta2 - c[j]);
    theta[j] = mu1[i]*theta1[j]+(1-mu1[i])*theta2[j]; 
    }
    theta1[J] = Phi_approx(eta1 - c[J-1]);
    theta2[J] = Phi_approx(eta2 - c[J-1]);
    theta[J] = mu1[i]*theta1[J]+(1-mu1[i])*theta2[J]; 
    y[i,1] ~ categorical(theta);
}  
}

