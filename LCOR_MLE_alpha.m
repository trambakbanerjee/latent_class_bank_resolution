function [ NegProbitMLE] = LCOR_MLE_alpha(theta,beta1,sig,X,Y,W)
%Negative Likelihood function of the latent class ordinal probit model

Alpha = theta';
Betac1 = beta1(:,1);
Betac2 = beta1(:,2);
sigma1 = sqrt(sig(1));
sigma2 = sqrt(sig(2));
r0 = length(Alpha);

mu0 = zeros(1,r0);
S0=3*eye(r0);


gamma = [-Inf,0,1,Inf];


    
    ProbS = (normcdf(W*Alpha,0,1));
    B1 = normcdf(gamma(Y+1)'-X*Betac1,0,sigma1)- normcdf(gamma(Y)'-X*Betac1,0,sigma1);
    B2 = normcdf(gamma(Y+1)'-X*Betac2,0,sigma2)- normcdf(gamma(Y)'-X*Betac2,0,sigma2);
    B = B1.*(1-ProbS)+B2.*ProbS;
    ProbitMLE = log(B);
  


NegProbitMLE = -sum(ProbitMLE)-log(mvnpdf(Alpha',mu0,S0));

end

