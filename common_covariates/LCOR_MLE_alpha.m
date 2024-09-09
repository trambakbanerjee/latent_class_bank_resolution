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

% Gradient and Hessian

% d0 = normpdf(W*Alpha)./B;
% d1 = repmat(d0,1,r0);
% Bd1 = repmat(B1,1,r0);
% Bd2 = repmat(B2,1,r0);
% g0 = d1.*(W.*Bd1 - W.*Bd2);
% g = sum(g0,1)' + (S0)\(Alpha - mu0);
% 
% p0 = W.*W;
% Wx = [W(:,2) W(:,1)];
% p1 = Wx.*W;
% q = zeros(r0,r0,n);
% q(1,1,:) = p0(:,1);
% q(1,2,:) = p1(:,1);
% q(2,1,:) = p1(:,1);
% q(2,2,:) = p0(:,2);
% 
% A0 = repmat(Alpha',n,1);
% t0 = W.*A0;
% t1 = Wx.*A0;
% v = zeros(r0,r0,n);
% v(1,1,:) = t0(:,1);
% v(1,2,:) = t1(:,2);
% v(2,1,:) = t1(:,1);
% v(2,2,:) = t0(:,2);
% 
% 
% u = zeros(r0,r0,n);
% for in = 1:n
% u(:,:,in) = v(:,:,in)*q(:,:,in);
% end
% 
% d4 = zeros(r0,r0,n);
% d4(1,1,:) = -d0;
% d4(1,2,:) = -d0;
% d4(2,1,:) = -d0;
% d4(2,2,:) = -d0;
% 
% Bd14 = zeros(r0,r0,n);
% Bd24 = zeros(r0,r0,n);
% 
% Bd14(1,1,:) = B1;
% Bd14(1,2,:) = B1;
% Bd14(2,1,:) = B1;
% Bd14(2,2,:) = B1;
% 
% Bd24(1,1,:) = B2;
% Bd24(1,2,:) = B2;
% Bd24(2,1,:) = B2;
% Bd24(2,2,:) = B2;
% 
% h1 = Bd14.*d4.*u - Bd24.*d4.*u; 
% e4 = d4.*d4;
% hg = g0.*g0;
% g00 = [g0(:,2) g0(:,1)];
% hg1 = g00.*g0;
% h02 = zeros(r0,r0,n);
% h02(1,1,:) = hg(:,1);
% h02(1,2,:) = hg1(:,1);
% h02(2,1,:) = hg1(:,1);
% h02(2,2,:) = hg(:,2);
% h2 = e4.*h02;
% H = sum(h1,3) + sum(h2,3) + inv(S0);

end

