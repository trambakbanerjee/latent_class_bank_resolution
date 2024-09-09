
function postSamp = GibbsAlgo_Sim(W,X,Y,n)
% MCMC procedure for the latent class model with ordinal response

p = size(X,2);
k = size(W,2);

% Priors

BetaVec1 = zeros(p,1);
BetaVec2 = zeros(p,1);
%% ******************************

% Beta prior hyperparameters

BetaCovar1 = 1*eye(p);
invBCovar = inv(BetaCovar1);

% Sigma prior hyperparameters

nu1=8.658;
d1=2.6;

nu2=8.658;
d2=2.6;

% number of iterations

T = 11000;

% Burn-in

burn = 1000;

% Initialize vectors and matrices

z = zeros(n,1);
BetaOProbit = zeros(p*2,T);
Atmp=zeros(T,k);	
alpha = zeros(1,k);
sigma = ones(T,2);

s1 = binornd(1,0.5,1,n);
s = s1+1;

%Cut-points

gamma(1,:) = [-Inf,0,1,Inf];
gamma(2,:) = [-Inf,0,1,Inf];

b = [-Inf,0,Inf];
   
sig = [1 1];    
beta1 = zeros(p,2);
beta1(:,2) = 1;
    
for it = 1:T-1
    s=updateS1(alpha,Y,W,X,beta1,sig,gamma);    

     cntS = tabulate(s); 

%% Conditional of z given beta, s

for is = 1:2

    z(s == is) = truncnormFuncVec(X(s==is,:)*beta1(:,is),sqrt(sig(is))*ones(cntS(is,2),1),gamma(is,Y(s==is))',gamma(is,Y(s==is)+1)');
    
end    
    
    X1 = X(s==1,:);
    X2 = X(s==2,:);
    
    z1 = z(s==1);
    z2 = z(s==2);
    
    BetaCovarHat1 = inv(invBCovar + (X1'*X1)/sig(1));
    BetaCovarHat2 = inv(invBCovar + (X2'*X2)/sig(2));

%% Draw from conditional of Beta given s, z
      
    BetaHat1 = BetaCovarHat1*(invBCovar*BetaVec1 + (X1'*z1)/sig(1));
    BetaHat2 = BetaCovarHat2*(invBCovar*BetaVec2 + (X2'*z2)/sig(2));

    beta1(:,1) = mvnrnd(BetaHat1,BetaCovarHat1);
    beta1(:,2) = mvnrnd(BetaHat2,BetaCovarHat2);
    
%% Draw from conditional of alpha given s    

% Independence MH - Collapsed Gibbs
 alpha=updateAlphaGibbs3(alpha,beta1,sig,X,Y,W);  
 
 %% Draw from conditional of sigma given y, beta, s, z and alpha
 
 C1 = (z1 - X1*beta1(:,1))'*(z1 - X1*beta1(:,1));
 C2 = (z2 - X2*beta1(:,2))'*(z2 - X2*beta1(:,2));

 
    dHat1 = 2/(d1+C1);
    nuHat1 = (nu1+cntS(1,2))/2;
    prec1 = gamrnd(nuHat1,dHat1);
    sig(1) = prec1^(-1);
     
    dHat2 = 2/(d2+C2);
    nuHat2 = (nu2+cntS(2,2))/2;
    prec2 = gamrnd(nuHat2,dHat2);
    sig(2) = prec2^(-1);
    
%% Updated draws

BetaOProbit(1:p,it+1) = beta1(:,1);
BetaOProbit(p+1:2*p,it+1) = beta1(:,2);
sigma(it+1,:) = sig;
Atmp(it+1,:)= alpha;	
end


postMean(1:2*p,:) = mean(BetaOProbit(:,burn:end),2);
postMean(2*p+1:2*p+k,:) = mean(Atmp(burn:end,:),1)';
postMean(2*p+1+k:2*p+2+k,:) = mean(sigma(burn:end,:),1)';

postSD(1:2*p,:) = std(BetaOProbit(:,burn:end),1,2);
postSD(2*p+1:2*p+k,:) = std(Atmp(burn:end,:),1,1)';
postSD(2*p+1+k:2*p+2+k,:) = std(sigma(burn:end,:),1,1)';

postSamp.BetaOProbit = BetaOProbit;
postSamp.sigma = sigma;
postSamp.Atmp = Atmp;
postSamp.postMean = postMean;
postSamp.postSD = postSD;
