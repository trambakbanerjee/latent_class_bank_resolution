
function simdata = GenDataStructureSet1_1_noint(rep,N)
% ***Simulate Data for Estimation of Latent Class Ordinal Probit models***
% Variables: W,X,Y
% Generated latent variables: TrueL,r,TrueW
% Change from previous version: 4 covariates instead of 2
% Y is the ordinal response, X is the matrix of covariates with
% coefficients TrueBetaX. TrueW = X*TrueBetaX
% r is the latent response variable determining class membership. W is the
% matrix of covariates determining class membership with coefficents
% TrueAlpha1. TrueL = W*TrueAlpha1
% Change from version 2: Include more counts from category 3 into latent
% class 1
% =======================================================================

%% Class membership model
% global W;

% Matrix of W

rng(rep);
W(:,1) = ones(N,1);
W(:,2) = normrnd(0,1,1,N);

% True BetaXs
TrueAlpha1 = [-0.3; 1.5;];

% Param. Spec 1

 TrueBetaX(:,1) = [0; -0.9; -0.8; -0.6];  % Beta for latent class 1
 TrueBetaX(:,2) = [0; 0.6; 0.6; 1];% Beta for latent class 2

% True Latent Variable L

TrueL = W*TrueAlpha1;

% True sigma^2

sig2_1 = 0.25; % sigma^2 for latent class 1
sig2_2 = 0.25; % sigma^2 for latent class 2


% Generate Vector of s

r = zeros(N,1);
pi = normcdf(TrueL,0,1);

r = binornd(1,pi);
r1 = r;


%% Ordinal model

% Matrix of X

X(:,1) = ones(N,1);
X(:,2) = normrnd(0.5,1,1,N);
X(:,3) = normrnd(0.5,1,1,N);
X(:,4) = normrnd(0,0.8,1,N);

% Generate Y for the two latent classes using true latent variable TrueW
% TrueW here corresponds to the true continuous variable underlying Y -
% referred to as "z" in the paper

indus1 = find(r==0);
indus2 = find(r==1);

TrueW = ones(N,1);

e = zeros(N,1);
N1 = length(indus1);
N2 = length(indus2);

% Draw normal random errors for the two latent classes

e(indus1) = normrnd(0,sqrt(sig2_1),N1,1);
e(indus2) = normrnd(0,sqrt(sig2_2),N2,1);

% Generate TrueW (or z in the paper) as X*beta + e 

TrueW(indus1) = X(indus1,:)*TrueBetaX(:,1)+ e(indus1);
TrueW(indus2) = X(indus2,:)*TrueBetaX(:,2) + e(indus2);

Y = zeros(N,1);

Y(TrueW<0) = 1;
Y(TrueW >= 0 & TrueW <1) = 2;
Y(TrueW >= 1) = 3;

% Write out true values, X, W and Y matrices for ue by stan later
truevals = [TrueBetaX(:,1);TrueBetaX(:,2);TrueAlpha1;sig2_1;sig2_2];
simdata.W = W;
simdata.X = X;
simdata.Y = Y;
simdata.truevals = truevals;
simdata.classmeans = [mean(X(indus1,:)*TrueBetaX(:,1));
    mean(X(indus2,:)*TrueBetaX(:,2))];

FileName = strcat('set_1_n',num2str(N),'_',num2str(rep));
save(FileName, 'simdata');


% writematrix(Y,'C:\Users\j1prs01\Dropbox (Fed Reserve KC)\bayes_latent_class\revision_0324\sim_data\y_spec2.csv');
% writematrix(X,'C:\Users\j1prs01\Dropbox (Fed Reserve KC)\bayes_latent_class\revision_0324\sim_data\X_spec2.csv');
% writematrix(W,'C:\Users\j1prs01\Dropbox (Fed Reserve KC)\bayes_latent_class\revision_0324\sim_data\W_spec2.csv');


