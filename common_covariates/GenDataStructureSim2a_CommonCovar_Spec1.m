
function simdata = GenDataStructureSim2a_CommonCovar_Spec1(rep,N,common)
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
W(:,2) = normrnd(0.5,1,1,N);
W(:,3) = normrnd(0.5,1,1,N);
W(:,4) = normrnd(0,1,1,N);

% True BetaXs
TrueAlpha1 = [-0.3; 0.7;0.1; 0.6];

% Param. Spec 1
TrueSig = 0.25;
TrueBetaX(:,1) = [0.6; -0.7; -0.6; 0.5];
TrueBetaX(:,2) = [0.1; 0.6; 0.2; 0.8];
TrueSig1 = TrueSig; TrueSig2 = TrueSig;
Trueval = [TrueBetaX(:,1);TrueBetaX(:,2);TrueAlpha1;TrueSig1;TrueSig2];

% True Latent Variable L

TrueL = W*TrueAlpha1;

% Generate Vector of s

r = zeros(N,1);
pi = normcdf(TrueL,0,1);

r = binornd(1,pi);
r1 = r;

%% Check histogram

%histogram(W(r==0,2));hold on;
%histogram(W(r==1,2));hold off;

%% Ordinal model

% Matrix of X

X(:,1) = ones(N,1);
X(:,2) = normrnd(0.5,1,1,N);
X(:,3) = normrnd(0.5,0.8,1,N);
X(:,4) = normrnd(0,0.8,1,N);

if common == 1
    X(:,2) = W(:,2);
elseif common == 2
    X(:,[2 3]) = W(:,[2 3]);
elseif common == 3
    X(:,[2 3 4]) = W(:,[2 3 4]);
end

% Generate Y for the two latent classes using true latent variable TrueW

indus1 = find(r==0);
indus2 = find(r==1);


TrueW = ones(N,1);
e = normrnd(0,sqrt(TrueSig),N,1);

TrueW(indus1) = X(indus1,:)*TrueBetaX(:,1)+ e(indus1);
TrueW(indus2) = X(indus2,:)*TrueBetaX(:,2) + e(indus2);

Y = zeros(N,1);

Y(TrueW<0) = 1;
Y(TrueW >= 0 & TrueW <1) = 2;
Y(TrueW >= 1) = 3;

simdata.W = W;
simdata.X = X;
simdata.Y = Y;
simdata.truevals = Trueval;
simdata.classmeans = [mean(X(indus1,:)*TrueBetaX(:,1));
    mean(X(indus2,:)*TrueBetaX(:,2))];

% %% Scatter plots
% 
% % X2
% 
% scatter(X(indus1,2),Y(indus1));hold on;
% scatter(X(indus2,2),Y(indus2));hold off;
% 
% % X3
% 
% scatter(X(indus1,3),Y(indus1));hold on;
% scatter(X(indus2,3),Y(indus2));hold off;
% 
% % X4
% 
% scatter(X(indus1,4),Y(indus1));hold on;
% scatter(X(indus2,4),Y(indus2));hold off;
% 
% 
% 
% histogram(X(r==0,:)*TrueBetaX(:,1));hold on;
% histogram(X(r==1,:)*TrueBetaX(:,2));hold off;
% 
% [f1,xi1] = ksdensity(X(r==0,:)*TrueBetaX(:,1));
% [f2,xi2] = ksdensity(X(r==1,:)*TrueBetaX(:,2));
% 
% 
% plot(xi1,f1,'r-.');hold on;
% plot(xi2,f2,'k');hold off;
% legend('Latent Class 1', 'Latent Class 2');

