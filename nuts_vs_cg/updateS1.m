function [s_updated] = updateS1(g,Y,W,X,beta1,sig,cpts)
% Update S (Class Indicators and Corresponding Sample Sizes)
 
 n=size(W,1);				% Number of subjects
 r=size(W,2);				% Number of class reg parms, incl intcpt			
 q=size(X,2);			  % Number of covariates 

%  g=reshape(g,[,r]);	% Convert g to matrix if vector
 K=size(g,1)+1;  			% Number of classes

 gamma = cpts;
 
 eta=[-W*g' W*g'];
 pi = zeros(n,2);
 
 for ik = 1:size(eta,2)
 pi(:,ik)=normcdf(eta(:,ik),0,1); % Class proportions, given W
 end   

 % Likelihood contributions, Classes 1-K
   l=zeros(n,K);		% Like contribution for subject i, class k
   pnum=zeros(n,K);		% Numerator for post prob of C_i=k
 
 for k = 1:K
   mu1=X*beta1(:,k);
   Sig1 = sqrt(sig(k));
   l(:,k)=normcdf(gamma(k,Y+1)',mu1,Sig1*ones(n,1))-normcdf(gamma(k,Y)',mu1,Sig1*ones(n,1));
   pnum(:,k)=pi(:,k).*l(:,k);
 end

 p = zeros(n,2);
 for ik = 1:size(pnum,2)
 p(:,ik)=pnum(:,ik)./(sum(pnum,2));
 end
 p(isnan(p))=1/K;


 s=binornd(1,p(:,1));
 s_updated(s==1) = 1;
 s_updated(s==0) = 2;


end

