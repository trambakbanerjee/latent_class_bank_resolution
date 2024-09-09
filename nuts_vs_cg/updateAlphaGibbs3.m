function[alpha] = updateAlphaGibbs3(a,beta,sig,X,Y,W)
 % Function Parameters:
  % a    = k x 1 matrix of class-membership regr. parms 
  % W    = n x k matrix of covariates including intercept
  % Returns k x 1 alpha matrix
  % Independence proposal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = [-Inf,0,1,Inf];

 n=size(W,1);				% Number of subjects
 k=size(W,2);				% Number of regression parms
 mu0=zeros(k,1);            % Prior mean of alpha (for all classes)   
%  S0=10*eye(r); 			    % Assume diffuse prior variance
 S0=3*eye(k);               % changing to ensure consistency with other simulations

 nu = 10;
 Q = zeros(n,2);

% Draw Candidate	
  options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'MaxFunEvals', 10000,'TolFun', 10^-40, 'TolX', 10^-40);

  g = @(theta)LCOR_MLE_alpha(theta,beta,sig,X,Y,W);

  [aStar,~,~,~,~,Hessian] = fminunc(g,a,options);

  Covar = inv(Hessian);
  [corr,sd] = corrcov(Covar); 
    
  anew = aStar + sd'.*mvtrnd(corr,nu);
    
  lpold=log(mvnpdf(a,mu0',S0));		% Old Log Prior
  lpnew=log(mvnpdf(anew,mu0',S0));     % New Log Prior  

  lpropNew = log(mvtpdf((anew-aStar)./sd',corr,nu)); % New Log proposal
  lpropOld = log(mvtpdf((a-aStar)./sd',corr,nu)); % Old Log proposal
  
  Q(:,1) = normcdf(gamma(Y+1)',X*beta(:,1),sqrt(sig(1)))-normcdf(gamma(Y)',X*beta(:,1),sqrt(sig(1)));
  Q(:,2) = normcdf(gamma(Y+1)',X*beta(:,2),sqrt(sig(2)))-normcdf(gamma(Y)',X*beta(:,2),sqrt(sig(2)));
 
 % Old and New Class Probabilities 
 
  etaold=[-W*a' W*a'];
  
  pold1 = zeros(n,2);
  
  for ik = 1:size(etaold,2)
  pold1(:,ik) =normcdf(etaold(:,ik),0,1);
  end
  
  pold = pold1.*Q;
    
  etanew=[-W*anew' W*anew'];
  
  pnew1 = zeros(n,2);
  
  for ik = 1:size(etanew,2)
  pnew1(:,ik) =normcdf(etanew(:,ik),0,1);  
  end
  
  pnew = pnew1.*Q;

 % Old and New log likelihoods

  lold =log(sum(pold,2));	
  lnew =log(sum(pnew,2));


 % MH Acceptance Ratio on Log Scale
   ratio=min(0,sum(lnew)+sum(lpnew)+lpropOld-(sum(lold)+sum(lpold)+lpropNew));
   if(log(unifrnd(0,1))<ratio); 
   a=anew;
   end
 
  alpha = a;


end

