function [ normrndvalue ] = truncnormFuncVec( mu, sigma0, A,B )
% truncated normal random numbers/vectors

dim = length(sigma0);
normrndvalue = zeros(dim,1);


U = unifrnd(0,1,dim,1);

normrndvalue = norminv(U.*(normcdf(B,mu,sigma0)-normcdf(A,mu,sigma0))+normcdf(A,mu,sigma0),mu,sigma0);
    
end


