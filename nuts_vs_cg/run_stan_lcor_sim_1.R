library(rstan)
library(ggmcmc)
library(R.matlab)
library(doParallel)

reps = 100
n = 1200


cl <- makeCluster(6)
registerDoParallel(cl)
results <- foreach(r = 1:reps, .packages = c("rstan","R.matlab"))%dopar%{

# read in data files
  
  dat = readMat(paste('set_1_n',n,'_',r,'.mat',sep=''))
  X1 <- dat$simdata[[2]]
  W1 <- dat$simdata[[1]]
  y <- dat$simdata[[3]]
  trueVal <- dat$simdata[[4]]
  p = ncol(X1);
  k = ncol(W1);
  trueValFlip_both = matrix(0,length(trueVal),1)
  trueValFlip_both[1:p] = trueVal[(p+1):(2*p)] 
  trueValFlip_both[(p+1):(2*p)] = trueVal[1:p] # flip across the latent classes
  trueValFlip_both[(2*p+1):(2*p+k)] = -1*trueVal[(2*p+1):(2*p+k)] 
  trueValFlip_both[(2*p+k+1):length(trueVal)] = trueVal[(2*p+k+1):length(trueVal)]
  
  trueValFlip_alpha = trueVal
  trueValFlip_alpha[(2*p+1):(2*p+k)] = -1*trueVal[(2*p+1):(2*p+k)]
  
  trueValFlip_betas = trueVal
  trueValFlip_betas[1:p] = trueVal[(p+1):(2*p)] 
  trueValFlip_betas[(p+1):(2*p)] = trueVal[1:p]
  
  truevals = c(trueVal[9:10],
               trueVal[2:4], trueVal[6:8])
  truevalsflip_both = c(trueValFlip_both[9:10],
                        trueValFlip_both[2:4], trueValFlip_both[6:8])
  truevalsflip_alpha = c(trueValFlip_alpha[9:10],
                         trueValFlip_alpha[2:4], trueValFlip_alpha[6:8])
  truevalsflip_betas = c(trueValFlip_betas[9:10],
                         trueValFlip_betas[2:4], trueValFlip_betas[6:8])
  
  cutoffs <- c(-Inf,0,1,Inf)
  
  X = X1[,2:4]; # in settings where there is only one covariate apart from the intercept
  W = W1[,2,drop = FALSE];
  N = nrow(y);
  p = ncol(X);
  k = ncol(W);
  J = 3;
  
  # Specify hyperparameters
  alpha0_mean = 0;
  alpha0_scale = 1;
  alpha1_mean = 0;
  alpha1_scale = 3;
  
  fit1 <- stan(
    file = "lcor.stan",  # Stan program
    data = list(N=N,p = p,k = k, y = y, X = X, W = W,J = J, alpha0_mean = alpha0_mean, alpha0_scale = alpha0_scale, alpha1_mean = alpha1_mean, alpha1_scale = alpha1_scale), # data elements
    chains = 1,             # number of Markov chains
    #diagnostic_file = "check_stan_run.csv", # write out diagnostic file
    #init = initf2,          # initial values
    warmup = 1000,          # number of warmup iterations per chain
    iter = 10000,            # total number of iterations per chain
    cores = 1,              # number of cores (could use one per chain)
    refresh = 1             # no progress shown
  )
  
  list_of_draws <- rstan::extract(fit1)
  print(names(list_of_draws))
  fit_summary <- summary(fit1)
  c2 = fit_summary$summary[10]
  c1 = fit_summary$summary[9]
  postmean = c(fit_summary$summary[1:2,1],
                   fit_summary$summary[3:8,1]/(c2-c1))
  postsd= c(fit_summary$summary[1:2,3],
                 fit_summary$summary[3:8,3]/(c2-c1))
  
  return(list('postmean'=postmean,'postsd'=postsd,'truevals'=truevals,
              'truevalsflip_both'=truevalsflip_both,
              'truevalsflip_alpha'=truevalsflip_alpha,
              'truevalsflip_betas'=truevalsflip_betas))
}

save.image('stan_lcor_sim_1.RData')
stopCluster(cl)
registerDoSEQ()

coverage_1sd = coverage_2sd = postmean = postsd = matrix(0,8,reps)
postmean_clean = postsd_clean = matrix(0,8,reps)

truevals = results[[1]]$truevals
truevalsflip_both = results[[1]]$truevalsflip_both
truevalsflip_alpha = results[[1]]$truevalsflip_alpha
truevalsflip_betas = results[[1]]$truevalsflip_betas

for(r in 1:reps){
  
  postmean[,r] = results[[r]]$postmean
  postsd[,r] = results[[r]]$postsd
  
  cmtrue = truevals[1:2] 
  cmmean = postmean[1:2,r]
  cmsd = postsd[1:2,r] 
  flipped_alpha = 0
  ii = sum(1*(-cmtrue > -2*cmsd+cmmean & -cmtrue < 2*cmsd+cmmean))
  flipped_alpha =  ii
  
  cmtrue = truevals[3:8] 
  cmmean = postmean[3:8,r]
  cmsd = postsd[3:8,r] 
  flipped_betas = 0
  ii = sum(1*(cmtrue[4:6] > -2*cmsd[1:3]+cmmean[1:3] & cmtrue[4:6] < 2*cmsd[1:3]+cmmean[1:3]))
  jj = sum(1*(cmtrue[1:3] > -2*cmsd[4:6]+cmmean[4:6] & cmtrue[1:3] < 2*cmsd[4:6]+cmmean[4:6]))
  flipped_betas = ii+jj
  
  flipped_both = 1*(flipped_alpha >= 1 & flipped_betas >= 1)
  
  if (flipped_both == 0){
    
    if(flipped_alpha >= 1){
      
      coverage_1sd[,r] = 1*(truevalsflip_alpha>(postmean[,r]-postsd[,r]) &
                              truevalsflip_alpha<(postmean[,r]+postsd[,r]))
      
      coverage_2sd[,r] = 1*(truevalsflip_alpha>(postmean[,r]-2*postsd[,r]) &
                              truevalsflip_alpha<(postmean[,r]+2*postsd[,r]))
      tt = postmean[1:2,r]
      postmean_clean[,r] = c(-tt,postmean[3:8,r])
      postsd_clean[,r] = postsd[,r]
      
    } else if(flipped_betas >= 1){
      
      coverage_1sd[,r] = 1*(truevalsflip_betas>(postmean[,r]-postsd[,r]) &
                              truevalsflip_betas<(postmean[,r]+postsd[,r]))
      
      coverage_2sd[,r] = 1*(truevalsflip_betas>(postmean[,r]-2*postsd[,r]) &
                              truevalsflip_betas<(postmean[,r]+2*postsd[,r]))
      tt = postmean[3:8,r]
      postmean_clean[,r] = c(postmean[1:2,r],tt[4:6],tt[1:3])
      tt = postsd[3:8,r]
      postsd_clean[,r] = c(postsd[1:2,r],tt[4:6],tt[1:3])
      
    } else{
      
      coverage_1sd[,r] = 1*(truevals>(postmean[,r]-postsd[,r]) &
                              truevals<(postmean[,r]+postsd[,r]))
      
      coverage_2sd[,r] = 1*(truevals>(postmean[,r]-2*postsd[,r]) &
                              truevals<(postmean[,r]+2*postsd[,r]))
      postmean_clean[,r] = postmean[,r]
      postsd_clean[,r] = postsd[,r]
      
    }
  } else if (flipped_both == 1){
    coverage_1sd[,r] = 1*(truevalsflip_both > (postmean[,r]-postsd[,r]) &
                            truevalsflip_both < (postmean[,r]+postsd[,r]))
    
    coverage_2sd[,r] = 1*(truevalsflip_both > (postmean[,r]-2*postsd[,r]) &
                            truevalsflip_both < (postmean[,r]+2*postsd[,r]))
    
    tt = postmean[,r]
    postmean_clean[,r] = c(-tt[1:2],tt[4:6],tt[1:3])
    tt = postsd[,r]
    postsd_clean[,r] = c(postsd[1:2,r],tt[4:6],tt[1:3])
  }
}
save.image('stan_lcor_sim_1.RData')
