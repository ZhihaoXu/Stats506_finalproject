## Final Project - R code
##
## Author(s): Zhihao Xu, xuzhihao@umich.edu
## Updated: December 4, 2020 

# 79: -------------------------------------------------------------------------

# libraries: ------------------------------------------------------------------
library(parallel)
library(MASS)

# monte carlo study
compute_fdr = function(rej, sig){
  # compute the fdr of rejection set
  # Inputs: 
  #   rej - rejection set
  #   sig - index of signal
  # Output: estimation of fdr
  if (length(rej)==0){
    fdr = 0
  } else{
    fdr = sum(rej %in% sig == FALSE)/length(rej)
  }
  return(fdr)
}

compute_power = function(rej, sig){
  # compute the power of rejection set
  # Inputs: 
  #   rej - rejection set
  #   sig - index of signal
  # Output: estimation of power
  if (length(rej)==0){
    power = 0
  } else{
    power = sum(rej %in% sig)/length(sig)
  }
  return(power)
}

## Independent
ind_sim = function(idx, mu, method="BH", q=0.05){
  # simulation of independent case
  # Inputs: 
  #   mu - mean of simulation
  #   method - method for p-value correction
  #   q - significant level
  # Output: rejection set
  d = matrix(rnorm(200*100, mean=mu,sd=1), nrow=100)
  p = 2*(1-pnorm(abs(apply(d, 1, mean)), mean=0, sd=1/sqrt(200)))
  rej = which(p.adjust(p, method)<q)
  return(rej)
}
result_ind_null_bh = mclapply(1:1000, ind_sim,mu=0, method="BH", 
                              q=0.05, mc.cores = 4)
result_ind_null_by = mclapply(1:1000, ind_sim,mu=0, method="BY", 
                              q=0.05, mc.cores = 4)
bh_ind_null = mean(unlist(lapply(result_ind_null_bh, compute_fdr, sig=0)))
by_ind_null = mean(unlist(lapply(result_ind_null_by, compute_fdr, sig=0)))
result_ind_null = rbind(bh_ind_null, by_ind_null)

result_ind_sig_bh = mclapply(1:1000, ind_sim,
                             mu=c(rep(seq(0.2,1,0.4), each=2), rep(0,94)),
                             method="BH", q=0.05, mc.cores = 4)
result_ind_sig_by = mclapply(1:1000, ind_sim,
                             mu=c(rep(seq(0.2,1,0.4), each=2), rep(0,94)),
                             method="BY", q=0.05, mc.cores = 4)

bh_ind_sig = c(mean(unlist(lapply(result_ind_sig_bh, compute_fdr, sig=1:6))),
               mean(unlist(lapply(result_ind_sig_bh, compute_power, sig=1:6))))
by_ind_sig = c(mean(unlist(lapply(result_ind_sig_by, compute_fdr, sig=1:6))),
               mean(unlist(lapply(result_ind_sig_by, compute_power, sig=1:6))))
result_ind_sig = rbind(bh_ind_sig, by_ind_sig)

## Positive Correlation
Sigma = matrix(rep(0.8,100*100),100)
post_cor_sim = function(idx, mu, method="BH", q=0.05){
  # simulation of positive correlation case
  # Inputs: 
  #   mu - mean of simulation
  #   method - method for p-value correction
  #   q - significant level
  # Output: rejection set
  d = mvrnorm(n=200, mu=mu, Sigma=Sigma)
  p = 2*(1-pnorm(abs(apply(d, 2, mean)), mean=0, sd=1/sqrt(200)))
  rej = which(p.adjust(p, method)<q)
  return(rej)
}

result_post_cor_null_bh = mclapply(1:1000, post_cor_sim,mu=rep(0,100), 
                                   method="BH", q=0.05, mc.cores = 4)
result_post_cor_null_by = mclapply(1:1000, post_cor_sim,mu=rep(0,100), 
                                   method="BY", q=0.05, mc.cores = 4)
bh_post_cor_null = mean(unlist(lapply(result_post_cor_null_bh, 
                                      compute_fdr, sig=0)))
by_post_cor_null = mean(unlist(lapply(result_post_cor_null_by, 
                                      compute_fdr, sig=0)))
result_post_cor_null = rbind(bh_post_cor_null, by_post_cor_null)

result_post_cor_sig_bh = mclapply(1:1000, post_cor_sim,
                                  mu=c(rep(seq(0.2,1,0.4), each=2), rep(0,94)),
                                  method="BH",q=0.05, mc.cores = 4)
result_post_cor_sig_by = mclapply(1:1000, post_cor_sim,
                                  mu=c(rep(seq(0.2,1,0.4), each=2), rep(0,94)),
                                  method="BY",q=0.05, mc.cores = 4)
bh_post_cor_sig = c(mean(unlist(lapply(result_post_cor_sig_bh, 
                                       compute_fdr, sig=1:6))),
                    mean(unlist(lapply(result_post_cor_sig_bh, 
                                       compute_power, sig=1:6))))
by_post_cor_sig = c(mean(unlist(lapply(result_post_cor_sig_by, 
                                       compute_fdr, sig=1:6))),
                   mean(unlist(lapply(result_post_cor_sig_by, 
                                      compute_power, sig=1:6))))
result_post_cor_sig = rbind(bh_post_cor_sig, by_post_cor_sig)

## Negative Correlation
Sigma = matrix(rep(0.8,50*50),50)
neg_cor_sim = function(idx, mu, method="BH", q=0.05){
  # simulation of negative correlation case
  # Inputs: 
  #   mu - mean of simulation
  #   method - method for p-value correction
  #   q - significant level
  # Output: rejection set
  d = mvrnorm(n=200, mu=mu, Sigma=Sigma)
  d = cbind(d,-d)
  p = 2*(1-pnorm(abs(apply(d, 2, mean)), mean=0, sd=1/sqrt(200)))
  rej = which(p.adjust(p,method)<q)
  return(rej)
}

result_neg_cor_null_bh = mclapply(1:1000, neg_cor_sim, mu=rep(0,50), 
                                  method="BH", q=0.05, mc.cores = 4)
result_neg_cor_null_by = mclapply(1:1000, neg_cor_sim, mu=rep(0,50), 
                                  method="BY", q=0.05, mc.cores = 4)
bh_neg_cor_null = mean(unlist(lapply(result_neg_cor_null_bh, 
                                     compute_fdr, sig=0)))
by_neg_cor_null = mean(unlist(lapply(result_neg_cor_null_by, 
                                     compute_fdr, sig=0)))
result_neg_cor_null = rbind(bh_neg_cor_null, by_neg_cor_null)

result_neg_cor_sig_bh = mclapply(1:1000, neg_cor_sim,
                                  mu=c(rep(seq(0.2,1,0.4),each=2), rep(0,44)),
                                  method="BH",q=0.05, mc.cores = 4)
result_neg_cor_sig_by = mclapply(1:1000, neg_cor_sim,
                                  mu=c(rep(seq(0.2,1,0.4),each=2), rep(0,44)),
                                  method="BY",q=0.05, mc.cores = 4)
bh_neg_cor_sig = c(mean(unlist(lapply(result_neg_cor_sig_bh, 
                                      compute_fdr, sig=c(1:6,51:56)))),
                   mean(unlist(lapply(result_neg_cor_sig_bh, 
                                      compute_power, sig=c(1:6,51:56)))))
by_neg_cor_sig = c(mean(unlist(lapply(result_neg_cor_sig_by, 
                                      compute_fdr, sig=c(1:6,51:56)))),
                   mean(unlist(lapply(result_neg_cor_sig_by, 
                                      compute_power, sig=c(1:6,51:56)))))
result_neg_cor_sig = rbind(bh_neg_cor_sig, by_neg_cor_sig)


result_null = rbind(result_ind_null, result_post_cor_null, result_neg_cor_null)
result_sig = rbind(result_ind_sig, result_post_cor_sig, result_neg_cor_sig)
