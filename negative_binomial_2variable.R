library(rstan)
wf <- read.csv('C:\\Users\\gas01500\\Desktop\\Stan Hurdle Model\\whitefly.csv')
data_list = list(N = nrow(wf), N_trt = length(unique(wf$trt)), N_block = length(unique(wf$block)), 
                 y = wf$imm, trt = wf$trt, blk = wf$block)

stan_code = '
  data {
  int<lower=0> N;
  int<lower=0> N_trt;
  int<lower =0> N_block;
  int<lower=0> y[N];
  int trt[N];
  int blk[N];
  }
parameters {
  real<lower=0, upper=1> theta;

  real <lower = 0> lambda[N_trt];
  real <lower = 0> beta_block[N_block];

  real trt_mean;
  real<lower=0> trt_sd;

  real block_mean;
  real<lower=0> block_sd;
}

model {

  for(j in 1:N_trt)
    lambda[j] ~ normal(trt_mean, trt_sd);
  for(j in 1:N_block)
    beta_block[j] ~ normal(block_mean, block_sd);
    
  trt_mean ~ normal(0, 10);
  trt_sd ~ cauchy(0,5);

  block_mean ~ normal(0, 10);
  block_sd ~ cauchy(0,10);

  for (i in 1:N) {
    if (y[i] == 0)
      target += log(theta);
    else
      target += log1m(theta) + poisson_log_lpmf(y[i]|lambda[trt[i]]+beta_block[blk[i]]) - poisson_lccdf(0|exp(lambda[trt[i]]+beta_block[blk[i]]));
}
}
'
stan_run <- stan(data = data_list, model_code = stan_code, control = list(adapt_delta = 0.80))

### > 200 divergent points which is no good 

beta_means = apply(extract(stan_run, pars = c('lambda'))$lambda,2,'mean')
beta_means2 = apply(extract(stan_run, pars = c('beta_block'))$betablock,2,'mean')

theta_mean = mean(extract(stan_run, pars = 'theta')$theta)
y_sim_mean = exp(beta_means[wf$trt])
rZIP = function(mean, theta) {
  pois = rpois(length(mean), mean)
  pois[runif(length(mean))<theta] = 0
  return(pois)
}
y_sim = rZIP(y_sim_mean, theta_mean)

hist(wf$imm, breaks = seq(0,max(wf$imm)))
hist(y_sim, breaks = seq(0,max(wf$imm)), 
     add = TRUE, col = rgb(0.75,0.75,0.75,0.4))
