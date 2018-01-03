library(rstan)

### Taken from Andrew Parnell 
stan_code = '
data {
int<lower=0> N;
int<lower=0> N_trt;
int<lower=0> y[N];
int trt[N];    
}
parameters {
real<lower=0, upper=1> q_0;
real beta_trt[N_trt];
real trt_mean;
real<lower=0> trt_sd;
}
model {
for(j in 1:N_trt)
beta_trt[j] ~ normal(trt_mean, trt_sd);
trt_mean ~ normal(0, 10);
trt_sd ~ cauchy(0, 5);

for (i in 1:N) {
if (y[i] == 0)
target += log(q_0);
else
target += log1m(q_0) + poisson_log_lpmf(y[i] | beta_trt[trt[i]])
- poisson_lccdf(0 | exp(beta_trt[trt[i]]));
}
}
'

stan_run = stan(data = list(N = nrow(wf), 
                            N_trt = length(unique(wf$trt)),
                            y = wf$imm,
                            trt = wf$trt),
                model_code = stan_code)

beta_means = apply(extract(stan_run, pars = 'beta_trt')$beta_trt,2,'mean')
q_0_mean = mean(extract(stan_run, pars = 'q_0')$q_0)
y_sim_mean = exp(beta_means[wf$trt])
rZIP = function(mean, q_0) {
  pois = rpois(length(mean), mean)
  pois[runif(length(mean))<q_0] = 0
  return(pois)
}
y_sim = rZIP(y_sim_mean, q_0_mean)