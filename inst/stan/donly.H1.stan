

data {
 int<lower = 0> N1;  // replicates condition A
 int<lower = 0> N2;  // replicates condition B
 real<lower = 0> f0; // common initial abundance
 real<lower = 0> t;  // time
 int<lower = 0> co1[N1];  // total counts condition A
 real<lower = 0> a1[N1];  // alpha for condition A
 real<lower = 0> b1[N1];  // beta for condition A
 int<lower = 0> co2[N2];  // total counts condition B
 real<lower = 0> a2[N2];  // alpha for condition B
 real<lower = 0> b2[N2];  // beta for condition B
 real<lower = 0> SP;  // prior standard deviation
}

parameters {
 real pl1[N1];
 real mu1;
 real<lower=0> sigma1;
 real pl2[N2];
 real mu2;
 real<lower=0> sigma2;
}

transformed parameters {
 real<lower=0,upper=1> p1[N1];
 real<lower=0,upper=1> p2[N2];
 real lfcd = mu1-mu2;
 for(i in 1:N1) {
  p1[i] = 1-exp(-t*2^(pl1[i]))*f0/co1[i];
 };
 for(i in 1:N2) {
  p2[i] = 1-exp(-t*2^(pl2[i]))*f0/co2[i];
 };
}

model {
 mu1 ~ normal(0,3);
 sigma1 ~ normal(0,SP);
 mu2 ~ normal(0,3);
 sigma2 ~ normal(0,SP);
 for (i in 1:N1) {
   pl1[i] ~ normal(mu1,sigma1);
   p1[i] ~ beta(a1[i],b1[i]);
 };
 for (i in 1:N2) {
   pl2[i] ~ normal(mu2,sigma2);
   p2[i] ~ beta(a2[i],b2[i]);
 };
}
