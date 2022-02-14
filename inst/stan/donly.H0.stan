

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
 real mu;
 real<lower=0> sigma;
 real pl2[N2];
}

transformed parameters {
 real<lower=0,upper=1> p1[N1];
 real<lower=0,upper=1> p2[N2];
 for(i in 1:N1) {
  p1[i] = 1-exp(-t*2^(pl1[i]))*f0/co1[i];
 };
 for(i in 1:N2) {
  p2[i] = 1-exp(-t*2^(pl2[i]))*f0/co2[i];
 };
}

model {
 mu ~ normal(0,3);
 sigma ~ normal(0,SP);
 for (i in 1:N1) {
   pl1[i] ~ normal(mu,sigma);
   p1[i] ~ beta(a1[i],b1[i]);
 };
 for (i in 1:N2) {
   pl2[i] ~ normal(mu,sigma);
   p2[i] ~ beta(a2[i],b2[i]);
 };
}
