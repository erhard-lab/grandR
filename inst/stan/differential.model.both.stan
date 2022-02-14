

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
 real log2d1[N1];
 real log2d2[N2];
 real mu1;
 real<lower=0> sigma1;
 real mu2;
 real<lower=0> sigma2;
 real smu1;
 real<lower=0> ssigma1;
 real smu2;
 real<lower=0> ssigma2;
}

transformed parameters {
 real log2s1[N1];
 real log2s2[N2];
 real<lower=0,upper=1> ntr1[N1];
 real<lower=0,upper=1> ntr2[N2];
 real lfcd = mu1-mu2;
 real lfcs = smu1-smu2;
 for(i in 1:N1) {
  ntr1[i] = 1-exp(-t*2^(log2d1[i]))*f0/co1[i];
  log2s1[i] = log2(ntr1[i]*co1[i]*2^log2d1[i]/(1-exp(-t*2^log2d1[i])));
 };
 for(i in 1:N2) {
  ntr2[i] = 1-exp(-t*2^(log2d2[i]))*f0/co2[i];
  log2s2[i] = log2(ntr2[i]*co2[i]*2^log2d2[i]/(1-exp(-t*2^log2d2[i])));
 };
}

model {
 mu1 ~ normal(0,3);
 sigma1 ~ normal(0,SP);
 mu2 ~ normal(0,3);
 sigma2 ~ normal(0,SP);
 smu1 ~ normal(log2(f0),4);
 ssigma1 ~ normal(0,SP);
 smu2 ~ normal(log2(f0),4);
 ssigma2 ~ normal(0,SP);
 for (i in 1:N1) {
   log2d1[i] ~ normal(mu1,sigma1);
   log2s1[i] ~ normal(smu1,ssigma1);
   //target += log( 1- exp(-t*2^log2d1[i])*t*2^log2d1[i]/(1-exp(-t*2^log2d1[i])));
   ntr1[i] ~ beta(a1[i],b1[i]);
   //target += log(f0/co1[i])+t*2^log2d1[i]+log(t)+log2d1[i]*log(2)+log(log(2));
 };
 for (i in 1:N2) {
   log2d2[i] ~ normal(mu2,sigma2);
   log2s2[i] ~ normal(smu2,ssigma2);
   //target += log( 1- exp(-t*2^log2d2[i])*t*2^log2d2[i]/(1-exp(-t*2^log2d2[i])));
   ntr2[i] ~ beta(a2[i],b2[i]);
   //target += log(f0/co2[i])+t*2^log2d2[i]+log(t)+log2d2[i]*log(2)+log(log(2));
 };
}
