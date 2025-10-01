functions {
  real gpareto_lpdf(vector y, real k, real sigma) {
    // generalised Pareto log pdf 
    int N = rows(y);
    real inv_k = inv(k);
    if (k<0 && max(y)/sigma > -inv_k)
      reject("k<0 and max(y)/sigma > -1/k; found k, sigma =", k, sigma);
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma);
    if (fabs(k) > 1e-15)
      return -(1+inv_k)*sum(log1p((y) * (k/sigma))) -N*log(sigma);
    else
      return -sum(y)/sigma -N*log(sigma); // limit k->0
  }

  real gpareto_lcdf(vector y, real k, real sigma) {
    // generalised Pareto log cdf
    real inv_k = inv(k);
    if (k<0 && max(y)/sigma > -inv_k)
      reject("k<0 and max(y)/sigma > -1/k; found k, sigma =", k, sigma);
    if (sigma<=0)
      reject("sigma<=0; found sigma =", sigma);
    if (fabs(k) > 1e-15)
      return sum(log1m_exp((-inv_k)*(log1p((y) * (k/sigma)))));
    else
      return sum(log1m_exp(-(y)/sigma)); // limit k->0
  }
}

data {
  // the input data
  int<lower = 1> n;               // total number of obs
  int<lower = 0> n_cens;          // number of censored obs
  vector<lower = 0>[n] helium;      // number of regular obs
  int<lower=1> cens_id[n_cens];
  int<lower=1> nocens_id[n-n_cens];
  
  // parameters for the prior
  real<lower = 0> a;
  real<lower = 0> b;
}

parameters {
  real<lower = 0> k;
  real<lower = 0> sigma;
}

model {
 // prior
 k ~ normal(0, b);
 sigma ~ normal(0,a);

 // likelihood
 target += gpareto_lcdf(helium[cens_id] | k, sigma);  
 target += gpareto_lpdf(helium[nocens_id] | k, sigma);
 
 }
