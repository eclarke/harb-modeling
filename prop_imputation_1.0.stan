data{
  int<lower=1> n_obs;
  int<lower=0> n_mis;
  int<lower=1> n;
  int<lower=1> n_subject_id;
  int<lower=1, upper=n_obs+n_mis> ii_obs[n_obs];
  int<lower=1, upper=n_obs+n_mis> ii_mis[n_mis];
  int specimen_id[n];
  int subject_id[n];
  real<lower=0> prop_obs[n_obs];
}
parameters {
  real a_subject[n_subject_id];
  real<lower=0> sigma;
  real prop_mis[n_mis];
}
transformed parameters {
  real mu[n];
  real prop[n];
  prop[ii_obs] = prop_obs;
  prop[ii_mis] = prop_mis;
  for (i in 1:n) {
    mu[i] = a_subject[subject_id[i]];
  }
}
model {
  a_subject ~ normal(0, 0.5);
  sigma ~ cauchy(0, 1);
  prop ~ normal(mu, sigma);
}
