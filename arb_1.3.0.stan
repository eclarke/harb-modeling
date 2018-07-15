data{
  int<lower=1> n;
  int<lower=1> n_specimen_id2;
  int read_count[n];
  int total_reads[n];
  int specimen_id2[n];
}
parameters{
  vector[n] alpha;
  real mu;
  real<lower=0> sigma;
}
model{
  mu ~ normal(-1, 1);
  sigma ~ normal(0, 1);
  alpha ~ normal(mu, sigma);
  read_count ~ binomial_logit(total_reads, alpha)
}
