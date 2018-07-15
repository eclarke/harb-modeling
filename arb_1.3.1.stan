data{
  int<lower=1> n;
  int<lower=1> n_specimen_id2;
  int read_count[n];
  int total_reads[n];
  int specimen_id2[n];
}
transformed data {
  vector[n] emp_prob;
  for (i in 1:n) {
    emp_prob[i] = (read_count[i] * 1.0)/total_reads[i];
  }
}
parameters{
  vector[n] alpha_std;
  real mu;
  real<lower=0> sigma;
  real b_lag;
}
model{
  mu ~ normal(-1, 1);
  sigma ~ normal(0, 1);
  alpha_std ~ normal(0, 1);
  b_lag ~ normal(0, 1);
  read_count[2:n] ~ binomial_logit(total_reads[2:n], (mu + alpha_std[2:n]*sigma) + b_lag * emp_prob[1:(n-1)]);
}
