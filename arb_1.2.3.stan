data{
  int<lower=1> n;
  int<lower=1> n_specimen_id2;
  int read_count[n];
  int total_reads[n];
  int specimen_id2[n];
  int<lower=0, upper=1> abx_yn[n];
  vector[n] emp_prob;  // the empirical, or observed, probabilities of the genus
}
parameters{
  real a_0;
  real b_lag;
  real b_abx;
  vector[n_specimen_id2] a_specimen;
  real<lower=0, upper=1> s_tilde;
}
transformed parameters {
  real sigma_specimen = tan(pi() * s_tilde * 0.5); // alternative parameterization for half-cauchy
}
model{
  vector[n] prob;
  // implicit uniform prior for s_tilde
  a_specimen ~ normal(0 , sigma_specimen);
  a_0 ~ normal(0 , 10);
  b_lag ~ normal(0, 1);
  b_abx ~ normal(0, 1);
  prob[1] = emp_prob[1];  // Set the first value of prob to be the empirical probability
  for (i in 2:n) {
    prob[i] = a_0 + a_specimen[specimen_id2[i]] + b_lag * prob[i-1] + b_abx * abx_yn[i-1];
  }
  read_count ~ binomial_logit(total_reads, prob);
}
