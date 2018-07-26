data{
  int<lower=1> n;
  int<lower=1> n_subject_id;
  int read_count[n];
  vector[n] lag_emp_prop;  // previous day's empirical proportions 
  int total_reads[n];
  int specimen_id2[n];
  int subject_id[n];
  vector[n] on_abx; // whether the subject is under the effects of abx on this timepoint
}
parameters{
  vector[n] a_spec_std;
  vector[n_subject_id] a_subj_std;
  vector[n_subject_id] b_abx_std;
  real mu;
  real<lower=0> sigma_spec;
  real<lower=0> sigma_subj;
  real<lower=0> sigma_abx;
  real b_lag;
}
// transformed parameters {
//   real prob[n];
//   real a_spec[n];
//   real a_subj[n_subject_id];
//   real b_abx[n_subject_id];
//   for (i in 1:n) {
//     a_spec[i] = sigma_spec * a_spec_std[specimen_id2[i]];
//     a_subj[subject_id[i]] = sigma_subj * a_subj_std[subject_id[i]];
//     b_abx[subject_id[i]] = sigma_abx * b_abx_std[subject_id[i]];
//     prob[i] = mu + a_spec[specimen_id2[i]] + a_subj[subject_id[i]] + 
//       b_lag * lag_emp_prop[i] + b_abx[subject_id[i]] * on_abx[i];
//   }
// }
transformed parameters{
  vector[n] prob;
  for (i in 1:n) {
    prob[i] = mu + sigma_spec * a_spec_std[specimen_id2[i]] + sigma_subj * a_subj_std[subject_id[i]] +
    b_lag * lag_emp_prop[i] + sigma_abx * b_abx_std[subject_id[i]] * on_abx[i];
  }
}
model{
  // vector[n] prob;
  // Global intercept term
  mu ~ normal(-1, 1);
  // Sigma and unit normal intercepts for specimen
  sigma_spec ~ cauchy(0, 1);
  a_spec_std ~ normal(0, 1);
  // Sigma and unit normal intercepts for subject
  sigma_subj ~ cauchy(0, 1);
  a_subj_std ~ normal(0, 1);
  // Normal prior for lag coefficient
  b_lag ~ normal(0, 1);
  // Normal priors for abx coefficient
  sigma_abx ~ cauchy(0, 1);
  b_abx_std ~ normal(0, 1);
  read_count ~ binomial_logit(total_reads, prob);
}
generated quantities {
  vector[n] a_spec = sigma_spec * a_spec_std;
  vector[n_subject_id] a_subj = sigma_subj * a_subj_std;
  vector[n_subject_id] b_abx = sigma_abx * b_abx_std;
}

