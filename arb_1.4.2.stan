data{
  int<lower=1> n;
  int<lower=1> n_subject_id;
  int n_obs_subject[n_subject_id];  // observations per subject
  int read_count[n];
  int total_reads[n];
  int specimen_id2[n];
  int subject_id[n];
  vector[n] abx_yn;
}
transformed data {
  vector[n] emp_prob;
  for (i in 1:n) {
    emp_prob[i] = (read_count[i] * 1.0)/total_reads[i];
  }
}
parameters{
  vector[n] a_spec_std;
  vector[n_subject_id] a_subj_std;
  real mu;
  real<lower=0> sigma_spec;
  real<lower=0> sigma_subj;
  real b_lag;
  // real b_abx;
}
model{
  int pos = 1;
  mu ~ normal(-1, 1);
  sigma_spec ~ exponential(1);
  a_spec_std ~ normal(0, 1);
 // mu_subj ~ normal(-1, 1);
  sigma_subj ~ exponential(1);
  a_subj_std ~ normal(0, 1);
  b_lag ~ normal(0, 1);
  
  for (subj in 1:n_subject_id) {
    // Read counts from 2:N, effectively (N being within-subject)
    segment(read_count, pos+1, n_obs_subject[subj]-1) ~ binomial_logit(
      // Total reads from 2:N, effectively
      total_reads[pos+1:(pos+n_obs_subject[subj])-1],
      // Subject- and specimen-level intercepts, specimen sliced as with total reads above
      mu + (sigma_spec * a_spec_std[pos+1:(pos+n_obs_subject[subj])-1] + sigma_subj * a_subj_std[subj]) +
      // Autoregression lag on effectively 1:(N-1)
      b_lag * emp_prob[pos:(pos+n_obs_subject[subj])-2]
    );
    pos = pos + n_obs_subject[subj];
  }
}
generated quantities {
  vector[n] a_spec = mu + sigma_spec * a_spec_std;
  vector[n_subject_id] a_subj = mu + sigma_subj * a_subj_std;
}

