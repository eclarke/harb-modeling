data{
  int<lower=1> n;
  int<lower=1> n_subject_id;
  int read_count[n];
  vector[n] lag_emp_prop;  // previous day's empirical proportions 
  int total_reads[n];
  int specimen_id2[n];
  int subject_id[n];
  vector[n] on_abx; // whether the subject is under the effects of abx on this timepoint
  vector[n] lag_nonzero; // if the counts in the previous day were zero
}

parameters{
  vector[n] a_spec_std;
  vector[n_subject_id] a_subj_std;
  vector[n_subject_id] b_abx1_std;
  vector[n_subject_id] b_abx2_std;
  real mu;
  real<lower=0> sigma_spec;
  real<lower=0> sigma_subj;
  real<lower=0> sigma_abx1;
  real<lower=0> sigma_abx2;
  real b_lag;
}
transformed parameters{
  vector[n] prob;
  vector[n] a_spec = sigma_spec * a_spec_std;
  vector[n_subject_id] a_subj = sigma_subj * a_subj_std;
  vector[n_subject_id] b_abx1 = sigma_abx1 * b_abx1_std;
  vector[n_subject_id] b_abx2 = sigma_abx2 * b_abx2_std;
  
  for (i in 1:n) {
    prob[i] = (
      // Global, specimen, and subject-level intercepts
      mu + a_spec[specimen_id2[i]] + a_subj[subject_id[i]] +
      // Global lag term
      b_lag * lag_emp_prop[i] + 
      b_abx1[subject_id[i]] * on_abx[i] +
      // Subject-level abx effect (dropped if previous day is zero, or no abx)
      b_abx2[subject_id[i]] * on_abx[i] * lag_nonzero[i]
    );
  }
}
model{
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
  // Normal priors for abx coefficients
  sigma_abx1 ~ cauchy(0, 1);
  b_abx1_std ~ normal(0, 1);
  sigma_abx2 ~ cauchy(0, 1);
  b_abx2_std ~ normal(0, 1);
  // Final model
  read_count ~ binomial_logit(total_reads, prob);
}


