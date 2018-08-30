data {
  int<lower=1> n_timepoints;        // Number of timepoints
  int<lower=1> n_subjects;          // Number of subjects
  int<lower=1> n_abx;               // Number of antibiotics
  int ii_subjects[n_timepoints];    // Index of subject for each timepoint
  matrix[n_timepoints, n_abx] abx;  // Antibiotics on/off for each timepoint/abx
  int reads[n_timepoints];          // Reads of target bacteria
  int total[n_timepoints];          // Total reads 
}
parameters {
  real a;                           // Global intercept
  real a_subj[n_subjects];          // Subject-specific intercepts
  real<lower=0> s_subj;             // Pooled subject variance
  matrix[n_subjects, n_abx] b_abx;  // Slopes for each subject + abx pair
  real<lower=0> s_abx[n_abx];       // Abx-specific variances
}
transformed parameters {
  real phi[n_timepoints];
  for (i in 1:n_timepoints) {
    phi[i] = a + a_subj[ii_subjects[i]];
    for (j in 1:n_abx) {
      phi[i] += b_abx[ii_subjects[i], j] * abx[i, j];
    }
  }
}
model {
  a ~ normal(0, 1);
  a_subj ~ normal(0, s_subj);
  s_subj ~ exponential(1);
  s_abx ~ exponential(1);
  for (i in 1:n_subjects) {
    b_abx[i, ] ~ normal(0, s_abx);
  }
  reads ~ binomial_logit(total, phi);
}
