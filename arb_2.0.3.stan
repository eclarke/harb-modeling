data {
  int<lower=1> n_timepoints;        // Number of timepoints
  int<lower=1> n_subjects;          // Number of subjects
  int<lower=1> n_abx;               // Number of antibiotics
  int ii_subjects[n_timepoints];    // Index of subject for each timepoint
  matrix[n_timepoints, n_abx] abx;  // Antibiotics on/off for each timepoint/abx
  int reads[n_timepoints];          // Reads of target bacteria
  int total[n_timepoints];          // Total reads 
  real prev[n_timepoints];          // Previous proportion
}
parameters {
  real a;                               // Global intercept
  real a_subj[n_subjects];              // Subject-specific intercepts
  real b_lag;                           // Global lag coefficient
  real<lower=0> s_subj;                 // Pooled subject variance

  matrix[n_subjects, n_abx] b_abx_nc;   // Non-centered slopes for each subject + abx
  matrix[n_subjects, n_abx] b_abxp_nc;  // NC slopes for each subject + (abx * prev prop)
  row_vector<lower=0>[n_abx] s_abx_nc;  // Abx-specific variances for NCP
  row_vector<lower=0>[n_abx] s_abxp_nc; // Abx-specific variances for b_abxp term in NCP
  row_vector[n_abx] a_abx_nc;           // Abx-specific intercepts for NCP
  row_vector[n_abx] a_abxp_nc;          // Abx-specific intercepts for b_abxp in NCP
}
transformed parameters {
  real phi[n_timepoints];
  matrix[n_subjects, n_abx] b_abx;
  matrix[n_subjects, n_abx] b_abxp;
  for (i in 1:n_timepoints) {
    phi[i] = a + a_subj[ii_subjects[i]] + b_lag * prev[i];
    // Add each antibiotic-specific coefficient separately
    for (j in 1:n_abx) {
      // Non-centered parameterization for b_abx and b_abxp:
      b_abxp[ii_subjects[i], j] = a_abxp_nc[j] + b_abxp_nc[ii_subjects[i], j] * s_abxp_nc[j];
      b_abx[ii_subjects[i], j] = a_abx_nc[j] + b_abx_nc[ii_subjects[i], j] * s_abx_nc[j];
      phi[i] += (
        b_abxp[ii_subjects[i], j] * abx[i, j] * prev[i] +
        b_abx[ii_subjects[i], j] * abx[i, j]
      );
    }
  }
}
model {
  a       ~ normal(0, 1);
  a_subj  ~ normal(0, s_subj);
  b_lag   ~ normal(0, 1);
  s_subj  ~ exponential(1);
  
  a_abx_nc  ~ normal(0, 1);
  a_abxp_nc ~ normal(0, 1);
  s_abx_nc  ~ exponential(1);
  s_abxp_nc ~ exponential(1);
  for (i in 1:n_subjects) {
    b_abx_nc[i, ] ~ normal(0, 1);
    b_abxp_nc[i,] ~ normal(0, 1);
  }
  
  reads ~ binomial_logit(total, phi);
}
