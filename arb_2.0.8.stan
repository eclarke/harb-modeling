data {
  int<lower=1> n_specimens;         // Number of timepoints
  int<lower=1> n_subjects;          // Number of subjects
  int<lower=1> n_abx;               // Number of antibiotics
  int subjects[n_specimens];        // Index of subject for each timepoint
  matrix[n_specimens, n_abx] abx;   // Antibiotics on/off for each timepoint/abx
  int reads[n_specimens];           // Reads of target bacteria
  int total[n_specimens];           // Total reads 
  real prev[n_specimens];           // Previous proportion
}

parameters {
  real a;     // Global intercept
  real b_lag; // Global lag coefficient

  // Non-centered parameterization for a_spec
  // real za_spec[n_specimens];  // Standardized intercepts for each specimen
  real<lower=0> s_spec;       // Pooled specimen variance
  
  // Non-centered parameterization for a_subj
  // real za_subj[n_subjects]; // Subject-specific intercepts
  real<lower=0> s_subj;     // Pooled subject variance

  real a_spec[n_specimens];
  real a_subj[n_subjects];

  matrix[n_subjects, n_abx] b_abx;
  matrix[n_subjects, n_abx] b_abxp;

  // Non-centered parameterization for b_abx and b_abxp
  // matrix[n_subjects, n_abx] zb_abx;   // Non-centered slopes for each subject + abx
  // matrix[n_subjects, n_abx] zb_abxp;  // NC slopes for each subject + (abx * prev prop)
  row_vector<lower=0>[n_abx] s_abx;   // Abx-specific variances for NCP
  row_vector<lower=0>[n_abx] s_abxp;  // Abx-specific variances for b_abxp term in NCP
  // row_vector[n_abx] za_abx;           // Abx-specific intercepts for NCP
  // row_vector[n_abx] za_abxp;          // Abx-specific intercepts for b_abxp in NCP
}
transformed parameters {
  real phi[n_specimens];
  // real a_spec[n_specimens];
  // real a_subj[n_subjects];

  // matrix[n_subjects, n_abx] b_abxp;
  for (i in 1:n_specimens) {
    int subject_id = subjects[i];
    // NCP for a_spec and a_subj
    // a_spec[i] = za_spec[i] * s_spec;
    // a_subj[subject_id] = za_subj[subject_id] * s_subj;
    
    phi[i] = a + a_subj[subject_id] + b_lag * prev[i] + a_spec[i];
    // Add each antibiotic-specific coefficient separately
    for (j in 1:n_abx) {
      // NCP for b_abx and b_abxp:
      // b_abxp[subject_id, j] = za_abxp[j] + zb_abxp[subject_id, j] * s_abxp[j];
      // b_abx[subject_id, j] = za_abx[j] + zb_abx[subject_id, j] * s_abx[j];
      phi[i] += (
        b_abxp[subject_id, j] * abx[i, j] * prev[i] +
        b_abx[subject_id, j] * abx[i, j]
      );
    }
  }
}
model {
  a ~ normal(0, 1);
  b_lag ~ normal(0, 1);
  
  a_subj ~ normal(0, s_subj);
  // za_subj ~ normal(0, 5);
  s_subj ~ cauchy(0, 5);
  
  a_spec ~ normal(0, s_spec);
  // za_spec ~ normal(0, 5);
  s_spec ~ cauchy(0, 5);
  
  // za_abx  ~ normal(0, 2);
  // za_abxp ~ normal(0, 2);
  s_abx  ~ cauchy(0, 2);
  s_abxp ~ cauchy(0, 2);

  for (i in 1:n_subjects) {
    b_abx[i, ] ~ normal(0, s_abx);
    b_abxp[i,] ~ normal(0, s_abxp);
  }
  
  reads ~ binomial_logit(total, phi);
}
