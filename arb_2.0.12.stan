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

  vector[n_abx] zb_abx1;
  real za_abx1;
  real<lower=0> s_abx1;
  
  matrix[n_subjects, n_abx] zb_abx2;
  vector[n_abx] za_abx2;
  vector<lower=0>[n_abx] s_abx2;

}
transformed parameters {
  real phi[n_specimens];
  // real a_spec[n_specimens];
  // real a_subj[n_subjects];

  vector[n_abx] b_abx1;
  matrix[n_subjects, n_abx] b_abx2;
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
      b_abx1[j] = za_abx1 + zb_abx1[j] * s_abx1;
      b_abx2[subject_id, j] = za_abx2[j] + zb_abx2[subject_id, j] * s_abx2[j];
      phi[i] += (
        b_abx1[j] * abx[i, j] +
        b_abx2[subject_id, j] * abx[i, j]
      );
    }
  }
}
model {
  a ~ normal(0, 1);
  b_lag ~ normal(0, 1);
  
  a_subj ~ normal(0, s_subj);
  s_subj ~ normal(0, 5);
  
  a_spec ~ normal(0, s_spec);
  s_spec ~ normal(0, 5);
  
  zb_abx1 ~ normal(0, 1);
  za_abx1 ~ normal(0, 1);
  s_abx1 ~ normal(0, 1);
  
  za_abx2 ~ normal(0, 1);
  s_abx2  ~ normal(0.5, 0.5);

  for (i in 1:n_subjects) {
    zb_abx2[i, ] ~ normal(0, 1);
  }
  
  reads ~ binomial_logit(total, phi);
}

