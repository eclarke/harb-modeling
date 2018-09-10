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

  real a_spec[n_specimens];
  real<lower=0> s_spec;       
  
  real a_subj[n_subjects];
  real<lower=0> s_subj;     
  
  real b_abx[n_abx];
  // real<lower=0> s_abx[n_abx];
  
  // matrix[n_subjects, n_abx] b_sabx;
  // real<lower=0> s_sabx[n_abx];

}
transformed parameters {
  real phi[n_specimens];
  for (i in 1:n_specimens) {
    int subject_id = subjects[i];
    phi[i] = a + a_subj[subject_id] + b_lag * prev[i] + a_spec[i];
    // Add each antibiotic-specific coefficient separately
    for (j in 1:n_abx) {
      phi[i] += (
        b_abx[j] * abx[i, j] // +
        // b_sabx[subject_id, j] * abx[i, j]
      );
    }
  }
}
model {
  a ~ normal(0, 1);
  b_lag ~ normal(0, 1);
  
  a_subj ~ normal(0, s_subj);
  s_subj ~ cauchy(0, 5);
  
  a_spec ~ normal(0, s_spec);
  s_spec ~ cauchy(0, 5);
  
  b_abx ~ normal(0, 5);
  // s_abx ~ normal(1, 1);
  
  // s_sabx  ~ normal(1, 1);
  // 
  // for (i in 1:n_subjects) {
  //   b_sabx[i, ] ~ normal(0, s_sabx);
  // }
  // 
  reads ~ binomial_logit(total, phi);
}
