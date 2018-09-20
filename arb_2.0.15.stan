/* 
* Autoregressive binomial model, v2.0.15
* Changes from v2.0.13:
*  - negate b_abx1 and b_abx2 when multiplied by prev day abundance
* 
* reads ~ Binomial(total, logit(phi))
* phi[i,j,k] = a + a_spec[i] + a_subj[j] - ((b_abx1[k] * abx[i,k]) + (b_abx2[j,k] * abx[i,k])) * prev[i] + (b_lag * prev[i])
* 
* a = global intercept
* a_spec[i] = specimen i intercept (partially pooled between specimens)
* a_subj[j] = subject j intercept (partially pooled btwn subjects)
* b_abx1[k] = antibiotic k slope (partially pooled btwn abx)
* abx[i,k] = indicator variable for antibiotic k at specimen i
* b_abx2[j,k] = subject j + antibiotic k slope (partially pooled between subjects)
* b_lag = coefficient for previous specimen's reads
* prev[i] = previous specimen's reads (as proportion of total)
*
*/
data {
  int<lower=1> n_specimens;         // Number of timepoints
  int<lower=1> n_subjects;          // Number of subjects
  int<lower=1> n_abx;               // Number of antibiotics
  int subjects[n_specimens];        // Index of subject for each specimen
  matrix[n_specimens, n_abx] abx;   // Antibiotics on/off for each specimen/abx
  int reads[n_specimens];           // Reads of target bacteria
  int total[n_specimens];           // Total reads 
  real prev[n_specimens];           // Previous proportion
}

parameters {
  real a;     // Global intercept
  real b_lag; // Global lag coefficient

  real a_spec[n_specimens]; // Specimen intercepts
  real<lower=0> s_spec;     // Pooled specimen variance
  
  real a_subj[n_subjects];  // Subject intercepts
  real<lower=0> s_subj;     // Pooled subject variance

  // Non-centered parameterization for b_abx1
  // b_abx1 = za_abx1 + zb_abx1 * s_abx1
  vector[n_abx] zb_abx1;
  real za_abx1;
  real<lower=0> s_abx1;
  
  // Non-centered parameterization for b_abx2
  // b_abx2 = za_abx2 + zb_abx2 * s_abx2
  matrix[n_subjects, n_abx] zb_abx2;
  vector[n_abx] za_abx2;
  vector<lower=0>[n_abx] s_abx2;
}
transformed parameters {
  real phi[n_specimens];

  vector[n_abx] b_abx1;
  matrix[n_subjects, n_abx] b_abx2;

  for (i in 1:n_specimens) {
    int subject_id = subjects[i];
    phi[i] = a + a_spec[i] + a_subj[subject_id] + b_lag * prev[i];
    // Add each antibiotic-specific coefficient separately
    for (j in 1:n_abx) {
      b_abx1[j] = za_abx1 + zb_abx1[j] * s_abx1;
      b_abx2[subject_id, j] = za_abx2[j] + zb_abx2[subject_id, j] * s_abx2[j];
      phi[i] -= (
        b_abx1[j] * abx[i, j] * prev[i] +
        b_abx2[subject_id, j] * abx[i, j] * prev[i]
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

