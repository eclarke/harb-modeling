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
    real<lower=0> sigma_specimen;
}
model{
    vector[n] prob;
    sigma_specimen ~ cauchy( 0 , 1 );
    a_specimen ~ normal( 0 , sigma_specimen );
    a_0 ~ normal( 0 , 10 );
    b_lag ~ normal( 0, 1 );
    prob[1] = emp_prob[1];  // Set the first value of prob to be the empirical probability
    read_count[1] ~ binomial_logit(total_reads[1], prob[1]);
    for (i in 2:n) {
      prob[i] = a_0 + a_specimen[specimen_id2[i]] + b_lag * prob[i-1] + b_abx * abx_yn[i-1];     
      read_count[i] ~ binomial_logit(total_reads[i], prob[i]);
    }
}

