data{
    int<lower=1> n;
    int<lower=1> n_specimen_id2;
    int<lower=0> read_count[n];
    int<lower=0> total_reads[n];
    int<lower=0, upper=1> abx_yn[n];
    int specimen_id2[n];
    // int<lower=0, upper=1> abx_yn[n];
    vector[n] emp_prob;  // the empirical, or observed, probabilities of the genus
}
parameters{
    real a_0;
    real b_lag;
    real b_abx;
    real s_a;
    real<lower=0> s_b;
    vector[n_specimen_id2] a_specimen;
    // real<lower=0> sigma_specimen;
}
transformed parameters {
    real<lower=0> sigma_specimen = s_a ./ sqrt(s_b);
    vector[n] prob;
    prob[1] = a_0 + a_specimen[specimen_id2[1]] + b_lag * emp_prob[1];
    for (i in 2:n) {
      prob[i] = a_0 + a_specimen[specimen_id2[i]] + b_lag * prob[i-1] + b_abx * abx_yn[i-1];
    }
}
model{
    s_a ~ normal(0,1);
    s_b ~ gamma(0.5, 0.5);
    a_specimen ~ normal( 0 , sigma_specimen );
    a_0 ~ normal( 0 , 10 );
    b_lag ~ normal( 0, 1 );
    b_abx ~ normal( 0, 1 );
    read_count ~ binomial_logit(total_reads, prob);
}

