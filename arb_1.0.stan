data{
    int<lower=1> N;
    int<lower=1> N_specimen_id2;
    int<lower=1> N_subject_id;
    int read_count[N];
    int total_reads[N];
    int specimen_id2[N];
    int subject_id[N];
}
parameters{
    real a_0;
    vector[N_specimen_id2] a_specimen;
    vector[N_subject_id] a_subject;
    real<lower=0> sigma_specimen;
    real<lower=0> sigma_subject;
}
model{
    vector[N] prob;
    sigma_subject ~ cauchy( 0 , 1 );
    sigma_specimen ~ cauchy( 0 , 1 );
    a_subject ~ normal( 0 , sigma_subject );
    a_specimen ~ normal( 0 , sigma_specimen );
    a_0 ~ normal( 0 , 10 );
    for ( i in 1:N ) {
        prob[i] = a_0 + a_specimen[specimen_id2[i]] + a_subject[subject_id[i]];
    }
    read_count ~ binomial_logit( total_reads , prob );
}
generated quantities{
    vector[N] prob;
    real dev;
    dev = 0;
    for ( i in 1:N ) {
        prob[i] = a_0 + a_specimen[specimen_id2[i]] + a_subject[subject_id[i]];
    }
    dev = dev + (-2)*binomial_logit_lpmf( read_count | total_reads , prob );
}