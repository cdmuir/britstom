data{
    int<lower=1> N;
    real Y[N];
    real X1[N];
    real X2[N];
}
parameters{
    real a;
    real b1;
    real b2;
    real lambda;
}
model{
    vector[N] beta;
    vector[N] alpha;
    vector[N] phi;
    lambda ~ pareto( 0.1 , 1.5 );
    a ~ normal(0, 10);
    b1 ~ normal(0, 10);
    b2 ~ normal(0, 10);
    for ( i in 1:N ) {
        phi[i] = a + b1 * X1[i] + b2 * X2[i];
    }
    for ( i in 1:N ) {
        beta[i] = lambda * (1 - phi[i]);
    }
    for ( i in 1:N ) {
        alpha[i] = lambda * (phi[i]);
    }
    Y ~ beta( alpha , beta );
}
generated quantities{
    vector[N] beta;
    vector[N] alpha;
    vector[N] phi;
    real dev;
    dev = 0;
    for ( i in 1:N ) {
        phi[i] = a + b1 * X1[i] + b2 * X2[i];
    }
    for ( i in 1:N ) {
        beta[i] = lambda * (1 - phi[i]);
    }
    for ( i in 1:N ) {
        alpha[i] = lambda * (phi[i]);
    }
    dev = dev + (-2)*beta_lpdf( Y | alpha , beta );
}
