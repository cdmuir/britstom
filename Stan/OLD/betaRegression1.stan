data{
    int<lower=1> N;
    real X[N];
}
parameters{
    real phi;
    real lambda;
}
model{
    vector[N] beta;
    vector[N] alpha;
    lambda ~ pareto( 0.1 , 1.5 );
    phi ~ beta( 1 , 1 );
    for ( i in 1:N ) {
        beta[i] = lambda * (1 - phi);
    }
    for ( i in 1:N ) {
        alpha[i] = lambda * phi;
    }
    X ~ beta( alpha , beta );
}
generated quantities{
    vector[N] beta;
    vector[N] alpha;
    real dev;
    dev = 0;
    for ( i in 1:N ) {
        beta[i] = lambda * (1 - phi);
    }
    for ( i in 1:N ) {
        alpha[i] = lambda * phi;
    }
    dev = dev + (-2)*beta_lpdf( X | alpha , beta );
}
