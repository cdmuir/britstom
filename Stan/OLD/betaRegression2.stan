data{
    int<lower=1> N;
    int<lower=1> N_X;
    real Y[N];
    int X[N];
}
parameters{
    real phi[N_X];
    real lambda;
}
model{
    vector[N] beta;
    vector[N] alpha;
    lambda ~ pareto( 0.1 , 1.5 );
    phi ~ beta( 1 , 1 );
    for ( i in 1:N ) {
        beta[i] = lambda * (1 - phi[X[i]]);
    }
    for ( i in 1:N ) {
        alpha[i] = lambda * phi[X[i]];
    }
    Y ~ beta( alpha , beta );
}
generated quantities{
    vector[N] beta;
    vector[N] alpha;
    real dev;
    dev = 0;
    for ( i in 1:N ) {
        beta[i] = lambda * (1 - phi[X[i]]);
    }
    for ( i in 1:N ) {
        alpha[i] = lambda * phi[X[i]];
    }
    dev = dev + (-2)*beta_lpdf( Y | alpha , beta );
}
