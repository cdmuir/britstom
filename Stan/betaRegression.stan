data {
	int<lower=0> N;	// number of data items
    int<lower=0> K;	// number of predictors
	matrix[N, K] x;	// predictor matrix
	vector[N] y;	// outcome vector
}
parameters {
	real a;					// intercept
	vector[K] b;			// coefficients for predictors
	real<lower=0> lambda;	// precision 
}
model {                      
	vector[N] phi;
	vector[N] alpha;
	vector[N] beta;
	lambda ~ pareto( 0.1 , 1.5 );
	a ~ normal(0, 10);
	b ~ normal(0, 10);
	phi = inv_logit(x * b + a);
	// a ~ normal(0.25, 0.01);
	// b ~ normal(0.5, 0.01);
	// phi = x * b + a;
	alpha = lambda * phi;
	beta = lambda * (1 - phi); 
	y ~ beta(alpha, beta);  // likelihood
}
generated quantities{
	vector[N] phi;
	vector[N] alpha;
	vector[N] beta;
    real dev;
	phi = inv_logit(x * b + a);
	alpha = lambda * phi;
	beta = lambda * (1 - phi); 
    dev = 0;
    dev = dev + (-2) * beta_lpdf( y | alpha, beta );
}
