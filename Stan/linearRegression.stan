data {
	int<lower=0> N;	// number of data items
    int<lower=0> K;	// number of predictors
	matrix[N, K] x;	// predictor matrix
	vector[N] y;	// outcome vector
}
parameters {
	real alpha;				// intercept
	vector[K] beta;			// coefficients for predictors
	real<lower=0> sigma;	// error scale 
}
model {                      
	vector[N] mu;
	mu = x * beta + alpha;
	y ~ normal(mu, sigma);  // likelihood
}
generated quantities{
    real dev;
    dev = 0;
    dev = dev + (-2) * normal_lpdf( y | x * beta + alpha, sigma );
}