data {
	int<lower=1> K;	// number of mixture components
	int<lower=1> N;	// number of data points
	real y[N];		// observations
}
parameters {
	simplex[K] theta;			// mixing proportions
	real phi[K];				// locations of mixture components
	real<lower=0> lambda[K];	// precisions of mixture components
}
transformed parameters {
	real alpha[K];
	real beta[K];
	for (k in 1:K) {
		alpha[k] = lambda[k] * phi[k];
		beta[k] = lambda[k] * (1 - phi[k]); 
	}
}
model {
	real ps[K];	// temp for log component densities
	lambda ~ pareto( 0.1 , 1.5 );
	phi ~ beta(1, 1);
	for (n in 1:N) {
		for (k in 1:K) {
			ps[k] = log(theta[k]) + beta_lpdf(y[n] | alpha[k], beta[k]);
		}
	target += log_sum_exp(ps);
	}
}