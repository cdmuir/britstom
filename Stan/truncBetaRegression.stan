data {
	int<lower=0> N_obs;			// number of observed data items
	int<lower=0> N_censL;		// number of data items below lower censor
	int<lower=0> N_censU;		// number of data items above upper censor
	vector[N_obs] y_obs;		// outcome vector
	real<upper=min(y_obs)> L;	// lower censor
	real<lower=max(y_obs)> U;	// upper censor
}
parameters {
	real phi;				// mean
	real<lower=0> lambda;	// precision 
}
model {                      
	real alpha;
	real beta;
	lambda ~ pareto( 0.1 , 1.5 );
	phi ~ beta( 1 , 1 );
	alpha = lambda * phi;
	beta = lambda * (1 - phi); 
	y_obs ~ beta(alpha, beta);
	target += N_censL * beta_lcdf(L | alpha, beta);
	target += N_censU * beta_lccdf(U | alpha, beta);
}
generated quantities{
	real alpha;
	real beta;
    real dev;
	alpha = lambda * phi;
	beta = lambda * (1 - phi); 
    dev = 0;
    dev = dev + (-2) * beta_lpdf( y_obs | alpha, beta );
}
