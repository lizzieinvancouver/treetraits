
// Traits by site, species, and individual
// Starting more simply, with just site and species

data {
  int<lower=0> N;
  int<lower=0> n_ind;
  int<lower=0> n_sp;
  int<lower=0> n_site;
  int<lower=1, upper=n_site> site[N];
  int<lower=1, upper=n_sp> sp[N];
  int<lower=1, upper=n_ind> ind[N];
  vector[N] sla;
 
}

parameters {
  vector[n_site] a_site;
  vector[n_sp] a_sp;
  vector[n_sp] b_warm;
 
  real mu_b_warm; 

  real<lower=0> sigma_b_warm;
  
  real<lower=0> sigma_y; 
  }


transformed parameters {
	vector[N] y_hat;
		
	for(i in 1:N){

		y_hat[i] <- a_site[site[i]] + a_sp[sp[i]] + b_warm[sp[i]] * warm[i];
		
		}
	
}

model {
	// Priors. Make them flat
	mu_b_warm ~ normal(0, 35); // 100 = 3 months on either side. Narrow down to 35

	sigma_b_warm ~ normal(0, 10); // Start big at 10, go smaller if introduces problems
	
	b_warm ~ normal(mu_b_warm, sigma_b_warm);

	lday ~ normal(y_hat, sigma_y);
}