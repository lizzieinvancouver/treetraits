

data {
  int<lower=0> N;
  int<lower=0> n_ind;
  int<lower=0> n_sp;
  int<lower=0> n_site;
  int<lower=1, upper=n_site> site[N];
  int<lower=1, upper=n_sp> sp[N];
  int<lower=1, upper=n_ind> ind[N];
  vector[N] y; // response
  vector[N] lat; // predictor

	int<lower=1> splookup[n_ind]; // will this work if unbalanced?
		
	}

parameters {

  real a_0;
  real b_lat_0;

  real mu_a_sp[n_sp]; 
  real mu_a_sp_ind[n_ind];  
  real mu_b_lat_sp[n_sp]; 

  real<lower=0> sigma_b_lat_sp; // 
  real<lower=0> sigma_b_lat_sp_ind; // 
    
  real<lower=0> sigma_y; // this only gets used in model block
  }


transformed parameters {
	real y_hat[N];
	real a_sp[n_sp]; // intercept for species
  	real b_lat_sp[n_sp]; // 

	real a_sp_ind[n_ind]; // intercept for individuals within species
	
			
	// Species level. Random slopes
	for (k in 1:n_sp) {
		a_sp[k] <- a_0 + mu_a_sp[k]; 
		b_lat_sp[k] <- b_lat_0 + mu_b_lat_sp[k];
				
		}
	
	// individual level	
	for (j in 1:n_ind){
	
		a_sp_ind[j] <- a_sp[splookup[j]] + mu_a_sp_ind[j];

	}
	
	// cutting level
	for(i in 1:N){

		y_hat[i] <- a_sp_ind[ind[i]] + 
					b_lat_sp[sp[i]] * lat[i]
					;		
		}
	
}

model {
	mu_b_lat_sp ~ normal(0, 35); // 100 = 3 months on either side. Narrow down to 35

	sigma_b_lat_sp ~ normal(0, 10); // Start big at 10, go smaller if introduces problems
			
	b_lat_sp ~ normal(mu_b_lat_sp, sigma_b_lat_sp);

	y ~ normal(y_hat, sigma_y);
}