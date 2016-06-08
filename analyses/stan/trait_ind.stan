// Trait analysis at individual level, using similar structure as the phenology model, lday_ind5.stan.

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

	int<lower=1> splookup[n_ind]; 
		
	}

parameters {

  real a_0;
  real b_lat_0;
  real b_site_0;
  
  real mu_a_sp[n_sp]; 
  real mu_a_sp_ind[n_ind];  
  real mu_b_lat_sp[n_sp]; 
  real mu_b_site_sp[n_sp]; 
  
  real<lower=0> sigma_a_sp;
  real<lower=0> sigma_a_sp_ind;
  real<lower=0> sigma_b_lat_sp; 
  real<lower=0> sigma_b_site_sp; 
    
  real<lower=0> sigma_y; // 
  }


transformed parameters {
	real y_hat[N];
	real a_sp[n_sp]; // intercept for species
  	real b_lat_sp[n_sp]; // 
  	real b_site_sp[n_sp]; //
  	
	real a_sp_ind[n_ind]; // intercept for individuals within species
	
			
	// Species level. Random slopes
	for (k in 1:n_sp) {
		a_sp[k] <- a_0 + mu_a_sp[k]; 
		b_lat_sp[k] <- b_lat_0 + mu_b_lat_sp[k];
		b_site_sp[k] <- b_site_0 + mu_b_site_sp[k];

				
		}
	
	// individual level	
	for (j in 1:n_ind){
	
		a_sp_ind[j] <- a_sp[splookup[j]] + mu_a_sp_ind[j];

	}
	
	// cutting level
	for(i in 1:N){

		y_hat[i] <- a_sp_ind[ind[i]] + 
					b_lat_sp[sp[i]] * lat[i] +
					b_site_sp[sp[i]] * site[i] 
					;		
		}
	
}

model {
	a_0 ~ normal(0, 100); // 
	b_lat_0 ~ normal(0, 100); // 
	b_site_0 ~ normal(0, 100); // 

	mu_a_sp ~ normal(0, 10); // 
	mu_a_sp_ind ~ normal(0, 10); // 
	mu_b_lat_sp ~ normal(0, 35); // 
	mu_b_site_sp ~ normal(0, 35); // 

	sigma_a_sp ~ normal(0, 10); //
	sigma_a_sp_ind ~ normal(0, 10); //
	sigma_b_lat_sp ~ normal(0, 10); //
	sigma_b_site_sp ~ normal(0, 10); //
	
	y ~ normal(y_hat, sigma_y);
	
}