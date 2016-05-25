// flynn@fas.harvard.edu
// 3 level model for leafout day (lday) as a function for now of just warming and photoperiod. 
// based off of previous model, pooling just a species level for warming, photoperiod, and chilling effects, and off of "Three-level nested linear model in Stan": 
// http://rstudio-pubs-static.s3.amazonaws.com/64315_bc3a395edd104095a8384db8d9952f43.html
// Difference in this model is that we are not just doing "beta_0" intercept as in the other model, but also slopes. I call intercept "a" and slopes "b".
// Top level: Species
// 2nd level: Individuals (6 per species per site, which is currently being ignored
// 3rd level: Cuttings (4-12 per individual, spread out across different treatments)


data {
  int<lower=0> N;
  int<lower=0> n_ind;
  int<lower=0> n_sp;
  int<lower=1, upper=n_sp> sp[N];
  int<lower=1, upper=n_ind> ind[N];
  vector[N] lday; // response
  vector[N] warm; // predictor
  vector[N] photo; // predictor

	int<lower=1> splookup[N]; // cluster id for species
	int<lower=1> indlookup[N]; // cluster id for individual
		
	}

parameters {
  real a_0; // overall intercept, same as "beta_0" in classroom example model
  real b_warm_0; // overall warming effect
  real b_photo_0; // overall photoperiod effect  

  real mu_a_sp[n_sp];
  real mu_b_warm_sp[n_sp]; 
  real mu_b_photo_sp[n_sp];

  real mu_a_sp_ind[n_ind]; 
  real mu_b_warm_sp_ind[n_ind]; 
  real mu_b_photo_sp_ind[n_ind];

  real<lower=0> sigma_b_warm_sp; 
  real<lower=0> sigma_b_photo_sp;

  real<lower=0> sigma_b_warm_sp_ind; 
  real<lower=0> sigma_b_photo_sp_ind;
  
  real<lower=0> sigma_y; // this only gets used in model block
  }


transformed parameters {
	real y_hat[N];
	real a_sp[n_sp]; // intercept for species
  	real b_warm_sp[n_sp]; // slope of warming effect at species level
	real b_photo_sp[n_sp]; // slope of photoperiod effect, at species level

	real a_sp_ind[n_sp]; // intercept for individuals within species
	real b_warm_sp_ind[n_sp]; // 
	real b_photo_sp_ind[n_sp]; // 
	
		
	// Species level. Random intercept (a) and slopes for warming and photoperiod
	for (k in 1:n_sp) {
		
		a_sp[k] <- a_0 + mu_a_sp[k];
				
		}
	
	// individual level. again both random intercepts and slopes
	
	for (j in 1:n_ind){
		a_sp_ind[j] <- a_sp[splookup[j]] + mu_a_sp_ind[j];
	
	}
	
	// cutting level. note that "warm" and "photo" are vectors of data, not parameters.
	for(i in 1:N){
		y_hat[i] <- a_sp_ind[indlookup[i]];
		}
	
}

model {

	mu_a_sp ~ normal(0, 35);
//	mu_b_warm_sp ~ normal(0, 35); // 100 = 3 months on either side. Narrow down to 35
//	mu_b_photo_sp ~ normal(0, 35);

	mu_a_sp_ind ~ normal(0,10)
//	mu_b_warm_sp_ind ~ normal(0, 10); // 10 d on either side at individual level
//	mu_b_photo_sp_ind ~ normal(0, 10);


//	sigma_b_warm_sp ~ normal(0, 10); // Start big at 10, go smaller if introduces problems
//	sigma_b_photo_sp ~ normal(0, 10); 
	
//	sigma_b_warm_ind ~ normal(0, 10); // Reduce sd of sigma at individual level? 
//	sigma_b_photo_ind ~ normal(0, 10); 

	
//	b_warm_sp ~ normal(mu_b_warm_sp, sigma_b_warm_sp);
//	b_photo_sp ~ normal(mu_b_photo_sp, sigma_b_photo_sp);

//	b_warm_sp_ind ~ normal(mu_b_warm_sp_ind, sigma_b_warm_sp_ind);
//	b_photo_sp_ind ~ normal(mu_b_photo_sp_ind, sigma_b_photo_sp_ind)

	lday ~ normal(y_hat, sigma_y);
}