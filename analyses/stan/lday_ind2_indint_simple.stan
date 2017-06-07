// flynn@fas.harvard.edu
// 3 level model for leafout day (lday) as a function of warming and photoperiod treatments, each 2 levels. 
// based off of previous model, pooling just a species level for warming, photoperiod, and chilling effects, and off of "Three-level nested linear model in Stan": 
// http://rstudio-pubs-static.s3.amazonaws.com/64315_bc3a395edd104095a8384db8d9952f43.html
// Difference in this model is that we are not just doing "beta_0" intercept as in the other model, but also slopes. I call intercept "a" and slopes "b".
// Top level: Species (20)
// 2nd level: Individual trees (10 per species).
// 3rd level: Cuttings (4 per individual tree, spread out across different treatments)

data {
  int<lower=1> N;
  int<lower=1> n_ind;
  int<lower=1> n_sp;
  int<lower=1, upper=n_sp> sp[N];
  int<lower=1, upper=n_ind> ind[N]; // this serves also as individual-level cluster id
  vector[N] y; // response, day of budburst or leafout
  vector[N] warm; // predictor
  vector[N] photo; // predictor

	int<lower=1> splookup[n_ind]; // cluster id for species
		
	}

parameters {
  real a_0; // overall intercept, same as "beta_0" in classroom example model
  real b_warm_0; // overall warming effect
  real b_photo_0; // overall photoperiod effect  

  real mu_a_sp[n_sp];
  real mu_b_warm_sp[n_sp]; 
  real mu_b_photo_sp[n_sp];

  real mu_a_sp_ind[n_ind]; 


  real<lower=0> sigma_a_sp; 
  real<lower=0> sigma_b_warm_sp; 
  real<lower=0> sigma_b_photo_sp;

  real<lower=0> sigma_a_sp_ind; 
  
  real<lower=0> sigma_y; 
  }


transformed parameters {
	real y_hat[N];
	real a_sp[n_sp]; // intercept for species
  	real b_warm_sp[n_sp]; // slope of warming effect at species level
	real b_photo_sp[n_sp]; // slope of photoperiod effect, at species level

	real a_sp_ind[n_ind]; // intercept for individuals within species
	
		
	// Species level. Random intercept (a) and slopes for warming and photoperiod
	for (k in 1:n_sp) {
		
		a_sp[k] <- a_0 + mu_a_sp[k]; // 
		b_warm_sp[k] <- b_warm_0 + mu_b_warm_sp[k];
		b_photo_sp[k] <- b_photo_0 + mu_b_photo_sp[k];
				
		}
	
	// individual level. Random intercepts only
	
	for (j in 1:n_ind){
		a_sp_ind[j] <- a_sp[splookup[j]] + mu_a_sp_ind[j];
	}
	
	// cutting (branch) level. 
	for(i in 1:N){

		y_hat[i] <- a_sp_ind[ind[i]] + 
					b_warm_sp[sp[i]] * warm[i] + 
					b_photo_sp[sp[i]] * photo[i]
					;
		}
}

model {
	a_0 ~ normal(35, 100); // possibly reduce sigmas down... to 1 (perhaps according to andrew). With fewer observations, consider stronger prior
	b_warm_0 ~ normal(0, 100); 
	b_photo_0 ~ normal(0, 100);
		
	mu_a_sp ~ normal(0, sigma_a_sp);
	mu_b_warm_sp ~ normal(0, sigma_b_warm_sp); // 
	mu_b_photo_sp ~ normal(0, sigma_b_photo_sp);

	mu_a_sp_ind ~ normal(0, sigma_a_sp_ind);

	sigma_a_sp ~ normal(0, 20); 
	sigma_b_warm_sp ~ normal(0, 20); 
	sigma_b_photo_sp ~ normal(0, 20); 

	sigma_a_sp_ind ~ normal(0, 10); 

	y ~ normal(y_hat, sigma_y);
}