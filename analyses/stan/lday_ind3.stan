// flynn@fas.harvard.edu
// 3 level model for leafout day (lday) as a function for now of just warming and photoperiod. 
// based off of previous model, pooling just a species level for warming, photoperiod, and chilling effects, and off of "Three-level nested linear model in Stan": 
// http://rstudio-pubs-static.s3.amazonaws.com/64315_bc3a395edd104095a8384db8d9952f43.html
// Difference in this model is that we are not just doing "beta_0" intercept as in the other model, but also slopes. I call intercept "a" and slopes "b".
// Top level: Species
// 2nd level: Individuals (6 per species per site, which is currently being ignored
// 3rd level: Cuttings (4-12 per individual, spread out across different treatments)
// adding site and interactions between warming, photoperiod, and site back in

data {
  int<lower=1> N;
  int<lower=1> n_ind;
  int<lower=1> n_sp;
  int<lower=1, upper=n_sp> sp[N];
  int<lower=1, upper=n_ind> ind[N]; // this serves also as individual-level cluster id
  vector[N] lday; // response
  vector[N] warm; // predictor
  vector[N] photo; // predictor
  vector[N] site; // predictor

	int<lower=1> splookup[n_ind]; // cluster id for species
		
	}

transformed data { 			// 3 interaction terms (no chilling yet)
  vector[N] inter_wp;           
  vector[N] inter_ws;           
  vector[N] inter_ps;           

  inter_wp    <- warm .* photo;  
  inter_ws    <- warm .* site;  
  inter_ps    <- photo .* site;  

}

parameters {
  real a_0; // overall intercept, same as "beta_0" in classroom example model
  real b_warm_0; // overall warming effect
  real b_photo_0; // overall photoperiod effect  
  real b_site_0; // overall site effect  
  real b_inter_wp_0; // overall warming x photoperiod effect  
  real b_inter_ws_0; // overall warming x site effect  
  real b_inter_ps_0; // overall photoperiod x site effect  

  real mu_a_sp[n_sp];
  real mu_b_warm_sp[n_sp]; 
  real mu_b_photo_sp[n_sp];
  real mu_b_site_sp[n_sp];
  real mu_b_inter_wp_sp[n_sp];
  real mu_b_inter_ws_sp[n_sp];
  real mu_b_inter_ps_sp[n_sp];
      
  real mu_a_sp_ind[n_ind]; 
  real mu_b_warm_sp_ind[n_ind]; 
  real mu_b_photo_sp_ind[n_ind];
  real mu_b_site_sp_ind[n_ind];
  real mu_b_inter_wp_sp_ind[n_ind];
  real mu_b_inter_ws_sp_ind[n_ind];
  real mu_b_inter_ps_sp_ind[n_ind];
  
  real<lower=0> sigma_a_sp; 
  real<lower=0> sigma_b_warm_sp; 
  real<lower=0> sigma_b_photo_sp;
  real<lower=0> sigma_b_site_sp;
  real<lower=0> sigma_b_inter_wp_sp;
  real<lower=0> sigma_b_inter_ws_sp;
  real<lower=0> sigma_b_inter_ps_sp;

  real<lower=0> sigma_a_sp_ind; 
  real<lower=0> sigma_b_warm_sp_ind; 
  real<lower=0> sigma_b_photo_sp_ind;
  real<lower=0> sigma_b_site_sp_ind;
  real<lower=0> sigma_b_inter_wp_sp_ind;
  real<lower=0> sigma_b_inter_ws_sp_ind;
  real<lower=0> sigma_b_inter_ps_sp_ind;
  
  real<lower=0> sigma_y; 
  }


transformed parameters {
	real y_hat[N];
	
	real a_sp[n_sp]; // intercept for species
  	real b_warm_sp[n_sp]; // slope of warming effect at species level
	real b_photo_sp[n_sp]; // slope of photoperiod effect, at species level
	real b_site_sp[n_sp]; // slope of photoperiod effect, at species level
	real b_inter_wp_sp[n_sp]; // slope of warming x photoperiod effect, at species level
	real b_inter_ws_sp[n_sp]; // slope of warming x photoperiod effect, at species level
	real b_inter_ps_sp[n_sp]; // slope of warming x photoperiod effect, at species level
		
	real a_sp_ind[n_ind]; // intercept for individuals within species
	real b_warm_sp_ind[n_ind]; // 
	real b_photo_sp_ind[n_ind]; // 
	real b_site_sp_ind[n_ind]; 
	real b_inter_wp_sp_ind[n_ind]; 
	real b_inter_ws_sp_ind[n_ind]; 
	real b_inter_ps_sp_ind[n_ind]; 
		
	// Species level. Random intercept (a) and slopes for warming and photoperiod
	for (k in 1:n_sp) {
		
		a_sp[k] <- a_0 + mu_a_sp[k];
		b_warm_sp[k] <- b_warm_0 + mu_b_warm_sp[k];
		b_photo_sp[k] <- b_photo_0 + mu_b_photo_sp[k];
		b_site_sp[k] <- b_site_0 + mu_b_site_sp[k];
		b_inter_wp_sp[k] <- b_inter_wp_0 + mu_b_inter_wp_sp[k];
		b_inter_ws_sp[k] <- b_inter_ws_0 + mu_b_inter_ws_sp[k];
		b_inter_ps_sp[k] <- b_inter_ps_0 + mu_b_inter_ps_sp[k];										
		}
	
	// individual level. again both random intercepts and slopes
	
	for (j in 1:n_ind){
		a_sp_ind[j] <- a_sp[splookup[j]] + mu_a_sp_ind[j];
		b_warm_sp_ind[j] <- b_warm_sp[splookup[j]] + mu_b_warm_sp_ind[j];
		b_photo_sp_ind[j] <- b_photo_sp[splookup[j]] + mu_b_photo_sp_ind[j];
		b_site_sp_ind[j] <- b_site_sp[splookup[j]] + mu_b_site_sp_ind[j];
		b_inter_wp_sp_ind[j] <- b_inter_wp_sp[splookup[j]] + mu_b_inter_wp_sp_ind[j];
		b_inter_ws_sp_ind[j] <- b_inter_ws_sp[splookup[j]] + mu_b_inter_ws_sp_ind[j];
		b_inter_ps_sp_ind[j] <- b_inter_ps_sp[splookup[j]] + mu_b_inter_ps_sp_ind[j];										

	}
	
	// cutting level. note that "warm" and "photo" are vectors of data, not parameters.
	for(i in 1:N){

		y_hat[i] <- a_sp_ind[ind[i]] + 
					b_warm_sp_ind[ind[i]] * warm[i] + 
					b_photo_sp_ind[ind[i]] * photo[i]+
					b_site_sp_ind[ind[i]] * site[i]+
					b_inter_wp_sp_ind[ind[i]] * inter_wp[i]+
					b_inter_ws_sp_ind[ind[i]] * inter_ws[i]+
					b_inter_ps_sp_ind[ind[i]] * inter_ps[i]
					;
		
		}
	
}

model {

	mu_a_sp ~ normal(0, sigma_a_sp); // 20 d on either site at sp level
	mu_b_warm_sp ~ normal(0, sigma_b_warm_sp); 
	mu_b_photo_sp ~ normal(0, sigma_b_photo_sp);
	mu_b_site_sp ~ normal(0, sigma_b_site_sp);
	mu_b_inter_wp_sp ~ normal(0, sigma_b_inter_wp_sp);
	mu_b_inter_ws_sp ~ normal(0, sigma_b_inter_ws_sp);
	mu_b_inter_ps_sp ~ normal(0, sigma_b_inter_ps_sp);
	
	mu_a_sp_ind ~ normal(0, sigma_a_sp_ind); // 10 d on either side at individual level
	mu_b_warm_sp_ind ~ normal(0, sigma_b_warm_sp_ind); 
	mu_b_photo_sp_ind ~ normal(0, sigma_b_photo_sp_ind);
	mu_b_site_sp_ind ~ normal(0, sigma_b_site_sp_ind);
	mu_b_inter_wp_sp_ind ~ normal(0, sigma_b_inter_wp_sp_ind);
	mu_b_inter_ws_sp_ind ~ normal(0, sigma_b_inter_ws_sp_ind);
	mu_b_inter_ps_sp_ind ~ normal(0, sigma_b_inter_ps_sp_ind);
	
	sigma_a_sp ~ normal(0, 20); // Start big at 10, go smaller if introduces problems
	sigma_b_warm_sp ~ normal(0, 20); 
	sigma_b_photo_sp ~ normal(0, 20); 
	sigma_b_site_sp ~ normal(0, 20); 
	sigma_b_inter_wp_sp ~ normal(0, 20); 
	sigma_b_inter_ws_sp ~ normal(0, 20); 
	sigma_b_inter_ps_sp ~ normal(0, 20); 
		
	sigma_a_sp_ind ~ normal(0, 10); // Reduce sd of sigma further at individual level? 
	sigma_b_warm_sp_ind ~ normal(0, 10); 
	sigma_b_photo_sp_ind ~ normal(0, 10); 
	sigma_b_site_sp_ind ~ normal(0, 20); 
	sigma_b_inter_wp_sp_ind ~ normal(0, 20); 
	sigma_b_inter_ws_sp_ind ~ normal(0, 20); 
	sigma_b_inter_ps_sp_ind ~ normal(0, 20); 
	
	lday ~ normal(y_hat, sigma_y);
}