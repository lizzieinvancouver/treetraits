// flynn@fas.harvard.edu
// 3 level model for leafout day (lday) as a function of warming (2 levels) and photoperiod (2), chilling treatment (2 levels of 2 treatments), site (2), and their two-way interactions. 
// based off of previous model, pooling just a species level for warming, photoperiod, and chilling effects, and off of "Three-level nested linear model in Stan": 
// http://rstudio-pubs-static.s3.amazonaws.com/64315_bc3a395edd104095a8384db8d9952f43.html
// Difference in this model is that we are not just doing "beta_0" intercept as in the other model, but also slopes. I call intercept "a" and slopes "b".
// Top level: Species
// 2nd level: Individuals (6 per species per site, which is currently being ignored
// 3rd level: Cuttings (4-12 per individual, spread out across different treatments)

data {
	int<lower=1> N;
	int<lower=1> n_ind;
	int<lower=1> n_sp;
	int<lower=1, upper=n_sp> sp[N];
	int<lower=1, upper=n_ind> ind[N]; // this serves also as individual-level cluster id
	vector[N] lday; 	// response
	vector[N] warm; 	// predictor
	vector[N] photo; 	// predictor
	vector[N] site; 	// predictor
	vector[N] chill1;	// predictor
	vector[N] chill2;	// predictor

	int<lower=1> splookup[n_ind]; // cluster id for species
		
	}

transformed data { 			// 9 interaction terms 
	vector[N] inter_wp;           
	vector[N] inter_ws;           
	vector[N] inter_ps;           
	vector[N] inter_wc1;           
	vector[N] inter_wc2;           
	vector[N] inter_pc1;           
	vector[N] inter_pc2;           
	vector[N] inter_sc1;           
	vector[N] inter_sc2;    

	inter_wp	<- warm .* photo;  
	inter_ws	<- warm .* site;  
	inter_ps	<- photo .* site;  
	inter_wc1   <- warm .* chill1;  
	inter_wc2   <- warm .* chill2;  
	inter_pc1   <- photo .* chill1;  
	inter_pc2   <- photo .* chill2;  
	inter_sc1   <- site .* chill1;  
	inter_sc2   <- site .* chill2;  

}

parameters {
  real a_0; // overall intercept, same as "beta_0" in classroom example model
  real b_warm_0; // overall warming effect
  real b_photo_0; // overall photoperiod effect  
  real b_site_0; // overall site effect  
  real b_chill1_0; // overall chill1 effect  
  real b_chill2_0; // overall chill2 effect  
  real b_inter_wp_0; // overall warming x photoperiod effect  
  real b_inter_ws_0; // overall warming x site effect  
  real b_inter_ps_0; // overall photoperiod x site effect  
  real b_inter_wc1_0; // overall warming x chill1 effect  
  real b_inter_wc2_0; // overall warming x chill2 effect  
  real b_inter_pc1_0; // overall photoperiod x chill1 effect  
  real b_inter_pc2_0; // overall photoperiod x chill2 effect  
  real b_inter_sc1_0; // overall site x chill1 effect  
  real b_inter_sc2_0; // overall site x chill2 effect  

  real mu_a_sp[n_sp];
  real mu_b_warm_sp[n_sp]; 
  real mu_b_photo_sp[n_sp];
  real mu_b_site_sp[n_sp];
  real mu_b_chill1_sp[n_sp];
  real mu_b_chill2_sp[n_sp];
  real mu_b_inter_wp_sp[n_sp];
  real mu_b_inter_ws_sp[n_sp];
  real mu_b_inter_ps_sp[n_sp];
  real mu_b_inter_wc1_sp[n_sp];
  real mu_b_inter_wc2_sp[n_sp];
  real mu_b_inter_pc1_sp[n_sp];
  real mu_b_inter_pc2_sp[n_sp];
  real mu_b_inter_sc1_sp[n_sp];
  real mu_b_inter_sc2_sp[n_sp];
      
  real mu_a_sp_ind[n_ind]; // indexed with individual. All others with species

  real<lower=0> sigma_a_sp; 
  real<lower=0> sigma_b_warm_sp; 
  real<lower=0> sigma_b_photo_sp;
  real<lower=0> sigma_b_site_sp;
  real<lower=0> sigma_b_chill1_sp;
  real<lower=0> sigma_b_chill2_sp;
  real<lower=0> sigma_b_inter_wp_sp;
  real<lower=0> sigma_b_inter_ws_sp;
  real<lower=0> sigma_b_inter_ps_sp;
  real<lower=0> sigma_b_inter_wc1_sp;
  real<lower=0> sigma_b_inter_wc2_sp;
  real<lower=0> sigma_b_inter_pc1_sp;
  real<lower=0> sigma_b_inter_pc2_sp;
  real<lower=0> sigma_b_inter_sc1_sp;
  real<lower=0> sigma_b_inter_sc2_sp;

  real<lower=0> sigma_a_sp_ind; 

  real<lower=0> sigma_y; 
  }

transformed parameters {
	real y_hat[N];
	
	real a_sp[n_sp]; // intercept for species
  	real b_warm_sp[n_sp]; // slope of warming effect at species level
	real b_photo_sp[n_sp]; // slope of photoperiod effect, at species level
	real b_site_sp[n_sp]; // slope of photoperiod effect, at species level
	real b_chill1_sp[n_sp]; // slope of chill1 effect, at species level
	real b_chill2_sp[n_sp]; // slope of chill2 effect, at species level
	real b_inter_wp_sp[n_sp]; // slope of warming x photoperiod effect, at species level
	real b_inter_ws_sp[n_sp]; // slope of warming x photoperiod effect, at species level
	real b_inter_ps_sp[n_sp]; // slope of warming x photoperiod effect, at species level
	real b_inter_wc1_sp[n_sp]; 
	real b_inter_wc2_sp[n_sp]; 
	real b_inter_pc1_sp[n_sp]; 
	real b_inter_pc2_sp[n_sp]; 
	real b_inter_sc1_sp[n_sp]; 
	real b_inter_sc2_sp[n_sp]; 
		
	real a_sp_ind[n_ind]; // intercept for individuals within species
		
	// Species level. Random intercept (a) and slopes (b) for warming, photoperiod, chilling, site, 2-way interax
	for (k in 1:n_sp) {
		
		a_sp[k] <- a_0 + mu_a_sp[k];
		b_warm_sp[k] <- b_warm_0 + mu_b_warm_sp[k];
		b_photo_sp[k] <- b_photo_0 + mu_b_photo_sp[k];
		b_site_sp[k] <- b_site_0 + mu_b_site_sp[k];
		b_chill1_sp[k] <- b_chill1_0 + mu_b_chill1_sp[k];
		b_chill2_sp[k] <- b_chill2_0 + mu_b_chill2_sp[k];
		b_inter_wp_sp[k] <- b_inter_wp_0 + mu_b_inter_wp_sp[k];
		b_inter_ws_sp[k] <- b_inter_ws_0 + mu_b_inter_ws_sp[k];
		b_inter_ps_sp[k] <- b_inter_ps_0 + mu_b_inter_ps_sp[k];										
		b_inter_wc1_sp[k] <- b_inter_wc1_0 + mu_b_inter_wc1_sp[k];
		b_inter_wc2_sp[k] <- b_inter_wc2_0 + mu_b_inter_wc2_sp[k];
		b_inter_pc1_sp[k] <- b_inter_pc1_0 + mu_b_inter_pc1_sp[k];										
		b_inter_pc2_sp[k] <- b_inter_pc2_0 + mu_b_inter_pc2_sp[k];
		b_inter_sc1_sp[k] <- b_inter_sc1_0 + mu_b_inter_sc1_sp[k];
		b_inter_sc2_sp[k] <- b_inter_sc2_0 + mu_b_inter_sc2_sp[k];										

		}
	
	// individual level. Random intercepts only
	for (j in 1:n_ind){
	
		a_sp_ind[j] <- a_sp[splookup[j]] + mu_a_sp_ind[j];
	
	}
	
	// cutting (branch) level. 
	for(i in 1:N){

		y_hat[i] <- a_sp_ind[ind[i]] + // indexed with individual
					b_warm_sp[sp[i]] * warm[i] + // indexed with species
					b_photo_sp[sp[i]] * photo[i]+
					b_site_sp[sp[i]] * site[i] +
					b_chill1_sp[sp[i]] * chill1[i] +
					b_chill2_sp[sp[i]] * chill2[i] +
					b_inter_wp_sp[sp[i]] * inter_wp[i]+
					b_inter_ws_sp[sp[i]] * inter_ws[i]+
					b_inter_ps_sp[sp[i]] * inter_ps[i]+
					b_inter_wc1_sp[sp[i]] * inter_wc2[i]+
					b_inter_wc2_sp[sp[i]] * inter_wc2[i]+
					b_inter_pc1_sp[sp[i]] * inter_pc1[i]+
					b_inter_pc2_sp[sp[i]] * inter_pc2[i]+
					b_inter_sc1_sp[sp[i]] * inter_sc1[i]+
					b_inter_sc2_sp[sp[i]] * inter_sc2[i]
					;
					}
	}

model {
	mu_a_sp ~ normal(0, sigma_a_sp); // 20 d on either site at sp level
	mu_b_warm_sp ~ normal(0, sigma_b_warm_sp); 
	mu_b_photo_sp ~ normal(0, sigma_b_photo_sp);
	mu_b_site_sp ~ normal(0, sigma_b_site_sp);
	mu_b_chill1_sp ~ normal(0, sigma_b_chill1_sp);
	mu_b_chill2_sp ~ normal(0, sigma_b_chill2_sp);
	mu_b_inter_wp_sp ~ normal(0, sigma_b_inter_wp_sp);
	mu_b_inter_ws_sp ~ normal(0, sigma_b_inter_ws_sp);
	mu_b_inter_ps_sp ~ normal(0, sigma_b_inter_ps_sp);
	mu_b_inter_wc1_sp ~ normal(0, sigma_b_inter_wc1_sp);
	mu_b_inter_wc2_sp ~ normal(0, sigma_b_inter_wc2_sp);
	mu_b_inter_pc1_sp ~ normal(0, sigma_b_inter_pc1_sp);
	mu_b_inter_pc2_sp ~ normal(0, sigma_b_inter_pc2_sp);
	mu_b_inter_sc1_sp ~ normal(0, sigma_b_inter_sc1_sp);
	mu_b_inter_sc2_sp ~ normal(0, sigma_b_inter_sc2_sp);
	
	mu_a_sp_ind ~ normal(0, sigma_a_sp_ind); // 10 d on either side at individual level
	
	sigma_a_sp ~ normal(0, 20); // Start big at 10, go smaller if introduces problems
	sigma_b_warm_sp ~ normal(0, 20); 
	sigma_b_photo_sp ~ normal(0, 20); 
	sigma_b_site_sp ~ normal(0, 20); 
	sigma_b_chill1_sp ~ normal(0, 20); 
	sigma_b_chill2_sp ~ normal(0, 20); 
	sigma_b_inter_wp_sp ~ normal(0, 20); 
	sigma_b_inter_ws_sp ~ normal(0, 20); 
	sigma_b_inter_ps_sp ~ normal(0, 20); 
	sigma_b_inter_wc1_sp ~ normal(0, 20);
	sigma_b_inter_wc2_sp ~ normal(0, 20);
	sigma_b_inter_pc1_sp ~ normal(0, 20);
	sigma_b_inter_pc2_sp ~ normal(0, 20);
	sigma_b_inter_sc1_sp ~ normal(0, 20);
	sigma_b_inter_sc2_sp ~ normal(0, 20);
		
	sigma_a_sp_ind ~ normal(0, 10); // individual level sigma
	
	lday ~ normal(y_hat, sigma_y);
}