

// Adding individuals to the random effects
// Now with chilling levels as dummy variables for each level, two levels. Omitting for simplicity now, also omitting interactions
// Including site as intercept
// Leafout day as a function of species and site of origin as modeled group level factors, and temperature, photoperiod, and chilling as unmodeled factors (experimental manipulation)

data {
  int<lower=0> N;
  int<lower=0> n_ind;
  int<lower=0> n_sp;
  int<lower=0> n_site;
  int<lower=1, upper=n_site> site[N];
  int<lower=1, upper=n_sp> sp[N];
  int<lower=1, upper=n_ind> ind[N];
  vector[N] lday;
  vector[N] warm;
  vector[N] photo;
//  vector[N] chill1; 
//  vector[N] chill2;
}

parameters {
  vector[n_site] a_site;
  vector[n_sp] a_sp;
  vector[n_sp] b_warm;
  vector[n_sp] b_photo;
//  vector[n_sp] b_chill1;
//  vector[n_sp] b_chill2;
 
  real mu_b_warm; 
//  real mu_b_chill1;
//  real mu_b_chill2;
  real mu_b_photo;

  real<lower=0> sigma_b_warm;
  real<lower=0> sigma_b_photo;
//  real<lower=0> sigma_b_chill1;
//  real<lower=0> sigma_b_chill2;
  
  real<lower=0> sigma_y; 
  }


transformed parameters {
	vector[N] y_hat;
	// Using ind within sp as in earnings_latin_square
		
	for(i in 1:N){

		y_hat[i] <- a_site[site[i]] + a_sp[sp[i],ind[i]] + b_warm[sp[i],ind[i]] * warm[i] + b_photo[sp[i],ind[i]] * photo[i];
		
		}
	
}

model {
	// Priors. Make them flat
	mu_b_warm ~ normal(0, 35); // 100 = 3 months on either side. Narrow down to 35
	mu_b_photo ~ normal(0, 35);
//	mu_b_chill1 ~ normal(0, 35);
//	mu_b_chill2 ~ normal(0, 35);

	sigma_b_warm ~ normal(0, 10); // Start big at 10, go smaller if introduces problems
	sigma_b_photo ~ normal(0, 10); 
//	sigma_b_chill1 ~ normal(0, 10);
//	sigma_b_chill2 ~ normal(0, 10);
	
	b_warm ~ normal(mu_b_warm, sigma_b_warm);
	b_photo ~ normal(mu_b_photo, sigma_b_photo);
//	b_chill1 ~ normal(mu_b_chill1, sigma_b_chill1);
//	b_chill2 ~ normal(mu_b_chill2, sigma_b_chill2);

	lday ~ normal(y_hat, sigma_y);
}