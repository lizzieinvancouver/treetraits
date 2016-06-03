# Fake data testing of pheno budburst experiment at individual level
dostan = FALSE

library(rstan)
library(ggplot2)
library(shinystan)

setwd("~/Documents/git/treetraits/analyses")
source('stan/savestan.R')
# get latest .Rdata file

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("Fake_Budburst_ind.RData") # not ind level, just testing

splookup <- unique(fake[c("ind","sp")])[,"sp"] #200 long = 10 ind x 20 sp 

# To Stan!

datalist.f <- with(fake, 
 list(lday = bb, 
        warm = as.numeric(warm), 
        photo = as.numeric(photo),
        chill1 = as.numeric(chill1), 
        chill2 = as.numeric(chill2),
        site = as.numeric(site),
        sp = as.numeric(sp),
        ind = as.numeric(ind),
        N = nrow(fake),
        splookup = splookup,
        n_sp = length(unique(sp)),
        n_ind = length(unique(ind))
        )
)
if(dostan){
  doym.f <- stan('stan/lday_ind5.stan', data = datalist.f, 
                 iter = 2002
                  ) 
  sf <- summary(doym.f)$summary
  
  ssm.f <- as.shinystan(doym.f)
  
  savestan("Fake ind")
}

# Load lastest fake data output. Grep for both Fake and Stan Output.
if(!exists('doym.f')){
  
  fakes <- na.omit(grep("Stan Output", dir())[match(grep("Fake", dir()), grep("Stan Output", dir()))])

  load(sort(dir()[fakes], T)[1])
}

launch_shinystan(ssm.f) 

sf[grep("^a_sp_ind\\[", rownames(sf)),] # Was centered on 35 in fake data, why now 65?
sf[grep("^a_sp\\[", rownames(sf)),] # Was centered on 35 in fake data, why now 65?
sf[grep("^a_0", rownames(sf)),] # 61 is center

sf[grep("^b_warm_sp\\[", rownames(sf)),] # -16.8 each sp
sf[grep("^b_warm_0", rownames(sf)),] # -1.6 overall. Sum with species level to re-capture -19
sf[grep("^mu_b_warm_sp\\[", rownames(sf)),] # 0

sf[grep("^b_photo_sp\\[", rownames(sf)),] # -17.1 each sp
sf[grep("^b_photo_0", rownames(sf)),] # -1.7 overall. Sum with species level, -19.8, but only put in -14!

sf[grep("^b_site_sp\\[", rownames(sf)),] # 1 each sp
sf[grep("^b_site_0", rownames(sf)),] # 0.8 overall. Sum with species level, 1.8, put in 2, ok

sf[grep("^b_chill1_sp\\[", rownames(sf)),] # -15.5 each sp
sf[grep("^b_chill1_0", rownames(sf)),] # -1.5 overall. Sum with species level, -17, should be -20


sf[grep("^b_chill2_sp\\[", rownames(sf)),] # -15.5 each sp
sf[grep("^b_chill2_0", rownames(sf)),] # -1.5 overall. Sum with species level, -17, should be -20

sf[grep("^b_inter_wp_sp\\[", rownames(sf)),] # 3.4, right on
sf[grep("^b_inter_wp_0", rownames(sf)),] # 3.4 also! overall. Don't sum with species level, -19.8, but only put in -14!


plot(sf['b_warm_0','mean'], sf['b_photo_0','mean'])

plotlet("b_warm_sp", "b_photo_sp", 
#         xlim = c(-24, -23),
#         ylim = c(-18, -17),
        #xlab = "Advance due to 30d 4° chilling", 
        #ylab = "Advance due to 30d 1.5° chilling", 
        dat = sf)

plot(density(sf[grep("^a_sp\\[", rownames(sf)),'mean']),
     main = "Species level intercepts",
     col = "midnightblue")

lines(density(sf[grep("^a_sp_ind\\[", rownames(sf)),'mean']),
     col = "darkred")

plotlet("mu_b_warm_sp", "mu_b_photo_sp", 
        dat = sf)

plot(density(sf[grep("^b_warm_sp\\[", rownames(sf)),'mean']),
     main = "Species level warming effects")

plot(density(sf[grep("^b_warm_sp_ind\\[", rownames(sf)),'mean']),
     main = "Individual level warming effects")

plotlet("b_warm_sp_ind", "b_photo_sp_ind", 
#         xlim = c(-30, -20),
#         ylim = c(-20, -15),
        dat = sf)

plotlet("mu_b_warm_sp_ind", "mu_b_photo_sp_ind", 
        dat = sf)


###### Posterior predictive checks
# pull out coefficients at each level

nsite = length(unique(fake$site)) # 2 sites
nsp = length(unique(fake$sp)) # 20 species
ntreat = with(fake, length(unique(paste(site, warm, photo))))
# all(tapply(fake$ind, fake$sp, length) == 80) # make sure all are the same
nind = mean(tapply(fake$ind, fake$sp, length) / ntreat)  #10


nwarm = length(unique(fake$warm)) # 2 temperature treatments
nphoto = length(unique(fake$photo)) # 2 photoperiod treatments
#nchill = 3

rep = 1 # only 1 cutting from an individual within each combination of treatments. 

(ntot = nsite*nwarm*nphoto) # 8 rows. But will be looping over individuals and species below

# Build up the data frame
site = gl(nsite, rep, length = ntot)

warm = gl(nwarm, rep*nsite, length = ntot)
photo = gl(nphoto, rep*nsite*nwarm, length = ntot)

treatcombo = paste(warm, photo, sep = "_")

(d <- data.frame(site, warm, photo, treatcombo))

# Extracting fitted values from the stan fit object
sitediff = sf[grep("^b_site_0", rownames(sf)),'mean']
warmdiff = sf[grep("^b_warm_0", rownames(sf)),'mean']
photodiff = sf[grep("^b_photo_0", rownames(sf)),'mean']

# interactions. 9 two-way interactions
sitewarm = sf[grep("^b_inter_ws_0", rownames(sf)),'mean']
sitephoto = sf[grep("^b_inter_ps_0", rownames(sf)),'mean']
warmphoto = sf[grep("^b_inter_wp_0", rownames(sf)),'mean']

######## SD for each treatment
sitediff.sd = sf[grep("^b_site_0", rownames(sf)),'sd'] 
warmdiff.sd = sf[grep("^b_warm_0", rownames(sf)),'sd']  
photodiff.sd = sf[grep("^b_photo_0", rownames(sf)),'sd'] 
sitewarm.sd = sf[grep("^b_inter_ws_0", rownames(sf)),'sd'] 
sitephoto.sd = sf[grep("^b_inter_ps_0", rownames(sf)),'sd'] 
warmphoto.sd = sf[grep("^b_inter_wp_0", rownames(sf)),'sd'] 

mm <- model.matrix(~(site+warm+photo)^2, data.frame(warm, photo))

############ !
spint <- sf[grep("^a_sp\\[", rownames(sf)),'mean'] # Was centered on 35 in fake data, why now 61?

poster <- vector() # to hold the posterior predictive checks

for(i in 1:nsp){ # loop over species, as these are the random effect modeled. 
  # Within species, have a loop for individuals
  
  indx <- which(splookup == i)
  
  coeff <- data.frame(
      sf[rownames(sf) %in% paste("a_sp_ind[", indx, "]", sep = ""),'mean'],
      sf[rownames(sf) %in% paste("b_site_sp_ind[", indx, "]", sep = ""),'mean'],
      sf[rownames(sf) %in% paste("b_warm_sp_ind[", indx, "]", sep = ""),'mean'],
      sf[rownames(sf) %in% paste("b_photo_sp_ind[", indx, "]", sep = ""),'mean'],
      sf[rownames(sf) %in% paste("b_inter_ws_sp_ind[", indx, "]", sep = ""),'mean'],
      sf[rownames(sf) %in% paste("b_inter_ps_sp_ind[", indx, "]", sep = ""),'mean'],
      sf[rownames(sf) %in% paste("b_inter_wp_sp_ind[", indx, "]", sep = ""),'mean']
  )
    
    for(j in 1:nind){ # simulate data for each individual, using these coefficients
      
      bb <- rnorm(n = length(warm), mean = mm %*% as.numeric(coeff[j,]), sd = 0.1)
    
      
      posterx <- data.frame(bb, sp = i, ind = paste(i, j, sep="_"),
                        site, warm, photo)
    
    poster <- rbind(poster, posterx)  
  }
}


plot(fake$bb, poster$bb,
     xlab = "Simulated data",
     ylab = "Posterior predictive check",
     xlim = c(-10, 90),
     ylim = c(-10, 90),
     pch = 16,
     col = alpha('midnightblue', 0.2)
     )
abline(a=0, b=1, lty=3)



##################################################################################################################################################################################################


plotlet <- function(x, y, xlab=NULL, ylab=NULL, dat, groups = NULL, ...){
  if(is.null(xlab)) xlab = x; if(is.null(ylab)) ylab = y
  if(is.null(groups)) { col.pch = "black"; col.lines = "grey50" }
  else {
    colz = c("brown", "blue3")
    ccolz = rep(colz[1], length(groups))
    ccolz[groups == 2] = colz[2]
    col.pch = ccolz
    col.lines = alpha(ccolz, 0.4)
  }
  
  plot(
    dat[grep(paste("^", x,"\\[",sep=""), rownames(dat)),1],
    dat[grep(paste("^", y,"\\[",sep=""), rownames(dat)),1],
    pch = "+",
    ylab = ylab,
    xlab = xlab,
    col = col.pch,
    ...
  )
  
  abline(h=0, lty = 3, col = "grey60")
  abline(v=0, lty = 3, col = "grey60")
  
  arrows(
    dat[grep(paste("^", x,"\\[",sep=""), rownames(dat)),"mean"],
    dat[grep(paste("^", y,"\\[",sep=""), rownames(dat)),"25%"],
    dat[grep(paste("^", x,"\\[",sep=""), rownames(dat)),"mean"],
    dat[grep(paste("^", y,"\\[",sep=""), rownames(dat)),"75%"],
    length = 0, col = col.lines)
  
  arrows(
    dat[grep(paste("^", x,"\\[",sep=""), rownames(dat)),"25%"],
    dat[grep(paste("^", y,"\\[",sep=""), rownames(dat)),"mean"],
    dat[grep(paste("^", x,"\\[",sep=""), rownames(dat)),"75%"],
    dat[grep(paste("^", y,"\\[",sep=""), rownames(dat)),"mean"],
    length = 0, col = col.lines)
  

}