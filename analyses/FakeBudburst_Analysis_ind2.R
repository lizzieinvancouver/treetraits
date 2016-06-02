# Fake data testing of pheno budburst experiment at individual level
dostan = T

library(rstan)
library(ggplot2)
library(shinystan)

setwd("~/Documents/git/treetraits/analyses")
source('stan/savestan.R')
# get latest .Rdata file

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("Fake_Budburst_ind2.RData") # not ind level, just testing

splookup <- unique(fake[c("ind","sp")])[,"sp"] #200 long = 10 ind x 20 sp 

# To Stan

datalist.f <- with(fake, 
 list(lday = bb, # budburst as respose 
                   warm = as.numeric(warm), 
                   site = as.numeric(site),
                   sp = as.numeric(sp), 
                   ind = as.numeric(ind),
                   photo = as.numeric(photo), 
                   N = nrow(fake), 
                   splookup = splookup,
                   n_sp = length(unique(sp)),
                   n_ind = length(unique(ind))
                   )
)
if(dostan){
  doym.f <- stan('stan/lday_ind2.stan', data = datalist.f, 
                 iter = 5005
                  ) 
  
  savestan("Fake ind interax 2")
}

# Load lastest fake data output. Grep for both Fake and Stan Output.
if(!exists('doym.f')){
  
  fakes <- na.omit(grep("Stan Output", dir())[match(grep("Fake", dir()), grep("Stan Output", dir()))])

  load(sort(dir()[fakes], T)[1])
}

sf <- summary(doym.f)$summary

ssm.f <- as.shinystan(doym.f)
launch_shinystan(ssm.f) 

plot(sf['b_warm_0','mean'], sf['b_photo_0','mean'])

plotlet("b_warm_sp", "b_photo_sp", 
#         xlim = c(-24, -23),
#         ylim = c(-18, -17),
        #xlab = "Advance due to 30d 4° chilling", 
        #ylab = "Advance due to 30d 1.5° chilling", 
        dat = sf)


plot(density(sf[grep("^a_sp\\[", rownames(sf)),'mean']),
     main = "Species level intercepts")

plot(density(sf[grep("^a_sp_ind\\[", rownames(sf)),'mean']),
     main = "Ind level intercepts")


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

nsp = length(unique(fake$sp)) # 20 species
ntreat = with(fake, length(unique(paste(warm, photo))))
# all(tapply(fake$ind, fake$sp, length) == 80) # make sure all are the same
nind = mean(tapply(fake$ind, fake$sp, length) / ntreat)  #10

nwarm = length(unique(fake$warm)) # 2 temperature treatments
nphoto = length(unique(fake$photo)) # 2 photoperiod treatments
#nchill = 3

rep = 1 # only 1 cutting from an individual within each combination of treatments. 

(ntot = nwarm*nphoto) # 4 rows. But will be looping over individuals and species below

# Build up the data frame

warm = gl(nwarm, rep, length = ntot)
photo = gl(nphoto, rep*nwarm, length = ntot)

treatcombo = paste(warm, photo, sep = "_")

(d <- data.frame(warm, photo, treatcombo))

# Extracting fitted values from the stan fit object
warmdiff = sf[grep("^b_warm_0", rownames(sf)),'mean']
photodiff = sf[grep("^b_photo_0", rownames(sf)),'mean']

######## SD for each treatment
warmdiff.sd = sf[grep("^b_warm_0", rownames(sf)),'sd']  
photodiff.sd = sf[grep("^b_photo_0", rownames(sf)),'sd'] 

mm <- model.matrix(~(warm+photo), data.frame(warm, photo))

############ !
spint <- sf[grep("^a_sp\\[", rownames(sf)),'mean'] # Was centered on 35 in fake data, why now 65??

poster <- vector() # to hold the posterior predictive checks

for(i in 1:nsp){ # loop over species, as these are the random effect modeled. 
  # Within species, have a loop for individuals
  
  indx <- which(splookup == i)
  
  coeff <- data.frame(
      sf[rownames(sf) %in% paste("a_sp_ind[", indx, "]", sep = ""),'mean'],
      sf[rownames(sf) %in% paste("b_warm_sp_ind[", indx, "]", sep = ""),'mean'],
      sf[rownames(sf) %in% paste("b_photo_sp_ind[", indx, "]", sep = ""),'mean']
  )
    
    for(j in 1:nind){ # simulate data for each individual, using these coefficients
      
      bb <- rnorm(n = length(warm), mean = mm %*% as.numeric(coeff[j,]), sd = 0.1)
    
      
      posterx <- data.frame(bb, sp = i, ind = paste(i, j, sep="_"),
                        warm, photo)
    
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