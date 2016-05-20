# Fake data for buburst stan work
library(dplyr)

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# Set up: same as experiment, with two sites, 28 species, two levels each of warming and photoperiod, and three levels of chilling. 2016-04-01 adding interactions. This ends up generating expected differences, but variation in effect sizes across species is minimal currently.
# 2016-05-16 simplifying a lot, but adding individuals
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

nsite = 2
nsp = 28
nind = 12
  
nwarm = 2
nphoto = 2
#nchill = 3

rep = 1 # only 1 individual within each combination of treatments. 

(ntot = nsite*nwarm*nphoto) # 8 rows. But will be looping over individuals and species below

# Build up the data frame
site = gl(nsite, rep, length = ntot)

warm = gl(nwarm, rep*nsite, length = ntot)
photo = gl(nphoto, rep*nsite*nwarm, length = ntot)

treatcombo = paste(warm, photo, sep = "_")

(d <- data_frame(site, warm, photo, treatcombo))

###### Set up differences for each level
sitediff = 2 
warmdiff = -20 # days earlier from 1 to 2
photodiff = -14

# interactions. 9 two-way interactions
sitewarm = 0
sitephoto = 0
warmphoto = 3.5 # positive 3.5. So at the warm level, the effect of longer days is muted by 3.5 days.

######## SD for each treatment
sitediff.sd = 1.5 
warmdiff.sd = 1 
photodiff.sd = 1
sitewarm.sd = 1
sitephoto.sd = 1
warmphoto.sd = 1
mm <- model.matrix(~(site+warm+photo)^2, data.frame(site, warm, photo))
colnames(mm)

coeff <- c(1, sitediff, warmdiff, photodiff, 
           sitewarm, sitephoto, 
           warmphoto
           )

bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 1) # should be able to do sd = mm %*% sd.coeff as well, with a different sd for each parameter.

(fake <- data_frame(bb, site, warm, photo))

##### Again, now with species now.

baseinter = 35 # baseline intercept across all species 
spint <- baseinter + c(1:nsp)-mean(1:nsp) # different intercepts by species

fake <- vector()

for(i in 1:nsp){ # loop over species, as these are the random effect modeled
  
  # Give species different difference values, drawn from normal. Could make dataframe of diffs and diff.sds, and use apply..
  coeff <- c(spint[i], 
             rnorm(1, sitediff, sitediff.sd),
             rnorm(1, warmdiff, warmdiff.sd),
             rnorm(1, photodiff, photodiff.sd), 
             rnorm(1, sitewarm, sitewarm.sd), 
             rnorm(1, sitephoto, sitephoto.sd),
             rnorm(1, warmphoto, warmphoto.sd)
             )
  
  bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 0.1)
  
  fakex <- data.frame(bb, sp = i, site, warm, photo)
      
  fake <- rbind(fake, fakex)  
  }

summary(lm(bb ~ (site+warm+photo)^2, data = fake)) # sanity check 

#summary(lmer(bb ~ (site|sp) + (warm|sp) + (photo|sp) + (chill1|sp) + (chill2|sp), data = fake)) # too hard for lmer.



# and again, with individuals

##### Again, now with species now.

baseinter = 35 # baseline intercept across all species 
spint <- baseinter + c(1:nsp)-mean(1:nsp) # different intercepts by species

fake <- vector()

for(i in 1:nsp){ # loop over species, as these are the random effect modeled. Within species, have a loop for individuals
  
  
  indint <- spint[i] + 1:nind-mean(1:nind)
  
  # Give species different difference values, drawn from normal.
  
  coeff <- c(spint[i], 
             rnorm(1, sitediff, sitediff.sd),
             rnorm(1, warmdiff, warmdiff.sd),
             rnorm(1, photodiff, photodiff.sd), 
             rnorm(1, sitewarm, sitewarm.sd), 
             rnorm(1, sitephoto, sitephoto.sd),
             rnorm(1, warmphoto, warmphoto.sd)
  )
  
  bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 0.1)
  
  fakex <- data.frame(bb, sp = i, site, warm, photo)
  
  fake <- rbind(fake, fakex)  
}

summary(lm(bb ~ (site+warm+photo)^2, data = fake)) # sanity check 

save(list=c("fake"), file = "Fake Budburst.RData")











save(list=c("fake"), file = "Fake Budburst.RData")



