# Fake data for buburst stan work

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# Set up: same as experiment, with two sites, 28 species, two levels each of warming and photoperiod, and three levels of chilling. 2016-04-01 adding interactions. This ends up generating expected differences, but variation in effect sizes across species is minimal currently.
# 2016-05-16 simplifying a lot, but adding individuals. Removed chilling and interactions for now.
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

nsp = 20
nind = 10
  
nwarm = 2
nphoto = 2

rep = 1 # only 1 individual within each combination of treatments. 

(ntot = nwarm*nphoto) # 8 rows. But will be looping over individuals and species below

# Build up the data frame

warm = gl(nwarm, rep, length = ntot)
photo = gl(nphoto, rep*nwarm, length = ntot)

treatcombo = paste(warm, photo, sep = "_")

(d <- data.frame(warm, photo, treatcombo))

###### Set up differences for each level
warmdiff = -20 # days earlier from 1 to 2
photodiff = -14

######## SD for each treatment
warmdiff.sd = 1 
photodiff.sd = 1

mm <- model.matrix(~(warm+photo), data.frame(warm, photo))

#  with individuals

baseinter = 35 # baseline intercept across all species 
spint <- baseinter + c(1:nsp)-mean(1:nsp) # different intercepts by species

fake <- vector()

for(i in 1:nsp){ # loop over species
  
  # Give variation for individual trees within species in the intercept
  indint <- spint[i] + 1:nind-mean(1:nind)
  
  for(j in 1:nind){ # Within species, loop for individuals
  # Give species different difference values, drawn from normal.
  
  coeff <- c(indint[j], 
             rnorm(1, warmdiff, warmdiff.sd),
             rnorm(1, photodiff, photodiff.sd)
  )
  
  bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 0.1)
  
  fakex <- data.frame(bb, sp = i, ind = paste(i, j, sep="_"),
                      warm, photo)
  
  fake <- rbind(fake, fakex)  
  }
  }

summary(lm(bb ~ (warm+photo)^2, data = fake)) # sanity check 

save(list=c("fake"), file = "Fake_Budburst_ind2.RData")




