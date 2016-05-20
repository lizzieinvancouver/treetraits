runstan = T

library(lme4)
library(sjPlot)
# Compiled tree traits


# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# Read in data

setwd("~/Documents/git/treetraits/analyses")
d <- read.csv("input/Summer 2015 Tree Traits.csv")

d$sla <- d$Leaf.area/d$Dry.mass
d$wd <- d$Stem.mass/d$Stem.volume

d$height = d$Height
d$height = as.numeric(sub("<1",0.5, as.character(d$height)))

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# prep for stan

if(runstan){ # things needed only if running the stan models
  library(rstan)
  library(shinystan) 
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  source('stan/savestan.R')
}


dS <- as.factor(d$Site)
levels(dS) = c(3, 1, 4, 2)
d$site = as.numeric(as.character(dS))

d$sp <- as.numeric(as.factor(d$Species))

d$ind <- as.numeric(as.factor(d$Individual))


# Individual level
m1 <- lmer(sla ~ wd + Latitude + (1|sp) + (1|site), data = d)

ranef(m1)
fixef(m1)

sjp.lmer(m1, type = "re") # with increasing latitude, decreasing SLA = thicker/denser at higher lat
sjp.lmer(m1, type = "fe", showIntercept = F)

m2 <- lmer(sla ~ wd + Latitude + (1|site/sp), data = d)

ranef(m2)
fixef(m2)

sjp.lmer(m2, type = "re")


# Can we include individual in this kind of model?

m3 <- lmer(sla ~ wd + Latitude + (1|site/sp/ind), data = d)

ranef(m3)
fixef(m3)

sjp.lmer(m3, type = 're') # too many here. 

# Stan versions
# Load fake data (for now just previous version, without ind level variation)


if(runstan){
  datalist.sla <- list(sla = d$sla, 
                     site = d$site, 
                     sp = d$sp, 
                     ind = d$ind,
                     lat = d$Latitude,
                     N = nrow(d), 
                     n_site = length(unique(d$site)), 
                     n_sp = length(unique(d$sp)),
                     n_ind = length(unique(d$ind))
                     )
  
  stan.sla <- stan('stan/Traits_ind.stan', 
                 data = datalist.sla, iter = 1001, chains = 4
                 ) 
  
}
sumerb <- summary(doym.b)$summary
sumerb[grep("mu_", rownames(sumerb)),]


# look at consistency of performance within individuals, across treatments, as measure of plasticity.

# Is consistency related to earlier leafout?

vars <- aggregate(lday ~ sp + site + ind + wd + sla + X.N + Pore.anatomy, FUN = function(x) sd(x, na.rm=T) 
                  / mean(x, na.rm=T)
                  , data = dxt)
