# Fake data testing of pheno budburst experiment at individual level

library(rstan)
library(ggplot2)
library(shinystan)

setwd("~/Documents/git/treetraits/analyses")
source('stan/savestan.R')
# get latest .Rdata file

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("Fake_Budburst_ind.RData") # not ind level, just testing

# load("Stan Output 2016-04-04 Fake Interax.RData")
# fakeout <- dir()[grep("Stan Output", dir())[is.na(match(grep("Stan Output", dir()), grep("Fake", dir())))]]
# load(sort(realout, T)[1])
# ls() 

## 
splookup <- unique(fake[c("ind","sp")])[,"sp"]

# To Stan!
datalist.f <- list(lday = fake$bb, # budburst as respose 
                   warm = as.numeric(fake$warm), 
                   site = as.numeric(fake$site), 
                   sp = as.numeric(fake$sp), 
                   ind = as.numeric(fake$ind),
                   photo = as.numeric(fake$photo), 
                   N = nrow(fake), 
                   splookup = splookup,
                   n_site = length(unique(fake$site)), 
                   n_sp = length(unique(fake$sp)),
                   n_ind = length(unique(fake$ind))
                  
                   )

doym.f <- stan('stan/lday_ind2.stan', data = datalist.f, 
               iter = 5005
#               , control = list(adapt_delta = 0.9,
#                               max_treedepth = 15)
                ) 

# lday_site_chill: < 120 seconds per chain, very fast
# lday_site_sp_chill: much slower.   
#doym.f <- stan('stan/lday0.stan', data = datalist.f, iter = 4000, chains = 4) 



sumer <- summary(doym.f)$summary
sumer[grep("mu_", rownames(sumer)),] # effects are perfectly captured now.

pairs(doym.f, pars = names(doym.f)[grep("mu_", names(doym.f))])

ssm.f <- as.shinystan(doym.f)
launch_shinystan(ssm.f) 

#setwd("~/Dropbox")

savestan("Fake ind test")

#setwd("~/Documents/git/buds/analyses")

# Load lastest fake data output. Grep for both Fake and Stan Output.
if(!exists('doym.f')){
  
  fakes <- na.omit(grep("Stan Output", dir())[match(grep("Fake", dir()), grep("Stan Output", dir()))])

  load(sort(dir()[fakes], T)[1])
}

sf <- summary(doym.f)$summary

plotlet("b_warm", "b_photo", 
        #xlab = "Advance due to 30d 4° chilling", 
        #ylab = "Advance due to 30d 1.5° chilling", 
        dat = sf)

plotlet("b_inter_wc2", "b_inter_wc1", 
        #xlab = "Advance due to 30d 4° chilling", 
        #ylab = "Advance due to 30d 1.5° chilling", 
        dat = sf)


di <- sf[grep("mu_b_inter", rownames(sf)),]

plot(seq(min(di[,"mean"]-di[,"sd"]*1.5), max(di[,"mean"]+di[,"sd"]*1.5), length.out = nrow(di)),
     1:nrow(di), type ="n")


# # now with fixed site and fixed sp
# 
# datalist.f <- list(lday = fake$bb, # budburst as respose 
#                    warm = as.numeric(fake$warm), 
#                    site = as.numeric(fake$site), 
#                    sp = as.numeric(fake$sp), 
#                    photo = as.numeric(fake$photo), 
#                    chill = as.numeric(fake$chill), 
#                    N = nrow(fake), 
#                    n_site = length(unique(fake$site)), 
#                    n_sp = length(unique(fake$sp))
# )
# 
# doym.f2 <- stan('stan/lday0_fixedsite_fixedsp.stan', data = datalist.f, iter = 4000, chains = 4) 
# 
# ssm.f <- as.shinystan(doym.f2)
# #launch_shinystan(ssm.f2) 
# 
# (sumer <- summary(doym.f)$summary)
# 
# setwd("~/Dropbox")
# 
# savestan()
# 
# setwd("~/Documents/git/buds/analyses")
# 
# # Tips for speeding up, from google groups
# set_cppo(mode = "fast")
# # For finding part of code that is slow
# dir(tempdir())



plotlet <- function(x, y, xlab=NULL, ylab=NULL, dat, groups = NULL){
  
  if(is.null(xlab)) xlab = x
  if(is.null(ylab)) ylab = y
  
  minmax = range(c(dat[grep(paste(x,"\\[",sep=""), rownames(dat)),1], dat[grep(paste(y,"\\[",sep=""), rownames(dat)),1]))
  
  if(is.null(groups)) { col.pch = "black"; col.lines = "grey50" }
  else {
    colz = c("midnightblue", "darkgreen")
    ccolz = rep(colz[1], length(groups))
    ccolz[groups == 2] = colz[2]
    col.pch = ccolz
    col.lines = alpha(ccolz, 0.4)
  }
  
    plot(
    dat[grep(paste(x,"\\[",sep=""), rownames(dat)),1],
    dat[grep(paste(y,"\\[",sep=""), rownames(dat)),1],
    pch = "+",
    xlim = c(floor(minmax)[1], ceiling(minmax)[2]),
    ylim = c(floor(minmax)[1], ceiling(minmax)[2]),
    ylab = ylab,
    xlab = xlab,
    col = col.pch
  )
  
  abline(h=0, lty = 3, col = "grey60")
  abline(v=0, lty = 3, col = "grey60")
  
  arrows(
    dat[grep(paste(x,"\\[",sep=""), rownames(dat)),"mean"],
    dat[grep(paste(y,"\\[",sep=""), rownames(dat)),"25%"],
    dat[grep(paste(x,"\\[",sep=""), rownames(dat)),"mean"],
    dat[grep(paste(y,"\\[",sep=""), rownames(dat)),"75%"],
    length = 0, col = col.lines)
  
  arrows(
    dat[grep(paste(x,"\\[",sep=""), rownames(dat)),"25%"],
    dat[grep(paste(y,"\\[",sep=""), rownames(dat)),"mean"],
    dat[grep(paste(x,"\\[",sep=""), rownames(dat)),"75%"],
    dat[grep(paste(y,"\\[",sep=""), rownames(dat)),"mean"],
    length = 0, col = col.lines)
  

}