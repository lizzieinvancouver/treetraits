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

# load("Stan Output 2016-05-25 Fake ind test 2.RData")
# fakeout <- dir()[grep("Stan Output", dir())]
# load(sort(fakeout, T)[1])
# ls() 

#lme1 <- lmer(bb ~ site + warm + photo + (warm|sp) + (photo|sp), data = fake)

#lme2 <- lmer(bb ~ site + warm + photo + (warm|sp/ind) + (photo|sp/ind), data = fake)

## 
splookup <- unique(fake[c("ind","sp")])[,"sp"] #200 long = 10 ind x 20 sp 

# To Stan!
datalist.f <- list(lday = fake$bb, # budburst as respose 
                   warm = as.numeric(fake$warm), 
                   site = as.numeric(fake$site),
                   sp = as.numeric(fake$sp), 
                   ind = as.numeric(fake$ind),
                   photo = as.numeric(fake$photo), 
                   N = nrow(fake), 
                   splookup = splookup,
                   n_sp = length(unique(fake$sp)),
                   n_ind = length(unique(fake$ind))
                   )

doym.f <- stan('stan/lday_ind3.stan', data = datalist.f, 
               iter = 5005
#               , control = list(adapt_delta = 0.9,
#                               max_treedepth = 15)
                ) 

savestan("Fake ind interax 1")

ssm.f <- as.shinystan(doym.f)
launch_shinystan(ssm.f) 


#setwd("~/Documents/git/buds/analyses")

# Load lastest fake data output. Grep for both Fake and Stan Output.
if(!exists('doym.f')){
  
  fakes <- na.omit(grep("Stan Output", dir())[match(grep("Fake", dir()), grep("Stan Output", dir()))])

  load(sort(dir()[fakes], T)[1])
}

sf <- summary(doym.f)$summary

plotlet("b_warm_sp", "b_photo_sp", 
        #xlab = "Advance due to 30d 4° chilling", 
        #ylab = "Advance due to 30d 1.5° chilling", 
        dat = sf)

plotlet("b_warm_sp_ind", "b_photo_sp_ind", 
        #xlab = "Advance due to 30d 4° chilling", 
        #ylab = "Advance due to 30d 1.5° chilling", 
        dat = sf)


# di <- sf[grep("mu_b_inter", rownames(sf)),]
# 
# plot(seq(min(di[,"mean"]-di[,"sd"]*1.5), max(di[,"mean"]+di[,"sd"]*1.5), length.out = nrow(di)),
#      1:nrow(di), type ="n")


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