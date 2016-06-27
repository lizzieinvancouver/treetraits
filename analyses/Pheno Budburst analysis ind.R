forlatex = F # set to F if just trying new figures, T if outputting for final
runstan = T # set to T to actually run stan models. F if loading from previous runs

# Analysis of budburst experiment 2015, at individual level. 

library(xtable)
library(ggplot2)

setwd("~/Documents/git/treetraits/analyses")

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# get latest .Rdata file
# To run from saved stan output (exclude Fake data output)
if(!runstan) {

  realout <- dir()[grep("Stan Output", dir())[is.na(match(grep("Stan Output", dir()), grep("Fake", dir())))]]
  if(!exists("doym.b")) load(sort(realout, T)[1]) # only run if stan output file is not yet in working memory.
  ls() 
  #launch_shinystan(doym.l)
}

if(runstan){ # things needed only if running the stan models
  library(rstan)
  library(shinystan) 
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  source('stan/savestan.R')
}

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> 
(toload <- sort(dir("./input")[grep("Budburst Data", dir('./input'))], T)[1])

load(file.path("input", toload))

if(forlatex) figpath = "../docs/ms/images" else figpath = "graphs"

# trait data

dt <- read.csv("./input/Summer 2015 Tree Traits.csv")

dt$sla <- with(dt, Leaf.area / Dry.mass)
dt$wd <- with(dt, Stem.mass / Stem.volume)

# Pheno Prep 
levels(dx$warm) = c(0,1); levels(dx$photo) = c(0, 1); levels(dx$site) = 1:2; 
levels(dx$chill) = 1:3
dx$warm <- as.numeric(dx$warm) - 1 # to set to 0 / 1
dx$photo <- as.numeric(dx$photo) - 1
dx$chill <- as.numeric(dx$chill)
dx$site <- as.numeric(dx$site) - 1 # 0 / 1 

dx <- dx[!is.na(dx$site),]

# Chill dummy variables
dx$chill1 = ifelse(dx$chill == 2, 1, 0) 
dx$chill2 = ifelse(dx$chill == 3, 1, 0) 

# with(dx, table(chill1, chill2)) # all three levels in here

dxb <- dx[!is.na(dx$bday),]
dxb$spn <- as.numeric(dxb$sp)
dxb$ind <- as.numeric(as.factor(as.character(dxb$ind)))

dxl <- dx[!is.na(dx$lday),]
dxl$spn <- as.numeric(dxl$sp)
dxl$ind <- as.numeric(as.factor(as.character(dxl$ind)))

# Trait prep
levels(dt$Site) = c(3, 1, 4, 2)

dt$site <- as.numeric(as.character(dt$Site))   # start at 1

dt <- dt[!is.na(dt$Latitude),]

dts <- dt[!is.na(dt$sla),]
dts$spn <- as.numeric(as.factor(as.character(dts$Species)))
dts$ind <- as.numeric(as.factor(as.character(dts$Individual)))

dtw <- dt[!is.na(dt$wd),]
dtw$spn <- as.numeric(as.factor(as.character(dtw$Species)))
dtw$ind <- as.numeric(as.factor(as.character(dtw$Individual)))

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

# Analyses:
# 1. Day of budburst by all factors, stan - using individual level model
# 2. Day of leaf out by all factors, stan - using individual level model
# 3. SLA by latitude
# 4. WD by latitude
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

# 1. Budburst day. 
#if(runstan){
  splookup <- unique(dxb[c("ind","spn")])[,"spn"] #271 long for budburst
  
  datalist.b <- with(dxb, list(lday = bday, # budburst as response 
                     warm = warm, 
                     site = site, 
                     sp = spn, 
                     ind = ind,
                     photo = photo, 
                    chill1 = chill1,
                    chill2 = chill2,
                     N = nrow(dxb), 
                     splookup = splookup,
                     n_site = length(unique(site)), 
                     n_sp = length(unique(spn)),
                     n_ind = length(unique(ind))
  ))

    doym.b <- stan('stan/lday_ind5.stan', 
                 data = datalist.b, iter = 5005, chains = 4
                  , control = list(adapt_delta = 0.9,
                                 max_treedepth = 15)
                      ) 

#}
  sumerb <- summary(doym.b)$summary
  sumerb[grep("mu_", rownames(sumerb)),]
  
  #ssm.b <- as.shinystan(doym.b)
  # launch_shinystan(ssm.b) 
  #yb = dxb$bday # for shinystan posterior checks

# plot effects
col4table <- c("mean","sd","25%","50%","75%","Rhat")

# manually to get right order
params <- c("b_warm_0","b_photo_0",
            "b_chill1_0","b_chill2_0",
            "b_site_0",
               "b_inter_wp_0",
                "b_inter_wc1_0","b_inter_wc2_0",
                "b_inter_pc1_0","b_inter_pc2_0",
               "b_inter_ws_0","b_inter_ps_0",
                "b_inter_sc1_0","b_inter_sc2"
)

meanzb <- sumerb[params,col4table]

rownames(meanzb) = c("Temperature",
                    "Photoperiod",
#                     "Chilling 4°",
#                     "Chilling 1.5°C",
                    "Site",
                    "Temperature x Photoperiod",
#                     "Temperature x Chilling 4°C",
#                     "Temperature x Chilling 1.5°C",
#                     "Photoperiod x Chilling 4°C",
#                     "Photoperiod x Chilling 1.5°C",
                    "Temperature x Site",
                    "Photoperiod x Site"
# ,
#                     "Site x Chilling 4°C",
#                     "Site x Chilling 1.5°C"
                    )

#if(runstan){
  splookup <- unique(dxl[c("ind","spn")])[,"spn"] #265 long
  
  datalist.l <- with(dxl, list(lday = lday, # budburst as response 
                               warm = warm, 
                               site = site, 
                               sp = spn, 
                               ind = ind,
                               photo = photo, 
                                chill1 = chill1,
                                chill2 = chill2,
                               N = nrow(dxl), 
                               splookup = splookup,
                               n_site = length(unique(site)), 
                               n_sp = length(unique(spn)),
                               n_ind = length(unique(ind))
  ))
  
  
  
  doym.l <- stan('stan/lday_ind5.stan', 
                 data = datalist.l, iter = 5005, chains = 4
                                   , control = list(adapt_delta = 0.9,
                                                  max_treedepth = 15)
  ) 

sumerl <- summary(doym.l)$summary

setwd("/Volumes/WeldShare/Wolkovich Lab/Dan")

savestan("Ind Models")



# ssm.l <- as.shinystan(doym.l)
# yl = dxl$lday # for shinystan posterior checks
# launch_shinystan(ssm.l) 

meanzl <- sumerl[mu_params,col4table]
rownames(meanzl) = rownames(meanzb)

gotchill <- tapply(dx$spn, dx$chill2, unique)$'1'
nochill <- unique(dx$spn)[is.na(match(unique(dx$spn), gotchill))]
sumerb[!is.na(match(rownames(sumerb), paste("b_chill1[", nochill, "]", sep=""))),] = NA
sumerb[!is.na(match(rownames(sumerb), paste("b_chill2[", nochill, "]", sep=""))),] = NA
sumerl[!is.na(match(rownames(sumerl), paste("b_chill1[", nochill, "]", sep=""))),] = NA
sumerl[!is.na(match(rownames(sumerl), paste("b_chill2[", nochill, "]", sep=""))),] = NA



###### Trait models

if(runstan){
  splookup <- unique(dts[c("ind","spn")])[,"spn"] #265 long
  
  datalist.s <- with(dts, list(y = sla, # SLA
                               lat = as.numeric(Latitude),
                               site = site, 
                               sp = spn, 
                               ind = ind,
                               N = nrow(dts), 
                               splookup = splookup,
                               n_site = length(unique(site)), 
                               n_sp = length(unique(spn)),
                               n_ind = length(unique(ind))
  ))
  
  
  
  latm.s <- stan('stan/trait_ind.stan', 
                 data = datalist.s, iter = 5005, chains = 4
                 #                  , control = list(adapt_delta = 0.9,
                 #                                 max_treedepth = 15)
  ) 
  
  splookup <- unique(dtw[c("ind","spn")])[,"spn"] #265 long
  
  datalist.w <- with(dtw, list(y = wd, # SLA
                               lat = as.numeric(Latitude),
                               site = site, 
                               sp = spn, 
                               ind = ind,
                               N = nrow(dtw), 
                               splookup = splookup,
                               n_site = length(unique(site)), 
                               n_sp = length(unique(spn)),
                               n_ind = length(unique(ind))
  ))
  
  
  
  latm.w <- stan('stan/trait_ind.stan', 
                 data = datalist.w, iter = 5005, chains = 4
                 #                  , control = list(adapt_delta = 0.9,
                 #                                 max_treedepth = 15)
  ) 

slat <- summary(latm.s)$summary
wlat <- summary(latm.w)$summary

savestan("Trait Models")
launch_shinystan(latm.s) # decreasing SLA with latitude; 

launch_shinystan(latm.w) # Very bad fit with wood density, why?




## Correlations between main effects and lo/bb
# warm, photo, chill1, chill2 vs. day of lo and day of bb

#bb, warm
blat.s <- slat[grep(paste("^b_lat_sp","\\[",sep=""), rownames(slat)),1]
blat.w <- wlat[grep(paste("^b_lat_sp","\\[",sep=""), rownames(wlat)),1]


pdf(file.path(figpath, "Sens_vs_day.pdf"), width = 9, height = 7)

plot(blat.s, blat.w)
  




on.exit(setwd("~/Documents/git/treetraits/docs/ms/"))




# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# Utility function to plot 'random effects' from stan output - used now mostly in Fig 3.
plotlet <- function(x, y, xlab=NULL, ylab=NULL, data, groups = NULL, ...){
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
    data[grep(paste(x,"\\[",sep=""), rownames(data)),1],
    data[grep(paste(y,"\\[",sep=""), rownames(data)),1],
    pch = "+",
    ylab = ylab,
    xlab = xlab,
    col = col.pch,
    ...
  )
  
  abline(h=0, lty = 3, col = "grey60")
  abline(v=0, lty = 3, col = "grey60")
  
  arrows(
    data[grep(paste(x,"\\[",sep=""), rownames(data)),"mean"],
    data[grep(paste(y,"\\[",sep=""), rownames(data)),"mean"]-data[grep(paste(y,"\\[",sep=""), rownames(data)),"se_mean"],
    data[grep(paste(x,"\\[",sep=""), rownames(data)),"mean"],
    data[grep(paste(y,"\\[",sep=""), rownames(data)),"mean"]+data[grep(paste(y,"\\[",sep=""), rownames(data)),"se_mean"],
    length = 0, col = col.lines)
  
  arrows(
    data[grep(paste(x,"\\[",sep=""), rownames(data)),"mean"]-data[grep(paste(x,"\\[",sep=""), rownames(data)),"se_mean"],
    data[grep(paste(y,"\\[",sep=""), rownames(data)),"mean"],
    data[grep(paste(x,"\\[",sep=""), rownames(data)),"mean"]+data[grep(paste(x,"\\[",sep=""), rownames(data)),"se_mean"],
    data[grep(paste(y,"\\[",sep=""), rownames(data)),"mean"],
    length = 0, col = col.lines)
  
  # match with species names
  text( data[grep(paste(x,"\\[",sep=""), rownames(data)),1],
        data[grep(paste(y,"\\[",sep=""), rownames(data)),1],
        sort(unique(dx$sp)),
        cex = 0.5, 
        pos = 3,
        col = col.pch)
}

# Groups
colz = c("brown", "blue3")

shrubs = c("VIBLAN","RHAFRA","RHOPRI","SPIALB","VACMYR","VIBCAS", "AROMEL","ILEMUC", "KALANG", "LONCAN", "LYOLIG")
trees = c("ACEPEN", "ACERUB", "ACESAC", "ALNINC", "BETALL", "BETLEN", "BETPAP", "CORCOR", "FAGGRA", "FRANIG", "HAMVIR", "NYSSYL", "POPGRA", "PRUPEN", "QUEALB" , "QUERUB", "QUEVEL")

treeshrub = levels(dx$sp)
treeshrub[treeshrub %in% shrubs] = 1
treeshrub[treeshrub %in% trees] = 2
treeshrub = as.numeric(treeshrub)