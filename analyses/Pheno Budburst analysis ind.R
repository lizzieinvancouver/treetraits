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
dx$warm <- as.numeric(dx$warm)
dx$photo <- as.numeric(dx$photo)
dx$chill <- as.numeric(dx$chill)
dx$site <- as.numeric(dx$site)

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

dt$site <- as.numeric(dt$Site)

dt <- dt[!is.na(dt$Latitude),]

dts <- dt[!is.na(dt$sla),]
dts$spn <- as.numeric(dts$Species)
dts$ind <- as.numeric(as.factor(as.character(dts$Individual)))

dtw <- dt[!is.na(dt$lday),]
dtw$spn <- as.numeric(dtw$Species)
dtw$ind <- as.numeric(as.factor(as.character(dtw$Individual)))

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
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

# Analyses:
# 1. Day of budburst by all factors, stan - using individual level model
# 2. Day of leaf out by all factors, stan - using individual level model

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

# 1. Budburst day. 
if(runstan){
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
                 data = datalist.b, iter = 4004, chains = 4
#                  , control = list(adapt_delta = 0.9,
#                                 max_treedepth = 15)
                      ) 

}
  sumerb <- summary(doym.b)$summary
  sumerb[grep("mu_", rownames(sumerb)),]
  
  #ssm.b <- as.shinystan(doym.b)
  # launch_shinystan(ssm.b) 
  #yb = dxb$bday # for shinystan posterior checks

# plot effects
col4table <- c("mean","sd","25%","50%","75%","Rhat")

# manually to get right order
params <- c("b_warm_0","b_photo_0",
            #"b_chill1_0","b_chill2_0",
            "b_site_0",
               "b_inter_wp_0",
#                "b_inter_wc1_0","b_inter_wc2_0",
#                "b_inter_pc1_0","b_inter_pc2_0",
               "b_inter_ws_0","b_inter_ps_0"
               # "b_inter_sc1_0","b_inter_sc2"
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

if(runstan){
  splookup <- unique(dxl[c("ind","spn")])[,"spn"] #265 long
  
  datalist.l <- with(dxl, list(lday = lday, # budburst as response 
                               warm = warm, 
                               site = site, 
                               sp = spn, 
                               ind = ind,
                               photo = photo, 
                               #                      chill1 = chill1,
                               #                      chill2 = chill2,
                               N = nrow(dxl), 
                               splookup = splookup,
                               n_site = length(unique(site)), 
                               n_sp = length(unique(spn)),
                               n_ind = length(unique(ind))
  ))
  
  
  
  doym.l <- stan('stan/lday_ind2.stan', 
                 data = datalist.l, iter = 5005, chains = 4
                 #                  , control = list(adapt_delta = 0.9,
                 #                                 max_treedepth = 15)
  ) 
}
sumerl <- summary(doym.l)$summary

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
                 data = datalist.s, iter = 2002, chains = 4
                 #                  , control = list(adapt_delta = 0.9,
                 #                                 max_treedepth = 15)
  ) 
}
slat <- summary(latm.s)$summary

launch_shinystan(latm.s)



## Correlations between main effects and lo/bb
# warm, photo, chill1, chill2 vs. day of lo and day of bb

#bb, warm
bwarm <- sumerb[grep(paste("b_warm","\\[",sep=""), rownames(sumerb)),1]
bphoto <- sumerb[grep(paste("b_photo","\\[",sep=""), rownames(sumerb)),1]
bchill1 <- sumerb[grep(paste("b_chill1","\\[",sep=""), rownames(sumerb)),1]

lwarm <- sumerl[grep(paste("b_warm","\\[",sep=""), rownames(sumerl)),1]
lphoto <- sumerl[grep(paste("b_photo","\\[",sep=""), rownames(sumerl)),1]
lchill1 <- sumerl[grep(paste("b_chill1","\\[",sep=""), rownames(sumerl)),1]

pdf(file.path(figpath, "Sens_vs_day.pdf"), width = 9, height = 7)

par(mfrow=c(2,3))
plot(adv$overallb, bwarm, ylab = "Warming sensitivity", pch = 16, cex = 2, col = alpha("grey20", 0.6), xlab = "Day of budburst")
legend("top", legend="Budburst", text.font=2, inset = 0.05, bty ="n", cex = 2)
plot(adv$overallb, bphoto, ylab = "Photoperiod sensitivity", pch = 16, cex = 2, col = alpha("grey20", 0.6), xlab = "Day of budburst")
plot(adv$overallb, bchill1, #ylim = c(-30, -10), 
     ylab = "Chilling sensitivity", pch = 16, cex = 2, col = alpha("grey20", 0.6), xlab = "Day of budburst")

plot(adv$overall, lwarm, ylab = "Warming sensitivity", pch = 16, cex = 2, col = alpha("grey20", 0.6), xlab = "Day of leafout")
legend("top", legend="Leafout", text.font=2, inset = 0.05, bty ="n", cex = 2)
plot(adv$overall, lphoto, ylab = "Photoperiod sensitivity", pch = 16, cex = 2, col = alpha("grey20", 0.6), xlab = "Day of leafout")
plot(adv$overall, lchill1, #  ylim = c(-30, -10), 
     ylab = "Chilling sensitivity", pch = 16, cex = 2, col = alpha("grey20", 0.6), xlab = "Day of leafout")

dev.off();#system(paste("open", file.path(figpath, "Sens_vs_day.pdf"), "-a /Applications/Preview.app"))

####### Trait pairs plot

panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, #col="darkblue",
       ...) }

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "complete.obs")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  rsig <- cor.test(x, y, use = "complete.obs")$p.value 
  
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  if(rsig <= 0.05) {    text(0.5, 0.5, txt, cex = 1, font=2)}
  else text(0.5, 0.5, txt, cex = 1, font=1)
}

pdf(file.path(figpath, "traitpairs.pdf"), width = 6, height = 6)
pairs(dxt[c("bday","lday","wd","sla","X.N","Pore.anatomy")],
      diag.panel = panel.hist, lower.panel = panel.cor,
      col = hsv(0.7,0.2,0.1,alpha = 0.1), pch = 16,
      labels = c("Budburst day","Leafout day","Stem density", "SLA", "Leaf N","Pore anatomy"),
      cex = 1.5,
      cex.labels = 1, oma = rep(2,4),
      font.labels = 2,
      gap = 0.5
)
dev.off() #; system(paste("open", file.path(figpath, "traitpairs.pdf"), "-a /Applications/Preview.app"))

on.exit(setwd("~/Documents/git/buds/docs/ms/"))