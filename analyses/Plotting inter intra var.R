# Plotting inter vs intraspecific variation using individual models

library(ggplot2)
library(rstan)
library(shinystan)
library(scales)

setwd("~/Documents/git/treetraits/analyses")

traitout <- dir()[grep("Stan Output", dir())[!is.na(match(grep("Stan Output", dir()), grep("Trait", dir())))]]

phenout <- dir()[grep("Stan Output", dir())[!is.na(match(grep("Stan Output", dir()), grep("Ind", dir())))]]


if(!exists("doym.s")) load(sort(traitout, T)[1]) # only run if stan output file is not yet in working memory.
if(!exists("doym.b")) load(sort(phenout, T)[1]) 

ls() 
#launch_shinystan(doym.l)


############### 2 panel figure for inter-intra
pdf("Inter Intra pheno trait plots.pdf", 
    width = 10, height = 14)
par(mfrow=c(2,1), mar = c(3,4,1,1))

################# Inter intra phenology plots


sum.s <- summary(doym.l)$summary
# intercepts
sp.a <- sum.s[grep("^a_sp\\[", rownames(sum.s)),]
# Intraspecific
ind.a <- sum.s[grep("^a_sp_ind\\[", rownames(sum.s)),]

denx <- density(sp.a[,"mean"])

par(mfrow=c(2,1), mar = c(3,4,1,1))

plot(denx,
     main = "Distribution of intercepts for Leafout model",
     ylim = c(0, 0.3)
)

# now make polygons for individuals for each species. 
splookup <- unique(dxl[c("ind","spn")])[,"spn"] 
ccolz <- heat.colors(length(unique(splookup)), alpha = 0.3)

for(i in unique(splookup)){
  dxx <- density(ind.a[splookup == i,"mean"])
  polygon(dxx$x, dxx$y, col=ccolz[i], 
          border = ccolz[i])
}

# Species-level distribution
polygon(denx$x, denx$y, col=alpha("purple", 0.5))

legend("topleft",
       bty="n",
       fill = c(alpha("purple",0.5), ccolz[32]),
       legend = c("Species-level intercepts", 
                  "Individual-level intercepts"))

# Plot inter vs intra.
plot(range(0:nrow(sp.a)),
     c(min(sp.a[,"mean"])*.95, max(sp.a[,"mean"])*1.05), type = "n",
     ylab= "Intercept value")

plotorder <- (1:max(splookup))[order(sp.a[,"mean"])]
counter = 1

for(i in plotorder){
  indx <- ind.a[splookup == i,"mean"]
  points(rep(counter, length(indx)),
         indx, pch = 16, col = alpha('grey50',0.4))
  points(counter, sp.a[i,"mean"], pch = 16)
  counter = counter + 1
}
title(main="Ordered by species")
legend("topleft",
      bty="n",
      col = c("grey80","black"),
      pch = 16,
      legend = c("Species-level intercepts", 
                  "Individual-level intercepts"))

# Interspecific. 
sum.s <- summary(latm.s)$summary

# intercepts
sp.a <- sum.s[grep("^a_sp\\[", rownames(sum.s)),]

# Intraspecific
ind.a <- sum.s[grep("^a_sp_ind\\[", rownames(sum.s)),]
dx <- density(sp.a[,"mean"])

plot(dx,
     main = "Distribution of intercepts for SLA model",
     ylim = c(0, 0.8)
     )

# now make polygons for individuals for each species

# lines(density(ind.a[,"mean"]))
splookup <- unique(dts[c("ind","spn")])[,"spn"] 
#lines(density(ind.a[,"mean"]))
ccolz <- heat.colors(length(unique(splookup)), alpha = 0.3)

for(i in unique(splookup)){
  dxx <- density(ind.a[splookup == i,"mean"])
  polygon(dxx$x, dxx$y, col=ccolz[i], 
          border = ccolz[i])
  
  }

# Species-level distribution
polygon(dx$x, dx$y, col=alpha("purple", 0.5))

legend("topleft",
       bty="n",
       fill = c(alpha("purple",0.5), ccolz[32]),
       legend = c("Species-level intercepts", 
                  "Individual-level intercepts"))


# Plot inter vs intra. 
plot(range(0:nrow(sp.a)),
  range(sp.a[,"mean"]), type = "n",
  ylab= "Intercept value")

plotorder <- (1:max(splookup))[order(sp.a[,"mean"])]
counter = 1

for(i in plotorder){
  indx <- ind.a[splookup == i,"mean"]
  points(rep(counter, length(indx)),
       indx, pch = 16, col = alpha('grey50',0.4))
  points(counter, sp.a[i,"mean"], pch = 16)
  counter = counter + 1
}
title(main="Ordered by species")

############### 
### Wood density
sum.s <- summary(latm.w)$summary

# intercepts
sp.a <- sum.s[grep("^a_sp\\[", rownames(sum.s)),]
# Intraspecific
ind.a <- sum.s[grep("^a_sp_ind\\[", rownames(sum.s)),]

dx <- density(sp.a[,"mean"])

par(mfrow=c(2,1), mar = c(3,4,1,1))

plot(dx,
     main = "Distribution of intercepts for Wood Density model",
     ylim = c(0, 1.5)
)

# now make polygons for individuals for each species
splookup <- unique(dtw[c("ind","spn")])[,"spn"] 
ccolz <- heat.colors(length(unique(splookup)), alpha = 0.3)

for(i in unique(splookup)){
  dxx <- density(ind.a[splookup == i,"mean"])
  polygon(dxx$x, dxx$y, col=ccolz[i], 
          border = ccolz[i])
  
}

# Species-level distribution
polygon(dx$x, dx$y, col=alpha("purple", 0.5))

legend("topleft",
       bty="n",
       fill = c(alpha("purple",0.5), ccolz[32]),
       legend = c("Species-level intercepts", 
                  "Individual-level intercepts"))


# Plot inter vs intra.
plot(range(0:nrow(sp.a)),
     c(min(sp.a[,"mean"])*0.95, max(sp.a[,"mean"])*1.05), type = "n",
     ylab= "Intercept value")

plotorder <- (1:max(splookup))[order(sp.a[,"mean"])]
counter = 1

for(i in plotorder){
  indx <- ind.a[splookup == i,"mean"]
  points(rep(counter, length(indx)),
         indx, pch = 16, col = alpha('grey50',0.4))
  points(counter, sp.a[i,"mean"], pch = 16)
  counter = counter + 1
}
title(main="Ordered by species")

################################
###### Height
sum.s <- summary(latm.h)$summary
# intercepts
sp.a <- sum.s[grep("^a_sp\\[", rownames(sum.s)),]
# Intraspecific
ind.a <- sum.s[grep("^a_sp_ind\\[", rownames(sum.s)),]

dx <- density(sp.a[,"mean"])

par(mfrow=c(2,1), mar = c(3,4,1,1))

plot(dx,
     main = "Distribution of intercepts for Height model",
     ylim = c(0, 1.5)
)

# now make polygons for individuals for each species. # !!!! change dts/dtw/dth/dtd as necessary for lookup
splookup <- unique(dth[c("ind","spn")])[,"spn"] 
ccolz <- heat.colors(length(unique(splookup)), alpha = 0.3)

for(i in unique(splookup)){
  dxx <- density(ind.a[splookup == i,"mean"])
  polygon(dxx$x, dxx$y, col=ccolz[i], 
          border = ccolz[i])
  }

# Species-level distribution
polygon(dx$x, dx$y, col=alpha("purple", 0.5))

legend("topleft",
       bty="n",
       fill = c(alpha("purple",0.5), ccolz[32]),
       legend = c("Species-level intercepts", 
                  "Individual-level intercepts"))

# Plot inter vs intra.
plot(range(0:nrow(sp.a)),
     c(min(sp.a[,"mean"])*0.95, max(sp.a[,"mean"])*1.05), type = "n",
     ylab= "Intercept value")

plotorder <- (1:max(splookup))[order(sp.a[,"mean"])]
counter = 1

for(i in plotorder){
  indx <- ind.a[splookup == i,"mean"]
  points(rep(counter, length(indx)),
         indx, pch = 16, col = alpha('grey50',0.4))
  points(counter, sp.a[i,"mean"], pch = 16)
  counter = counter + 1
}
title(main="Ordered by species")

################################
###### DBH
sum.s <- summary(latm.d)$summary
# intercepts
sp.a <- sum.s[grep("^a_sp\\[", rownames(sum.s)),]
# Intraspecific
ind.a <- sum.s[grep("^a_sp_ind\\[", rownames(sum.s)),]

dx <- density(sp.a[,"mean"])

par(mfrow=c(2,1), mar = c(3,4,1,1))

plot(dx,
     main = "Distribution of intercepts for DBH model",
     ylim = c(0, 1.5)
)

# now make polygons for individuals for each species. # !!!! change dts/dtw/dth/dtd as necessary for lookup
splookup <- unique(dtd[c("ind","spn")])[,"spn"] 
ccolz <- heat.colors(length(unique(splookup)), alpha = 0.3)

for(i in unique(splookup)){
  dxx <- density(ind.a[splookup == i,"mean"])
  polygon(dxx$x, dxx$y, col=ccolz[i], 
          border = ccolz[i])
}

# Species-level distribution
polygon(dx$x, dx$y, col=alpha("purple", 0.5))

legend("topleft",
       bty="n",
       fill = c(alpha("purple",0.5), ccolz[32]),
       legend = c("Species-level intercepts", 
                  "Individual-level intercepts"))

# Plot inter vs intra.
plot(range(0:nrow(sp.a)),
     c(min(sp.a[,"mean"])*1.05, max(sp.a[,"mean"])*1.05), type = "n",
     ylab= "Intercept value")

plotorder <- (1:max(splookup))[order(sp.a[,"mean"])]
counter = 1

for(i in plotorder){
  indx <- ind.a[splookup == i,"mean"]
  points(rep(counter, length(indx)),
         indx, pch = 16, col = alpha('grey50',0.4))
  points(counter, sp.a[i,"mean"], pch = 16)
  counter = counter + 1
}
title(main="Ordered by species")

dev.off()



