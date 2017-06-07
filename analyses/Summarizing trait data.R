# Compiled tree traits

library(ggplot2)
library(xtable)

setwd("~/Documents/git/treetraits/")
d <- read.csv("analyses/input/Summer 2015 Tree Traits.csv")

length(unique(d$Species))


dt <- with(d, table(Species, Site))

# Species we have at only one site

atonesite <- apply(dt, 1, function(x) length(x[x==0])==3)

atonesite[atonesite==TRUE]

# At all four sites

atfoursite <- apply(dt, 1, function(x) length(x[x==0])==0)

length(atfoursite[atfoursite==TRUE])

# Calcs

d$sla <- d$Leaf.area/d$Dry.mass
d$wd <- d$Stem.mass/d$Stem.volume

d$height = d$Height
d$height = as.numeric(sub("<1",0.5, as.character(d$height)))

# Clean DBH, to get a sum DBH across all stems

dbh <- vector() # i = 56 and 57 test
for(i in 1:nrow(d)){
  dx <- d[i, grep("DBH", names(d))]
  dsplit <- unlist(apply(dx, 1, function(x) strsplit(x, ",")))
  dbh <- c(dbh, sum(as.numeric(sub("<1", 0.5, dsplit))))
  }

d$DBH = dbh

# Plotting

plot(sla ~ height, data = d)

# Four sites

ggplot(d, aes(DBH, height, color = Site)) + geom_point()

ggplot(d, aes(wd, sla, color = Site)) + geom_point()

ggplot(d, aes(log(wd), log(sla), color = Site)) + geom_point() +  geom_smooth(method = "lm") 

ggplot(d, aes(wd, sla, color = Site)) + geom_point() +  geom_smooth(method = "lm") 

# Remove some outliers
mean(d$wd, na.rm=T)+sd(d$wd, na.rm=T)*2
mean(d$sla, na.rm=T)+sd(d$sla, na.rm=T)*2

d = d[d$sla <= 500 ,]
d = d[d$wd <= 1,]

ggplot(d, aes(wd, sla, color = Site)) + geom_point() +  geom_smooth(method = "lm") 


##### Summary tables

# Samples by species and site
xtable(dt)

# Total samples per species
xtable(data.frame(sort(rowSums(dt),T)))

# Total samples per site
xtable(data.frame(sort(colSums(dt),T)))

## from buds

(leafcount <- with(d, tapply(Leaf.area, list(Species, Site), function(x) length(x[!is.na(x)]))))

rowSums(leafcount, na.rm=T)
colSums(leafcount, na.rm=T)

