# Compiled tree traits
setwd("~/Documents/git/treetraits/")
d <- read.csv("Summer 2015 Tree Traits.csv")

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

# Plotting

plot(sla ~ height, data = d)

# Four sites

library(ggplot2)

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

<<<<<<< Updated upstream
# Species per site
xtable(apply(dt, 2, function(x) length(x[x>0])))
=======
## from buds

(leafcount <- with(d, tapply(Leaf.area, list(Species, Site), function(x) length(x[!is.na(x)]))))

rowSums(leafcount, na.rm=T)
colSums(leafcount, na.rm=T)

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# 3. Individual level
# look at consistency of performance within individuals, across treatments, as measure of plasticity.

# Is consistency related to earlier leafout?

vars <- aggregate(lday ~ sp + site + ind + wd + sla + X.N + Pore.anatomy, FUN = function(x) sd(x, na.rm=T) 
                  / mean(x, na.rm=T)
                  , data = dxt)

# remove extreme values 

vars$day = lday.agg[match(vars$sp, lday.agg$sp),"lday"]
vars$site = as.factor(vars$site); levels(vars$site) = c("HF","SH")
summary(lmer(lday ~ day + (1|sp), data = vars))

pdf(file.path(figpath, "indvar.pdf"), width = 5, height = 5)

ggplot(vars, aes(day, lday, group = site)) + geom_point(aes(col=site)) + geom_smooth(method = "lm") + 
  xlab("Day of leafout in Warm/Long") + ylab("CV of leafout across treatments within individuals") 

dev.off();system(paste("open", file.path(figpath, "indvar.pdf"), "-a /Applications/Preview.app"))
>>>>>>> Stashed changes
