### Started 2 February 2016 ###
### By Lizzie ### 

## Quick look at responses of species we collected tissue for genotyping from ##

## Updated 13 July 2016 by Lizzie to make it run! Before I realized there is another
# file in the budgenetics repo by the same name that runs fine  ...

## So, important note, SEE ALSO file here:
## https://github.com/lizzieinvancouver/budgenetics/tree/master/analyses ##
## before editing or worrying about this! ##

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## packages
library(plyr)
library(dplyr)

setwd("~/Documents/git/projects/treegarden/budburstexp2015/analyses")

# get latest data
print(toload <- sort(dir("./input")[grep("Budburst Data", dir('./input'))], T)[1])

load(file.path("input", toload))

source("source/simpleplot.R")

dx <- dx[!is.na(dx$site),] # one Betpap entry has no site, need to check

# Analysis of where the leafout cuttings were
lx <- dx[dx$nl == 1,]

summary(lx)

##
setwd("~/Documents/git/projects/treegarden/genetics/analyses")


## get the tissue individuals
hftissue <- read.csv("input/HF_TISSUE_data - Tissue individuals, HF.csv", header=TRUE)
shtissue <- read.csv("input/St.Hip_TISSUE_data - Tissue individuals.csv", header=TRUE)

hftissvector <- paste(hftissue$Individual, "HF", sep="_")
shtissvector <- paste(shtissue$Individual, "SH", sep="_")

tissues <- c(hftissvector, shtissvector)

budtiss <- dx[which(dx$ind %in% tissues),]

howmanychill <- aggregate(budtiss["bday"], budtiss[c("chill","sp" )], FUN=length) # ACEPEN, FAGGRA, POPGRA, QUERUB, VIBLAN also have full chilling


foo <- lm(bday~chill+warm+photo, data=subset(budtiss, ind=="ACEPEN01_HF"))

## summarize the budburst data
budsummary.chill <-
      ddply(budtiss, c("ind", "chill"), summarise,
      mean = mean(bday, na.rm=TRUE), sd = sd(bday, na.rm=TRUE),
      sem = sd(bday, na.rm=TRUE)/sqrt(length(bday)))

budsummary.warm <-
      ddply(budtiss, c("ind", "warm"), summarise,
      mean = mean(bday, na.rm=TRUE), sd = sd(bday, na.rm=TRUE),
      sem = sd(bday, na.rm=TRUE)/sqrt(length(bday)))

budsummary.photo <-
      ddply(budtiss, c("ind", "photo"), summarise,
      mean = mean(bday, na.rm=TRUE), sd = sd(bday, na.rm=TRUE),
      sem = sd(bday, na.rm=TRUE)/sqrt(length(bday)))

names(budsummary.chill)[names(budsummary.chill)=="chill"] <- "main.effect"
names(budsummary.warm)[names(budsummary.warm)=="warm"] <- "main.effect"
names(budsummary.photo)[names(budsummary.photo)=="photo"] <- "main.effect"

budsummary <- rbind(budsummary.chill, budsummary.warm, budsummary.photo)

## summarize site effects
whichspp <- unique(budtiss$sp)
siteffects <- c()
for (i in c(1:length(whichspp))){
    sp.model <- lm(bday~site*photo*warm, data=subset(budtiss, sp==whichspp[i]))
    siteffects[i] <- summary(sp.model)$coef[2,1] # site effect
}
# SPIALB has strong site effect
siteffects.out <- data.frame(spp=unique(budtiss$sp), siteeffect=siteffects)

write.csv(budsummary, "output/budsummary.csv", row.names=FALSE)
write.csv(siteffects.out, "output/siteffects.csv", row.names=FALSE)
