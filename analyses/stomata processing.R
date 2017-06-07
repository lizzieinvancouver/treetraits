# Processing stomata data for common garden

# For each individual, up to 3 replicates.
# Each image has five measured random stomates in focus with lengths. This is the straight line distance across the stomate.

# Area column: this 
# Scale: pixels per micrometer
# Stom_number: number of stomates in this field of view

setwd("~/Documents/git/treetraits/analyses/data")

d <- read.csv("Common Garden Stomata.csv")

# Goals: get mean values individual; mean length and stomatal density. Density is number per mm^2. Normal values ~  50 - 500.
# Area: 585632.274 um^2 = 0.585632274 mm^2

d$stom_dens <- d$Stom_number / d$Area * 1000000 

avglen <- apply(d[,2:6], 1, mean)

id <- substr(d$Plant.ID, 1, 11)

dx <- data.frame(stom_leng = tapply(avglen, id, mean))

dx$stom_dens <- tapply(d$stom_dens, id, mean)

write.csv(dx, "Stomata_Processed.csv", row.names = T)
